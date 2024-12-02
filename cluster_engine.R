run_batch <- function(
  experiments,           # Which experiments to run
  batch_replications,    # How many replications per experiment?
  job_function,          # This function is run on the workers after replication setup
  worker_setup_function = default_worker_setup(), # This function is run on the workers on worker startup
  job_wrapper           = default_job_wrapper,    # This function is run on the workers to set up each replication
  test_single_job_index = integer(),      # Use this to test single jobs (e.g. if job_id 11 crashes, set this to 11)
  test_single_job_per_experiment = FALSE, # Use to test 1 replication per experiment
  test_single_job_through_clustermq = FALSE,
  test_single_job_externally = FALSE,     # By default test locally, setting this to TRUE tests on a worker
  job_resources = list( # Each worker will request these resources
    walltime = 180, # Wall time in minutes
    memory = 3072,  # Memory in megabytes
    cores = 1       # Number of CPU cores per worker (R is single-threaded by default, so default to 1)
  ),
  n_jobs = 200,          # The maximum simultaneous jobs 
  timeout = 120,         # Job timeout in seconds
  log_worker = FALSE,    # Set this to true to create logs 
  max_calls_worker = Inf # This option can limit the number of jobs run per worker (useful for debugging random segfaults)
)
{
  if (test_single_job_per_experiment) 
  {
    assert_that(length(test_single_job_index) == 0)
  }
  
  options(clustermq.data.warning = 999)
  
  cat(paste0('Clustermq run - enumerating clustermq experiments...\n'))
  
  clustermq_args <- list()
  
  replication_index <- integer()
  experiment_vec <- character()
  seed_vec <- list()
  
  max_job_id <- 0L
  for (experiment_name in names(experiments))
  {
    experiment_args <- experiments[[experiment_name]]
    
    clustermq_args[[experiment_name]] <- experiment_args
    
    experiment_seeds <- setup_lecuyer_seeds(seed = experiment_args[["seed"]], num = batch_replications)
    
    if (test_single_job_per_experiment) {
      test_single_job_index <- c(test_single_job_index, max_job_id + 1)
    }
    
    max_job_id <- max_job_id + batch_replications
    replication_index <- c(replication_index, seq_len(batch_replications))
    experiment_vec <- c(experiment_vec, rep(experiment_name, batch_replications))
    
    seed_vec <- c(seed_vec, experiment_seeds)
    
  }
  
  clustermq_jobs <- tibble(
    job_id = seq_len(max_job_id), 
    replication_index = replication_index,
    experiment = experiment_vec,
    seed = seed_vec
  )
  
  Q_const_arg <- list(
    "job_function"     = job_function, 
    "setup_function"   = worker_setup_function,
    "num_cores"        = job_resources$cores,
    "experiment_args"  = clustermq_args
  )
  
  cat(paste0('Clustermq run - submitting and running experiments...\n'))
  cat(paste0("Clustermq run - submitting jobs at ", Sys.time(), "\n"))
  
  if (length(test_single_job_index) > 0) { # Set this for local testing
    outputs <- list()
    
    current_scheduler <- getOption("clustermq.scheduler")
    all_jobs_finished <- TRUE
    
    assert_that(
      all(test_single_job_index > 0 & test_single_job_index <= nrow(clustermq_jobs)),
      msg="Test job indices out of range"
    )
    
    worker_setup_complete <- FALSE
    
    for (job_index in test_single_job_index) {
      current_job <- clustermq_jobs[job_index, ]
      
      
      # current_job <- as.list(current_job)
      # current_job$seed <- current_job$seed[[1]]
      # 
      # # worker_setup_complete <- FALSE
      # outputs <- c(
      #   outputs,
      #   do.call("job_wrapper", args = c(current_job, Q_const_arg, list("capture_stack" = FALSE )))
      # 
      # 
      # )
      if (!test_single_job_through_clustermq)
      {
        current_job <- as.list(current_job)
        current_job$seed <- current_job$seed[[1]]
        
        outputs[[job_index]] <- do.call("job_wrapper", args = c(current_job, Q_const_arg, list("capture_stack" = FALSE )))
      } else {
        
      
      tryCatch(
        expr = {
          # if (!test_single_job_through_clustermq)
          # {
          #   current_job <- as.list(current_job)
          #   current_job$seed <- current_job$seed[[1]]
          # 
          #   outputs[[job_index]] <- do.call("job_wrapper", args = c(current_job, Q_const_arg, list("capture_stack" = FALSE )))
          # } else {
            if (!test_single_job_externally) {
              options("clustermq.scheduler" = "local")
            }
            outputs[[job_index]] <- 
              Q_rows(
                current_job,
                fun = job_wrapper,
                const = Q_const_arg,
                template = job_resources,
                n_jobs = 1,
                max_calls_worker = max_calls_worker,
                fail_on_error = TRUE,
                timeout = timeout,
                log_worker = log_worker
              )[[1]]
          # }
        },
        error=function(cond) {
          message(paste0("Job #", job_index, " failed."))
          message("Here's the original error message:")
          message(cond)
          # Choose a return value in case of error

          stop("At least one job failed.")
        },
        finally = {
          options("clustermq.scheduler" = current_scheduler)
        }
      )
        }
      
      errored_jobs <- sapply(outputs, \(x) {is(x, "simpleError")})
      if (any(errored_jobs)) {
        
        index <- which(errored_jobs)
        
        message("Error in job_id ", clustermq_jobs[job_index, ]$job_id, ": ", outputs[[index]]$message)
        message("Call stack:")
        for (i in seq_along(outputs[[index]]$calls)) {
          message(i, ": ", outputs[[index]]$calls[i])
        }
        cat("\n")
        
        stop("At least one job failed.")
      }
      
      cat("\n")
    }
    
    
  } else {
    n_jobs_final <- max(min(nrow(clustermq_jobs), n_jobs), ceiling(nrow(clustermq_jobs) / max_calls_worker))
    outputs <- Q_rows(
      clustermq_jobs,
      fun = job_wrapper,
      const = Q_const_arg,
      template = job_resources,
      n_jobs = n_jobs_final,
      max_calls_worker = max_calls_worker,
      fail_on_error = TRUE,
      timeout = timeout,
      log_worker = log_worker
    )
    
    errored_jobs <- sapply(outputs, \(x) {is(x, "simpleError")}) |
                    sapply(outputs, \(x) {is(x, "error")})
    if (any(errored_jobs)) {
      errored_job_indices <- which(errored_jobs)
      for (index in errored_job_indices) {
        message("Error in job_id ", clustermq_jobs[index, ]$job_id, " with call stack" )
        for (i in seq_along(outputs[[index]]$calls)) {
          message(i, ": ", outputs[[index]]$calls[i])
        }
        cat("\n")
      }
      browser()
      stop("At least one job failed.")
    }
  }
  
  cat(paste0("Clustermq run - finished waiting for jobs at ", Sys.time(), "\n"))
    
  cat(paste0('Clustermq run - reducing results...  \n'))
  cat(paste0("Clustermq run - starting reduce at ", Sys.time(), "\n"))
  
  
  
  trial_tibble <- clustermq_jobs %>% inner_join(
    bind_rows(lapply(outputs, function(x) x[["trial_tibble"]])),
    by = "job_id"
  )

  dose_tibble <- clustermq_jobs %>% inner_join(
    bind_rows(lapply(outputs, function(x) x[["dose_tibble"]])),
    by = "job_id"
  )
  
  result <- list(
    outputs = outputs, 
    trial_tibble = trial_tibble,
    dose_tibble = dose_tibble,
    outcome_overall = reduce_results_per_scenario(trial_tibble)
    # ,MTD_overview = reduce_results_to_MTD_overview(trial_tibble)
  )
  
  cat(paste0("Clustermq run - finished reduce at ", Sys.time(), "\n"))
  cat(paste0('Clustermq run - reducing results... done\n'))
  return(result)
}


default_worker_setup <- function(
  working_dir = getwd(),                   # Working directory
  additional_source_files = character(0L), # Source these additional files
  lib_paths = .libPaths(),                 # Use these .libPaths() - by default same as on master node
  transfer_options = FALSE,                # Transfer options() from master node to workers
  source_files = function(working_dir) {
    return(
      c(
        file.path(working_dir, "load_packages.R"),
        file.path(working_dir, "load_functions.R"),
        file.path(working_dir, "cluster_engine.R")
      ) 
    )
  },
  additional_setup = function(working_dir) {}
) {
  # Ensure arguments are evaluated so the function that this function returns has these arguments available
  working_dir
  lib_paths
  transfer_options
  source_files
  additional_source_files
  additional_setup
  
  if (transfer_options) {
    master_options <- options()
  }
  
  f <- function () {
    # This is the function that will get run on worker startup
    if (transfer_options) {
      options(master_options)
    }
    
    if (!is.null(lib_paths)) {
      .libPaths(lib_paths)
    }
    
    all_source_files <- c(
      source_files(working_dir),
      additional_source_files
    )
    
    if (length(all_source_files) > 0) {
      for (file_name in all_source_files) {
        source(file_name)
      }
    }
    
    if (transfer_options) {
      options(master_options)
    }
    
    additional_setup(working_dir = working_dir)
  }
  
  return(f)
}

  


tryStack <- function(
  expr,
  silent=FALSE
)
{
  tryenv <- new.env()
  out <- try(withCallingHandlers(expr, error=function(e)
  {
    stack <- sys.calls()
    stack <- stack[-(2:7)]
    stack <- head(stack, -2)
    stack <- sapply(stack, deparse)
    if(!silent && isTRUE(getOption("show.error.messages"))) 
      cat("This is the error stack: ", stack, sep="\n")
    assign("stackmsg", value=paste(stack,collapse="\n"), envir=tryenv)
  }), silent=silent)
  if(inherits(out, "try-error")) out[2] <- tryenv$stackmsg
  out
}

default_job_wrapper <- function(
  setup_function,
  job_function,
  job_id,
  replication_index,
  experiment,
  experiment_args,
  seed,
  num_cores,
  capture_stack = TRUE
) 
{
  cat("Current job: ", job_id, "\n")
  
  if(!exists('worker_setup_complete', mode='logical')) {
    cat("Calling setup function:\n")
    setup_function()

    worker_setup_complete <<- TRUE
  } else {
    cat("Calling GC:\n")
    print(gc())
  }
  
  RNGkind("L'Ecuyer-CMRG")
  print(seed)
  .Random.seed <<- seed
  
  # Use as many cores as specified
  options(mc.cores = num_cores)
  
  if (capture_stack) {
    return(
      evaluate::try_capture_stack({
        job_function(
          experiment = experiment,
          experiment_args = experiment_args[[experiment]],
          job_id = job_id,
          replication_index = replication_index
        )
      }, env=environment())
    )
  } else {
    job_function(
      experiment = experiment,
      experiment_args = experiment_args[[experiment]],
      job_id = job_id,
      replication_index = replication_index
    )
  }

}


setup_lecuyer_seeds <- function(seed = NULL, lecuyer_seed = NULL, num) {
  RNGkind("L'Ecuyer-CMRG")
  
  assert_that(is.null(lecuyer_seed) + is.null(seed) == 1)
  
  if (!is.null(seed)) {
    set.seed(seed)
    lecuyer_seed <- .Random.seed
  }
  
  job_seeds <- list()
  job_seeds[[1]] <- parallel::nextRNGStream(lecuyer_seed)
  i <- 2
  while(i < num+1) {
    job_seeds[[i]] <- parallel::nextRNGStream(job_seeds[[i-1]])
    i  <- i + 1
  }
  job_seeds
}