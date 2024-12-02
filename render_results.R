source("load_packages.R")

if (!exists("time_stamp")) {
  quarto_time_stamp <- "2024-11-04-16-04-42.547127"
} else {
  quarto_time_stamp <- time_stamp
}

quarto_render(
  "render_results.qmd", 
  execute_params = list("path.export" = file.path("results", quarto_time_stamp)),
  output_file = paste0("results-", quarto_time_stamp, ".html")
)
