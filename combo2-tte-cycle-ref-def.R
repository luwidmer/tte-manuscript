## We define time as follows here:
##
## - time-origin is the first patient, first visit date. Is defined within simulation code (not needed here yet)
## - calendar time is measured in units of a cycle length and is wrt to the time-origin.
## - for the analysis data set, the clock is always in relation to the accrual time of
##   the patient, which is used as date for the treatment start.
##

## cycle length in weeks

cycle_length <- weeks(6)
cycle_length_d <- as.numeric(cycle_length, "days")
cycle_length_w <- as.numeric(cycle_length, "weeks")
cycle_length_h <- as.numeric(cycle_length, "hours")

## pseudo date used to mark distant / unkown future time-point
inf_date      <- as_date(Inf)
inf_datetime  <- as_datetime(Inf)


## Reference doses
drugA_dose_ref <- 1
drugB_dose_ref <- 160