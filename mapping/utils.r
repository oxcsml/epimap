library(gsubfn)


process_dates_modelled = function(
  alldates,
  first_day_modelled = NULL,
  last_day_modelled = NULL,
  days_ignored = 7,
  weeks_modelled = NULL,
  days_per_step = 7
) {
  list_args = paste(
    ", first_day_modelled=", first_day_modelled,
    ", weeks_modelled=", weeks_modelled,
    ", last_day_modelled=", last_day_modelled,
    ", days_ignored=", days_ignored,
    ", days_per_step=", days_per_step,
    sep=""
  )
  Tall = length(alldates)
  first_day_of_data = alldates[1]
  last_day_of_data = alldates[Tall]
  if (!is.null(last_day_modelled)) {
    last_day_modelled = as.Date(last_day_modelled)
    if (!is.null(days_ignored) &&
        as.double(last_day_of_data-last_day_modelled)<days_ignored) {
      stop(paste(
        "Last day modelled is less than ", days_ignored, 
        " before last day of data", list_args,
        sep=""
      ))
    }
  } else if (!is.null(days_ignored)) {
    last_day_modelled = last_day_of_data - days_ignored
  } else {
    last_day_modelled = NULL
  }

  if (!is.null(first_day_modelled)) {
    first_day_modelled = as.Date(first_day_modelled)
    if (is.null(last_day_modelled)) {
      if (is.null(weeks_modelled)) {
        stop(paste("Cannot work out last day modelled", list_args, sep=""))
      } else {
        last_day_modelled = first_day_modelled + 7*opt$weeks_modelled - 1
      }
    }
  } else if (!is.null(last_day_modelled) && !is.null(opt$weeks_modelled)) {
    first_day_modelled = last_day_modelled - 7*opt$weeks_modelled + 1
    if (first_day_modelled < first_day_of_data) {
      stop(paste(
          "First day modelled is before first data of data", list_args,
          sep=""
      ))
    }
  } else {
    stop(paste(
      "Cannot work out first day modelled", list_args,
      sep=""
    ))
  }
  days_modelled = as.double(last_day_modelled-first_day_modelled)+1
  Tstep <- days_per_step
  Tlik = floor(days_modelled/Tstep)*Tstep
  if (Tlik <= 0) {
    stop(paste(
      "Number of days modelled is <= 0", list_args,
      sep=""
    ))
  }
  Nstep = Tlik / Tstep
  Tcond = sum(alldates<first_day_modelled)
  stopifnot(first_day_modelled==alldates[Tcond+1])
  Tcur = Tcond + Tlik
  Tignore = Tall - Tcur
  message(paste("Tstep =",Tstep))
  message(paste("Nstep =",Nstep))
  message(paste("Days modelled =",Tlik))
  message(paste("Days ignored =",Tignore))
  message(paste("First day modelled =",alldates[Tcond+1]))
  message(paste("Last day modelled =",alldates[Tcur]))

  return(list(
    Nstep = Nstep,
    Tstep = Tstep,
    Tcond = Tcond,
    Tlik = Tlik,
    Tcur = Tcur,
    Tignore = Tignore
  ))
    
}
