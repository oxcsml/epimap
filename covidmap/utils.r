library(gsubfn)


#' Compute date range used for modelling from specification.
#'
#' Requires specifying at least two of (\code{first_day_modelled}, 
#' \code{weeks_modelled}, and one of (\code{last_day_modelled}, \code{days_ignored})).
#'
#' @param alldates List of dates for which data is available.
#' @param first_day_modelled First day modelled (defaults to NULL)
#' @param last_day_modelled Last day modelled (defaults to NULL)
#' @param days_ignored Days ignored (defaults to 7)
#' @param weeks_modelled Weeks modelled (defaults to NULL)
#' @param days_per_step Days per step (defaults to 7)
#' @return Named list consisting of (Nstep, Tstep, Tcond, Tcur, Tignore)
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
  #print(paste("last_day_of_data",last_day_of_data))
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
    #print(paste("last_day_modelled",last_day_modelled))
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
    } else {
      weeks = floor((last_day_modelled-first_day_modelled)/7)
      if (!is.null(weeks_modelled) && weeks != weeks_modelled) {
        stop(paste("First/last day modelled and weeks modelled inconsistent", list_args, sep=""))
      }
      first_day_modelled = last_day_modelled - 7*weeks + 1
      #print(paste("first_day_modelled",first_day_modelled))
    }
  } else if (!is.null(last_day_modelled) && !is.null(opt$weeks_modelled)) {
    first_day_modelled = last_day_modelled - 7*opt$weeks_modelled + 1
    #print(paste("first_day_modelled",first_day_modelled))
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
  #Tlik = floor(days_modelled/Tstep)*Tstep
  Tlik = floor(days_modelled/7)*7 # model in units of weeks
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

time_distances = function(
    time_steps,
    Tstep,
    days,
    lockdown_day=as.Date("2020-03-23")
) {
    # time steps
    times = 1:(time_steps)
    timedist = matrix(0, time_steps, time_steps)
    for (i in 1:(time_steps)) {
      for (j in 1:(time_steps)) {
        timedist[i, j] = abs(times[i] - times[j]) * Tstep
      }
    }

    # precompute lockdown cutoff kernel
    lockdown_day = as.Date("2020-03-23")
    days_period_start = days[seq(1, length(days), Tstep)]
    days_period_start = vapply(days_period_start, (function (day) as.Date(day, format="%Y-%m-%d")), double(1))
    day_pre_lockdown = vapply(days_period_start, (function (day) day < lockdown_day), logical(1))
    
    time_corellation_cutoff = matrix(0,time_steps,time_steps)
    for (i in 1:(time_steps)) {
      for (j in 1:(time_steps)) {
        time_corellation_cutoff[i, j] = !xor(day_pre_lockdown[i], day_pre_lockdown[j])
      }
    }

    list(
        time_distances = timedist,
        time_correlations = time_corellation_cutoff
    )
}