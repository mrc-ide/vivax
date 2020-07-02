#' @title
#' Fitting eir.
#' @description
#' Provides basic functionality for fitting the EIR_equil parameter in a site file
#' to specified data (usually prevalence or incidence) at given time point(s).

#' @param variable The name of the output variable to be fitted (eg: "prev_2_10" or "clin_inc_all"). This must be
#' a variable name in the output.
#' @param target The target value(s) of the specified variable. These may be a single value or vector.
#' @param rows The output row(s)/position(s). These may be a single value or vector.
#' @param tolerance The desired tolerance of the fit. The fit being succesful when the \code{abs(variable - target) < tolerance)}.
#' @param interval The interval to be searched (see \link[stats]{uniroot} for more information, note: \code{extendInt = 'upX'}).
#' @param maxiter The maximum number of iterations (model runs) to search for. Importantn for keeping timescales reasonable
#' when fitting many sites.
#' @param verbose Output is printed during each iteration, including the current total_M and the output target value.
#' @param extendInt Extend the search interval, see \link[stats]{uniroot} for more details
#' @param check.conv Error or warning for convergence failure, see \link[stats]{uniroot} for more details
#' @param ... Additional arguments to \code{fit_v()}
#' @inheritParams run_simulation
#'
#' @export
fit_mv <- function (variable = "Pv_clin_2_10", target, rows, model = NULL, 
                    tolerance = 0.05, interval = c(1e-04, 10), maxiter = 10, 
                    verbose = FALSE, extendInt = "no", check.conv = FALSE, 
                    ...) 
{
  stopifnot(all(target >= 0), tolerance > 0, maxiter > 0)
  if (verbose) {
    message("~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~")
    message("Fitting ", variable, " = ", paste(target, 
                                               collapse = " "), " in row(s) ", paste(rows, 
                                                                                     collapse = " "), ".")
  }
  if (all(target == 0)) {
    new_M <- 0
  }
  else {
    fit <- stats::uniroot(fit_v, interval = interval, target = target, 
                          variable = variable, model = model, rows = rows, 
                          tolerance = tolerance, extendInt = extendInt, check.conv = check.conv, 
                          maxiter = maxiter, verbose = verbose, ...)
    new_M <- round(fit$root, 5)
  }
  message("Fitting, complete. Total_M = ", new_M, " .")
  return(new_M)
}

#' @title
#' Internal: Fitting EI.
#' @description
#' Internal function provided to \link[stats]{uniroot} for fitting eir
#'
#' @param total_M The eir for the current model run
#' @inheritParams fit_mv
#' @param ... Additional arguments to \code{fit_v()}
fit_v <- function (total_M, target, variable, rows, tolerance, verbose, 
                   model, ...) 
{
  model <- c(model, list(EIR_equil = total_M / 365))
  temp <- suppressMessages(run_simulation(model = model, ...))
  op <- tapply(temp[[variable]] / temp$N_pop, rep(1:(length(temp$time) / 365), each = 365), sum)
  if (verbose) {
    message("Current total_m = ", round(total_M, 8))
    message("Current iteration ", variable, " = ", 
            paste(round(op[rows], 4), collapse = " "), 
            ".")
  }
  if (sum(abs(op[rows] - target)) < tolerance & 
      all(op[rows] > 0)) {
    return(0)
  }
  else {
    if (any((op[rows] == 0) & (target != 
                               0))) {
      message("Zero transmission")
      return((total_M - 100) * 100)
    }
    message("Current absolute difference = ", sum(abs(op[rows] - target)))
    return(sum(op[rows] - target))
  }
  if (verbose) {
    message("~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~")
  }
}