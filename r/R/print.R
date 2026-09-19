#' @export
print.di_estimate <- function(x, ...) {
  cat("<debiasedInference estimate>\n")
  cat(" method:   ", x$method, "\n", sep = "")
  cat(" bandwidth:", format(x$bandwidth), "\n")
  cat(" tau:       ", format(x$tau), "\n", sep = "")
  cat(" grid size: ", length(x$estimate), "\n", sep = "")
  invisible(x)
}

#' @export
print.di_band <- function(x, ...) {
  cat("<debiasedInference simultaneous confidence band>\n")
  cat(" method:    ", x$method, "\n", sep = "")
  cat(" confidence:", format(x$confidence), "\n")
  cat(" bootstrap: ", x$n_boot, " replicates\n", sep = "")
  cat(" bandwidth: ", format(x$bandwidth), "\n", sep = "")
  invisible(x)
}

#' @export
print.di_set <- function(x, ...) {
  cat("<debiasedInference level set>\n")
  cat(" method: ", x$method, "\n", sep = "")
  cat(" level:  ", format(x$level), "\n", sep = "")
  root_count <- if (is.null(dim(x$roots))) length(x$roots) else nrow(x$roots)
  cat(" roots:  ", root_count, "\n", sep = "")
  invisible(x)
}
