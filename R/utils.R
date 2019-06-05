#' Check the category of fault data
#'
#' This is a function to decide whether a given fault data is time data or group data.
#'
#' @param data An object of faultdata
#' @return A string of "time", "group" and "unknown"

check.faultdata <- function(data) {
  if (all(data$fault == 0) && any(data$type == 1)) {
    "time"
  } else if (any(data$fault > 0) && all(data$type == 0)) {
    "group"
  } else {
    "unknown"
  }
}
