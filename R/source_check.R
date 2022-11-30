check_data <- function(x, nm="data") {

  if (!is.data.frame(x)) stop(paste0("ERROR: ", nm, " must be a data frame"))
  if (!nrow(x)) stop(paste0("ERROR: ", nm, " has 0 rows"))
  if (!ncol(x)) stop(paste0("ERROR: ", nm, " has 0 columns"))
  NULL
}

check_variable <- function(x, nm, valid) {

  if (!isString(x)) stop(paste0("ERROR: ", nm, " must be a column name in the input data"))
  if (!(x %in% valid)) stop(paste0("ERROR: ", nm, " must be a column name in the input data"))
  NULL
}

isString <- function(x) {
  is.character(x) && (length(x) == 1)
}

check_number <- function(x, nm, min=NULL, max=NULL) {

  if (!is.numeric(x)) stop(paste0("ERROR: ", nm, " must be a single numeric value"))
  if (length(x) != 1) stop(paste0("ERROR: ", nm, " must be a single numeric value"))
  if (length(min) && (x < min)) stop(paste0("ERROR: ", nm, " must be >= ", min))
  if (length(max) && (x > max)) stop(paste0("ERROR: ", nm, " must be <= ", max))
  NULL
}

check_tl_tu <- function(tl, tu) {

  tl.nm <- "tl"
  tu.nm <- "tu"
  check_number(tl, tl.nm, min=0, max=NULL) 
  check_number(tu, tu.nm, min=0, max=NULL) 
  if (tl >= tu) stop(paste0("ERROR: ", tl.nm, " >= ", tu.nm))
  NULL
}

check_dataAndVars <- function(data, obs.time, time.to.event, event.status, trt.group) {

  check_data(data, nm="data")
  valid <- colnames(data)
  check_variable(obs.time, "obs.time", valid)
  check_variable(time.to.event, "time.to.event", valid)
  check_variable(event.status, "event.status", valid)
  check_variable(trt.group, "trt.group", valid)

  # Remove missing and non-finite data
  vv  <- c(obs.time, time.to.event, event.status, trt.group)
  tmp <- rep(TRUE, nrow(data))
  for (v in vv) {
    vec <- data[, v, drop=TRUE]
    tmp <- tmp & is.finite(vec)
  }
  data <- data[tmp, , drop=FALSE]
  if (!nrow(data)) stop("ERROR: no observations left after removing missing/non-finite values")

  tmp <- data[, obs.time, drop=TRUE] < 0
  if (any(tmp)) stop("ERROR: obs.time contains negative values")
  tmp <- data[, time.to.event, drop=TRUE] < 0
  if (any(tmp)) stop("ERROR: time.to.event contains negative values")
  tmp <- !(data[, event.status, drop=TRUE] %in% 0:1)
  if (any(tmp)) stop("ERROR: event.status must be binary coded as 0 for non-events, 1 for events")
  tmp <- !(data[, trt.group, drop=TRUE] %in% 0:1)
  if (any(tmp)) stop("ERROR: trt.group must be binary coded as 0 for controls or 1 for treated group")

  data
}
