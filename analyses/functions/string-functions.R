# based on: https://statisticsglobe.com/format-number-as-percentage-in-r
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(x * 100, format = format, digits = digits, ...), "%")
}