library(tidyverse)
library(fastverse)
set_collapse(
  mask = c("helper", "special", "fast-fun"),
  na.rm = T,
  nthreads = 48
)
library(vroom)

cusDetect <- function(x, pattern) {
  grepl(
    pattern = pattern,
    x = x,
    ignore.case = TRUE,
    perl = TRUE
  )
}

cusDiscard <- function(x, pattern) {
  grep(
    pattern = pattern,
    x = x,
    ignore.case = TRUE,
    perl = TRUE,
    value = TRUE,
    invert = TRUE
  )
}

setwd("~/volumes/ultra") # setup your working directory
