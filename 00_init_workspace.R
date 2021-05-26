#*******************************************************************************
#
# Project: Medical University Innsbruck - Intrauterine insemination (IUI)
# Date:    2020-04-15
# Author:  Patrick Rockenschaub
# Purpose: Load common packages and set paths
#
#******************************************************************************* 

# Load packages -----------------------------------------------------------

library(tibble)
library(magrittr)
library(tidyverse)
library(lubridate)
library(glue)
library(knitr)


# Set relative paths ------------------------------------------------------

.dir_data <- "data"
.dir_raw <- glue("{.dir_data}/raw")
.dir_der <- glue("{.dir_data}/derived")
.dir_res <- "results"
.dir_src <- "src"