## A script to bring all of the different plots together into their separate
## figures and into one core document

#### Set Up                                                                 ####

## Libraries to load
require(dplyr)
require(ggplot2)
require(qpdf)

## Set up where to save
## Date of the run
init_date <- as.character(snakemake@params[["date"]])

## Fitness function
init_fit_func <- as.character(snakemake@params[["fit_func"]])

## Mutation rate
init_mu <- as.numeric(snakemake@params[["mu"]])

## Create the filepath where things will be saved
filepath <- paste0("runs/ddt_", init_fit_func, "_mu_", init_mu, "_", init_date, "/")

## Recursively list .rds files under /res/
pdf_files <- list.files(
  path = paste0(filepath, "res/"),
  pattern = "^fig[0-9]\\.pdf$",
  full.names = TRUE,
  recursive = TRUE
)

pdf_combine(input = pdf_files, output = paste0(filepath, "res/all_figs.pdf"))
