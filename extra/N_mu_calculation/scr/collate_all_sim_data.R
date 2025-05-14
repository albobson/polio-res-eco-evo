## This script takes in the data generated in the sims and
## collates them into one csv file.

# List all .csv files in the directory starting with 'bo'
files <- list.files(path = "sim_data/", pattern = "*.csv", full.names = TRUE)

## Read in data
csvs <- do.call(rbind,
                   lapply(files,
                          function(x) read.csv(x, stringsAsFactors = FALSE)))

## Write data
write.csv(csvs, paste0("sim_data/all_stoch-v2_var.csv"), fileEncoding = "UTF-8")