# Edit these two paths, or set ESCC_DATA / ESCC_RESULTS in your environment.
data_dir <- Sys.getenv("ESCC_DATA", unset = "data")
results_dir <- Sys.getenv("ESCC_RESULTS", unset = "results_core")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
input <- function(name) file.path(data_dir, name)
output <- function(name) file.path(results_dir, name)
read_tab <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
write_tab <- function(x, name) write.table(x, output(name), sep = "\t", row.names = FALSE, quote = FALSE)
