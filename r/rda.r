# Runs vegan function RDA on QIIME distance file.
#
# Usage:
# R --slave --args -d unifrac.txt < rda.r
# 
# Print help string:
# R --slave --args -h < rda.r
#
# Requires environment variable QIIME_DIR pointing to top-level QIIME directory.

# Load libraries and source files.
library('optparse',warn.conflicts=FALSE,quietly=TRUE)
library('vegan',warn.conflicts=FALSE,quietly=TRUE)
envvars <- as.list(Sys.getenv())
if(is.element('QIIME_DIR', names(envvars))){
    qiimedir <- envvars[['QIIME_DIR']]
    source(sprintf('%s/qiime/support_files/R/loaddata.r',qiimedir))
} else {
    stop("Please add QIIME_DIR environment variable pointing to the top-level QIIME directory.")
}

# Make option list and parse command line.
option_list <- list(
    make_option(c("-d", "--distmat"), type="character",
        help="Input distance matrix [required].")
)
opts <- parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=TRUE))

# File requirements.
if(is.null(opts$distmat)) stop('Please supply a distance matrix.')

# Load distance matrix.
distmat <- load.qiime.distance.matrix(opts$distmat)

# Run rda and create a plot of the results.
pdf("./rda_plot.pdf")
results <- rda(distmat)
plot(results)
