# Runs vegan function RDA on QIIME distance file.
#
# Usage:
# R --slave --args -d distance_matrix.txt -m mapping_file.txt -c Treatment < r/rda.r
# 
# Print help string:
# R --slave --args -h < r/rda.r
#
# Requires environment variable QIIME_DIR pointing to top-level QIIME directory.

# Load libraries and source files.
library('optparse', warn.conflicts=FALSE, quietly=TRUE)

# TODO - remove me
source('source/R-2.14.2/src/library/stats/R/cmdscale.R')

envvars <- as.list(Sys.getenv())
if (is.element('QIIME_DIR', names(envvars))) {
    qiimedir <- envvars[['QIIME_DIR']]
    source(sprintf('%s/qiime/support_files/R/loaddata.r', qiimedir))
} else {
    stop("Please add QIIME_DIR environment variable pointing to the top-level QIIME directory.")
}

# Make option list and parse command line.
option_list <- list(
    make_option(c("-d", "--distmat"), type="character",
        help="Input distance matrix [required]."),
    make_option(c("-o", "--outdir"), type="character", default='.',
        help="Output directory [default %default]")
)
opts <- parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=TRUE))

# Make sure we have our required files.
if (is.null(opts$distmat)) stop('Please supply a distance matrix.')

# Create output directory if needed.
if (opts$outdir != ".") dir.create(opts$outdir, showWarnings=FALSE, recursive=TRUE)

# Load data.
distmat <- load.qiime.distance.matrix(opts$distmat)
qiime.data <- remove.nonoverlapping.samples(distmat=distmat)

# Run cmdscale.
results <- cmdscale(as.dist(qiime.data$distmat), k=8, eig=TRUE)
#print(results)
x <- results$points[, 1]
y <- -results$points[, 2] # reflect so North is at the top
## note asp = 1, to ensure Euclidean distances are represented correctly
plot(x, y, type = "n", xlab = "", ylab = "", asp = 1, main = "cmdscale(distmat)")
text(x, y, rownames(results$points), cex = 0.6)
