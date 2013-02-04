# Runs vegan function mantel.partial on QIIME distance matrices
#
# Usage:
# R --slave --args -d1 mat1.txt -d2 mat2.txt -d3 control_matrix.txt < r/pmantel.r
# 
# Print help string:
# R --slave --args -h < r/pmantel.r
#
# Requires environment variable QIIME_DIR pointing to top-level QIIME directory.

# Load libraries and source files.
library('optparse', warn.conflicts=FALSE, quietly=TRUE)
library('vegan', warn.conflicts=FALSE, quietly=TRUE)
envvars <- as.list(Sys.getenv())
if (is.element('QIIME_DIR', names(envvars))) {
    qiimedir <- envvars[['QIIME_DIR']]
    source(sprintf('%s/qiime/support_files/R/loaddata.r', qiimedir))
} else {
    stop("Please add QIIME_DIR environment variable pointing to the top-level QIIME directory.")
}

# Make option list and parse command line.
option_list1 <- list(
    make_option(c("-a", "--distmat1"), type="character",
        help="Input distance matrix [required]."),
    make_option(c("-b", "--distmat2"), type="character",
        help="Input distance matrix [required]."),
    make_option(c("-c", "--distmat3"), type="character",
        help="Input distance matrix [required]."),
    make_option(c("-m", "--mapfile"), type="character",
        help="Input metadata mapping file [required]."),
    make_option(c("-o", "--outdir"), type="character", default='.',
        help="Output directory [default %default]")
)
opts <- parse_args(OptionParser(option_list = option_list1), args = commandArgs(trailingOnly=TRUE))

# Make sure we have our required files.
if (is.null(opts$distmat1)) stop('Please supply all three distance matrices.')
if (is.null(opts$distmat2)) stop('Please supply all three distance matrices.')
if (is.null(opts$distmat3)) stop('Please supply all three distance matrices.')

# Create output directory if needed.
if (opts$outdir != ".") dir.create(opts$outdir, showWarnings=FALSE, recursive=TRUE)

# Load data.
distmat1 <- load.qiime.distance.matrix(opts$distmat1)
distmat2 <- load.qiime.distance.matrix(opts$distmat2)
distmat3 <- load.qiime.distance.matrix(opts$distmat3)

print(mantel.partial(distmat1, distmat2, distmat3, method="pearson", permutations=999))
