# Runs the Mantel correlogram method on a QIIME distance matrix and geographic
# distance matrix.
#
# Usage:
# R --slave --args -d dist_mat.txt -s geo_dist_mat.txt < r/mantel_correlogram.r
# 
# Print help string:
# R --slave --args -h < r/mantel_correlogram.r
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
option_list <- list(
    make_option(c("-d", "--distmat"), type="character",
        help="Input distance matrix [required]."),
    make_option(c("-s", "--geodistmat"), type="character",
        help="Input geographic distance matrix [required]."),
    make_option(c("-o", "--outdir"), type="character", default='.',
        help="Output directory [default %default]")
)
opts <- parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=TRUE))

# Make sure we have our required files.
if (is.null(opts$distmat)) stop('Please supply a distance matrix.')
if (is.null(opts$geodistmat)) stop('Please supply a geographic distance matrix.')

# Create output directory if needed.
if (opts$outdir != ".") dir.create(opts$outdir, showWarnings=FALSE, recursive=TRUE)

# Load data.
distmat <- load.qiime.distance.matrix(opts$distmat)
geodistmat <- load.qiime.distance.matrix(opts$geodistmat)

# Run Mantel correlogram, write results summary to a file, and create a plot of
# the correlogram (save as a PDF).
results <- mantel.correlog(as.dist(distmat), D.geo=as.dist(geodistmat))
results.filepath <- sprintf('%s/mantel_correlogram_results.txt', opts$outdir)
sink(results.filepath)
print(results)
sink(NULL)

plot.filepath <- sprintf('%s/mantel_correlogram_plot.pdf', opts$outdir)
pdf(plot.filepath)
plot(results)
