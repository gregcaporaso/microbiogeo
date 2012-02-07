# Runs a repeated-measures PERMANOVA over a QIIME distance matrix using the
# given category representing a time series to control the permutation model.
#
# Much of this code was taken from http://thebiobucket.blogspot.com/2011/04/repeat-measure-adonis-lately-i-had-to.html#more
# and modified to work with R version 2.13.1 and vegan version 2.0-2.
#
# Usage:
# R --slave --args -d distance_matrix.txt -m mapping_file.txt -c Time < r/repeated_measures_permanova.r
# 
# Print help string:
# R --slave --args -h < r/repeated_measures_permanova.r
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
    make_option(c("-m", "--mapfile"), type="character",
        help="Input metadata mapping file [required]."),
    make_option(c("-c", "--category"), type="character",
        help="Metadata column header representing time series. Must be numeric [required]"),
    make_option(c("-o", "--outdir"), type="character", default='.',
        help="Output directory [default %default]")
)
opts <- parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=TRUE))

# Make sure we have our required files.
if (is.null(opts$mapfile)) stop('Please supply a mapping file.')
if (is.null(opts$category)) stop('Please supply a mapping file category.')
if (is.null(opts$distmat)) stop('Please supply a distance matrix.')

# Create output directory if needed.
if (opts$outdir != ".") dir.create(opts$outdir, showWarnings=FALSE, recursive=TRUE)

# Load data.
map <- load.qiime.mapping.file(opts$mapfile)
distmat <- load.qiime.distance.matrix(opts$distmat)
qiime.data <- remove.nonoverlapping.samples(map=map, distmat=distmat)

# More error checking.
if (nrow(qiime.data$map) == 0)
    stop('\n\nMapping file and distance matrix have no samples in common.\n\n')
if (!is.element(opts$category, colnames(qiime.data$map)))
    stop(sprintf('\n\nHeader %s not found in mapping file.\n\n', opts$category))

time = as.factor((qiime.data$map[[opts$category]]))
time.frame <- data.frame(time)

# Compute the true R2 value.
fit <- adonis(as.dist(qiime.data$distmat) ~ time, time.frame, permutations=1)
print("True R2 value:")
print(fit)

# The number of permutations (this should be user-configurable eventually).
num.perms <- 1999

# Set up frame which will be populated by random R2 values. The first entry will
# be the true R2.
pop <- rep(NA, num.perms + 1)
pop[1] <- fit$aov.tab[1, 5]

# Set up a "permControl" object to make the permuations respect the time series
# variable. Turn off mirroring as time should only flow in one direction.
ctrl <- permControl(within = Within(type = "series", mirror = FALSE))

# Number of samples:
num.samples <- nrow(qiime.data$distmat)

# In adonis(...) you need to put permutations = 1, otherwise adonis will not
# run.
for(i in 2:(num.perms + 1)) {
     idx <- shuffle(num.samples, control = ctrl)
     fit.rand <- adonis(as.dist(qiime.data$distmat) ~ time[idx], time.frame, permutations = 1)
     pop[i] <- fit.rand$aov.tab[1,5]
}

# Get the p-value.
pval <- sum(pop >= pop[1]) / (num.perms + 1)
print("p-value:")
print(pval)
