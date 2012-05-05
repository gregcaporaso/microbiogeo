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
library('vegan', warn.conflicts=FALSE, quietly=TRUE)

# TODO - remove me
source('source/vegan/R/capscale.R')
source('source/vegan/R/ordiGetData.R')
source('source/vegan/R/ordiParseFormula.R')

envvars <- as.list(Sys.getenv())
if (is.element('QIIME_DIR', names(envvars))) {
    qiimedir <- envvars[['QIIME_DIR']]
    source(sprintf('%s/qiime/support_files/R/loaddata.r', qiimedir))
} else {
    stop("Please add QIIME_DIR environment variable pointing to the top-level QIIME directory.")
}

# Make option list and parse command line.
option_list <- list(
    make_option(c("-i", "--distmat"), type="character",
        help="Input distance matrix [required]."),
    make_option(c("-m", "--mapfile"), type="character",
        help="Input metadata mapping file [required]."),
    make_option(c("-c", "--category"), type="character",
        help="Metadata column header giving cluster IDs [required]"),
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

# Run DB-RDA and create a plot of the results.
factor = as.factor((qiime.data$map[[opts$category]]))
factors.frame <- data.frame(factor)
capscale.results <- capscale(as.dist(qiime.data$distmat) ~ factor, factors.frame)
capscale.results.filepath <- sprintf('%s/rda_results.txt', opts$outdir)
sink(capscale.results.filepath)
print(capscale.results)
sink(NULL)

plot.filepath <- sprintf('%s/rda_plot.pdf', opts$outdir)
pdf(plot.filepath)
plot(capscale.results, display=c("wa", "bp"))
