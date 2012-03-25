# Runs a stepwise linear DFA on a QIIME OTU table and mapping file.
#
# Usage:
# R --slave --args -i otu_table.txt -m mapping_file.txt -c Treatment < r/dfa.r
# 
# Print help string:
# R --slave --args -h < r/dfa.r
#
# Requires environment variable QIIME_DIR pointing to top-level QIIME directory.

# Load libraries and source files.
library('optparse', warn.conflicts=FALSE, quietly=TRUE)
library('klaR', warn.conflicts=FALSE, quietly=TRUE)

envvars <- as.list(Sys.getenv())
if (is.element('QIIME_DIR', names(envvars))) {
    qiimedir <- envvars[['QIIME_DIR']]
    source(sprintf('%s/qiime/support_files/R/loaddata.r', qiimedir))
} else {
    stop("Please add QIIME_DIR environment variable pointing to the top-level QIIME directory.")
}

# Make option list and parse command line.
option_list <- list(
    make_option(c("-i", "--otutable"), type="character",
        help="Input OTU table [required]."),
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
if (is.null(opts$otutable)) stop('Please supply an OTU table.')

# Create output directory if needed.
if (opts$outdir != ".") dir.create(opts$outdir, showWarnings=FALSE, recursive=TRUE)

# Load data.
map <- load.qiime.mapping.file(opts$mapfile)
otutable <- load.qiime.otu.table(opts$otutable)
qiime.data <- remove.nonoverlapping.samples(map=map, otus=otutable)

# More error checking.
if (nrow(qiime.data$map) == 0)
    stop('\n\nMapping file and OTU table have no samples in common.\n\n')
if (!is.element(opts$category, colnames(qiime.data$map)))
    stop(sprintf('\n\nHeader %s not found in mapping file.\n\n', opts$category))

# Run stepwise linear DFA function and save results.
results <- stepclass(qiime.data$otus, qiime.data$map[[opts$category]], 'lda')
#results <- lda(qiime.data$otus, qiime.data$map[[opts$category]])
results.filepath <- sprintf('%s/dfa_results.txt', opts$outdir)
sink(results.filepath)
print(results)
sink(NULL)

plot.filepath <- sprintf('%s/dfa_plot.pdf', opts$outdir)
pdf(plot.filepath)
plot(results)
