# Runs vegan function BETADISPER on QIIME distance file.
# usage:
# R --slave --args -d unifrac.txt -m Fasting_Map.txt -c Treatment -o betadisper < betadisper.r
# 
# print help string:
# R --slave --args -h < betadisper.r
#
# Requires environment variable QIIME_DIR pointing to  top-level QIIME directory.

# load libraries and source files
library('optparse',warn.conflicts=FALSE,quietly=TRUE)
library('vegan',warn.conflicts=FALSE,quietly=TRUE)
envvars <- as.list(Sys.getenv())
if(is.element('QIIME_DIR', names(envvars))){
    qiimedir <- envvars[['QIIME_DIR']]
    source(sprintf('%s/qiime/support_files/R/loaddata.r',qiimedir))
} else {
    stop("Please add QIIME_DIR environment variable pointing to the top-level QIIME directory.")
}

# make option list and parse command line
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

# File requirements
if(is.null(opts$mapfile)) stop('Please supply a mapping file.')
if(is.null(opts$category)) stop('Please supply a mapping file header.')
if(is.null(opts$distmat)) stop('Please supply a distance matrix.')

# create output directory if needed
if(opts$outdir != ".") dir.create(opts$outdir,showWarnings=FALSE, recursive=TRUE)

# load qiime data
map <- load.qiime.mapping.file(opts$mapfile)
distmat <- load.qiime.distance.matrix(opts$distmat)
qiime.data <- remove.nonoverlapping.samples(map=map, distmat=distmat)

# begin permdisp/betadisper
dis <- vegdist(distmat)

# loads categories
attach(map)

# Calculate multivariate dispersions
mod <- betadisper(dis, qiime.data$map[[opts$category]])

# Perform test
results1 <- anova(mod)

# Permutation test for F
results2 <- permutest(mod, pairwise = TRUE)

# write output file
filepath <- sprintf('%s/betadisper_results.txt',opts$outdir)
sink(filepath)
print(results1)
print(results2)
sink(NULL)
