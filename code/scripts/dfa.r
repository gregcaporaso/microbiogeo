# Runs a stepwise linear DFA on a QIIME OTU table and mapping file.
# usage:
# R --slave --args --source_dir $QIIME_HOME/qiime/support_files/R/ -i otu_table.txt -m Fasting_Map.txt -c Treatment -o dfa < dfa.r
#
# print help string:
# R --slave --args -h --source_dir $QIIME_HOME/qiime/support_files/R/ < dfa.r
#
# Requires command-line param --source_dir pointing to QIIME R source dir

# load libraries and source files
args <- commandArgs(trailingOnly=TRUE)
if(!is.element('--source_dir', args)){
    stop("\n\nPlease use '--source_dir' to specify the R source code directory.\n\n")
}
sourcedir <- args[which(args == '--source_dir') + 1]
source(sprintf('%s/loaddata.r',sourcedir))
source(sprintf('%s/util.r',sourcedir))
load.library('optparse')
load.library('klaR')

# make option list and parse command line
option_list <- list(
    make_option(c("--source_dir"), type="character",
        help="Path to R source directory [required]."),
    make_option(c("-i", "--otutable"), type="character",
        help="Input OTU table [required]."),
    make_option(c("-m", "--mapfile"), type="character",
        help="Input metadata mapping file [required]."),
    make_option(c("-c", "--category"), type="character",
        help="Metadata column header giving cluster IDs [required]."),
    make_option(c("-o", "--outdir"), type="character", default='.',
        help="Output directory [default %default].")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# make sure we have our required files
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
results.filepath <- sprintf('%s/dfa_results.txt', opts$outdir)
sink(results.filepath)
print(results)
sink(NULL)

plot.filepath <- sprintf('%s/dfa_plot.pdf', opts$outdir)
pdf(plot.filepath)
plot(results)
