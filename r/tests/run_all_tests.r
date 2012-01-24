# R --slave --args --qiime_dir $QIIME_HOME < run_all_tests.r
# 
# print help string:
# R --slave --args -h < run_all_tests.r
#
args <- commandArgs(trailingOnly=TRUE)
if(!is.element('--source_dir', args)){
    stop("\n\nPlease use '--source_dir' to specify the R source code directory.\n\n")
}
sourcedir <- args[which(args == '--source_dir') + 1]

# load libraries and source files
source(sprintf('%s/loaddata.r',sourcedir))
source(sprintf('%s/util.r',sourcedir))
load.library('optparse')
load.library('RUnit')


# make option list and parse command line
option_list <- list(
    make_option(c("--source_dir"), type="character",
        help="Path to R source directory [required]."),
    make_option(c("--test_dir"), type="character", default='.',
        help="Path to R source directory [default %default].")
)
opts <- parse_args(OptionParser(option_list = option_list), args = args)

ts <- defineTestSuite('QIIME R Tests', opts$test_dir)
tr <- runTestSuite(ts, verbose=0)
printTextProtocol(tr)