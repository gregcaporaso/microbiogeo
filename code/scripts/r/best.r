# Runs vegan function BIOENV on QIIME distance file.
#
# Usage:
# R --slave --args -c unweighted_unifrac_dm.txt -e < r/best.r > best_result.txt
# 
# Print help string:
# R --slave --args -h < r/best.r
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
    make_option(c("-c", "--comm_dm"), type="character",
        help="Input distance matrix [required]."),
    make_option(c("-e", "--env_vars"), type="character",
        help="Environmental variable table(mapping file)")#,
    #make_option(c("-o", "--outdir"), type="character", default='.',
    #    help="Output directory [default %default]")
)
opts <- parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=TRUE))

# Make sure we have our required files.
if (is.null(opts$comm_dm)) stop('Please supply a mapping file.')
if (is.null(opts$env_vars)) stop('Please supply a distance matrix.')

# Create output directory if needed.
#if (opts$outdir != ".") dir.create(opts$outdir, showWarnings=FALSE, recursive=TRUE)

# Load data.
cdm <- load.qiime.distance.matrix(opts$comm_dm)
edm <- load.qiime.mapping.file(opts$env_vars)
#edm <- load.qiime.distance.matrix(opts$env_vars)
#qiime.data <- remove.nonoverlapping.samples(comm_dm=comm_dm, env_vars=env_vars)

#print(qiime.data$distmat)
#data(varespec)
#print(varespec)

#source('bioenv2.default.R')
#source('bioenv2.formula.R')
#source('bioenv2.R')
#print(scale(edm)
#print(ncol(cdm))
#cdm2 = vegdist(cdm, method="bray")
#print(ncol(cdm2))

out <- bioenv(comm = cdm, env = edm, method="spearman", index="euclidean")#qiime.data$map[[opts$category]])
print(out)
summary(out)

#source('bvstep.r')
#out <- bv.step(cdm,edm ,  
#        fix.dist.method="pearson", var.dist.method="pearson", 
#        scale.fix=FALSE, scale.var=FALSE,  
#        max.rho=0.95, min.delta.rho=0.001, 
#        random.selection=TRUE, 
#        prop.selected.var=0.3, 
#        num.restarts=50, 
#        output.best=10, 
#        var.always.include=NULL
#        ) 

