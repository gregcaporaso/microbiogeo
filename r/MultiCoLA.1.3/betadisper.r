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
    make_option(c("-s", "--samples"), type="integer",
        help="Input distance matrix [required]."),
    make_option(c("-l", "--levels"), type="integer",
        help="Input metadata mapping file [required]."),
    make_option(c("-w", "--workingdir"), type="character",
        help="Metadata column header giving cluster IDs [required]"),
    make_option(c("-m", "--map"), type="character",
        help="Input metadata mapping file [required]."),
    make_option(c("-o", "--otu"), type="character",
        help="Metadata column header giving cluster IDs [required]")
)
opts <- parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=TRUE))

# File requirements
setwd(opts$workingdir)

M<-read.table((opts$otu),header=TRUE,row.names=1)
source("../taxa.pooler.1.3.r")
source("../COtables.1.3.r")
source("../cutoff.impact.1.3.r")
source("../cutoff.impact.fig.1.3.r")
ENV<-read.table(opts$map,header=TRUE,row.names=1)
source("../VP.COL.1.3.r")
source("../corrcoeff.ENV.1.3.r")
source("../signif.1.4.r")

all_taxa_pooled<-taxa.pooler(M,(opts$samples),(opts$levels))
truncated.DS.phylum<-COtables(all_taxa_pooled[[1]], Type="ADS",typem="dominant")
truncated.DS.class<-COtables(all_taxa_pooled[[2]], Type="ADS",typem="dominant")
truncated.DS.order<-COtables(all_taxa_pooled[[3]], Type="ADS",typem="dominant")
truncated.DS.family<-COtables(all_taxa_pooled[[4]], Type="ADS",typem="dominant")
truncated.DS.genus<-COtables(all_taxa_pooled[[5]], Type="ADS",typem="dominant")
truncated.DS.OTUcompleteDS<-COtables(all_taxa_pooled[[6]], Type="ADS",typem="dominant")
truncated.DS.OTUwholeDS<-COtables(all_taxa_pooled[[7]], Type="ADS",typem="dominant")
corr.all<-cutoff.impact(all_taxa_pooled,"y",Type="ADS",corcoef="spearman",typem="dominant")
output.all<-cutoff.impact.fig(corr.all)
VP.1.taxa<-VP.COL(all_taxa_pooled,ENV,Type="ADS")
corrcoeff.table.ADS<-matrix(NA,21,5)
row.names(corrcoeff.table.ADS)<-c(paste("CO_",c(0.01,seq(0.05,0.95,by=0.05),0.99),sep=""))
colnames(corrcoeff.table.ADS)<-c("Sum",paste("RDA1.",colnames(ENV),sep=""))
OTU.ADS<-VP.1.taxa[[c(7,3)]]
SPE<-OTU.ADS[[1]];corrcoeff.table.ADS[1,]<-corrcoeff(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[2]];corrcoeff.table.ADS[2,]<-corrcoeff(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[3]];corrcoeff.table.ADS[3,]<-corrcoeff(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[4]];corrcoeff.table.ADS[4,]<-corrcoeff(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[5]];corrcoeff.table.ADS[5,]<-corrcoeff(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[6]];corrcoeff.table.ADS[6,]<-corrcoeff(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[6]];corrcoeff.table.ADS[7,]<-corrcoeff(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[7]];corrcoeff.table.ADS[8,]<-corrcoeff(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[8]];corrcoeff.table.ADS[9,]<-corrcoeff(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[9]];corrcoeff.table.ADS[10,]<-corrcoeff(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[10]];corrcoeff.table.ADS[11,]<-corrcoeff(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[11]];corrcoeff.table.ADS[12,]<-corrcoeff(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[13]];corrcoeff.table.ADS[13,]<-corrcoeff(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[14]];corrcoeff.table.ADS[14,]<-corrcoeff(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[15]];corrcoeff.table.ADS[15,]<-corrcoeff(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[16]];corrcoeff.table.ADS[16,]<-corrcoeff(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[17]];corrcoeff.table.ADS[17,]<-corrcoeff(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[18]];corrcoeff.table.ADS[18,]<-corrcoeff(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[19]];corrcoeff.table.ADS[19,]<-corrcoeff(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[20]];corrcoeff.table.ADS[20,]<-corrcoeff(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[21]];corrcoeff.table.ADS[21,]<-corrcoeff(SPE,ENV);rm(SPE)
SPE<-all_taxa_pooled[[7]]
corrcoeff.table.ADS.orig<-corrcoeff(SPE,ENV)
row.names(corrcoeff.table.ADS.orig)<-c("CO_1")
corrcoeff.table.ADS<-rbind(corrcoeff.table.ADS,corrcoeff.table.ADS.orig)
write.table(corrcoeff.table.ADS,"corrcoeff.table.ADS.txt",quote=FALSE)
signif.table.ADS<-matrix(NA,21,5)
row.names(signif.table.ADS)<-c(paste("CO_",c(0.01,seq(0.05,0.95,by=0.05),0.99),sep=""))
colnames(signif.table.ADS)<- c("whole.sig","ENV1.sig","ENV2.sig","ENV3.sig","ENV4.sig")
OTU.ADS<-VP.1.taxa[[c(7,3)]]
SPE<-OTU.ADS[[1]];signif.table.ADS[1,]<-signif(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[2]];signif.table.ADS[2,]<-signif(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[3]];signif.table.ADS[3,]<-signif(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[4]];signif.table.ADS[4,]<-signif(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[5]];signif.table.ADS[5,]<-signif(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[6]];signif.table.ADS[6,]<-signif(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[7]];signif.table.ADS[7,]<-signif(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[8]];signif.table.ADS[8,]<-signif(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[9]];signif.table.ADS[9,]<-signif(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[10]];signif.table.ADS[10,]<-signif(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[11]];signif.table.ADS[11,]<-signif(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[12]];signif.table.ADS[12,]<-signif(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[13]];signif.table.ADS[13,]<-signif(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[14]];signif.table.ADS[14,]<-signif(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[15]];signif.table.ADS[15,]<-signif(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[16]];signif.table.ADS[16,]<-signif(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[17]];signif.table.ADS[17,]<-signif(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[18]];signif.table.ADS[18,]<-signif(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[19]];signif.table.ADS[19,]<-signif(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[20]];signif.table.ADS[20,]<-signif(SPE,ENV);rm(SPE)
SPE<-OTU.ADS[[21]];signif.table.ADS[21,]<-signif(SPE,ENV);rm(SPE)
SPE<-all_taxa_pooled[[7]]
signif.table.ADS.orig<-signif(SPE,ENV)
row.names(signif.table.ADS.orig)<-c("CO_1")
signif.table.ADS<-rbind(signif.table.ADS, signif.table.ADS.orig)
write.table(signif.table.ADS,"signif.table.ADS.txt",quote=FALSE)


