"load.qiime.mapping.file" <- function(filepath){
    map <- read.table(filepath,sep='\t',comment='',head=T,row.names=1,check=FALSE)
    return(map)
}

"load.qiime.otu.table" <- function(filepath,include.lineages=FALSE){
    otus <- read.table(filepath,sep='\t',comment='',head=T,row.names=1,check=FALSE, skip=1)
    # drop "Consensus Lineage" column if present
    C = ncol(otus)
    if(colnames(otus)[C]=='Consensus Lineage'){
        lineages <- otus[,C]
        otus <- otus[,-C]
    } else {
        lineages <- NULL
    }
    otus <- as.matrix(t(otus))
    
    if(include.lineages){
        return(list(otus=otus,lineages=lineages))
    } else {
        return(otus=otus)
    }
}

"load.qiime.taxon.table" <- function(filepath){
    taxa <- as.matrix(t(read.table(filepath,sep='\t',head=T,row.names=1,check=FALSE)))
    return(taxa)
}

"load.qiime.distance.matrix" <- function(filepath){
    d <- as.matrix(read.table(filepath,sep='\t',head=T,row.names=1,check=FALSE))
    return(d)
}

# ensure map, data table, etc., contain the same samples in the same order
"remove.nonoverlapping.samples" <- function(map=NULL,otus=NULL,taxa=NULL,distmat=NULL){
    IDs <- NULL
    objects <- list(map=map,otus=otus,taxa=taxa,distmat=distmat)

    # find overlapping samples in all tables
    for(obj in objects){
        if(!is.null(obj)) {
            if(is.null(IDs)){
                IDs <- rownames(obj)
            } else {
                IDs <- intersect(rownames(obj), IDs)
            }
        }
    }
    
    # drop non-overlapping samples 
    for(i in 1:length(objects)){
        if(!is.null(objects[[i]])) {
            objects[[i]] <- objects[[i]][IDs,]
        }
    }
    
    return(objects)
}
