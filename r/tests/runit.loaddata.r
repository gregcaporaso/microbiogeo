# Creates temporary mapping, otu table, distance matrix, and taxon files
".setUp" <- function(){

    # create the 8 types of mapping file
    sampleids <- sprintf('sample%d',1:4)
    headers <- sprintf('header%d',1:4)
    map.columns <- list(factor(c('a','b','c','a')), 1:4, 1:4/4,c(3,NA,4,5))
    fullmap <- data.frame(map.columns)
    rownames(fullmap) <- sampleids
    colnames(fullmap) <- headers
    
    map.fnames <- replicate(8,tempfile())
    # type1: no comments before header, no comments after header, header w/#
    sink(map.fnames[[1]])
    cat('#SampleID\t')
    write.table(fullmap,quote=F,sep='\t')
    sink(NULL)
    
    # type2: comments before header, no comments after header, header w/#
    sink(map.fnames[[2]])
    cat('# comment 1\n# comment 2\n');
    cat('#SampleID\t')
    write.table(fullmap,quote=F,sep='\t')
    sink(NULL)
    
    # type3: no comments before header, comments after header, header w/#
    sink(map.fnames[[3]])
    cat('#SampleID\t', paste(colnames(fullmap),collapse='\t'),'\n', sep='')
    cat('# comment 1\n# comment 2\n');
    write.table(fullmap,quote=F,sep='\t',col.names=FALSE)
    sink(NULL)

    # type4: comments before header, comments after header, header w/#
    sink(map.fnames[[4]])
    cat('# comment 1\n# comment 2\n');
    cat('#SampleID\t', paste(colnames(fullmap),collapse='\t'),'\n', sep='')
    cat('# comment 3\n# comment 4\n');
    write.table(fullmap,quote=F,sep='\t',col.names=FALSE)
    sink(NULL)

    # type5: no comments before header, no comments after header, header w/o #
    sink(map.fnames[[5]])
    cat('SampleID\t')
    write.table(fullmap,quote=F,sep='\t')
    sink(NULL)

    # type6: comments before header, no comments after header, header w/o #
    sink(map.fnames[[6]])
    cat('# comment 1\n# comment 2\n');
    cat('SampleID\t')
    write.table(fullmap,quote=F,sep='\t')
    sink(NULL)

    # type7: no comments before header, comments after header, header w/o #
    sink(map.fnames[[7]])
    cat('#SampleID\t', paste(colnames(fullmap),collapse='\t'),'\n', sep='')
    cat('# comment 1\n# comment 2\n');
    write.table(fullmap,quote=F,sep='\t',col.names=FALSE)
    sink(NULL)

    # type8: comments before header, comments after header, header w/o #
    sink(map.fnames[[8]])
    cat('# comment 1\n# comment 2\n');
    cat('#SampleID\t', paste(colnames(fullmap),collapse='\t'),'\n', sep='')
    cat('# comment 3\n# comment 4\n');
    write.table(fullmap,quote=F,sep='\t',col.names=FALSE)
    sink(NULL)

    # create 2 types of OTU table, with/without lineage
    otuids <- sprintf('otu%d',1:3)
    otus <- as.data.frame(matrix(1:9,3,3), row.names=otuids)
    colnames(otus) <- sampleids[-4]
    otu.fnames <- replicate(2,tempfile())

    # without lineage
    sink(otu.fnames[[1]])
    cat('# comment\n#OTU ID\t')
    write.table(otus, quote=F, sep='\t')
    sink(NULL)
    
    # with lineage
    lineages <- sprintf('lineage%d',1:3)
    otus[['Consensus Lineage']] <- lineages
    sink(otu.fnames[[2]])
    cat('# comment\n#OTU ID\t')
    write.table(otus, quote=F, sep='\t')
    sink(NULL)
    
    # create taxon table
    taxonomies <- sprintf('taxonomy%d',1:3)
    taxa <- as.data.frame(matrix(1:9,3,3), row.names=taxonomies)
    colnames(taxa) <- sampleids[-4]
    taxon.fname <- tempfile()
    sink(taxon.fname)
    cat('Taxon\t')
    write.table(taxa, quote=F, sep='\t')
    sink(NULL)

    # create distance matrix
    dist.mat <- matrix(1:16,4,4)
    rownames(dist.mat) <- sampleids
    colnames(dist.mat) <- sampleids
    dist.fname <- tempfile()
    sink(dist.fname)
    cat('\t')
    write.table(dist.mat, quote=F, sep='\t')
    sink(NULL)


    # store these as global variables
    .GlobalEnv$true.map <- fullmap
    .GlobalEnv$true.otus <- as.matrix(t(otus[,-4])) # no lineage
    .GlobalEnv$true.lineages <- lineages
    .GlobalEnv$true.taxa <- as.matrix(t(taxa))
    .GlobalEnv$true.dist.mat <- dist.mat
    .GlobalEnv$map.fnames <- map.fnames
    .GlobalEnv$otu.fnames <- otu.fnames
    .GlobalEnv$taxon.fname <- taxon.fname
    .GlobalEnv$dist.fname <- dist.fname
}

# Removes all temporary files
".tearDown" <- function(){
    filelist <- c(map.fnames, otu.fnames, taxon.fname, dist.fname)
    sapply(filelist, unlink)
}

# Ensure all mapping file formats are read correctly
"test.load.qiime.mapping.file" <- function(){
    for(i in seq_along(map.fnames)){
        map <- load.qiime.mapping.file(map.fnames[[i]])
        checkEquals(true.map, map)
    }
}

# Ensure all otu file formats are read correctly
"test.load.qiime.otu.table" <- function(){
    for(i in seq_along(otu.fnames)){
        otus <- load.qiime.otu.table(otu.fnames[[i]])
        checkEquals(true.otus, otus)
    }
    # check correct reading of lineages
    otu.res <- load.qiime.otu.table(otu.fnames[[i]], include.lineages=TRUE)
    checkEquals(true.otus, otu.res$otus)
    checkEquals(true.lineages, otu.res$lineages)
}

# Ensure metadata checking works correctly
"test.otu.table.has.metadata" <- function(){
    headers <- sprintf('lineage%d', 1:3)
    result <- otu.table.has.metadata(headers)
    checkTrue(!result)
    result <- otu.table.has.metadata(c(headers,'Consensus Lineage'))
    checkTrue(result)
    result <- otu.table.has.metadata(c(headers,'OTU Metadata'))
    checkTrue(result)
}

# Ensure header is correctly located in all file formats
"test.get.header.index" <- function(){
    exp <- c(1,3,1,3,1,3,1,3)
    for(i in seq_along(map.fnames)){
        header.ix <- get.header.index(map.fnames[[i]])
        checkEquals(exp[i], header.ix)
    }
}

# Ensure taxon table is read correctly
"test.load.qiime.taxon.table" <- function(){
    result <- load.qiime.taxon.table(taxon.fname)
    checkEquals(true.taxa, result)
}

# Ensure distance matrix is read correctly
"load.qiime.distance.matrix" <- function(){
    result <- load.qiime.distance.matrix(dist.fname)
    checkEquals(true.dist.mat, result)
}

# Ensure nonoverlapping samples are correctly removed
"test.remove.nonoverlapping.samples" <- function(){
    result <- remove.nonoverlapping.samples(map=true.map, otus=true.otus, distmat=true.dist.mat, taxa=true.taxa)
    all.sampleids <- NULL # will hold vector of sample ids in each column
    for(i in seq_along(result)){
        all.sampleids <- cbind(rownames(result[[i]]))
        checkEquals(3, nrow(result[[i]]))
    }
    
    num.unique.ids <- apply(all.sampleids, 1, function(xx) length(unique(xx)))
    checkEquals(c(1,1,1), num.unique.ids)
}
