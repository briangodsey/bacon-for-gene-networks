## load data
datdir <- "./data/DREAM4/"

numgenes <- 10
##numgenes <- 100

for( datasetnum in 1:5 ){
  ##datasetnum <- 1

  goldstdname <- paste(datdir,
                       "DREAM4_Challenge2_GoldStandards/Size ",
                       numgenes,
                       "/DREAM4_GoldStandard_InSilico_Size",
                       numgenes,"_",
                       datasetnum,".tsv",sep="")

  goldstd <- read.table(goldstdname,
                        stringsAsFactors=FALSE)


  ## create a gold standard matrix
  gstdmat <- matrix(0,numgenes,numgenes)
  rownames(gstdmat) <- paste("G",1:numgenes,sep="")
  colnames(gstdmat) <- paste("G",1:numgenes,sep="")

  for( i in 1:nrow(goldstd) ){
    gstdmat[ goldstd[i,1],goldstd[i,2] ] <- goldstd[i,3]
  }

  write.table(gstdmat,
              file=paste("GoldStandardMat_",numgenes,"_",
                datasetnum,".tsv",sep=""),
              row.names=FALSE, col.names=FALSE,quote=FALSE)
}
