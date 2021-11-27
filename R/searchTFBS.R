tfbs.enrichmentTest.multiple<-function( tfbs, file.twoBit, enh.up.bed, enh.unc.bed, mTH, min.size, run.repeats, ncores ){

  message(paste("repeating how many times:", run.repeats, sep=" "))
  message(paste("TREs changed/background: ", NROW(enh.up.bed), NROW(enh.unc.bed)))

  motifs_list<-list()

  for (i in 1: run.repeats){
    message(paste("LOOP:", i, sep=" "))
    set.seed(i)
    t.comp <- tfbs.enrichmentTest(
          tfbs,
          file.twoBit,
          enh.up.bed,
          negative.bed= enh.unc.bed,
          gc.correction=TRUE,
          gc.min.sample= min.size,
          threshold = mTH,
          pv.adj = "fdr",
          ncores = ncores,
          use.cluster=FALSE);
      res<-t.comp$result
      motifs_list[[i]]<-res
    }

  return(motifs_list)

}

searchTFBS <- function(tfTar, file.tfs, file.twoBit, pval.cutoff.up=0.01, pval.cutoff.down=0.1, half.size=150, mTH=7, min.size=150, run.repeats=2, ncores=1 ){
  
  if(class(tfTar)!="tfTarget")
     stop("The first parameter is not a 'tfTarget' object!");

  load(file.tfs)
  if(class(tfs)!="tfbs")
     stop("The second parameter is not a 'tfbs' object!");

  if(!all(file.exists(file.twoBit)))
     stop("The 2bit file is not found!");

  deseq.table.TRE <-tfTar$tab.dif.tres
  #ncores <- tfTar$ncores;
  
  deseq.table.sig <- center.bed(deseq.table.TRE[!is.na(deseq.table.TRE$padj) & deseq.table.TRE$padj <pval.cutoff.up,], half.size, half.size)
  enh.unc.bed     <- center.bed(deseq.table.TRE[!is.na(deseq.table.TRE$padj) & deseq.table.TRE$padj> pval.cutoff.down,], half.size, half.size)
  enh.up.bed      <- deseq.table.sig[deseq.table.sig$log2FoldChange>0,]
  enh.down.bed    <- deseq.table.sig[deseq.table.sig$log2FoldChange<0,]

  # remove negative coordinates
  enh.unc.bed[,2]<-sapply(enh.unc.bed[,2],function(x)max(x,0))
  enh.up.bed[,2]<-sapply(enh.up.bed[,2],function(x)max(x,0))
  enh.down.bed[,2]<-sapply(enh.down.bed[,2],function(x)max(x,0))

  motif.list.up <- tfbs.enrichmentTest.multiple(tfs, file.twoBit,  enh.up.bed, enh.unc.bed, mTH, min.size, run.repeats, ncores);
  motif.list.down <- tfbs.enrichmentTest.multiple(tfs, file.twoBit, enh.down.bed, enh.unc.bed, mTH, min.size, run.repeats, ncores);

  tfTar$pval.cutoff.up <- pval.cutoff.up;
  tfTar$pval.cutoff.down <- pval.cutoff.down;
  tfTar$half.size <- half.size;
  tfTar$mTH <- mTH;
  tfTar$min.size <- min.size;
  tfTar$run.repeats <- run.repeats;
  tfTar$tfs <- tfs;
  tfTar$file.twoBit <- file.twoBit;
  
  tfTar$enh.up.bed <- enh.up.bed;
  tfTar$enh.down.bed <- enh.down.bed;
  tfTar$enh.unc.bed <- enh.unc.bed;
  tfTar$motif.list.up <- motif.list.up;
  tfTar$motif.list.down <- motif.list.down;

  return(tfTar);

}