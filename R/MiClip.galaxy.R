MiClip.galaxy <-
  function(file="",control=NULL,mut.type="T2C",step=5,max.hmm=100,paired=F,suffix=NULL,
           empirical="auto",model.cut=0.2,max.iterats=20,conver.cut=0.01)
{
  mic=MiClip(file=file,mut.type=mut.type,step=step,max.hmm=max.hmm,
    paired=paired,suffix=suffix,empirical=empirical,model.cut=model.cut,
    max.iterats=max.iterats,conver.cut=conver.cut)
  mic=MiClip.read(mic)
  mic=MiClip.enriched(mic,quiet=F)
  mic=MiClip.binding(mic,quiet=F)
      
  if (! is.null(control)) 
  {
    mic=MiClip.snp(mic,file=control,mut.type=mut.type,paired=paired,suffix=suffix)
    write.table(file="snps.csv",mic$snps,quote=F,row.names=F,sep=",")
  }
  
  write.table(file="enriched.csv",mic$enriched,quote=F,row.names=F,sep=",")
  write.table(file="sites.csv",mic$sites,quote=F,row.names=F,sep=",")
  write.table(file="clusters.csv",mic$clusters,quote=F,row.names=F,sep=",")
}

