MiClip.snp <-
function(mic,file="",mut.type="T2C",paired=F,suffix=NULL)
{
  ##########  read snp data  ###################

  if (file=="") {stop("No control alignment file provided!")}
  if (paired==T & is.null(suffix)) {suffix=c("1","2")}
  
  control=MiClip(file=file,mut.type=mut.type,paired=paired,suffix=suffix)
  control=MiClip.read(control)   #  read control experiments
  control=control$raw         
  
  ts=as.vector(unlist(t(control[,5:9])))   # expand bins to single bps
  ms=as.vector(unlist(t(control[,10:14])))
  snps=data.frame(tag=ts,mutant=ms)
  snps$chr=rep(control$V2,each=5)
  snps$strand=rep(control$V4,each=5)
  snps$pos=0
  snps[1+5*(0:(dim(control)[1]-1)),"pos"]=control$V3
  
  for (i in 2:5)
  {
    snps[i+5*(0:(dim(control)[1]-1)),"pos"]=snps[i-1+5*(0:(dim(control)[1]-1)),"pos"]+1
  }
  
  ##############  test null hypothesis  #################
  
  snps=snps[snps$mutant>0 & snps$tag>=3,] # test null hypothesis of being SNPs
  p0=sum(snps$mutant)/sum(snps$tag)
  snps$snp=FALSE
  options(warn=-1)
  for (i in 1:dim(snps)[1])
  {
    test=prop.test(x=snps$mutant[i],n=snps$tag[i],p=p0,alternative="greater")
    if (test$p.value<=0.1) {snps$snp[i]=TRUE}
  }
  options(warn=0)
  snps=snps[snps$snp,]   
  snps$string=paste(snps$chr,snps$strand,snps$pos,sep="_")
  
  ##########  prepare binding site matrix  #################
  
  mic$sites$SNP=FALSE      # step up the columns to store SNP information
  mic$clusters$SNP=FALSE
  sites=mic$sites   
  sites$row=seq(dim(sites)[1])
  sites=sites[sites$sites,c("row","region_id","chr","strand","pos")]
  sites$string=paste(sites$chr,sites$strand,sites$pos,sep="_")
  sites=sites[sites$string %in% snps$string,]
  
  ##########  mark SNPs  #######################
  
  for (i in 1:dim(sites)[1])
  {
    row=sites$row[i]
    region_id=sites$region_id[i]
    mic$sites$SNP[row]=TRUE
    mic$clusters$SNP[region_id]=TRUE
  }
  
  ###############  return results  #############
  
  mic$snps=snps[,c("chr","strand","pos","tag","mutant")]
  class(mic)="MiClip"
  return(mic)
}
