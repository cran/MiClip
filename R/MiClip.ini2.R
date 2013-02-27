MiClip.ini2 <-
function(mic)
{
  ################  read data  ##############################
  
  raw=mic$raw
  max.hmm=mic$max.hmm
  fit1=mic$fit[1]
  model.cut=mic$model.cut
  max.iterats=mic$max.iterats
  conver.cut=mic$conver.cut
  enriched=mic$enriched
  clusters=mic$regions
  
  clusters=clusters[,c("mark","chr","strand","cluster_start","cluster_end")]
  colnames(clusters)=c("region_id","chr","strand","start","end")
  
  bins=raw
  step=bins[2,3]-bins[1,3]
  bins=cbind(bins,enriched$enriched)
  colnames(bins)=c("region_id","chr","start","strand",paste("t",1:step,sep=""),paste("m",1:step,sep=""),"state")
  bins=bins[bins$state==2,]
  
  ###############  fit model  ##############################
  
  ts=unlist(bins[,paste("t",1:step,sep="")])
  ms=unlist(bins[,paste("m",1:step,sep="")])
  ts[ts>max.hmm]=max.hmm
  ms[ms>max.hmm]=max.hmm
  ps=ms/ts
  
  ts=ts[! is.na(ps)]
  ms=ms[! is.na(ps)]
  ps=ps[! is.na(ps)]
  
  ts1=ts[ps<model.cut] # for zero-inflated binomial 
  ms1=ms[ps<model.cut]
  ts2=ts[ps>=model.cut] # for binomial
  ms2=ms[ps>=model.cut] 
  
  # fit zero-inflated binomial
  
  options(warn=-1)
  fit=try(vglm(cbind(ms1,ts1-ms1)~1,family=zibinomial,maxit=50),silent=TRUE)
  options(warn=0)
  
  if (length(fit)==2)
  {
    pstr0=Coef(fit)[1]
    prob=Coef(fit)[2]
  }else
  {
    fit=MiClip.grid(ms1,ts1-ms1,model.cut)
    pstr0=fit[1]
    prob=fit[2]
  }
  
  # fit binomial
  
  fit=glm(cbind(ms2,ts2-ms2)~1,binomial(link=logit))
  intercept=VGAM::logit(as.vector(fit$coef),inverse=T)
  
  if (pstr0>=1 || pstr0<=0 || prob>=1 || prob<=0 || intercept>=1 || intercept<=0) {stop("Model fails to fit for the second HMM!")}
  
  #################  expand bins  ###############################
  
  ts=as.vector(unlist(t(bins[,paste("t",1:step,sep="")])))
  ms=as.vector(unlist(t(bins[,paste("m",1:step,sep="")])))
  
  bases=data.frame(tag=ts,mutant=ms)
  bases$region_id=rep(bins$region_id,each=step)
  bases$chr=rep(bins$chr,each=step)
  bases$strand=rep(bins$strand,each=step)
  bases$pos=0
  bases[1+step*(0:(dim(bins)[1]-1)),"pos"]=bins$start
  
  for (i in 2:step)
  {
    bases[i+step*(0:(dim(bins)[1]-1)),"pos"]=bases[i-1+step*(0:(dim(bins)[1]-1)),"pos"]+1
  }
  
  len=dim(bases)[1]
  marks=rep(0,len)
  marks[bases$tag<fit1]=-1
  mark=0
  
  jump_id=bases$region_id[2:len]!=bases$region_id[1:(len-1)]
  jump_pos=bases$pos[2:len]!=bases$pos[1:(len-1)]+1
  jump_tag=bases$tag[2:len]>=fit1 & bases$tag[1:(len-1)]<fit1
  jump=c(1,1*(jump_id | jump_pos | jump_tag))
  
  for (i in 1:len)
  {
    if (marks[i]>=0)
    {
      mark=mark+jump[i]
      marks[i]=mark
    }
  }
  
  bases$mark=marks
  
  bases=bases[bases$mark>0,]
  len=dim(bases)[1]
  
  #################  calculate emission table  ##################
  
  # 1 means in state 2 (not binding), 2 means in state 2 (binding)
  
  dist_table1=matrix(data=1,ncol=max.hmm+1,nrow=max.hmm+1) # i is tag count, j is mutant count
  dist_table2=matrix(data=1,ncol=max.hmm+1,nrow=max.hmm+1)
  
  for (i in 1:(max.hmm+1))
  {
    for (j in 1:i)
    {
      dist_table1[i,j]=choose(i-1,j-1)*prob^(j-1)*(1-prob)^(i-j)
      dist_table2[i,j]=choose(i-1,j-1)*intercept^(j-1)*(1-intercept)^(i-j)
    }
  }
  
  dist_table1=(1-pstr0)*dist_table1
  dist_table1[,1]=dist_table1[,1]+pstr0
  
  if (range(dist_table1)[1]<1e-300 || range(dist_table2)[1]<1e-300)
  {
    cat("emission table is out of the dynamic range of R!\n")
    stop("Set max # of tags allowed at each bin to a smaller value!")
  }
  
  emission=matrix(data=0,ncol=2,nrow=len) 
  # first column is state 2 (not binding), second column is state 2 (binding)
  
  for (i in 1:len)
  {
    tag=bases$tag[i]
    mutant=bases$mutant[i]
    if (tag>max.hmm) {tag=max.hmm}
    if (mutant>max.hmm) {mutant=max.hmm}
    emission[i,]=c(dist_table1[tag+1,mutant+1],dist_table2[tag+1,mutant+1])
  }
  
  #################  calculate initial gamma table  #############
  
  gammas=matrix(data=0,ncol=2,nrow=len)
  colnames(gammas)=1:2
  rownames(gammas)=1:len
  
  gammas[,1]=0.1+0.8*(bases$mutant==0)
  gammas[,2]=0.1+0.8*(bases$mutant>0)
  
  ##################  region table  ##############################
  
  regions=data.frame(mark=unique(bases$mark))
  regions$start=c(1,1+which(bases$mark[2:len]!=bases$mark[1:(len-1)]))
  regions$end=c(which(bases$mark[2:len]!=bases$mark[1:(len-1)]),len)
  
  ##################  initialize other parameters  ###############
  
  xis=as.data.frame(matrix(data=0,ncol=4,nrow=len))
  colnames(xis)=c("11","12","21","22")
  rownames(xis)=c(1:len)
  
  for (i in 1:2)
  {
    for (j in 1:2)
    {
      xis[1:(len-1),paste(i,j,sep="")]=gammas[1:(len-1),i]*gammas[2:len,j]
    }
  }
  
  a=matrix(data=0,ncol=2,nrow=2)
  rownames(a)=c("1","2")
  colnames(a)=c("1","2") # a[i,j] means the probability from i to j
  
  fb=as.data.frame(matrix(data=0,ncol=4,nrow=len))
  colnames(fb)=c("f1","f2","b1","b2")
  rownames(fb)=c(1:len)
  
  ##########################  save data  #############################
  
  emission=as.matrix(emission)
  xis=as.matrix(xis)
  fb=as.matrix(fb)
  
  result=list(max.iterats=max.iterats,conver.cut=conver.cut,clusters=clusters,
              enriched=enriched,gammas=gammas,regions=regions,emission=emission,
              xis=xis,a=a,fb=fb,len=len,bases=bases)
  class(result)="MiClip"
  
  return(result)
}
