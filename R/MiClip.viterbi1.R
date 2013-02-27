MiClip.viterbi1 <-
function(mic)
{
  raw=mic$raw
  max.hmm=mic$max.hmm
  fit=mic$fit
  model.cut=mic$model.cut
  max.iterats=mic$max.iterats
  conver.cut=mic$conver.cut
  gammas=mic$gammas
  regions=mic$regions
  emission=mic$emission
  a=mic$a
  len=mic$len
  bins=mic$bins
  rows_back=mic$rows_back
  rows_forward=mic$rows_forward
  
  #################  initialize  ####################################
  
  delta=matrix(data=0,ncol=2,nrow=len)
  colnames(delta)=c(1:2)
  rownames(delta)=1:len
  
  phi=matrix(data=0,ncol=2,nrow=len)
  colnames(phi)=c(1:2)
  rownames(phi)=1:len
  
  delta[regions$start,1]=gammas[regions$start,1]*emission[regions$start,1]
  delta[regions$start,2]=gammas[regions$start,2]*emission[regions$start,2]
  
  #################  viterbi  ##########################################
  
  for (i in rows_forward)
  {
    temp=a*delta[i-1,]
    delta[i,]=apply(temp,2,max)*emission[i,]
    delta[i,]=delta[i,]/sum(delta[i,])
    phi[i,]=apply(temp,2,function(x) return(which(x==max(x))[1]))
  }
  
  states=matrix(data=0,ncol=2,nrow=len)
  colnames(states)=c("state","probability")
  rownames(states)=1:len
  
  states[regions$end,"state"]=apply(delta[regions$end,],1,function(x) return(which(x==max(x))[1]))
  
  for (i in rows_back)
  {
    states[i,"state"]=phi[i+1,states[i+1,"state"]]
  }
  
  for (i in 1:len)
  {
    states[i,"probability"]=gammas[i,states[i,"state"]]
  }
  
  ##############  format output  ###########################
  
  output=cbind(bins[,c("mark","chr","start","end","strand","tag")],states)
  colnames(output)=c("region_id","chr","start","end","strand","tag","enriched","probability")
  
  result=list(raw=raw,max.hmm=max.hmm,fit=fit,regions=regions,
    model.cut=model.cut,max.iterats=max.iterats,conver.cut=conver.cut,
    enriched=output)
  class(result)="MiClip"
  
  return(result)
}
