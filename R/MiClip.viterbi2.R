MiClip.viterbi2 <-
function(mic)
{
  enriched=mic$enriched
  gammas=mic$gammas
  regions=mic$regions
  emission=mic$emission
  a=mic$a
  len=mic$len
  bases=mic$bases
  rows_forward=mic$rows_forward
  rows_back=mic$rows_back
  clusters=mic$clusters
  
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
  
  sites=cbind(bases[,c("region_id","mark","strand","chr","pos","tag","mutant")],states)
  colnames(sites)=c("region_id","sub_region_id","strand","chr","pos","tag","mutant","sites","probability")
  
  sites$probability=round(sites$probability,digits=3)
  enriched$probability=round(enriched$probability,digits=3)
  
  result=list(enriched=enriched,sites=sites,clusters=clusters)
  class(result)="MiClip"
  
  return(result)
}
