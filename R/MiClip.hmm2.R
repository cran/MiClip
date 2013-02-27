MiClip.hmm2 <-
function(mic,quiet)
{
  max.iterats=mic$max.iterats
  conver.cut=mic$conver.cut
  enriched=mic$enriched
  gammas=mic$gammas
  regions=mic$regions
  emission=mic$emission
  xis=mic$xis
  a=mic$a
  fb=mic$fb
  len=mic$len
  bases=mic$bases
  clusters=mic$clusters
  
  ##############  initialize  ###############################
  
  iterat=0
  scales=rep(0,len)
  
  rows_back=which(!(rownames(xis) %in% regions[,"end"]))
  rows_back=rows_back[order(rows_back,decreasing=T)]
  rows_forward=which(!(rownames(xis) %in% regions[,"start"]))
  
  a_conver=a
  a_conver[]=-1
  conver=1
  
  ######################  HMM  ############################
  
  while (conver>conver.cut & iterat<max.iterats)
  { 
    # part I
    
    for (i in 1:2)
    {
      for (j in 1:2)
      {
        temp1=sum(gammas[rows_back,i])
        temp2=sum(xis[rows_back,paste(i,j,sep="")])
        a[i,j]=temp2/temp1
      }
    }
    
#    a[2,1]=1
#    a[2,2]=0
    
    # part II
    
    fb[regions$start,"f1"]=gammas[regions$start,1]*emission[regions$start,1]
    fb[regions$start,"f2"]=gammas[regions$start,2]*emission[regions$start,2]
    
    scales[regions$start]=1/apply(fb[regions$start,c("f1","f2")],1,sum)
    fb[regions$start,c("f1","f2")]=fb[regions$start,c("f1","f2")]*scales[regions$start]
    
    for (i in rows_forward)
    {
      fb[i,c("f1","f2")]=(fb[i-1,c("f1","f2")] %*% a) * emission[i,]  
      scales[i]=1/sum(fb[i,c("f1","f2")])
      fb[i,c("f1","f2")]=fb[i,c("f1","f2")]*scales[i]
    }
    
    fb[regions$end,c("b1","b2")]=1
    fb[regions$end,c("b1","b2")]=fb[regions$end,c("b1","b2")]*scales[regions$end]
    
    for (i in rows_back)
    {
      fb[i,c("b1","b2")]=(fb[i+1,c("b1","b2")] * emission[i+1,]) %*% t(a)
      fb[i,c("b1","b2")]=fb[i,c("b1","b2")]*scales[i]
    }
    
    # part III
    
    gammas[,1]=fb[,"f1"]*fb[,"b1"]
    gammas[,2]=fb[,"f2"]*fb[,"b2"]
    gammas=gammas/scales
    
    for (i in 1:2)
    {
      for (j in 1:2)
      {
        xis[1:(len-1),paste(i,j,sep="")]=a[i,j]*emission[2:len,j]*fb[1:(len-1),paste("f",i,sep="")]*fb[2:len,paste("b",j,sep="")]
      }
    }
    
    # check convergence
    
    iterat=iterat+1
    
    conver=max(abs(a-a_conver))
    a_conver=a
    
    if (quiet==FALSE) {cat(">")}
  }
  
  if (quiet==FALSE) {cat("\n")}
  
  #################  save data  ##################################
  
  result=list(enriched=enriched,gammas=gammas,regions=regions,clusters=clusters,
    emission=emission,a=a,len=len,bases=bases,rows_forward=rows_forward,
    rows_back=rows_back)
  class(result)="MiClip"
  
  return(result)
}
