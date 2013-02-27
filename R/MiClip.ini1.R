MiClip.ini1 <-
function(mic)
{
  ################  read data  ##############################

  raw=mic$raw
  max.hmm=mic$max.hmm
  empirical=mic$empirical
  model.cut=mic$model.cut
  max.iterats=mic$max.iterats
  conver.cut=mic$conver.cut
  step=mic$step
  
  bins=raw
  step=bins[2,3]-bins[1,3]
  
  bins[,5]=apply(bins[,5:(4+step)],1,sum)
  bins=bins[,1:5]

  colnames(bins)=c("mark","chr","start","strand","tag")
  bins$end=bins$start+step-1
  
  #############  check bin number per sequence  ###############
  
  n_bins=dim(raw)[1]/length(unique(raw[,1]))
    
  if (n_bins>15)
  {
    cat("The step size might be too small. ")
    cat(paste("Suggested step size",round(step*n_bins/10),"\n"))
  }
  
  #############  generate regions table  ######################

  len=dim(bins)[1]

  regions=as.data.frame(matrix(data=0,ncol=3,nrow=length(unique(bins$mark))))
  colnames(regions)=c("mark","start","end")

  regions$mark=unique(bins$mark)
  regions$start=c(1,which(bins$mark[1:(len-1)]!=bins$mark[2:len])+1)
  regions$end=c(which(bins$mark[1:(len-1)]!=bins$mark[2:len]),len)
  regions$sum=aggregate(bins$tag,by=list(bins$mark),sum)[,2]
  colnames(regions)=c("mark","start","end","sum")

  regions$chr=bins$chr[regions$start]
  regions$strand=bins$strand[regions$start]
  regions$cluster_start=bins$start[regions$start]
  regions$cluster_end=bins$start[regions$end]+step-1
  
  ###############  define fit function  #####################

  fit2pois=function(x)
  {
    temp=all.moments(x,order.max=3)
  
    a=temp[2]*temp[2]+temp[2]-temp[3]
    if (a==0) {return (-1)}
    b=temp[2]*temp[2]-temp[2]*temp[3]+2*temp[2]-3*temp[3]+temp[4]
    c=temp[3]*temp[3]-temp[2]*temp[2]+temp[2]*temp[3]-temp[2]*temp[4]
  
    if(b*b-4*a*c<=0) {return (-2)}
    theta1=(b+sqrt(b*b-4*a*c))/(-2*a)
    theta2=(b-sqrt(b*b-4*a*c))/(-2*a)
    alpha=1-(temp[2]-theta2)/(theta1-theta2)
    if (alpha<=0 || alpha>=1) {return (-3)}
    if (theta1<=0 && theta2<=0) {return (-4)}
    if (theta1<=0) {theta1=min(0.1,theta2/2)}
    if (theta2<=0) {theta2=min(0.1,theta1/2)}
    
    if (theta1<theta2) 
    {
      return (c(theta1,theta2))
    }else 
    {
      return (c(theta2,theta1))
    }
  }

  #################  fit model  #################################
  
  if (empirical=="auto") 
  {
    empirical=quantile(bins$tag,0.995)
  }else
  {
    empirical=as.numeric(empirical)
    empirical=2*empirical*step
  }

  fit=fit2pois(bins$tag[bins$tag<empirical])
  if (length(fit)==1) {stop("Model fitting fails for 1st HMM")}
  
  fit=fit/step
  bins$tag=round(bins$tag/step)
  
  bins_raw=bins
  bins$tag[bins$tag>max.hmm]=max.hmm
  
  #################  calculate emission table  ##################

  # 1 means in state 0, 2 means in state 1

  pois_table=matrix(data=0,nrow=max.hmm+1,ncol=2)
  pois_table[1,1]=exp(-fit[1])
  pois_table[1,2]=exp(-fit[2])

  for (i in 2:(max.hmm+1))
  {
    pois_table[i,1]=pois_table[i-1,1]*fit[1]/(i-1)
    pois_table[i,2]=pois_table[i-1,2]*fit[2]/(i-1)
  }

  if (range(pois_table)[1]<1e-300)
  {
    cat("emission table is out of the dynamic range of R!\n")
    stop("Set max # of tags allowed at each bin to a smaller value!\n")
  }

  test=min(round(1.5*fit[2]),max.hmm)
  if (pois_table[test,1]>=pois_table[test,2])
  {
    cat("Warning: model fitting may not be appropriate!\n")
    cat("Set max # of tags allowed at each bin to a smaller value!\n")
  }

  s1=pois_table[bins$tag+1,1]
  s2=pois_table[bins$tag+1,2]

  emission=data.frame(s1=s1,s2=s2) # emission table
  colnames(emission)=c(1,2)
  rownames(emission)=1:len

  #################  calculate initial gamma table  #############

  gammas=matrix(data=0,ncol=2,nrow=len)
  colnames(gammas)=1:2
  rownames(gammas)=1:len

  regions$sum=regions$sum/(regions$end-regions$start+1)
  ave=rep(regions$sum,regions$end-regions$start+1)

  gammas[,1]=0.1+0.8*(bins$tag<fit[2])
  gammas[,2]=0.1+0.8*(bins$tag>=fit[2])

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
  
  result=list(raw=raw,max.hmm=max.hmm,fit=fit,
              model.cut=model.cut,max.iterats=max.iterats,
              conver.cut=conver.cut,gammas=gammas,regions=regions,
              emission=emission,xis=xis,a=a,fb=fb,len=len,bins=bins_raw)
  class(result)="MiClip"
  
  return(result)
}
