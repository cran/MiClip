MiClip <-
function(file="",mut.type="T2C",step=5,max.hmm=100,paired=F,suffix=NULL,
  empirical="auto",model.cut=0.2,max.iterats=20,conver.cut=0.01,background=NULL)
{
  if (file=="") {stop("No input alignment file provided!")}
  if (paired==T & is.null(suffix)) {suffix=c("1","2")}
  
  result=list(file=file,mut.type=mut.type,step=step,paired=paired,
    max.hmm=max.hmm,empirical=empirical,suffix=suffix,background=background,
    model.cut=model.cut,max.iterats=max.iterats,conver.cut=conver.cut)
  class(result)="MiClip"
  
  return(result)
}
