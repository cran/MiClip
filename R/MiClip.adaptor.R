MiClip.adaptor<-
function(file="",format="fastq",adaptor="",min=15,mismatch=0.4)
{
  if (file=="" || adaptor=="") {stop("File or adaptor sequence is empty!")}
  
  ##############  check perl  ###############################
  
  status=system("perl -e \"1+1\"",ignore.stdout=T)
  if (status!=0) {stop("It seems Perl is not intalled!\n")}
  perlpath=system.file("exec",package="MiClip")
  
  ###############  trim adaptor  #############################
  
  command=paste("perl \"",perlpath,"/remove_adaptor.pl\" ",format," \"",file,
  "\" ",adaptor," ",min," ",mismatch,sep="")
  status=system(command,intern=FALSE) 
  if (status!=0) {stop("Adaptor removing failed! Possible corrupted sequencing file!\n")}
}
