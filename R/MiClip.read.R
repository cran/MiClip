MiClip.read <-
function(mic)
{
  file=mic$file
  mut.type=mic$mut.type
  step=mic$step
  max.hmm=mic$max.hmm
  empirical=mic$empirical
  model.cut=mic$model.cut
  max.iterats=mic$max.iterats
  conver.cut=mic$conver.cut
  paired=mic$paired
  suffix=mic$suffix
  background=mic$background
  
  ##############  check perl  ###############################
  
  status=system("perl -e \"1+1\"",ignore.stdout=T)
  if (status!=0) {stop("It seems Perl is not intalled!\n")}
  perlpath=system.file("exec",package="MiClip")
  tmpdir=tempdir()
  
  ##############  execute perl scripts  ###############################
  
  if (paired==F)
  {  
    #########  form clusters  ######################################
    
    cluster_file=tempfile(pattern="cluster.",tmpdir,fileext=".bed")
    cluster_file=sub('\\\\','/',cluster_file)
    command=paste("perl \"",perlpath,"/cluster.pl\" \"",file,"\" \"",cluster_file,"\"",sep="")
    status=system(command,intern=FALSE) # form clusters
    if (status!=0) {stop("Forming cluster file has non-zero exit status!\n")}
    cat("Identifying clusters finished!\n")
    
    #########  read tags and mutants  ###############################
    
    bin_file=tempfile(pattern="bin.",tmpdir,fileext=".bed") # find paths
    bin_file=sub('\\\\','/',bin_file)
    command=paste("perl \"",perlpath,"/preprocess.pl\" ",step," \"",cluster_file,
      "\" \"",mut.type,"\" \"",bin_file,"\"",sep="")
    status=system(command,intern=FALSE) # count base coverage
    if (status!=0) {stop("Generating bin file has non-zero exit status!\n")}
    cat("Generating bin file finished!\n")
  }else
  {
    #########  merge paired-end data  ##############################
    
    merge_file=tempfile(pattern="merge.",tmpdir,fileext=".bed")
    merge_file=sub('\\\\','/',merge_file)
    command=paste("perl \"",perlpath,"/merge_pair.pl\" ",suffix[1]," ",suffix[2],
                  " \"",file,"\" \"",merge_file,"\"",sep="")
    status=system(command,intern=FALSE) # merge paired-end read
    if (status!=0) {stop("Merging paired-end sequencing file has non-zero exit status!\n")}
    cat("Merging mate pairs finished!\n")
    
    #########  form clusters  ######################################
    
    cluster_file=tempfile(pattern="cluster.",tmpdir,fileext=".bed")
    cluster_file=sub('\\\\','/',cluster_file)
    command=paste("perl \"",perlpath,"/cluster_p.pl\" \"",merge_file,"\" \"",cluster_file,"\"",sep="")
    status=system(command,intern=FALSE) # form clusters
    if (status!=0) {stop("Forming cluster file has non-zero exit status!\n")}
    cat("Identifying clusters finished!\n")
    
    #########  read tags and mutants  ###############################
    
    bin_file=tempfile(pattern="bin.",tmpdir,fileext=".bed") # find paths
    bin_file=sub('\\\\','/',bin_file)
    command=paste("perl \"",perlpath,"/preprocess_p.pl\" ",step," \"",cluster_file,
                  "\" \"",mut.type,"\" \"",bin_file,"\"",sep="")
    status=system(command,intern=FALSE) # count base coverage
    if (status!=0) {stop("Generating bin file has non-zero exit status!\n")}
    cat("Generating bin file finished!\n")
  }
  
  ##############  read in data  ##################################
  
  raw=read.table(bin_file)
  unlink(cluster_file) 
  unlink(bin_file) # read and delete intermediate files
  
  result=list(raw=raw,max.hmm=max.hmm,empirical=empirical,model.cut=model.cut,
    max.iterats=max.iterats,conver.cut=conver.cut,step=step,scale=scale,
    background=background)
  class(result)="MiClip"
  
  return(result)
}
