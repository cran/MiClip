MiClip.binding <-
function(mic,quiet=FALSE)
{
  mic=MiClip.ini2(mic)
  if (quiet==FALSE) {cat("Initialization of the second HMM finished!\n")}
  
  mic=MiClip.hmm2(mic,quiet)
  if (quiet==FALSE) {cat("Iterations of the second HMM finished!\n")}
  
  mic=MiClip.viterbi2(mic)
  if (quiet==FALSE) {cat("Viterbi algorithm of the second HMM finished!\n")}
  
  mic$enriched$enriched=(mic$enriched$enriched==2)
  mic$sites$sites=(mic$sites$sites==2)
  
  mic$clusters$enriched=F
  mic$clusters$sites=F
  mic$clusters$enriched[unique(mic$enriched$region_id[mic$enriched$enriched])]=T
  mic$clusters$sites[unique(mic$sites$region_id[mic$sites$sites])]=T
  
  return(mic)
}
