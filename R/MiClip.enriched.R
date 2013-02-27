MiClip.enriched <-
function(mic,quiet=FALSE)
{
  mic=MiClip.ini1(mic)
  if (quiet==FALSE) {cat("Initialization of the first HMM finished!\n")}
  
  mic=MiClip.hmm1(mic,quiet)
  if (quiet==FALSE) {cat("Iterations of the first HMM finished!\n")}
  
  mic=MiClip.viterbi1(mic)
  if (quiet==FALSE) {cat("Viterbi algorithm of the first HMM finished!\n")}
  
  return(mic)
}
