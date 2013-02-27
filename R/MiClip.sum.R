MiClip.sum <-
function(mic,...)
{
  enriched=mic$enriched
  sites=mic$sites
  if (is.null(sites)) {return(cat("Run MiClip.binding() first!"))}
    
  summary.enriched=c()
  summary.sites=c()
    
  summary.enriched$num.cluster=length(unique(enriched$region_id))
  summary.enriched$num.enriched.cluster=length(unique(enriched$region_id[enriched$enriched]))
  summary.enriched$num.bins=dim(enriched)[1]
  summary.enriched$state.count=table(enriched$enriched)
  summary.enriched$prob.stat=summary(enriched$prob,...)
  summary.enriched$mean.tag.enriched=round(mean(enriched$tag[enriched$enriched]),...)
  summary.enriched$mean.tag.notenriched=round(mean(enriched$tag[! enriched$enriched]),...)
      
  cat("For identifying enriched regions\n")
  cat(paste("# of clusters:",summary.enriched$num.cluster,"\n"))
  cat(paste("# of identified enriched clusters:",summary.enriched$num.enriched.cluster,"\n"))
  cat(paste("# of bins:",summary.enriched$num.bins,"\n"))
  cat("# of bins in each state:")
  print(summary.enriched$state.count)
  cat("Statistics of probability\n")
  print(summary.enriched$prob.stat)
  cat(paste("Average tag count of enriched bin:",summary.enriched$mean.tag.enriched,"\n"))
  cat(paste("Average tag count of not enriched bin:",summary.enriched$mean.tag.notenriched,"\n"))
    
  cat("\n\n")
    
  summary.sites$num.cluster=length(unique(sites$region_id))
  summary.sites$num.sub_regions=length(unique(sites$sub_region_id))
  summary.sites$num.site.cluster=length(unique(sites$region_id[sites$sites]))
  summary.sites$num.bases=dim(sites)[1]
  summary.sites$state.count=table(sites$sites)
  summary.sites$prob.stat=summary(sites$prob,...)
  summary.sites$mean.tag.binding=round(mean(sites$tag[sites$sites]),...)
  summary.sites$mean.mutant.binding=round(mean(sites$mutant[sites$sites]),...)
  summary.sites$mean.tag.notbinding=round(mean(sites$tag[! sites$sites]),...)
  summary.sites$mean.mutant.notbinding=round(mean(sites$mutant[! sites$sites]),...)
      
  cat("For identifying binding sites\n")
  cat(paste("# of enriched clusters:",summary.sites$num.cluster,"\n"))
  cat(paste("# of sub enriched clusters:",summary.sites$num.sub_regions,"\n"))
  cat(paste("# of enriched clusters with identified binding sites:",summary.sites$num.site.cluster,"\n"))
  cat(paste("# of bases:",summary.sites$num.bases,"\n"))
  cat("# of bases in each state:")
  print(summary.sites$state.count)
  cat("Statistics of probability\n")
  print(summary.sites$prob.stat)
  cat(paste("Average tag count of binding site:",summary.sites$mean.tag.binding,"\n"))
  cat(paste("Average mutant count of binding site:",summary.sites$mean.mutant.binding,"\n"))
  cat(paste("Average tag count of not binding site:",summary.sites$mean.tag.notbinding,"\n"))
  cat(paste("Average mutant count of not binding site:",summary.sites$mean.mutant.notbinding,"\n"))
    
  cat("\n\n")
}