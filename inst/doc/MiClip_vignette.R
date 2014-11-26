### R code from vignette source 'MiClip_vignette.Rnw'

###################################################
### code chunk number 1: MiClip_vignette.Rnw:63-67
###################################################
library("MiClip")

MiClip.adaptor(file=system.file("extdata/test.fastq",package="MiClip"),
                adaptor="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC")


###################################################
### code chunk number 2: MiClip_vignette.Rnw:120-126
###################################################
library("MiClip")

test=MiClip(file=system.file("extdata/test.sam",package="MiClip"),mut.type="Del")

# for paired-end data
# test=MiClip(file="test.sam",paired=TRUE,suffix=c("F3","F5-RNA"))


###################################################
### code chunk number 3: MiClip_vignette.Rnw:137-138
###################################################
test=MiClip.read(test) # read raw data


###################################################
### code chunk number 4: MiClip_vignette.Rnw:147-148
###################################################
test=MiClip.enriched(test,quiet=FALSE) # identify enriched regions


###################################################
### code chunk number 5: MiClip_vignette.Rnw:159-160
###################################################
test=MiClip.binding(test,quiet=FALSE) # identify binding sites


###################################################
### code chunk number 6: MiClip_vignette.Rnw:171-172
###################################################
test=MiClip.snp(test,file=system.file("extdata/snp.sam",package="MiClip"),mut.type="Del")


###################################################
### code chunk number 7: MiClip_vignette.Rnw:185-197
###################################################
enriched=test$enriched # test will contain at least three data frames
sites=test$sites
clusters=test$clusters

head(enriched) # view these data frames
head(sites)
head(clusters)

head(enriched[enriched$enriched,]) # view enriched bins
head(sites[sites$sites,]) # view binding sites
head(clusters[clusters$enriched,]) # view clusters with enriched bins
head(clusters[clusters$sites,]) # view clusters with binding sites


###################################################
### code chunk number 8: MiClip_vignette.Rnw:212-217
###################################################
snps=test$snps
head(snps) # Inferred possible SNP sites are contained in this data frame

head(sites[sites$SNP,]) # In this dataset, three possible SNPs are found
head(clusters[clusters$SNP,])


###################################################
### code chunk number 9: MiClip_vignette.Rnw:226-227
###################################################
MiClip.sum(test)


###################################################
### code chunk number 10: session
###################################################
sessionInfo()


