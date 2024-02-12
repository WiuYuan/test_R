library(qrqc)
fqfiles <- c(none="untreated1_chr4.fq", sickle="untreated1_chr4_sickle.fq",
             trimfq="untreated1_chr4_trimfq.fq")
setwd("~/desktop/test_bio/bdsfile/chapter-10-sequence")
seq_info <- lapply(fqfiles, function(file) {
  readSeqFile(file, hash=FALSE, kmer=FALSE) })
quals <- mapply(function(sfq, name) {
  qs <- getQual(sfq)
  qs$trimmer <- name
  qs
}, seq_info, names(fqfiles), SIMPLIFY=FALSE)

d <- do.call(rbind, quals)
# Visualize qualities
p1 <- ggplot(d) + geom_line(aes(x=position, y=mean, linetype=trimmer))
p1 <- p1 + ylab("mean quality (sanger)") + theme_bw()
print(p1)
# Use qrqc's qualPlot with list produces panel plots
# Only shows 10% to 90% quantiles and lowess curve
p2 <- qualPlot(seq_info, quartile.color=NULL, mean.color=NULL) + theme_bw() 
p2 <- p2 + scale_y_continuous("quality (sanger)")
print(p2)