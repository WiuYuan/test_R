setwd("~/Desktop/test_bio/bdsfile/chapter-08-r")
d = read.csv("Dataset_S1.txt")
mean(d$depth)
summary(d$depth)
d[, "start", drop=FALSE]
d$Pi
which(d$Pi > 10)[1:4]
which.min(d$total.Bases)
library(ggplot2)
d$position <- (d$end + d$start) / 2
d$cent <- d$start >= 25800000 & d$end <= 29700000
d$diversity <- d$Pi / (10*1000)
colnames(d)[12] <- "percent.GC"
ggplot(d) + geom_point(aes(x=position, y=diversity))
ggplot(d) + geom_point(aes(x=position, y=diversity), alpha=0.01)
ggplot(d) + geom_density(aes(x=diversity), fill="black")
ggplot(d) + geom_density(aes(x=diversity, fill=cent), alpha=0.4)
ggplot(d, aes(x=depth, y=total.SNPs)) + geom_point() + geom_smooth()
ggplot(d, aes(x=percent.GC, y=depth)) + geom_point() + geom_smooth()

d$GC.binned <- cut(d$percent.GC, 5)
d$GC.binned
ggplot(d) + geom_bar(aes(x=GC.binned))
ggplot(d) + geom_bar(aes(x=percent.GC))
ggplot(d) + geom_density(aes(x=depth, linetype=GC.binned), alpha=0.5)
ggplot(d) + geom_bar(aes(x=Pi), binwidth=1) + scale_x_continuous(limits=c(0.01, 80))


reps <- read.delim("chrX_rmsk.txt.gz", header=TRUE)
head(reps, 3)
common_repclass <- c("SINE", "LINE", "LTR", "DNA", "Simple_repeat")
mtfs <- read.delim("motif_recombrates.txt", header=TRUE)
rpts <- read.delim("motif_repeats.txt", header=TRUE)

mtfs$pos <- paste(mtfs$chr, mtfs$motif_start, sep="-")
rpts$pos <- paste(rpts$chr, rpts$motif_start, sep="-")
i <- match(mtfs$pos, rpts$pos)
mtfs$repeat_name <- rpts$name[i]
mtfs_inner <- mtfs[!is.na(mtfs $repeat_name), ]
recm <- merge(mtfs, rpts, by.x="pos", by.y="pos")
recm <- merge(mtfs, rpts, by.x="pos", by.y="pos", all.x=TRUE)
p <- ggplot(mtfs, aes(x=dist, y=recom)) + geom_point(size=1)
p <- p + geom_smooth(method="loess", se=FALSE, span=1/10)
print(p)

ggplot(mtfs, aes(x=dist, y=recom)) + geom_point(size=1) +
  geom_smooth(aes(color=motif), method="loess", se=FALSE, span=1/10)

p <- ggplot(mtfs, aes(x=dist, y=recom)) + geom_point(size=1, color="grey")
p <- p + geom_smooth(method='loess', se=FALSE, span=1/16)
p <- p + facet_wrap( ~ motif, scales="free_y")
print(p)




d_split <- split(d$depth, d$GC.binned)
dpth_summ <- lapply(split(d$depth, d$GC.binned), summary)


install.packages("dplyr")
library(dplyr)
d_df <- tbl_df(d)
d_df
d
dd = select(d, start, end, Pi, Recombination, depth)
select(d_df, -(start:cent))
filter(d_df, Pi > 16, percent.GC > 80)
arrange(d_df, depth)
d_df %>% mutate(GC.scaled = scale(percent.GC)) %>% filter(GC.scaled > 4, depth > 4) %>%
  select(start, end, depth, GC.scaled, percent.GC) %>%
  arrange(desc(depth))
dd = d_df %>% mutate(GC.scaled = scale(percent.GC))
mtfs_df <- tbl_df(mtfs)
dd = mtfs_df %>% group_by(chr)

chrs <- c("chrom6", "chr2 2", "chr6", "chr4", "chr1", "chr16", " chrom8")
regexpr("[^\\d]6", chrs, perl=TRUE)
pos <- regexpr("\\d+", chrs, perl=TRUE)
pos
sub(pattern="1", replacement="a", x="1231")
sub("gene=(\\w)(\\d)", "\\1 \\2", "gene=a1", perl=TRUE)
sub(">[^ ]+ *(.*)", "\\1", ">1  length=301354135 type=dna")
sub("(.*)(\\d+)", "\\1 \\2", "chr19chr2", perl=TRUE)
paste("chr", c(1:22, "X", "Y"), sep="_")
leafy <- "gene=LEAFY;locus=2159208;gene_model=AT5G61850.1"
dd = strsplit(leafy, ";")
class(dd)
getwd()
list.files("hotspots", pattern="hotspots.*\\.bed")
hs_files <- list.files("hotspots", pattern="hotspots.*\\.bed", full.names=TRUE)
bedcols <- c("chr", "start", "end")
loadFile <- function(x) read.delim(x, header=FALSE, col.names=bedcols)
hs <- lapply(hs_files, loadFile)
names(hs) <- list.files("hotspots", pattern="hotspots.*\\.bed")
hsd <- do.call(rbind, hs)


setwd("~/desktop/test_bio/test_R")
write.table(mtfs, file="hotspot_motifs.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
hs_gzf <- gzfile("hotspot_motifs.txt.gz")
write.table(mtfs, file=hs_gzf, quote=FALSE, sep="\t", row.names=FALSE,
            col.names=TRUE)
ls()
