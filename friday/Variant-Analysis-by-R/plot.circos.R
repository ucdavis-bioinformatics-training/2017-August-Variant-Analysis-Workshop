# simple pie chart
  ## First, input the number of each type of variants after filtering using quality score of 60.
slices <- c(8539874, 137816, 407312, 470210, 324538)  
labels <- c("SNP", "MNP", "INS", "DEL", "COMPLEX")
pct <- round(slices/sum(slices)*100, digits=2)
labels <- paste(labels, pct)
labels <- paste(labels, "%", sep="")
pie(slices, labels=labels, col=rainbow(length(labels)), main="Different types of variants")

# histogram plot of the same data
variants <- round(slices/sum(slices), digits=2)
names(variants) <- c("SNP", "MNP", "INS", "DEL", "COMPLEX")
barplot(variants, col="red", main="Different types of variants", xlab="Type of variant", ylab="Percentage")


# circos plot
#source("http://bioconductor.org/biocLite.R")
#biocLite("RCircos")
#biocLite("IdeoViz")
library(IdeoViz)
library(RCircos)

## download cytoband ideogram data from UCSC using package IdeoViz
ideo <- getIdeo("equCab2")


## set up RCircos core components
chr.exclude <- NULL
cyto.info <- ideo
### may plot both to the inside and outside of the ideogram track
tracks.inside <- 10
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

## plot ideogram
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

gene.label.data <- read.table(file="gene.label", sep="\t", header=F)
colnames(gene.label.data) <- c("Chromosome", "chromStart", "chromEnd", "Gene")

## subset gene list, first remove genes that do not have gene symbols, because ENSEMBL IDs are too long in character.
tmp.list <- gene.label.data[-grep("ENSECAT", gene.label.data$Gene),]
idx <- sample(1:dim(tmp.list)[1], 30, replace=F)
gene.list <- tmp.list[idx,]
name.col <- 4
side <- "in"
track.num <- 1
RCircos.Gene.Connector.Plot(gene.list, track.num, side)

track.num <- 2
RCircos.Gene.Name.Plot(gene.list, name.col, track.num, side)

## check plot parameters
RCircos.Get.Gene.Name.Plot.Parameters()


## add extra tracks -- Histogram plot of CNV data
gene.list$CNV <- floor(runif(30, 1, 7))
data.col <- 5
track.num <- 5
side <- "in"
#RCircos.Histogram.Plot(gene.list, data.col, track.num, side, is.sorted=FALSE, min.value=-2)
my.RCircos.Histogram.Plot(gene.list, data.col, track.num, side, is.sorted=FALSE, max.value=max(gene.list$CNV), min.value=-4, colors=TRUE)


## add extra tracks -- Scatter plot of RNASeq results
gene.expr <- gene.list
colnames(gene.expr)
colnames(gene.expr) <- c("chromosome", "start", "stop", "gene.name", "CNV")
### generate random logFC data of 30 in length between -3 and 3
gene.expr$logFC <- runif(30, -3, 3)
### first three columns are necessary and one column that corresponds to the data to plot
head(gene.expr)

data.col <- 6
track.num <- 6
side <- "in"
### by.fold is a zero or positive number. If it's positive, then any data point with a value >= by.fold will be plotted as red color; any data point with a value <= -by.fold will be plotted as blue color; otherwise, data point will be plotted in black color.
by.fold <- 1.5
### plot scatter plot
RCircos.Scatter.Plot(gene.expr, data.col, track.num, side, by.fold, is.sorted=FALSE)

## add extra tracks -- Line plot of Coverage results
genome.cov <- gene.label.data[-grep("ENSECAT", gene.label.data$Gene),]
colnames(genome.cov) <- c("chromosome", "start", "stop", "gene.name")
genome.cov$logCOV <- rnorm(dim(genome.cov)[1], 0, 0.9)
data.col <- 5
track.num <- 7
side <- "in"
RCircos.Line.Plot(genome.cov, data.col, track.num, side, is.sorted=FALSE)

## add most inside track -- link lines and ribbons
### generate random translocation variants data
link.data <- data.frame(Chromosome=character(), chromStart=integer(), chromEnd=integer(), Chromosome.1=character(), chromStart.1=integer(), chromEnd.1=integer(), stringsAsFactors=F)
for (i in 1:15) {
  n.rand <- floor(runif(1, 1,34))
  chrom <- ideo$chrom[n.rand]
  str <- floor(runif(1, ideo$chromStart[n.rand], ideo$chromEnd[n.rand]))
  ed <- floor(runif(1, ideo$chromStart[n.rand], ideo$chromEnd[n.rand]))
  n.rand <- floor(runif(1, 1, 34))
  chrom.1 <- ideo$chrom[n.rand]
  str.1 <- floor(runif(1, ideo$chromStart[n.rand], ideo$chromEnd[n.rand]))
  ed.1 <- str.1
  if (ed < str) {
    tmp <- ed
    ed <- str
    str <- tmp
  }
  link.data <- rbind(link.data, data.frame(Chromsome=chrom, chromStart=str, chromEnd=ed, Chromosome.1=chrom.1, chromStart.1=str.1, chromEnd.1=ed.1))
  i <- i + 1
}

track.num <- 9
### plot link lines
RCircos.Link.Plot(link.data, track.num, TRUE, is.sorted=FALSE)

### ribbon data
ribbon.data <- link.data
colnames(ribbon.data) <- c("chromA", "chromStartA", "chromEndA", "chromB", "chromStartB", "chromEndB")
RCircos.Ribbon.Plot(ribbon.data=ribbon.data, track.num=9, by.chromosome=FALSE, twist=FALSE, is.sorted=FALSE)

## GWAS data plotting
#biocLite("devtools")
library(devtools)
install_github("stephenturner/qqman")
library(qqman)

# using the example data from qqman package
head(gwasResults)
manhattan(gwasResults)
# use colors for chromosomes
manhattan(gwasResults, col=c("red", "blue", "green"))
# change default horizontal line position
manhattan(gwasResults, suggestiveline=-log10(1e-6), genomewideline=-log10(1e-8), col=c("red", "blue", "green"))

