library(stringr)
library(ggplot2)
setwd("~/CaptiveBreeding/Rscripts")
print("xxxxxx")
args = commandArgs(trailingOnly=TRUE)
print(args)
Chr1L <-34295999
Chr2L <-23919585
Chr3L <-26478413
Chr4L <-26766062
Chr5L <-33503191
Chr6L <-32490106
Chr7L <-30731630
Chr8L <-30684215
Chr9L <-28654017
Chr10L <-26527961
Chr11L <-24262129
Chr12L <-24153179
Chr13L <-24061126
Chr14L <-23679165


fst2019<-read.delim(args[1], header=TRUE)
fst2020<-read.delim(args[2], header=TRUE)
fst2019<-na.omit(fst2019)
fst2020<-na.omit(fst2020)
addition_vector_2019<-c(rep(0, nrow(subset(fst2019,CHROM=="chr1"))),
                   rep(Chr1L, nrow(subset(fst2019,CHROM=="chr2"))),
                   rep(Chr1L+Chr2L, nrow(subset(fst2019,CHROM=="chr3"))),
                   rep(Chr1L+Chr2L + Chr3L , nrow(subset(fst2019,CHROM=="chr4"))),
                   rep(Chr1L+Chr2L + Chr3L + Chr4L, nrow(subset(fst2019,CHROM=="chr5"))),
                   rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L, nrow(subset(fst2019,CHROM=="chr6"))),
                   rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L, nrow(subset(fst2019,CHROM=="chr7"))),
                   rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L, nrow(subset(fst2019,CHROM=="chr8"))),
                   rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L, nrow(subset(fst2019,CHROM=="chr9"))),
                   rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L, nrow(subset(fst2019,CHROM=="chr10"))),
                   rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L + Chr10L, nrow(subset(fst2019,CHROM=="chr11"))),
                   rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L + Chr10L + Chr11L, nrow(subset(fst2019,CHROM=="chr12"))),
                   rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L + Chr10L + Chr11L + Chr12L, nrow(subset(fst2019,CHROM=="chr13"))),
                   rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L + Chr10L + Chr11L + Chr12L + Chr13L, nrow(subset(fst2019,CHROM=="chr14"))))
print("additioning 2019 done")
addition_vector_2020<-c(rep(0, nrow(subset(fst2020,CHROM=="chr1"))),
                   rep(Chr1L, nrow(subset(fst2020,CHROM=="chr2"))),
                   rep(Chr1L+Chr2L, nrow(subset(fst2020,CHROM=="chr3"))),
                   rep(Chr1L+Chr2L + Chr3L , nrow(subset(fst2020,CHROM=="chr4"))),
                   rep(Chr1L+Chr2L + Chr3L + Chr4L, nrow(subset(fst2020,CHROM=="chr5"))),
                   rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L, nrow(subset(fst2020,CHROM=="chr6"))),
                   rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L, nrow(subset(fst2020,CHROM=="chr7"))),
                   rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L, nrow(subset(fst2020,CHROM=="chr8"))),
                   rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L, nrow(subset(fst2020,CHROM=="chr9"))),
                   rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L, nrow(subset(fst2020,CHROM=="chr10"))),
                   rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L + Chr10L, nrow(subset(fst2020,CHROM=="chr11"))),
                   rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L + Chr10L + Chr11L, nrow(subset(fst2020,CHROM=="chr12"))),
                   rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L + Chr10L + Chr11L + Chr12L, nrow(subset(fst2020,CHROM=="chr13"))),
                   rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L + Chr10L + Chr11L + Chr12L + Chr13L, nrow(subset(fst2020,CHROM=="chr14"))))
print("additioning 2020 done")
chr_pos_2019<-fst2019$POS+addition_vector_2019
chr_pos_2020<-fst2020$POS+addition_vector_2020

fst2019$chr_pos<-chr_pos_2019
fst2020$chr_pos<-chr_pos_2020
ggsave(filename = "fstplot2020_20210917_rosalindf.png", ggplot(fst2020, aes(x = chr_pos_2020, y = WEIR_AND_COCKERHAM_FST)) +
        geom_point(size = 1,color="blue",alpha=0.4)+
        geom_point(data=fst2019,size=1,alpha=0.4,aes(chr_pos_2019, WEIR_AND_COCKERHAM_FST, color="red"))+
         theme_bw(),
       width = 8, height = 4, dpi = 300, units = "in", device='png')
print("saved 2020 plot")


