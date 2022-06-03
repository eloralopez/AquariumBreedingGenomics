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
chr_pos=NULL
cond<-fst2019$CHROM
chr_pos_func<-function(fst2019) {
  for (i in 1:nrow(fst2019)) {
    if (fst2019$CHROM[i] == "chr2") {
      chr_pos[i] = fst2019$POS[i] + Chr1L
    } else if (fst2019$CHROM[i] == "chr3"){
      chr_pos[i] = fst2019$POS[i]+Chr1L+Chr2L
    } else if (fst2019$CHROM[i] == "chr4"){
      chr_pos[i] = fst2019$POS[i]+Chr1L+Chr2L +Chr3L
    } else if (fst2019$CHROM[i] == "chr5"){
      chr_pos[i] = fst2019$POS[i]+Chr1L+Chr2L+ Chr3L + Chr4L 
    } else if (fst2019$CHROM[i] == "chr6"){
      chr_pos[i] = fst2019$POS[i]+Chr1L+Chr2L+ Chr3L + Chr4L + Chr5L
    } else if (fst2019$CHROM[i] == "chr7"){
      chr_pos[i] = fst2019$POS[i]+Chr1L+Chr2L+ Chr3L + Chr4L + Chr5L + Chr6L
    } else if (fst2019$CHROM[i] == "chr8"){
      chr_pos[i] = fst2019$POS[i]+Chr1L+Chr2L+ Chr3L + Chr4L + Chr5L + Chr6L + Chr7L 
    } else if (fst2019$CHROM[i] == "chr9"){
      chr_pos[i] = fst2019$POS[i]+Chr1L+Chr2L+ Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L 
    } else if (fst2019$CHROM[i] == "chr10"){
      chr_pos[i] = fst2019$POS[i]+Chr1L+Chr2L+ Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L 
    } else if (fst2019$CHROM[i] == "chr11"){
      chr_pos[i] = fst2019$POS[i]+Chr1L+Chr2L+ Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L + Chr10L
    } else if (fst2019$CHROM[i] == "chr12"){
      chr_pos[i] = fst2019$POS[i]+Chr1L+Chr2L+ Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L + Chr10L + Chr11L 
    } else if (fst2019$CHROM[i] == "chr13"){
      chr_pos[i] = fst2019$POS[i]+Chr1L+Chr2L+ Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L + Chr10L + Chr11L + Chr12L
    } else if (fst2019$CHROM[i] == "chr14"){
      chr_pos[i] = fst2019$POS[i]+Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L + Chr10L + Chr11L + Chr12L + Chr13L
    } else {
      chr_pos[i] = fst2019$POS[i]
    }
  }  
  return(chr_pos)
}  
system.time(chr_pos_func(fst2019))
chr_pos=NULL
pos
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

cutoff2019<-quantile(fst2019$WEIR_AND_COCKERHAM_FST, 0.99)
cutoff2020<-quantile(fst2020$WEIR_AND_COCKERHAM_FST, 0.99)
testfunction<-function(test,cutoff2019) {
  #cutoff2019<-0.431722
  #cutoff2019<-0.208314 
  if (test>= cutoff2019) {
    g = "outlier"
  } else {
    g = "background"
  }
  return(g)
}
#fst2019$outliers<-sapply(fst2019$WEIR_AND_COCKERHAM_FST, testfunction)
#fst2020$outliers<-sapply(fst2020$WEIR_AND_COCKERHAM_FST, testfunction)

fst2019$outliers<-mapply(testfunction, fst2019$WEIR_AND_COCKERHAM_FST, cutoff2019)
fst2020$outliers<-mapply(testfunction, fst2020$WEIR_AND_COCKERHAM_FST, cutoff2020)
fst2020outliers<-subset(fst2020, outliers=="outlier")
fst2019outliers<-subset(fst2019, outliers=="outlier")
sharedoutliers<-rbind(na.omit(fst2019outliers[match(fst2020outliers$POS, fst2019outliers$POS),]),
                      na.omit(fst2020outliers[match(fst2019outliers$POS, fst2020outliers$POS),]))
print(sharedoutliers)
write.table(sharedoutliers,"sharedoutlier.txt",sep="\t",row.names=FALSE)
ggplot(fst2020, aes(x = POS, y = WEIR_AND_COCKERHAM_FST)) +
  #geom_point(size = 3,color="blue",alpha=0.4)+
  #geom_point(data=fst2019,size=3,alpha=0.4,aes(POS, WEIR_AND_COCKERHAM_FST, color="red"))+
  geom_point(data=sharedoutliers,size=3,alpha=0.4,aes(POS, WEIR_AND_COCKERHAM_FST, color="purple"))+
  theme_bw()
ggsave(filename="fstplot_sharedoutliers20210917_personalcomp.png", ggplot(fst2020, aes(x = POS, y = WEIR_AND_COCKERHAM_FST)) +
         geom_point(size = 3,color="blue",alpha=0.4)+
         geom_point(data=fst2019,size=3,alpha=0.4,aes(POS, WEIR_AND_COCKERHAM_FST, color="red"))+
         geom_point(data=sharedoutliers,size=3,alpha=0.4,aes(POS, WEIR_AND_COCKERHAM_FST, color="purple"))+
         theme_bw(), width=8, height=4,dpi=300,units="in",device='png')
ggsave("plot.pdf", device=cairo_pdf)
ggsave(filename = "fstplot.png", ggplot(fst2020, aes(x = chr_pos_2020, y = WEIR_AND_COCKERHAM_FST)) +
         geom_point(size = 1,color="blue",alpha=0.4)+
         geom_point(data=fst2019,size=1,alpha=0.4,aes(chr_pos_2019, WEIR_AND_COCKERHAM_FST, color="red"))+
         theme_bw(),
       width = 8, height = 4, dpi = 300, units = "in", device='png')
year<-c(rep("2019",nrow(fst2019)),rep("2020",nrow(fst2020)))
fstcomb<-data.frame(year, "POS"=c(fst2019$POS,fst2020$POS), "FST"=c(fst2019$WEIR_AND_COCKERHAM_FST, fst2020$WEIR_AND_COCKERHAM_FST))
ggplot(fstcomb, aes(x=FST, color=year))+
  geom_density(alpha=0.5)+theme_bw()
theme_bw()
#plot(difference)
comparefst<-data.frame("fst2019"=fst2019$WEIR_AND_COCKERHAM_FST,"fst2020"=fst2020$WEIR_AND_COCKERHAM_FST, "POS"=fst2020$POS)
compareplot<-ggplot(comparefst, aes(x = fst2019, y = fst2020)) +
  geom_point(size = 3,color="black",alpha=0.4)+theme_bw()
subset(comparefst, fst2019>0.4 & fst2020>0.28)

differenceplot | fstplot | compareplot

