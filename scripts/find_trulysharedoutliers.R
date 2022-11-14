library(data.table)
library(ggplot2)
library(patchwork)

print("xxxxxx")
args = commandArgs(trailingOnly=TRUE)
print(args)

fst2019<-fread(args[1], header=TRUE)
fst2020<-fread(args[2], header=TRUE)
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
differences<-setdiff(fst2019$chr_pos, fst2020$chr_pos)
differences2020<-setdiff(fst2020$chr_pos, fst2019$chr_pos)

fst2019cleaned <- fst2019[!fst2019$chr_pos %in% differences, ]
fst2020cleaned <- fst2020[!fst2020$chr_pos %in% differences2020, ]

cutoff2019<-quantile(fst2019cleaned$WEIR_AND_COCKERHAM_FST, 0.99, na.rm=TRUE)
cutoff2020<-quantile(fst2020cleaned$WEIR_AND_COCKERHAM_FST, 0.99, na.rm=TRUE)
# density2019<-ggplot(data=fst2019cleaned, aes(x= WEIR_AND_COCKERHAM_FST, ..density..)) +
#   geom_density(alpha=0.5, fill="blue")+
#   ylab("Probability Density Function")+ xlab("FST")+
#   geom_vline(data=fst2019cleaned, aes(xintercept=cutoff2019, color="blue"),size=1)+
#   geom_vline(data=fst2019cleaned, aes(xintercept=mean(na.omit(WEIR_AND_COCKERHAM_FST)), color="black"),size=1)+
#   
#   scale_color_manual(values = c("blue","black"))+
#   xlim(-0.332,1)+
#   geom_text(aes(x=(mean(na.omit(WEIR_AND_COCKERHAM_FST))+0.05), y=8.0), label="mean FST")+
#   geom_text(aes(x=(cutoff2019+0.05),y=8), label = "99th percentile")+
#   theme_bw() +
#   theme(axis.text=element_text(size=15), 
#         axis.text.x=element_text(size=15, angle = 0),
#         axis.title=element_text(size=10))+
#   theme(legend.text=element_blank(), legend.title = element_blank())

# ggsave(filename = "density2019.png",ggplot(data=fst2019cleaned, aes(x= WEIR_AND_COCKERHAM_FST, ..density..)) +
#          geom_density(alpha=0.5, fill="blue")+
#          ylab("Probability Density Function")+ xlab("FST")+
#          geom_vline(aes(xintercept=cutoff2019, color="blue"),size=0.5)+
#          geom_vline(aes(xintercept=mean(na.omit(WEIR_AND_COCKERHAM_FST)), color="blue"),size=0.5)+
#          
#          scale_color_manual(values = c("blue","blue"))+
#          xlim(-0.332,1)+
#          geom_text(aes(x=(mean(na.omit(WEIR_AND_COCKERHAM_FST))+0.15), y=12.0), label="mean FST")+
#          geom_text(aes(x=(cutoff2019+0.15),y=12), label = "99th percentile")+
#                      theme_bw(),
#        width = 8, height = 4, dpi = 300, units = "in", device='png')
# ggsave(filename = "density2020.png",ggplot(data=fst2020cleaned, aes(x= WEIR_AND_COCKERHAM_FST, ..density..)) +
#          geom_density(alpha=0.5, fill="red")+
#          ylab("Probability Density Function")+ xlab("FST")+
#          geom_vline(aes(xintercept=cutoff2020, color="red"),size=0.5)+
#          geom_vline(aes(xintercept=mean(na.omit(WEIR_AND_COCKERHAM_FST)), color="red"),size=0.5)+
#          
#          scale_color_manual(values = c("red","red"))+
#          xlim(-0.332,1)+
#          #geom_text(aes(x=(mean(na.omit(WEIR_AND_COCKERHAM_FST))+0.15), y=12.0), label="mean FST")+
#          #geom_text(aes(x=(cutoff2020+0.15),y=12), label = "99th percentile")+
#          theme_bw(),
#        width = 8, height = 4, dpi = 300, units = "in", device='png')
density2019<-ggplot(data=fst2019cleaned, aes(x= WEIR_AND_COCKERHAM_FST, ..density..)) +
  geom_density(alpha=0.5, fill="blue")+
  ylab("Probability Density Function")+ xlab("FST")+
  geom_vline(aes(xintercept=cutoff2019, color="blue"),size=0.5)+
  geom_vline(aes(xintercept=mean(na.omit(WEIR_AND_COCKERHAM_FST)), color="blue"),size=0.5)+
  
  scale_color_manual(values = c("blue","blue"))+
  xlim(-0.332,1)+
  ylim(0,30)+
  theme_bw()+
  geom_text(aes(x=(mean(na.omit(WEIR_AND_COCKERHAM_FST))+0.06), y=20.0), label="mean = 0.02")+
  geom_text(aes(x=(cutoff2019+0.06),y=20), label = "99th percentile \n = 0.47")+
  
  
  theme(legend.position="none")+labs(tag = "A")
density2020<-ggplot(data=fst2020cleaned, aes(x= WEIR_AND_COCKERHAM_FST, ..density..)) +
  geom_density(alpha=0.5, fill="red")+
  ylab("Probability Density Function")+ xlab("FST")+
  geom_vline(aes(xintercept=cutoff2020, color="red"),size=0.5)+
  geom_vline(aes(xintercept=mean(na.omit(WEIR_AND_COCKERHAM_FST)), color="red"),size=0.5)+
  
  scale_color_manual(values = c("red","red"))+
  xlim(-0.332,1)+
  ylim(0,30)+
  geom_text(aes(x=(mean(na.omit(WEIR_AND_COCKERHAM_FST))+0.06), y=20.0), label="mean = 0.01")+
  geom_text(aes(x=(cutoff2020+0.061),y=20), label = "99th percentile \n = 0.20")+
  theme_bw()+
  theme(legend.position="none")+labs(tag = "B")
ggsave(filename = "densities20221109.png", density2019 / density2020,
       width = 13.33, height = 7.5, dpi = 300, units = "in", device='png')

testfunction<-function(test,cutoff2019) {
  if (test=="NaN") {
    g = "background"
  } else if (test>= cutoff2019) {
    g = "outlier"
  } else {
    g = "background"
  }
  return(g)
}

fst2019cleaned$outliers<-mapply(testfunction, fst2019cleaned$WEIR_AND_COCKERHAM_FST, cutoff2019)
fst2020cleaned$outliers<-mapply(testfunction, fst2020cleaned$WEIR_AND_COCKERHAM_FST, cutoff2020)
fst2019cleaned$cohort<-rep("2019",nrow(fst2019cleaned))
fst2020cleaned$cohort<-rep("2020",nrow(fst2020cleaned))
fst2020cleanedoutliers<-subset(fst2020cleaned, outliers=="outlier")
fst2019cleanedoutliers<-subset(fst2019cleaned, outliers=="outlier")
fst2019cleanedoutliers$CHROM_POS<-paste(fst2019cleanedoutliers$CHROM,fst2019cleanedoutliers$POS,sep=".")
fst2020cleanedoutliers$CHROM_POS<-paste(fst2020cleanedoutliers$CHROM,fst2020cleanedoutliers$POS,sep=".")

sharedoutliers<-rbind(na.omit(fst2019cleanedoutliers[match(fst2020cleanedoutliers$POS, fst2019cleanedoutliers$POS),]),
                      na.omit(fst2020cleanedoutliers[match(fst2019cleanedoutliers$POS, fst2020cleanedoutliers$POS),]))
sharedoutliers$CHROM_POS<-paste(sharedoutliers$CHROM,sharedoutliers$POS,sep=".")
#sharedoutliers$cohort<-c(rep("2019outlier",nrow(na.omit(fst2019cleanedoutliers[match(fst2020cleanedoutliers$POS, fst2019cleanedoutliers$POS),]))),
#                         rep("2020outlier", nrow(na.omit(fst2020cleanedoutliers[match(fst2019cleanedoutliers$POS, fst2020cleanedoutliers$POS),]))))
##
outlierfunction<-function(outlierlist2019,outlierlist2020) {
  if (outlierlist2019== "outlier" & outlierlist2020=="outlier") {
    g = "sharedoutlier"
  } else if  (outlierlist2019== "background" & outlierlist2020=="outlier") {
    g = "2020 outlier"
  } else if (outlierlist2019== "outlier" & outlierlist2020=="background") {
    g = "2019 outlier"
  } else {
    g = "background"
  }
  return(g)
}
tst<-mapply(outlierfunction, fst2019cleaned$outliers, fst2020cleaned$outliers)

fst_2019<-fst2019cleaned$WEIR_AND_COCKERHAM_FST
fst_2020<-fst2020cleaned$WEIR_AND_COCKERHAM_FST
outliers_2019<-fst2019cleaned$outliers
outliers_2020<-fst2020cleaned$outliers
cohort<-tst
freq_F0_2019<-fread(args[3],header=TRUE) #need to manually change #CHROM to CHROM in the afreq files before reading them in
freq_F1_2019<-fread(args[4],header=TRUE)
freq_F0_2020<-fread(args[5],header=TRUE)
freq_F1_2020<-fread(args[6],header=TRUE)
differencefunction<-function(dataframe) {
  addition_vector_2019<-c(rep(0, nrow(subset(dataframe,CHROM=="1"))),
                          rep(Chr1L, nrow(subset(dataframe,CHROM=="2"))),
                          rep(Chr1L+Chr2L, nrow(subset(dataframe,CHROM=="3"))),
                          rep(Chr1L+Chr2L + Chr3L , nrow(subset(dataframe,CHROM=="4"))),
                          rep(Chr1L+Chr2L + Chr3L + Chr4L, nrow(subset(dataframe,CHROM=="5"))),
                          rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L, nrow(subset(dataframe,CHROM=="6"))),
                          rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L, nrow(subset(dataframe,CHROM=="7"))),
                          rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L, nrow(subset(dataframe,CHROM=="8"))),
                          rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L, nrow(subset(dataframe,CHROM=="9"))),
                          rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L, nrow(subset(dataframe,CHROM=="10"))),
                          rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L + Chr10L, nrow(subset(dataframe,CHROM=="11"))),
                          rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L + Chr10L + Chr11L, nrow(subset(dataframe,CHROM=="12"))),
                          rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L + Chr10L + Chr11L + Chr12L, nrow(subset(dataframe,CHROM=="13"))),
                          rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L + Chr10L + Chr11L + Chr12L + Chr13L, nrow(subset(dataframe,CHROM=="14"))))
  
  chr_pos_2019<-dataframe$POS+addition_vector_2019
  dataframe$chr_pos<-chr_pos_2019
  differences<-setdiff(dataframe$chr_pos, fst2020cleaned$chr_pos)
  g<-dataframe[!dataframe$chr_pos %in% differences, ]
  return(g)
}
freq_F0_2019<-differencefunction(freq_F0_2019)
freq_F1_2019<-differencefunction(freq_F1_2019)
freq_F0_2020<-differencefunction(freq_F0_2020)
freq_F1_2020<-differencefunction(freq_F1_2020)

F0_2019_afreq<-freq_F0_2019$ALT_FREQS
F1_2019_afreq<-freq_F1_2019$ALT_FREQS
F0_2020_afreq<-freq_F0_2020$ALT_FREQS
F1_2020_afreq<-freq_F1_2020$ALT_FREQS

wholeframe<-data.frame("CHROM"=fst2019cleaned$CHROM, "POS"=fst2019cleaned$POS, "chr_pos"=fst2019cleaned$chr_pos, fst_2019, fst_2020,
                       outliers_2019, outliers_2020, cohort, F0_2019_afreq, F1_2019_afreq,
                       F0_2020_afreq, F1_2020_afreq)
wholeframe_sharedoutliers<-subset(wholeframe,cohort=="sharedoutlier")
wholeframe_sharedoutliers<-rbind(subset(wholeframe_sharedoutliers, F1_2019_afreq<F0_2019_afreq & F1_2020_afreq< F0_2020_afreq),
                                 subset(wholeframe_sharedoutliers, F1_2019_afreq>F0_2019_afreq & F1_2020_afreq> F0_2020_afreq)
                                 
)
write.table(wholeframe_sharedoutliers, file="trulysharedoutliers20210928.txt",quote=FALSE,row.names=FALSE,sep = "\t")
# wholeframe$diff2019<-wholeframe$F1_2019_afreq-wholeframe$F0_2019_afreq
# wholeframe$diff2020<-wholeframe$F1_2020_afreq-wholeframe$F0_2020_afreq
# wholeframe$diffoutliers2019_high<-mapply(testfunction, wholeframe$diff2019, cutoff2019_afreq_high)
# wholeframe$diffoutliers2019_low<-mapply(testfunction, wholeframe$diff2019, cutoff2019_afreq_low)
# 
# wholeframe$diffoutliers2020_high<-mapply(testfunction, wholeframe$diff2020, cutoff2020_afreq_high)
# wholeframe$diffoutliers2020_low<-mapply(testfunction, wholeframe$diff2020, cutoff2020_afreq_low)
# 
# wholeframe_sharedoutliers$diff2019<-wholeframe_sharedoutliers$F1_2019_afreq-wholeframe_sharedoutliers$F0_2019_afreq
# wholeframe_sharedoutliers$diff2020<-wholeframe_sharedoutliers$F1_2020_afreq-wholeframe_sharedoutliers$F0_2020_afreq
# wholeframe$diff2019abs<-abs(wholeframe$diff2019)
# wholeframe$diff2020abs<-abs(wholeframe$diff2020)
# wholeframe_sharedoutliers$diff2019abs<-abs(wholeframe_sharedoutliers$diff2019)
