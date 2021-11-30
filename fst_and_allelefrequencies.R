setwd("~/Documents/CaptiveBreeding/plots")
library(stringr)
library(ggplot2)
library(kinship2)
library(patchwork)
library(ggforce)
library(rehh)
library(data.table)
library(ggpubr)
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
fst2019<-read.delim("F1_vs_F0_2019.weir.fst", header=TRUE)
fst2020<-read.delim("F1_vs_F0_2020.weir.fst", header=TRUE)
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
fst2019cleaned <- fst2019[!fst2019$chr_pos %in% differences, ]
fst2020cleaned <- fst2020[!fst2020$chr_pos %in% differences, ]

differencefunction<-function(dataframe) {
  addition_vector_2019<-c(rep(0, nrow(subset(dataframe,X.CHROM=="1"))),
                          rep(Chr1L, nrow(subset(dataframe,X.CHROM=="2"))),
                          rep(Chr1L+Chr2L, nrow(subset(dataframe,X.CHROM=="3"))),
                          rep(Chr1L+Chr2L + Chr3L , nrow(subset(dataframe,X.CHROM=="4"))),
                          rep(Chr1L+Chr2L + Chr3L + Chr4L, nrow(subset(dataframe,X.CHROM=="5"))),
                          rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L, nrow(subset(dataframe,X.CHROM=="6"))),
                          rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L, nrow(subset(dataframe,X.CHROM=="7"))),
                          rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L, nrow(subset(dataframe,X.CHROM=="8"))),
                          rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L, nrow(subset(dataframe,X.CHROM=="9"))),
                          rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L, nrow(subset(dataframe,X.CHROM=="10"))),
                          rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L + Chr10L, nrow(subset(dataframe,X.CHROM=="11"))),
                          rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L + Chr10L + Chr11L, nrow(subset(dataframe,X.CHROM=="12"))),
                          rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L + Chr10L + Chr11L + Chr12L, nrow(subset(dataframe,X.CHROM=="13"))),
                          rep(Chr1L+Chr2L + Chr3L + Chr4L + Chr5L + Chr6L + Chr7L + Chr8L + Chr9L + Chr10L + Chr11L + Chr12L + Chr13L, nrow(subset(dataframe,X.CHROM=="14"))))
  
  chr_pos_2019<-dataframe$POS+addition_vector_2019
  dataframe$chr_pos<-chr_pos_2019
  differences<-setdiff(dataframe$chr_pos, fst2020cleaned$chr_pos)
  g<-dataframe[!dataframe$chr_pos %in% differences, ]
  return(g)
}
fst2019<-subset(fst2019, chr_pos != differences[1] & chr_pos != differences[2] &
                  chr_pos != differences[3] & chr_pos != differences[4] &
                  chr_pos != differences[5] & chr_pos != differences[6] &
                  chr_pos != differences[7] & chr_pos != differences[8] &
                  chr_pos != differences[9] & chr_pos != differences[10] &
                  chr_pos != differences[11] & chr_pos != differences[12] &
                  chr_pos != differences[13] & chr_pos != differences[14] &
                  chr_pos != differences[15] & chr_pos != differences[16] &
                  chr_pos != differences[17] & chr_pos != differences[18])
cutoff2019<-quantile(fst2019$WEIR_AND_COCKERHAM_FST, 0.99, na.rm=TRUE)
cutoff2020<-quantile(fst2020$WEIR_AND_COCKERHAM_FST, 0.99, na.rm=TRUE)

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

fst2019$outliers<-mapply(testfunction, fst2019$WEIR_AND_COCKERHAM_FST, cutoff2019)
fst2020$outliers<-mapply(testfunction, fst2020$WEIR_AND_COCKERHAM_FST, cutoff2020)
fst2019$cohort<-rep("2019",nrow(fst2019))
fst2020$cohort<-rep("2020",nrow(fst2020))
fst2020outliers<-subset(fst2020, outliers=="outlier")
fst2019outliers<-subset(fst2019, outliers=="outlier")
fst2019outliers$CHROM_POS<-paste(fst2019outliers$CHROM,fst2019outliers$POS,sep=".")
fst2020outliers$CHROM_POS<-paste(fst2020outliers$CHROM,fst2020outliers$POS,sep=".")

sharedoutliers<-rbind(na.omit(fst2019outliers[match(fst2020outliers$POS, fst2019outliers$POS),]),
                      na.omit(fst2020outliers[match(fst2019outliers$POS, fst2020outliers$POS),]))
sharedoutliers$CHROM_POS<-paste(sharedoutliers$CHROM,sharedoutliers$POS,sep=".")
sharedoutliers$cohort<-c(rep("2019outlier",nrow(na.omit(fst2019outliers[match(fst2020outliers$POS, fst2019outliers$POS),]))),
                         rep("2020outlier", nrow(na.omit(fst2020outliers[match(fst2019outliers$POS, fst2020outliers$POS),]))))
ggsave(filename = "fst_withtrulysharedoutliers20210929.png", ggplot(fst2020, aes(x = chr_pos, y = WEIR_AND_COCKERHAM_FST,color=cohort)) +
  scale_color_manual(values=c("skyblue","blue3","tomato","red3"))+
  geom_point(size = 3,alpha=0.1)+
  geom_point(data=fst2019,size=3,alpha=0.1,aes(chr_pos, WEIR_AND_COCKERHAM_FST, color=cohort))+
  geom_point(data=wholeframe_sharedoutliers,size=3,alpha=0.7,aes(chr_pos, fst_2019))+
  geom_point(data=wholeframe_sharedoutliers,size=3,alpha=0.7,aes(chr_pos, fst_2020))+
    
    #geom_hline(yintercept=mean(na.omit(fst2019$WEIR_AND_COCKERHAM_FST)))+
  #geom_hline(yintercept=mean(na.omit(fst2020$WEIR_AND_COCKERHAM_FST)))+
  theme_bw(),
  width = 8, height = 4, dpi = 300, units = "in", device='png')

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
tst<-mapply(outlierfunction, fst2019$outliers, fst2020$outliers)

fst_2019<-fst2019$WEIR_AND_COCKERHAM_FST
fst_2020<-fst2020$WEIR_AND_COCKERHAM_FST
outliers_2019<-fst2019$outliers
outliers_2020<-fst2020$outliers
cohort<-tst
freq_F0_2019<-read.delim("F0_2019.afreq",header=TRUE)
freq_F1_2019<-read.delim("F1_2019.afreq",header=TRUE)
freq_F0_2020<-read.delim("F0_2020.afreq",header=TRUE)
freq_F1_2020<-read.delim("F1_2020.afreq",header=TRUE)
freq_F0_2019<-differencefunction(freq_F0_2019)
freq_F1_2019<-differencefunction(freq_F1_2019)
freq_F0_2020<-differencefunction(freq_F0_2020)
freq_F1_2020<-differencefunction(freq_F1_2020)

F0_2019_afreq<-freq_F0_2019$ALT_FREQS
F1_2019_afreq<-freq_F1_2019$ALT_FREQS
F0_2020_afreq<-freq_F0_2020$ALT_FREQS
F1_2020_afreq<-freq_F1_2020$ALT_FREQS

wholeframe<-data.frame("CHROM"=fst2019$CHROM, "POS"=fst2019$POS, "chr_pos"=fst2019$chr_pos, fst_2019, fst_2020,
                       outliers_2019, outliers_2020, cohort, F0_2019_afreq, F1_2019_afreq,
                       F0_2020_afreq, F1_2020_afreq)
wholeframe_sharedoutliers<-subset(wholeframe,cohort=="sharedoutlier")
wholeframe_sharedoutliers<-rbind(subset(wholeframe_sharedoutliers, F1_2019_afreq<F0_2019_afreq & F1_2020_afreq< F0_2020_afreq),
                                 subset(wholeframe_sharedoutliers, F1_2019_afreq>F0_2019_afreq & F1_2020_afreq> F0_2020_afreq)
                                 
)
#write.table(wholeframe_sharedoutliers, file="trulysharedoutliers20210928.txt",quote=FALSE,
            row.names=FALSE,sep = "\t")
wholeframe_sharedoutliers_2<-data.frame("CHROM"=rep(wholeframe_sharedoutliers$CHROM,4), "POS"=rep(wholeframe_sharedoutliers$POS,4),
                                        "chr_pos"=rep(wholeframe_sharedoutliers$chr_pos,4), "afreqs"=c(wholeframe_sharedoutliers$F0_2019_afreq, wholeframe_sharedoutliers$F1_2019_afreq,wholeframe_sharedoutliers$F0_2020_afreq, wholeframe_sharedoutliers$F1_2020_afreq),
                                        "cohort"=c(rep("F0_2019",length(wholeframe_sharedoutliers$F0_2019_afreq)), rep("F1_2019",length(wholeframe_sharedoutliers$F1_2019_afreq)), rep("F0_2020",length(wholeframe_sharedoutliers$F0_2020_afreq)), rep("F1_2020",length(wholeframe_sharedoutliers$F1_2020_afreq))))

#ggplot(wholeframe, aes(x=chr_pos,y=fst_2019))+
  geom_point(size=3)
ggsave(filename = "AF_22trulysharedoutliers20210929.png",ggplot(wholeframe_sharedoutliers_2, aes(x=as.factor(chr_pos),y=afreqs,fill=cohort))+
  geom_bar(position="dodge", stat="identity") +
  scale_fill_brewer(palette="Paired")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=15),
        axis.title=element_blank(),
        axis.text.x = element_text(angle = 90)),
  width = 8, height = 4, dpi = 300, units = "in", device='png')

##just the 2019 big family:
freq_2019bigfamily_F0<-read.delim("CA2019parentsandF1_GATHERED_biSNPs_F0only_AF.afreq",header=TRUE)
freq_2019bigfamily_F1<-read.delim("CA2019parentsandF1_GATHERED_biSNPs_F1only_AF.afreq",header=TRUE)
f1_genocounts<-read.delim("CA2019parentsandF1_GATHERED_biSNPs_F1only_genocount.gcount",header=TRUE)
#freq_2019bigfamily_F0<-differencefunction(freq_2019bigfamily_F0)
#freq_2019bigfamily_F1<-differencefunction(freq_2019bigfamily_F1)
F0freq<-freq_2019bigfamily_F0$ALT_FREQS
F1freq<-freq_2019bigfamily_F1$ALT_FREQS
f1_homreffreq<-f1_genocounts$HOM_REF_CT/(f1_genocounts$HOM_REF_CT+f1_genocounts$HET_REF_ALT_CTS+f1_genocounts$TWO_ALT_GENO_CTS)
f1_hetfreq<-f1_genocounts$HET_REF_ALT_CTS/(f1_genocounts$HOM_REF_CT+f1_genocounts$HET_REF_ALT_CTS+f1_genocounts$TWO_ALT_GENO_CTS)
f1_homaltfreq<-f1_genocounts$TWO_ALT_GENO_CTS/(f1_genocounts$HOM_REF_CT+f1_genocounts$HET_REF_ALT_CTS+f1_genocounts$TWO_ALT_GENO_CTS)

freqcomp<-data.frame(F0freq, F1freq, "chrom"=freq_2019bigfamily_F0$X.CHROM,"pos"=freq_2019bigfamily_F0$POS,
                     f1_homreffreq, f1_hetfreq, f1_homaltfreq)
ggsave(filename = "CA2019bigfamilyAF.png", ggplot(freqcomp, aes(x = F1freq, y = )) +
  #scale_color_manual(values=c("skyblue","blue3","tomato","red3"))+
  geom_point(size = 3,alpha=0.1)+theme_bw(),
width = 8, height = 8, dpi = 300, units = "in", device='png')

ggsave(filename = "test.png",ggplot(freqcomp, aes(x = F1freq, y = f1_hetfreq)) +
  geom_point(),
  #scale_color_manual(values=c("skyblue","blue3","tomato","red3"))+
  #geom_point(size = 3,alpha=0.1)+theme_bw(),
  width = 8, height = 8, dpi = 300, units = "in", device='png')
 
ggsave(filename = "test.png",ggplot(freqcomp, aes(x = F1freq, y = f1_hetfreq)) +
  geom_point(alpha=0.1,color="blue")+
  geom_point(data=freqcomp, alpha=0.1,aes(F1freq, f1_homreffreq,color="red"))+
  geom_point(data=freqcomp, alpha=0.1,aes(F1freq, f1_homaltfreq,color="green"))+
  theme_bw(),
  width = 8, height = 5, dpi = 300, units = "in", device='png')

#for R on alice:
library(data.table)
library(ggplot2)
library(patchwork)
setwd("~/CaptiveBreeding/analyses_dp10")
fst2019<-fread("F1_vs_F0_2019.weir.fst", header=TRUE)
fst2020<-fread("F1_vs_F0_2020.weir.fst", header=TRUE)
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
density2019<-ggplot(data=fst2019cleaned, aes(x= WEIR_AND_COCKERHAM_FST, ..density..)) +
  geom_density(alpha=0.5, fill="blue")+
  ylab("Probability Density Function")+ xlab("FST")+
  geom_vline(data=fst2019cleaned, aes(xintercept=cutoff2019, color="blue"),size=1)+
  geom_vline(data=fst2019cleaned, aes(xintercept=mean(na.omit(WEIR_AND_COCKERHAM_FST)), color="black"),size=1)+
  
  scale_color_manual(values = c("blue","black"))+
  xlim(-0.332,1)+
  geom_text(aes(x=(mean(na.omit(WEIR_AND_COCKERHAM_FST))+0.05), y=8.0), label="mean FST")+
  geom_text(aes(x=(cutoff2019+0.05),y=8), label = "99th percentile")+
  theme_bw() +
  theme(axis.text=element_text(size=15), 
        axis.text.x=element_text(size=15, angle = 0),
        axis.title=element_text(size=10))+
  theme(legend.text=element_blank(), legend.title = element_blank())
  
ggsave(filename = "density2019.png",ggplot(data=fst2019cleaned, aes(x= WEIR_AND_COCKERHAM_FST, ..density..)) +
         geom_density(alpha=0.5, fill="blue")+
         ylab("Probability Density Function")+ xlab("FST")+
         geom_vline(aes(xintercept=cutoff2019, color="blue"),size=0.5)+
         geom_vline(aes(xintercept=mean(na.omit(WEIR_AND_COCKERHAM_FST)), color="blue"),size=0.5)+
         
         scale_color_manual(values = c("blue","blue"))+
         xlim(-0.332,1)+
         geom_text(aes(x=(mean(na.omit(WEIR_AND_COCKERHAM_FST))+0.15), y=12.0), label="mean FST")+
         geom_text(aes(x=(cutoff2019+0.15),y=12), label = "99th percentile")+
                     theme_bw(),
       width = 8, height = 4, dpi = 300, units = "in", device='png')
ggsave(filename = "density2020.png",ggplot(data=fst2020cleaned, aes(x= WEIR_AND_COCKERHAM_FST, ..density..)) +
         geom_density(alpha=0.5, fill="red")+
         ylab("Probability Density Function")+ xlab("FST")+
         geom_vline(aes(xintercept=cutoff2020, color="red"),size=0.5)+
         geom_vline(aes(xintercept=mean(na.omit(WEIR_AND_COCKERHAM_FST)), color="red"),size=0.5)+
         
         scale_color_manual(values = c("red","red"))+
         xlim(-0.332,1)+
         #geom_text(aes(x=(mean(na.omit(WEIR_AND_COCKERHAM_FST))+0.15), y=12.0), label="mean FST")+
         #geom_text(aes(x=(cutoff2020+0.15),y=12), label = "99th percentile")+
         theme_bw(),
       width = 8, height = 4, dpi = 300, units = "in", device='png')
density2019<-ggplot(data=fst2019cleaned, aes(x= WEIR_AND_COCKERHAM_FST, ..density..)) +
  geom_density(alpha=0.5, fill="blue")+
  ylab("Probability Density Function")+ xlab("FST")+
  geom_vline(aes(xintercept=cutoff2019, color="blue"),size=0.5)+
  geom_vline(aes(xintercept=mean(na.omit(WEIR_AND_COCKERHAM_FST)), color="blue"),size=0.5)+
  
  scale_color_manual(values = c("blue","blue"))+
  xlim(-0.332,1)+
  #geom_text(aes(x=(mean(na.omit(WEIR_AND_COCKERHAM_FST))+0.15), y=12.0), label="mean FST")+
  #geom_text(aes(x=(cutoff2019+0.15),y=12), label = "99th percentile")+
  theme_bw()
density2020<-ggplot(data=fst2020cleaned, aes(x= WEIR_AND_COCKERHAM_FST, ..density..)) +
  geom_density(alpha=0.5, fill="red")+
  ylab("Probability Density Function")+ xlab("FST")+
  geom_vline(aes(xintercept=cutoff2020, color="red"),size=0.5)+
  geom_vline(aes(xintercept=mean(na.omit(WEIR_AND_COCKERHAM_FST)), color="red"),size=0.5)+
  
  scale_color_manual(values = c("red","red"))+
  xlim(-0.332,1)+
  #geom_text(aes(x=(mean(na.omit(WEIR_AND_COCKERHAM_FST))+0.15), y=12.0), label="mean FST")+
  #geom_text(aes(x=(cutoff2020+0.15),y=12), label = "99th percentile")+
  theme_bw()
ggsave(filename = "densities.png", density2019 / density2020,
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
sharedoutliers$cohort<-c(rep("2019outlier",nrow(na.omit(fst2019cleanedoutliers[match(fst2020cleanedoutliers$POS, fst2019cleanedoutliers$POS),]))),
                         rep("2020outlier", nrow(na.omit(fst2020cleanedoutliers[match(fst2019cleanedoutliers$POS, fst2020cleanedoutliers$POS),]))))
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
freq_F0_2019<-fread("F0_2019.afreq",header=TRUE) #need to manually change #CHROM to CHROM in the afreq files before reading them in
freq_F1_2019<-fread("F1_2019.afreq",header=TRUE)
freq_F0_2020<-fread("F0_2020.afreq",header=TRUE)
freq_F1_2020<-fread("F1_2020.afreq",header=TRUE)
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
wholeframe$diff2019<-wholeframe$F1_2019_afreq-wholeframe$F0_2019_afreq
wholeframe$diff2020<-wholeframe$F1_2020_afreq-wholeframe$F0_2020_afreq

wholeframe_sharedoutliers<-subset(wholeframe,cohort=="sharedoutlier")
wholeframe_sharedoutliers$diff2019<-wholeframe_sharedoutliers$F1_2019_afreq-wholeframe_sharedoutliers$F0_2019_afreq
wholeframe_sharedoutliers$diff2020<-wholeframe_sharedoutliers$F1_2020_afreq-wholeframe_sharedoutliers$F0_2020_afreq

all<-ggplot(wholeframe, aes(x=diff2019,y=diff2020,color=cohort))+
  geom_point(size=1,alpha=0.5)+
  scale_color_manual(values=c("firebrick1","dodgerblue1","azure3","darkslateblue"))+
  geom_hline(yintercept=0)+
  geom_vline(xintercept = 0)+
  xlab("F1 allele freq minus F0 allele freq (2019)")+
  ylab("F1 allele freq minus F0 allele freq (2020)")+
  xlim(-0.8,0.8)+ylim(-0.8,0.8)+
  theme_bw()
nobg<-ggplot(subset(wholeframe,cohort!="background"), aes(x=diff2019,y=diff2020,color=cohort))+
  geom_point(size=1,alpha=0.5)+
  scale_color_manual(values=c("firebrick1","dodgerblue1","darkslateblue"))+
  geom_hline(yintercept=0)+
  geom_vline(xintercept = 0)+
  xlab("F1 allele freq minus F0 allele freq (2019)")+
  ylab("F1 allele freq minus F0 allele freq (2020)")+
  xlim(-0.8,0.8)+ylim(-0.8,0.8)+
  theme_bw()
shared<-ggplot(subset(wholeframe,cohort=="sharedoutlier"), aes(x=diff2019,y=diff2020,color=cohort))+
  geom_point(size=1,alpha=0.5)+
  scale_color_manual(values=c("darkslateblue"))+
  geom_hline(yintercept=0)+
  geom_vline(xintercept = 0)+
  xlab("F1 allele freq minus F0 allele freq (2019)")+
  ylab("F1 allele freq minus F0 allele freq (2020)")+
  xlim(-0.8,0.8)+ylim(-0.8,0.8)+
  theme_bw()
ggsave(filename="allelefreqplots.png",all | nobg | shared,
       width = 16, height = 7.5, dpi = 300, units = "in", device='png')

wholeframe_sharedoutliers<-rbind(subset(wholeframe_sharedoutliers, F1_2019_afreq<F0_2019_afreq & F1_2020_afreq< F0_2020_afreq),
                                 subset(wholeframe_sharedoutliers, F1_2019_afreq>F0_2019_afreq & F1_2020_afreq> F0_2020_afreq)
                                 
)
wholeframe_sharedoutliers_increasefreq<-subset(wholeframe_sharedoutliers, F1_2019_afreq>F0_2019_afreq & F1_2020_afreq> F0_2020_afreq)
wholeframe_sharedoutliers_decreasefreq<-subset(wholeframe_sharedoutliers, F1_2019_afreq<F0_2019_afreq & F1_2020_afreq< F0_2020_afreq)

#write.table(wholeframe_sharedoutliers, file="trulysharedoutliers20211018.txt",quote=FALSE,
#row.names=FALSE,sep = "\t")
write.table(wholeframe_sharedoutliers_increasefreq, file="trulysharedoutliers_increase20211018.txt",quote=FALSE,row.names=FALSE,sep = "\t")
write.table(wholeframe_sharedoutliers_decreasefreq, file="trulysharedoutliers_decrease20211018.txt",quote=FALSE,row.names=FALSE,sep = "\t")
wholeframe_sharedoutliers_2<-data.frame("CHROM"=rep(wholeframe_sharedoutliers$CHROM,4), "POS"=rep(wholeframe_sharedoutliers$POS,4),
                                        "chr_pos"=rep(wholeframe_sharedoutliers$chr_pos,4), "afreqs"=c(wholeframe_sharedoutliers$F0_2019_afreq, wholeframe_sharedoutliers$F1_2019_afreq,wholeframe_sharedoutliers$F0_2020_afreq, wholeframe_sharedoutliers$F1_2020_afreq),
                                        "cohort"=c(rep("F0_2019",length(wholeframe_sharedoutliers$F0_2019_afreq)), rep("F1_2019",length(wholeframe_sharedoutliers$F1_2019_afreq)), rep("F0_2020",length(wholeframe_sharedoutliers$F0_2020_afreq)), rep("F1_2020",length(wholeframe_sharedoutliers$F1_2020_afreq))))
#ggsave(filename = "trulysharedoutliers20211018.png",ggplot(wholeframe_sharedoutliers_2, aes(x=as.factor(chr_pos),y=afreqs,fill=cohort))+
         geom_bar(position="dodge", stat="identity") +
         scale_fill_brewer(palette="Paired")+
         theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
         theme(axis.text=element_text(size=15),
               axis.title=element_blank(),
               axis.text.x = element_text(angle = 90)),
       width = 8, height = 4, dpi = 300, units = "in", device='png')
wholeframe_sharedoutliers_increasefreq_2<-data.frame("CHROM"=rep(wholeframe_sharedoutliers_increasefreq$CHROM,4), "POS"=rep(wholeframe_sharedoutliers_increasefreq$POS,4),
                                        "chr_pos"=rep(wholeframe_sharedoutliers_increasefreq$chr_pos,4), "afreqs"=c(wholeframe_sharedoutliers_increasefreq$F0_2019_afreq, wholeframe_sharedoutliers_increasefreq$F1_2019_afreq,wholeframe_sharedoutliers_increasefreq$F0_2020_afreq, wholeframe_sharedoutliers_increasefreq$F1_2020_afreq),
                                        "cohort"=c(rep("F0_2019",length(wholeframe_sharedoutliers_increasefreq$F0_2019_afreq)), rep("F1_2019",length(wholeframe_sharedoutliers_increasefreq$F1_2019_afreq)), rep("F0_2020",length(wholeframe_sharedoutliers_increasefreq$F0_2020_afreq)), rep("F1_2020",length(wholeframe_sharedoutliers_increasefreq$F1_2020_afreq))))
ggsave(filename = "trulysharedoutliers_increasefreq20211014.png",ggplot(wholeframe_sharedoutliers_increasefreq_2, aes(x=as.factor(chr_pos),y=afreqs,fill=cohort))+
         geom_bar(position="dodge", stat="identity") +
         scale_fill_brewer(palette="Paired")+
         theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
         theme(axis.text=element_text(size=15),
               axis.title=element_blank(),
               axis.text.x = element_text(angle = 90)),
       width = 8, height = 4, dpi = 300, units = "in", device='png')
ggsave(filename = "fst_withtrulysharedoutliers20211018.png", ggplot(fst2020cleaned, aes(x = chr_pos, y = WEIR_AND_COCKERHAM_FST,color=cohort.outlier)) +
         scale_color_manual(values=c("skyblue","blue3","tomato","red3"))+
         geom_point(size = 3,alpha=0.1)+
         geom_point(data=fst2019cleaned,size=3,alpha=0.1,aes(chr_pos, WEIR_AND_COCKERHAM_FST, color=cohort.outlier))+
         #geom_point(data=wholeframe_sharedoutliers,size=3,alpha=0.7,aes(chr_pos, fst_2019))+
         #geom_point(data=wholeframe_sharedoutliers,size=3,alpha=0.7,aes(chr_pos, fst_2020))+
         
         #geom_hline(yintercept=mean(na.omit(fst2019cleaned$WEIR_AND_COCKERHAM_FST)))+
         #geom_hline(yintercept=mean(na.omit(fst2020cleaned$WEIR_AND_COCKERHAM_FST)))+
         theme_bw(),
       width = 8, height = 6, dpi = 300, units = "in", device='png')
fst2019cleaned$cohort.outlier<-paste(fst2019cleaned$cohort,fst2019cleaned$outliers,sep=".")
fst2020cleaned$cohort.outlier<-paste(fst2020cleaned$cohort,fst2020cleaned$outliers,sep=".")

ggsave(filename="trulyoutliers",ggplot(fst2020cleaned, aes(x = chr_pos, y = WEIR_AND_COCKERHAM_FST)) +
  scale_color_manual(values=c("blue3","tomato"))+
  geom_point(data=wholeframe_sharedoutliers,size=3,alpha=0.7,aes(chr_pos, fst_2019,color="blue3"))+
  geom_point(data=wholeframe_sharedoutliers,size=3,alpha=0.7,aes(chr_pos, fst_2020,color="tomato"))+
  
   theme_bw())
types<-c("nonsynonymous","synonymous","upstream","downstream","other")
trulysharedcounts<-c(15,26,173,122,(652-15-26-173-122))#,195,115,1,5)
fullvcfcounts<-c((681058+15878),570759,4882356,4874418, 
                 (16943427-681058-15878-570759-4882356-4874418))
trulyshared_increasecounts<-c(3,1,20,19,0)
trulyshared_decreasecounts<-c(12,25,153,103,0)#,165,112,1)

full_v_truly<-data.frame(fullvcfcounts,trulysharedcounts)
truly_v_decrease<-data.frame(trulysharedcounts, trulyshared_decreasecounts)
truly_v_increase<-data.frame(trulysharedcounts, trulyshared_increasecounts)
decrease_v_increase<-data.frame(trulyshared_decreasecounts, trulyshared_increasecounts)

chisq.test(full_v_truly)
chisq.test(truly_v_decrease)
chisq.test(truly_v_increase)
chisq.test(decrease_v_increase)
fullvcfpercent<-fullvcfcounts/16943427
trulysharedpercent<-trulysharedcounts/652
trulyshared_decreasepercent<-trulyshared_decreasecounts/576
trulyshared_increasepercent<-trulyshared_increasecounts/76

chiframe<-data.frame("percent"=100*(c(fullvcfpercent,trulysharedpercent, trulyshared_decreasepercent, trulyshared_increasepercent)),"types"=rep(types,4),
                     "dataset"=c(rep("all SNPs",length(types)),rep("shared outliers",length(types)),
                                 rep("shared outliers DECREASE FREQ",length(types)),rep("shared outliers INCREASE FREQ",length(types))))
chiframe2<-data.frame("counts"=c(trulysharedcounts, trulyshared_decreasecounts,trulyshared_increasecounts),"types"=rep(types,3),
                     "dataset"=c(rep("all shared outliers",length(types)),rep("shared outliers DECREASE FREQ",length(types)),rep("shared outliers INCREASE FREQ",length(types))))
chiframecounts<-data.frame("counts"=c(fullvcfcounts,trulysharedcounts, trulyshared_decreasecounts, trulyshared_increasecounts),"types"=rep(types,4),
                                     "dataset"=c(rep("all SNPs",length(types)),rep("shared outliers",length(types)),
                                                 rep("shared outliers DECREASE FREQ",length(types)),rep("shared outliers INCREASE FREQ",length(types))))
chiframecounts<-chiframecounts[1:10,]
chiframecounts$dataset<-factor(chiframecounts$dataset, levels = c("all SNPs","shared outliers"))
chiframecounts$types<-factor(chiframecounts$types, levels = types)
chiframe<-chiframe[1:10,]
chiframe$dataset<-factor(chiframe$dataset, levels = c("all SNPs","shared outliers"))
chiframe$types<-factor(chiframe$types, levels = types)

ggsave(file="SNPtypes20211019.png", ggplot(chiframecounts,aes(x=dataset,y=counts,color=types))+
  geom_bar(stat="identity",aes(fill=types),position="dodge")+
  scale_y_continuous(trans='log10')+
  theme_bw(),
  width = 8, height = 6, dpi = 300, units = "in", device='png')
  
ggsave(file="SNPtypespercents20211019.png",ggplot(chiframe,aes(x=dataset,y=percent,color=types))+
  geom_bar(stat="identity",aes(fill=types),position="dodge")+
  theme_bw(),
  width = 8, height = 6, dpi = 300, units = "in", device='png')

#plink heatmap
cluster1 <- read.csv("plink.genome.editednames.csv", header = TRUE)
pairs_c1 <- cluster1[,c(1,3,10)]
colnames(pairs_c1)[1] <- "ind1"
colnames(pairs_c1)[2] <- "ind2"
orderedpairs_c1 <- pairs_c1[with(pairs_c1, order(ind1, ind2)), ]
orderedpairs_c1$ind1<-factor(orderedpairs_c1$ind1, levels=c("19-1","19-10","19-11","19-12","19-13","19-14","19-15","19-16","19-17","19-18","19-19","19-2","19-20","19-21",
                                                     "19-22","19-24","19-3","19-4","19-5","19-6","19-7","19-8","19-9","20-1","20-10","20-2","20-4","20-5",
                                                     "20-6","20-7","20-8","20-9","57","58","62","70","71","72","73","74","76","77","79","80","56","60"))
orderedpairs_c1$ind2<-factor(orderedpairs_c1$ind2, levels=c("19-10","19-11","19-12","19-13","19-14","19-15","19-16","19-17","19-18","19-19","19-2","19-20","19-21","19-22",
                                                            "19-24","19-3","19-4","19-5","19-6","19-7","19-8","19-9","20-1","20-10","20-2","20-4","20-5","20-6",
                                                            "20-7","20-8","20-9","57","58","62","70","71","72","73","74",
                                                            "76","77","79","80","56","60","65"))

ggsave(file="plinkheatmap2021102.png",ggplot(data=orderedpairs_c1,aes(x=ind1, y= ind2, fill=PI_HAT)) + geom_tile(color="white")+
  scale_fill_gradient2(low="blue", mid="white",midpoint=max(orderedpairs_c1$PI_HAT)/2,high="red",limit=c(0,max(orderedpairs_c1$PI_HAT)),
                      space="Lab", name="Pi_hat")+
  
  theme(axis.text.x = element_text(angle = 90), axis.title = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
  width = 8, height = 8, dpi = 300, units = "in", device='png')
#ggsave(file="plinkheatmap20211018.png",ggplot(data=orderedpairs_c1,aes(x=ind1, y= ind2, fill=PI_HAT)) + geom_tile(color="white")+
  #scale_fill_gradient(low="blue",high="red",limit=c(0,max(orderedpairs_c1$PI_HAT)),
                       #space="Lab", name="Pi_hat")+
  #theme(axis.text.x = element_text(angle = 90), axis.title = element_blank()),
 # width = 8, height = 8, dpi = 300, units = "in", device='png')
fdc<-data.frame("chrpos"=rep(wholeframe_sharedoutliers$chr_pos,2),"fst"=c(wholeframe_sharedoutliers$fst_2019, 
              wholeframe_sharedoutliers$fst_2020), "year"=c(rep("2019",nrow(wholeframe_sharedoutliers)),
                                                                rep("2020",nrow(wholeframe_sharedoutliers))))

ggsave(file="justtrueoutliers.png",ggplot(fdc, aes(chrpos, fst, color=year))+
         scale_color_manual(values=c("blue3","tomato"))+
  geom_point(alpha=0.7)+theme_bw(),
  width = 8, height = 4, dpi = 300, units = "in", device='png')

##plink LD heatmap F0_2020:
#plink heatmap
cluster1 <- read.csv("/Users/elopez-nandam/Documents/CaptiveBreeding/data/just652outliers_F0_2020_no707273.ld.csv", header = TRUE)
chr1<-subset(cluster1, startsWith(cluster1$SNP_A, 'chr1:')==TRUE & startsWith(cluster1$SNP_B, 'chr1:')==TRUE)
ordered <- chr1[with(chr1, order(BP_A, BP_B)), ]
chr2<-subset(cluster1, startsWith(cluster1$SNP_A, 'chr2')==TRUE & startsWith(cluster1$SNP_B, 'chr2')==TRUE)
ordered <- chr2[with(chr2, order(BP_A, BP_B)), ]
chr3<-subset(cluster1, startsWith(cluster1$SNP_A, 'chr3')==TRUE & startsWith(cluster1$SNP_B, 'chr3')==TRUE)
ordered <- chr3[with(chr3, order(BP_A, BP_B)), ]
ordered$SNP_A<-factor(ordered$SNP_A, levels=unique(ordered$SNP_A))
ordered$SNP_B<-factor(ordered$SNP_B, levels=unique(ordered$SNP_B))

pairs_c1 <- cluster1[,c(3,6,7)]
colnames(pairs_c1)[1] <- "snp1"
colnames(pairs_c1)[2] <- "snp2"
pairs_c1$snp1<-factor(pairs_c1$snp1, levels=unique(pairs_c1$snp1))
pairs_c1$snp2<-factor(pairs_c1$snp2, levels=unique(pairs_c1$snp2))

#orderedpairs_c1 <- pairs_c1[with(pairs_c1, order(snp1, snp2)), ]
pairs_chr1<-subset(pairs_c1, startsWith(pairs_c1$snp1, 'chr1:')==TRUE & startsWith(pairs_c1$snp2, 'chr1:')==TRUE)
colorvector<-rep(c("red","blue","black"),length(unique(ordered$SNP_A))/3)
                 
ggsave(file="LDheatmap20211130_chr3_no707273.png",ggplot(data=ordered,aes(x=SNP_A, y= SNP_B, fill=R2)) + geom_tile(color="white")+
  scale_fill_gradient2(low="blue", mid="white",midpoint=max(ordered$R2)/2,high="red",limit=c(0,max(ordered$R2)),
                       space="Lab", name="Pi_hat")+
  
  theme(axis.text.x = element_text(angle = 90), axis.title = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
  width = 8, height = 8, dpi = 300, units = "in", device='png')

ggsave(file="LDheatmap20211130_allchrs_no707273.png",ggplot(data=pairs_c1,aes(x=snp1, y= snp2, fill=R2)) + geom_tile(color="white")+
  scale_fill_gradient2(low="blue", mid="white",midpoint=max(pairs_c1$R2)/2,high="red",limit=c(0,max(pairs_c1$R2)),
                       space="Lab", name="Pi_hat")+
  
  theme(axis.text.x = element_text(angle = 90), axis.title = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
  width = 8, height = 8, dpi = 900, units = "in", device='png')

