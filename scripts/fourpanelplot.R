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
wholeframe$diff2019abs<-abs(wholeframe$diff2019)
wholeframe$diff2020abs<-abs(wholeframe$diff2020)
wholeframe_sharedoutliers$diff2019abs<-abs(wholeframe_sharedoutliers$diff2019)
wholeframe_sharedoutliers$diff2020abs<-abs(wholeframe_sharedoutliers$diff2020)

all<-ggplot(wholeframe, aes(x=diff2019,y=diff2020,color=cohort))+
  geom_point(size=1,alpha=0.1, position=position_jitter(width=0.01, height=0.01))+
  scale_color_manual(values=c("firebrick1","dodgerblue1","azure3","darkslateblue"))+
  geom_hline(yintercept=0)+
  geom_vline(xintercept = 0)+
  xlab("F1 allele freq minus F0 allele freq (2019)")+
  ylab("F1 allele freq minus F0 allele freq (2020)")+
  xlim(-0.8,0.8)+ylim(-0.8,0.8)+
  theme_bw()
nobg<-ggplot(subset(wholeframe,cohort!="background"), aes(x=diff2019,y=diff2020,color=cohort))+
  geom_point(size=1,alpha=0.1, position=position_jitter(width=0.01, height=0.01))+
  scale_color_manual(values=c("firebrick1","dodgerblue1","darkslateblue"))+
  geom_hline(yintercept=0)+
  geom_vline(xintercept = 0)+
  xlab("F1 allele freq minus F0 allele freq (2019)")+
  ylab("F1 allele freq minus F0 allele freq (2020)")+
  xlim(-0.8,0.8)+ylim(-0.8,0.8)+
  theme_bw()
shared<-ggplot(subset(wholeframe,cohort=="sharedoutlier"), aes(x=diff2019,y=diff2020,color=cohort))+
  geom_point(size=1,alpha=0.1, position=position_jitter(width=0.01, height=0.01))+
  scale_color_manual(values=c("darkslateblue"))+
  geom_hline(yintercept=0)+
  geom_vline(xintercept = 0)+
  xlab("F1 allele freq minus F0 allele freq (2019)")+
  ylab("F1 allele freq minus F0 allele freq (2020)")+
  xlim(-0.8,0.8)+ylim(-0.8,0.8)+
  theme_bw()
r=0.3
xc=0.35
yc=0.5
r2=0.3
xc2=-0.35
yc2=-0.4
shared_ann<-shared+annotate("path",
                     x=xc+r*cos(seq(0,2*pi,length.out=100)),
                     y=yc+r*sin(seq(0,2*pi,length.out=100)),color="purple")+
  annotate("path",
                    x=xc2+r2*cos(seq(0,2*pi,length.out=100)),
                    y=yc2+r2*sin(seq(0,2*pi,length.out=100)),color="purple")
wholeframe_sharedoutliers<-subset(wholeframe,cohort=="sharedoutlier")
wholeframe_sharedoutliers$diff2019<-wholeframe_sharedoutliers$F1_2019_afreq-wholeframe_sharedoutliers$F0_2019_afreq
wholeframe_sharedoutliers$diff2020<-wholeframe_sharedoutliers$F1_2020_afreq-wholeframe_sharedoutliers$F0_2020_afreq
wholeframe$diff2019abs<-abs(wholeframe$diff2019)
wholeframe$diff2020abs<-abs(wholeframe$diff2020)
wholeframe_sharedoutliers$diff2019abs<-abs(wholeframe_sharedoutliers$diff2019)
wholeframe_sharedoutliers$diff2020abs<-abs(wholeframe_sharedoutliers$diff2020)

mean2019abs<-mean(wholeframe$diff2019abs)
mean2020abs<-mean(wholeframe$diff2020abs)
sd2019abs<-sd(wholeframe$diff2019abs)
sd2020abs<-sd(wholeframe$diff2020abs)
median2019abs<-median(wholeframe$diff2019abs)
median2020abs<-median(wholeframe$diff2020abs)

mean2019abs_shared<-mean(wholeframe_sharedoutliers$diff2019abs)
mean2020abs_shared<-mean(wholeframe_sharedoutliers$diff2020abs)
sd2019abs_shared<-sd(wholeframe_sharedoutliers$diff2019abs)
sd2020abs_shared<-sd(wholeframe_sharedoutliers$diff2020abs)
median2019abs_shared<-median(wholeframe_sharedoutliers$diff2019abs)
median2020abs_shared<-median(wholeframe_sharedoutliers$diff2020abs)


absdiffsframe<-data.frame(means=c(mean2019abs, mean2020abs, mean2019abs_shared, mean2020abs_shared),
                          sds=c(sd2019abs, sd2020abs, sd2019abs_shared, sd2020abs_shared),
                          year=c("2019","2020","2019","2020"),
                          type=c("All SNPs","All SNPs","652 outliers","652 outliers"))

absdiff<-ggplot(absdiffsframe, aes(x=year, y=means, color=type))+
  geom_point(position=position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=means-(sds*2), ymax=means+(sds*2)),
                width=0.1, position=position_dodge(width=0.3))+
  scale_color_manual(values=c("black","red"))+ylab("Absolute change in allele frequency from spawners to offspring")+
  #scale_y_continuous(expand = c(0,0.00))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
ggsave(filename="allelefreqplots_20220519.png", (all | nobg) / (shared_ann | absdiff),
       width = 16, height = 7.5, dpi = 300, units = "in", device='png')

