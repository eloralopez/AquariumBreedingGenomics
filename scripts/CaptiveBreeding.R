setwd("~/Documents/CaptiveBreeding/plots")
setwd("~/Documents/CaptiveBreeding/data")

library(stringr)
library(ggplot2)
library(kinship2)
library(patchwork)
library(ggforce)
library(rehh)
library(data.table)
library(ggpubr)
##PCA##
#snps<-read.delim('/Users/elopez-nandam/Documents/CaptiveBreeding/data/CSLparentsandF1_chr1_no2020-3_allsites_genotyped_nomissing_dp10gq30_biSNPs.012', header=F)
snps<-fread('/Users/elopez-nandam/Documents/CaptiveBreeding/data/CSLparentsandF1_chr1_no2020-3_allsites_genotyped_nomissing_dp10gq30_biSNPs.012', header=F)

names<-fread('/Users/elopez-nandam/Documents/CaptiveBreeding/data/CSLparentsandF1_chr1_no2020-3_allsites_genotyped_nomissing_dp10gq30_biSNPs.012.indv',header=F)
#snps<-fread('/Users/elopez-nandam/Documents/CaptiveBreeding/plots/CA2020-8shareduniques.012', header=F)
#names<-fread('/Users/elopez-nandam/Documents/CaptiveBreeding/plots/CA2020-8shareduniques.012.indv', header=F)
shortnames<-c("19-10","19-11","19-12","19-13","19-14","19-15","19-16","19-17","19-18","19-19","19-1","19-20","19-21",
  "19-22","19-24","19-2", "19-3","19-4","19-5","19-6","19-7","19-8","19-9","20-10","20-1","20-2","20-4","20-5",
  "20-6","20-7","20-8","20-9","57","58","62","70","71","72","73","76","77","79","80","56","60","65","74")
rownames(snps)<-shortnames
row.names(snps)
#snps<-as.matrix(snps[,-1])
pc.out<-prcomp(snps)

plot(jitter(pc.out$x[,1]),jitter(pc.out$x[,2]),cex=0, xlab="PC1", ylab="PC2")#,xlim=c(20,60),ylim=c(-10,10))
text(pc.out$x[,1], pc.out$x[,2],label=rownames(snps), cex=0.5)
plot(hclust(dist(snps)))
df_out <- as.data.frame(pc.out$x)#, pc.out$PC2, names)
names2<-shortnames
ggplot(df_out,aes(x=PC1,y=PC2))+
  geom_text(aes(label=shortnames),hjust=2, vjust=2, size=1)+theme_bw()
  theme_bw(), width = 8, height = 4, dpi = 300, units = "in", device='png')
pc.out
summary(pc.out)
E<-subset(snps, startsWith(row.names(snps), "CA70")==FALSE & startsWith(row.names(snps), "CA72")==FALSE & startsWith(row.names(snps), "CA73")==FALSE)
E<-subset(snps, startsWith(row.names(snps), "CA2019")==TRUE | startsWith(row.names(snps), "CAP")==TRUE)
E<-subset(snps, startsWith(row.names(snps), "CA2019")==FALSE & startsWith(row.names(snps), "CA2020")==FALSE &
            startsWith(row.names(snps), "CA70")==FALSE & startsWith(row.names(snps), "CA72")==FALSE & startsWith(row.names(snps), "CA73")==FALSE)

E<-subset(snps, startsWith(row.names(snps), "CAP")==FALSE & startsWith(row.names(snps), "CA2019")==FALSE & startsWith(row.names(snps), "CA70")==FALSE & startsWith(row.names(snps), "CA72")==FALSE & startsWith(row.names(snps), "CA73")==FALSE)# &
          startsWith(row.names(snps), "CA2020-1")==FALSE &
            startsWith(row.names(snps), "CA2020-2")==FALSE &
            startsWith(row.names(snps), "CA2020-3")==FALSE &
            startsWith(row.names(snps), "CA2020-4")==FALSE &
            startsWith(row.names(snps), "CA2020-5")==FALSE &
            startsWith(row.names(snps), "CA2020-6")==FALSE &
            startsWith(row.names(snps), "CA2020-9")==FALSE &
            startsWith(row.names(snps), "CA2020-10")==FALSE)

Eout<-prcomp(E)
df_out <- as.data.frame(Eout$x)#, pc.out$PC2, names)
#names2<-names$V1[c(30:31,33:35,37,40:43,47)]
names2<-names$V1[c(24:35,37,40:43,47)]
rownames(E)<-names2
ggplot(df_out,aes(x=PC1,y=PC2))+
  geom_point()+geom_text(aes(label=names2),hjust=0, vjust=0, size=2)#+
theme_bw(), width = 8, height = 4, dpi = 300, units = "in", device='png')

plot(jitter(Eout$x[,1]),jitter(Eout$x[,2]),cex=0, xlab="PC1", ylab="PC2")#,xlim=c(-10,10),ylim=c(20,35))
text(Eout$x[,1], Eout$x[,2],label=rownames(E), cex=0.5)
plot(hclust(dist(E)))
Eout
summary(Eout)

freqF0<-read.delim("freqchr1F0only.frq",header=TRUE)
split2F0<-str_split_fixed(freqF0$ALLELE_FREQ2, ":", 2)
minor_alleleF0<-as.numeric(split2F0[,2])


freqF1<-read.delim("freqchr1F1only.frq",header=TRUE)
#split1<-str_split_fixed(freqs$ALLELE_FREQ1, ":", 2)
split2F1<-str_split_fixed(freqF1$ALLELE_FREQ2, ":", 2)
minor_alleleF1<-as.numeric(split2F1[,2])
difference<-minor_alleleF1 - minor_alleleF0
#plot(freqF0$POS, difference)
freqs<-data.frame(freqF0$POS, minor_alleleF0, minor_alleleF1, difference)
freqs2<-data.frame("maf"=c(minor_alleleF0,minor_alleleF1),"f"=c(rep("F0",length(minor_alleleF0)),rep("F1",length(minor_alleleF1))))
ggplot(freqs2, aes(x=maf,color=f))+
  geom_density(alpha=0.5)  +theme_bw()
subset(freqs, difference>0.4)
freqF1$diff<-difference
ggplot(freqs, aes(x = minor_alleleF0, y = minor_alleleF1)) +
  geom_point(size = 3)
differenceplot<-ggplot(freqs, aes(x = freqF0$POS, y = difference)) +
  geom_point(size = 3, alpha=0.4)+ theme_bw()

fstall<-read.delim("F1_vs_F0.weir.fst",header=TRUE)
fst2019<-read.delim("F1_vs_F0_2019.weir.fst", header=TRUE)
#fst2019<-na.omit(fst2019)
fst2020<-read.delim("F1_vs_F0_2020.weir.fst", header=TRUE)
#fst2020<-na.omit(fst2020)

plot(difference, fstall$WEIR_AND_COCKERHAM_FST)
ggplot(fst2019, aes(x = POS, y = WEIR_AND_COCKERHAM_FST)) +
  geom_point(size = 3)

fstplot<-ggplot(fst2020, aes(x = POS, y = WEIR_AND_COCKERHAM_FST)) +
  geom_point(size = 3,color="blue",alpha=0.4)+
  geom_point(data=fst2019,size=3,alpha=0.4,aes(POS, WEIR_AND_COCKERHAM_FST, color="red"))+
  theme_bw()
 #fst distribs
#  hist(fst2020$WEIR_AND_COCKERHAM_FST)
#hist(fst2019$WEIR_AND_COCKERHAM_FST)
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
fstplot | compareplot
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
  differences<-setdiff(dataframe$chr_pos, fst2020$chr_pos)

  g<-subset(dataframe, chr_pos != differences[1] & chr_pos != differences[2] &
         chr_pos != differences[3] & chr_pos != differences[4] &
         chr_pos != differences[5] & chr_pos != differences[6] &
         chr_pos != differences[7] & chr_pos != differences[8] &
         chr_pos != differences[9] & chr_pos != differences[10] &
         chr_pos != differences[11] & chr_pos != differences[12] &
         chr_pos != differences[13] & chr_pos != differences[14] &
         chr_pos != differences[15] & chr_pos != differences[16] &
         chr_pos != differences[17] & chr_pos != differences[18])
  return(g)
}

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
#fst2019$outliers<-sapply(fst2019$WEIR_AND_COCKERHAM_FST, testfunction)
#fst2020$outliers<-sapply(fst2020$WEIR_AND_COCKERHAM_FST, testfunction)

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
write.table(sharedoutliers, file="sharedoutliers20210928.txt",quote=FALSE,
            row.names=FALSE,sep = "\t")
write.table(fst2019outliers, file="fst2019outliers20210928.txt",quote=FALSE,
            row.names=FALSE,sep = "\t")
write.table(fst2020outliers, file="fst2020outliers20210928.txt",quote=FALSE,
            row.names=FALSE,sep = "\t")
ggsave(filename = "fstplot20210924.png", ggplot(fst2020, aes(x = chr_pos, y = WEIR_AND_COCKERHAM_FST,color=cohort)) +
  scale_color_manual(values=c("skyblue","blue3","tomato","red3"))+
  geom_point(size = 3,alpha=0.1)+
  geom_point(data=fst2019,size=3,alpha=0.1,aes(chr_pos, WEIR_AND_COCKERHAM_FST, color=cohort))+
  geom_point(data=sharedoutliers,size=3,alpha=0.7,aes(chr_pos, WEIR_AND_COCKERHAM_FST, color=cohort))+
  geom_hline(yintercept=mean(fst2019$WEIR_AND_COCKERHAM_FST))+
    geom_hline(yintercept=mean(fst2020$WEIR_AND_COCKERHAM_FST))+theme_bw(),
  width = 8, height = 4, dpi = 300, units = "in", device='png')
ggsave(filename = "fstplot20210924_sharedoutliers.png", ggplot(fst2020, aes(x = chr_pos, y = WEIR_AND_COCKERHAM_FST)) +
         #geom_point(size = 3,alpha=0.1)+
         #geom_point(data=fst2019,size=3,alpha=0.1,aes(chr_pos, WEIR_AND_COCKERHAM_FST, color=cohort))+
         scale_color_manual(values=c("blue3","red3"))+
         geom_point(data=sharedoutliers,size=3,alpha=0.7,aes(chr_pos, WEIR_AND_COCKERHAM_FST, color=cohort))+ theme_bw(),
       width = 8, height = 4, dpi = 300, units = "in", device='png')
#ggsave(filename = "fstplot.png", ggplot(fst2020, aes(x = chr_pos_2020, y = WEIR_AND_COCKERHAM_FST)) +
         geom_point(size = 1,color="blue",alpha=0.4)+
         geom_point(data=fst2019,size=1,alpha=0.4,aes(chr_pos_2019, WEIR_AND_COCKERHAM_FST, color="red"))+
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
fst_2020<-fst2019$WEIR_AND_COCKERHAM_FST
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
write.table(wholeframe_sharedoutliers, file="trulysharedoutliers20210928.txt",quote=FALSE,
            row.names=FALSE,sep = "\t")
wholeframe_sharedoutliers_2<-data.frame("CHROM"=rep(wholeframe_sharedoutliers$CHROM,4), "POS"=rep(wholeframe_sharedoutliers$POS,4),
                                        "chr_pos"=rep(wholeframe_sharedoutliers$chr_pos,4), "afreqs"=c(wholeframe_sharedoutliers$F0_2019_afreq, wholeframe_sharedoutliers$F1_2019_afreq,wholeframe_sharedoutliers$F0_2020_afreq, wholeframe_sharedoutliers$F1_2020_afreq),
                                        "cohort"=c(rep("F0_2019",length(wholeframe_sharedoutliers$F0_2019_afreq)), rep("F1_2019",length(wholeframe_sharedoutliers$F1_2019_afreq)), rep("F0_2020",length(wholeframe_sharedoutliers$F0_2020_afreq)), rep("F1_2020",length(wholeframe_sharedoutliers$F1_2020_afreq))))
ggplot(wholeframe, aes(x=chr_pos,y=fst_2019))+
         geom_point(size=3)
ggplot(wholeframe_sharedoutliers_2, aes(x=as.factor(chr_pos),y=afreqs,fill=cohort))+
  geom_bar(position="dodge", stat="identity") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=15),
        axis.title=element_blank(),
        axis.text.x = element_text(angle = 90))
  geom_col(position="dodge",aes(x=as.factor(chr_pos),y=F1_2019_afreq, color="red"))+
  ylim(0.05,1)+theme_bw()
ndF0<-read.delim("chr1nucdivF0.sites.pi")
ndF1<-read.delim("chr1nucdiv.sites.pi")
f<-c(rep("F0",nrow(ndF0)),rep("F1",nrow(ndF1)))
nd<-data.frame(f, "POS"=c(ndF0$POS,ndF1$POS), "PI"=c(ndF0$PI, ndF1$PI))

ggplot(nd, aes(x=f, y=PI))+
  geom_violin()+theme_bw()+
  geom_sina(data=nd, aes(f,PI,color="red"))

ggplot(nd, aes(x=PI, color=f))+
  geom_density(alpha=0.5)+theme_bw()
  
#rehh##
head(readLines("chr1shapeittest.haps"), sep = "\n")
cat(readLines("plink.map"), sep = "\n")
haps<-read.delim("chr1shapeittest.haps",header=FALSE,  sep = "\n")
hh <- data2haplohh(hap_file ="chr1shapeittest.haps",
                   map_file = "plink.map",
                   allele_coding = "01")

library(fsthet)
gfile<-system.file("CSLparentsandF1_chr1_no2020-3_allsites_genotyped_nomissing_dp10gq30_biSNPs.genepop",package = 'fsthet')
gpop<-my.read.genepop(gfile)

#chr1germline<-fread("CA2019chr1germlinemuts.txt",header=F)
chr1germline<-fread("germlinemutationcountschr8-14.txt",header=F)

chr1germline$meanvaf<-as.numeric(as.character(chr1germline$V4))
chr1germline<-na.omit(chr1germline)
ggsave(file="mutationschr8-14.png",ggplot(chr1germline, aes(V2, meanvaf))+
  geom_point()+    geom_smooth(color="red", fill="blue", se=TRUE, level=0.95) + xlim(2,22)+
  #stat_cor(label.x=3,label.y=0.55, aes(group=1), color="red") +
  #stat_regline_equation(label.x = 3, label.y = 0.6,color="red")+xlim(2,22)+
  theme_bw(),
  width = 8, height = 4, dpi = 300, units = "in", device='png')
subset(chr1germline,V2==22)
ggplot(chr1germline, aes(V2))+
  geom_histogram(binwidth = 1)+xlim(0,22)+geom_vline(xintercept=1,color="red")+
  theme_bw()
nrow(chr1germline
     )
depth<-read.delim("CaptiveBreedingDepths.txt",header=TRUE)
qpois(0.0001,depth$AverageDepth,lower.tail=T) 
qpois(0.9999,depth$AverageDepth,lower.tail=T) 
depths_func<-function(samplenames, avgdepths) {
  upperbounds<-qpois(0.9999, avgdepths,lower.tail=T)  
  lowerbounds<-qpois(0.0001,avgdepths,lower.tail=T)  
  depthdf<-data.frame("samplename"=samplenames, "avgdepth"=avgdepths,"upper"=upperbounds,"lower"=lowerbounds)
  #return(depthdf)
  for (i in 1:length(samplenames)) {
    #sc<-print(paste0("vc.getGenotype(", samplenames[i],".getDP()<", upperbounds[i], 
    #"&& vc.getGenotype(", samplenames[i],").getDP()>",lowerbounds[i], "&&"))
    pc<-cat("vc.getGenotype(", "\"\'", samplenames[i],"\'\"", ")",".getDP()<", upperbounds[i], 
            " && vc.getGenotype(", "\"\'", samplenames[i],"\'\"", ").getDP()>",lowerbounds[i], " &&", "\n", sep = "")
    pc
  } 
  return(pc)
}  
depths_func(depth$sample, depth$AverageDepth)

agrilife<-subset(depth, startsWith(sample, "CAP")==FALSE)
CAPs<-subset(depth, startsWith(sample, "CAP"))

