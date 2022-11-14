#usage example: Rscript af_and_hw.R ../F0_just652_AF.afreq ../F1_just652_AF.afreq ../F1_just652.recode.vcf_genocounts.gcount ../foursibs2020_just652_F0_AF.afreq ../foursibs2020_just652_F1_AF.afreq ../foursibs2020_just652_F1.recode.vcf_genocounts.gcount
#packages installed on rosalindf as of 4/14/22
library(patchwork)
library(ggplot2)
library(ggforce)

#read in arguments:
args = commandArgs(trailingOnly=TRUE)

#read in files:
afreq_1<-read.delim(args[1],header=TRUE)

afreq_2<-read.delim(args[2],header=TRUE)

f1_genocounts<-read.delim(args[3],header=TRUE)

afreq_1_652<-read.delim(args[4],header=TRUE)

afreq_2_652<-read.delim(args[5],header=TRUE)

f1_genocounts_652<-read.delim(args[6],header=TRUE)

plotname<-args[7]

freq_func<-function(afreq_1, afreq_2, f1_genocounts) {
  
  F0freq<-afreq_1$ALT_FREQS
  
  F1freq<-afreq_2$ALT_FREQS
  
  f1_homreffreq<-f1_genocounts$HOM_REF_CT/(f1_genocounts$HOM_REF_CT+f1_genocounts$HET_REF_ALT_CTS+f1_genocounts$TWO_ALT_GENO_CTS)
  
  f1_hetfreq<-f1_genocounts$HET_REF_ALT_CTS/(f1_genocounts$HOM_REF_CT+f1_genocounts$HET_REF_ALT_CTS+f1_genocounts$TWO_ALT_GENO_CTS)
  
  f1_homaltfreq<-f1_genocounts$TWO_ALT_GENO_CTS/(f1_genocounts$HOM_REF_CT+f1_genocounts$HET_REF_ALT_CTS+f1_genocounts$TWO_ALT_GENO_CTS)
  F1_alt_counts<-round(afreq_2$ALT_FREQS*afreq_2$OBS_CT,digits=0)
  F1_ref_counts<-afreq_2$OBS_CT-round(afreq_2$ALT_FREQS*afreq_2$OBS_CT,digits=0)
  
  freqcomp<-data.frame(F0freq, F1freq, "chrom"=afreq_1$X.CHROM,"pos"=afreq_1$POS,
                       f1_homreffreq, f1_hetfreq, f1_homaltfreq, F1_alt_counts, F1_ref_counts)
  
  trimmed<-subset(freqcomp, 0 < F0freq & F0freq < 1)
  chi_null=rep("A",nrow(trimmed))
  for (i in 1:nrow(trimmed)) {
    obs_vec<-c(trimmed[i,]$F1_alt_counts, trimmed[i,]$F1_ref_counts)
    exp_vec<-c(trimmed[i,]$expected_alt_count, trimmed[i,]$expected_ref_count)
    if (trimmed$F0freq[i] == 0.25) {
      chi_p <- chisq.test(obs_vec, p=c(1/4,3/4)) } else if (trimmed$F0freq[i] == 0.5) {
        chi_p <- chisq.test(obs_vec, p=c(1/2, 1/2)) } else if (trimmed$F0freq[i]==0.75) {
          chi_p <- chisq.test(obs_vec, p=c(3/4, 1/4)) }  
    
    chi_null[i]<-chi_p$p.value
  }  
  
  trimmed$chi_null<-as.numeric(as.character(chi_null))
  print(nrow(subset(trimmed, chi_null>=0.05)))
  return(freqcomp)
  #return(trimmed)
}

freq_1<-freq_func(afreq_1, afreq_2, f1_genocounts)
freq_2<-freq_func(afreq_1_652, afreq_2_652, f1_genocounts_652)
freq_3<-freq_func(read.delim("bigfamily2020_allsnps_F0.recode.vcf.afreq",header=TRUE), read.delim("bigfamily2020_allsnps_F1.recode.vcf.afreq"), read.delim("bigfamily2020_allsnps_F1.recode.vcf.gcount"))
#freq_3<-freq_func(read.delim("../foursibs2020_just652_F0_AF.afreq"), read.delim("../foursibs2020_just652_F1_AF.afreq"), read.delim("../foursibs2020_just652_F1.recode.vcf_genocounts.gcount"))

freq_4<-freq_func(read.delim("bigfamily2020_just887_F0.recode.vcf.afreq"), read.delim("bigfamily2020_just887_F1.recode.vcf.afreq"), read.delim("bigfamily2020_just887_F1.recode.vcf.gcount"))


comb<-rbind(freq_1,freq_2)
comb$Dataset<-c(rep("All SNPs",nrow(freq_1)),rep("Shared Outliers",nrow(freq_2)))

bigfam0.5<-ggplot(subset(comb,F0freq==0.5), aes(x=F1freq,fill=Dataset))+
  geom_density(alpha=0.4,adjust=2,bw=0.003)+ylim(0,45)+xlab("2019 full siblings ALT allele frequency")+ylab("Probability Density")+
  scale_fill_manual(values=c("red","blue"))+scale_color_manual(values=c("red","blue"))+ theme_classic()+labs(tag = "B.")+
  ggtitle("Parents: CA56 x CA60, ALT allele frequency = 0.5")+ theme(plot.title = element_text(hjust = 0.5))

bigfam0.25<-ggplot(subset(comb,F0freq==0.25), aes(x=F1freq,fill=Dataset))+
  geom_density(alpha=0.4,adjust=2,bw=0.003)+ylim(0,45)+xlab("2019 full siblings ALT allele frequency")+ylab("Probability Density")+
  scale_fill_manual(values=c("red","blue"))+scale_color_manual(values=c("red","blue"))+ theme_classic()+labs(tag = "A.")+
  ggtitle("Parents: CA56 x CA60, ALT allele frequency = 0.25")+ theme(plot.title = element_text(hjust = 0.5))

bigfam0.75<-ggplot(subset(comb,F0freq==0.75), aes(x=F1freq,fill=Dataset))+
  geom_density(alpha=0.4,adjust=2,bw=0.003)+ylim(0,45)+xlab("2019 full siblings ALT allele frequency")+ylab("Probability Density")+
  scale_fill_manual(values=c("red","blue"))+scale_color_manual(values=c("red","blue"))+ theme_classic()+labs(tag = "C.")+
  ggtitle("Parents: CA56 x CA60, ALT allele frequency = 0.75")+ theme(plot.title = element_text(hjust = 0.5))

comb2<-rbind(freq_3,freq_4)
comb2$Dataset<-c(rep("All SNPs",nrow(freq_3)),rep("Shared Outliers",nrow(freq_4)))

bigfam20.5<-ggplot(subset(comb2,F0freq==0.5), aes(x=F1freq,fill=Dataset))+
  geom_density(alpha=0.4,adjust=2,bw=0.003)+ylim(0,45)+xlab("2020 full siblings ALT allele frequency")+ylab("Probability Density")+
  scale_fill_manual(values=c("red","blue"))+scale_color_manual(values=c("red","blue"))+ theme_classic()+labs(tag = "E.")+
  ggtitle("Parents: CA74 x CA80, ALT allele frequency = 0.5")+ theme(plot.title = element_text(hjust = 0.5))
bigfam20.25<-ggplot(subset(comb2,F0freq==0.25), aes(x=F1freq,fill=Dataset))+
  geom_density(alpha=0.4,adjust=2,bw=0.003)+ylim(0,45)+xlab("2020 full siblings ALT allele frequency")+ylab("Probability Density")+
  scale_fill_manual(values=c("red","blue"))+scale_color_manual(values=c("red","blue"))+ theme_classic()+labs(tag = "D.")+
  ggtitle("Parents: CA74 x CA80, ALT allele frequency = 0.25")+ theme(plot.title = element_text(hjust = 0.5))
bigfam20.75<-ggplot(subset(comb2,F0freq==0.75), aes(x=F1freq,fill=Dataset))+
  geom_density(alpha=0.4,adjust=2,bw=0.003)+ylim(0,45)+xlab("2020 full siblings ALT allele frequency")+ylab("Probability Density")+
  scale_fill_manual(values=c("red","blue"))+scale_color_manual(values=c("red","blue"))+ theme_classic()+labs(tag = "F.")+ 
  ggtitle("Parents: CA74 x CA80, ALT allele frequency = 0.75")+ theme(plot.title = element_text(hjust = 0.5))

bothfams<- bigfam0.25/ bigfam0.5 / bigfam0.75 / bigfam20.25 / bigfam20.5 / bigfam20.75
mendel_plotname<-paste0("mendel_",plotname,".png")


ggsave(mendel_plotname, bothfams,
       width = 8, height = 12, dpi = 300, units = "in", device='png')
