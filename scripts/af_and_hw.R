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
  #return(freqcomp)
  return(trimmed)
}

freq_1<-freq_func(afreq_1, afreq_2, f1_genocounts)
freq_2<-freq_func(afreq_1_652, afreq_2_652, f1_genocounts_652)
freq_3<-freq_func(read.delim("../foursibs2020_F0_AF.afreq",header=TRUE), read.delim("../foursibs2020_F1_AF.afreq"), read.delim("../foursibs2020_F1.recode.vcf_genocounts.gcount"))
freq_4<-freq_func(read.delim("../foursibs2020_just652_F0_AF.afreq"), read.delim("../foursibs2020_just652_F1_AF.afreq"), read.delim("../foursibs2020_just652_F1.recode.vcf_genocounts.gcount"))

#when the "return" is "trimmed":
mendel_filename1<-paste0("mendel_all_",plotname,".txt")
write.table(freq_1,file=mendel_filename1,quote=FALSE, row.names=FALSE, sep = "\t")

mendel_filename2<-paste0("mendel_652_",plotname,".txt")
write.table(freq_2,file=mendel_filename2,quote=FALSE, row.names=FALSE, sep = "\t")

af_plot_freq1<-ggplot(subset(freq_1, F0freq>0 & F0freq<1), aes(x = F0freq))+
#  geom_violin(aes(fill=as.factor(F0freq),y=F1freq),width=1)+
#  geom_sina(alpha=0.2,aes(y = F1freq,group=as.factor(F0freq)))+
  geom_jitter(alpha=0.2, aes(y = F1freq),width=0.01,height=0.01)+
  theme_bw()
af_plot_freq2<-ggplot(subset(freq_2, F0freq>0 & F0freq<1), aes(x = F0freq))+
#  geom_violin(aes(fill=as.factor(F0freq),y=F1freq),width=1)+
#  geom_sina(alpha=0.2,aes(y = F1freq,group=as.factor(F0freq)))+
  geom_jitter(alpha=0.2, aes(y = F1freq),width=0.02,height=0.02)+
  theme_bw()

hw_plot_freq1<-ggplot(freq_1) +
  geom_jitter(data=freq_1, alpha=0.1,aes(x = F1freq, y = f1_hetfreq),color="green",width=0.01,height=0.01)+
  geom_jitter(data=freq_1, alpha=0.1,aes(F1freq, f1_homreffreq),color="blue",width=0.01,height=0.01)+
  geom_jitter(data=freq_1, alpha=0.1,aes(F1freq, f1_homaltfreq),color="red",width=0.01,height=0.01)+
  xlab("alt allele freq")+ylab("genotype freq")+
  theme_bw()

hw_plot_freq2<-ggplot(freq_2) +
  geom_jitter(data=freq_2, alpha=0.1,aes(x = F1freq, y = f1_hetfreq),color="green",width=0.01,height=0.01)+
  geom_jitter(data=freq_2, alpha=0.1,aes(F1freq, f1_homreffreq),color="blue",width=0.01,height=0.01)+
  geom_jitter(data=freq_2, alpha=0.1,aes(F1freq, f1_homaltfreq),color="red",width=0.01,height=0.01)+
  xlab("alt allele freq")+ylab("genotype freq")+
  theme_bw()

comb<-rbind(freq_1,freq_2)
comb$label<-c(rep("allsnps",nrow(freq_1)),rep("just652",nrow(freq_2)))

bigfam0.5<-ggplot(subset(comb,F0freq==0.5), aes(x=F1freq,fill=label))+
  geom_density(alpha=0.4,adjust=2,bw=0.003)+ylim(0,40)+
  scale_fill_manual(values=c("red","blue"))+scale_color_manual(values=c("red","blue"))+ theme_classic()
bigfam0.25<-ggplot(subset(comb,F0freq==0.25), aes(x=F1freq,fill=label))+
  geom_density(alpha=0.4,adjust=2,bw=0.003)+ylim(0,40)+
  scale_fill_manual(values=c("red","blue"))+scale_color_manual(values=c("red","blue"))+ theme_classic()
bigfam0.75<-ggplot(subset(comb,F0freq==0.75), aes(x=F1freq,fill=label))+
  geom_density(alpha=0.4,adjust=2,bw=0.003)+ylim(0,40)+
  scale_fill_manual(values=c("red","blue"))+scale_color_manual(values=c("red","blue"))+ theme_classic()

#bigfam_mendel<-bigfam0.25 / bigfam0.5 / bigfam0.75
#mendel_plotname<-paste0("mendel_",plotname,".png")

comb2<-rbind(freq_1,freq_2)
comb2$label<-c(rep("allsnps",nrow(freq_1)),rep("just652",nrow(freq_2)))

bigfam20.5<-ggplot(subset(comb2,F0freq==0.5), aes(x=F1freq,fill=label))+
  geom_density(alpha=0.4,adjust=2,bw=0.003)+ylim(0,40)+
  scale_fill_manual(values=c("red","blue"))+scale_color_manual(values=c("red","blue"))+ theme_classic()
bigfam20.25<-ggplot(subset(comb2,F0freq==0.25), aes(x=F1freq,fill=label))+
  geom_density(alpha=0.4,adjust=2,bw=0.003)+ylim(0,40)+
  scale_fill_manual(values=c("red","blue"))+scale_color_manual(values=c("red","blue"))+ theme_classic()
bigfam20.75<-ggplot(subset(comb2,F0freq==0.75), aes(x=F1freq,fill=label))+
  geom_density(alpha=0.4,adjust=2,bw=0.003)+ylim(0,40)+
  scale_fill_manual(values=c("red","blue"))+scale_color_manual(values=c("red","blue"))+ theme_classic()

#bigfam2_mendel<-bigfam20.25 / bigfam20.5 / bigfam20.75
bothfams<- bigfam0.25+labs(tag = "A.") / bigfam0.5 +labs(tag = "B.")/ bigfam0.75 +labs(tag = "C.")/ bigfam20.25+labs(tag = "D.") / bigfam20.5+labs(tag = "E.") / bigfam20.75+labs(tag = "E.")
mendel_plotname<-paste0("mendel_",plotname,".png")


ggsave(bothfams, bigfam_mendel,
       width = 8, height = 12, dpi = 300, units = "in", device='png')

af_combined<- af_plot_freq1 | af_plot_freq2
hw_combined<- hw_plot_freq1 | hw_plot_freq2

af_plotname<-paste0("af_",plotname,".png")
hw_plotname<-paste0("hw_",plotname,".png")
#ggsave(af_plotname,af_combined,
#       width = 8, height = 5, dpi = 300, units = "in", device='png')

#ggsave(hw_plotname,hw_combined,
#       width = 8, height = 5, dpi = 300, units = "in", device='png')
x <- 0:1
dat <- data.frame(x, y=2*x*(1-x))
f <- function(x) y=2*x*(1-x)
dat2<-data.frame(x,y=x^2)
f2<- function(x) x^2
dat3<-data.frame(x, y=(1-x)^2)
f3<-function(x) (1-x)^2
hw_wlines_1<-ggplot(dat, aes(x,y)) + 

  geom_point(data=freq_1, alpha=0.1,aes(x = F1freq, y = f1_hetfreq),color="darkgoldenrod",position=position_jitter(width=0.01, height=0.01),size=2)+
  geom_point(data=freq_1, alpha=0.1,aes(F1freq, f1_homreffreq),color="blue4",position=position_jitter(width=0.01, height=0.01),size=2)+
  geom_point(data=freq_1, alpha=0.1,aes(F1freq, f1_homaltfreq),color="red",position=position_jitter(width=0.01, height=0.01),size=2)+
  stat_function(fun=f, colour="darkgoldenrod", size=0.5)+
  stat_function(fun=f2,colour="red", size=0.5)+
  stat_function(fun=f3, colour="blue4", size=0.5)+
  xlab("alt allele frequency")+ylab("genotype frequency")+
  theme_bw()
  #geom_hline(yintercept=0)+geom_vline(xintercept=0)+
#  scale_y_continuous(expand = c(0,0.00))+scale_x_continuous(expand=c(0,0))+
#  theme(axis.line = element_line(colour = "black"),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        panel.border = element_blank(),
#        panel.background = element_blank()) 

hw_wlines_2<-ggplot(dat, aes(x,y)) + 

  geom_point(data=freq_2, alpha=0.1,aes(x = F1freq, y = f1_hetfreq),color="darkgoldenrod",position=position_jitter(width=0.01, height=0.01),size=2)+
  geom_point(data=freq_2, alpha=0.1,aes(F1freq, f1_homreffreq),color="blue4",position=position_jitter(width=0.01, height=0.01),size=2)+
  geom_point(data=freq_2, alpha=0.1,aes(F1freq, f1_homaltfreq),color="red",position=position_jitter(width=0.01, height=0.01),size=2)+
  stat_function(fun=f, colour="darkgoldenrod", size=0.5)+
  stat_function(fun=f2,colour="red", size=0.5)+
  stat_function(fun=f3, colour="blue4", size=0.5)+
  xlab("alt allele frequency")+ylab("genotype frequency")+
  theme_bw() 
  #geom_hline(yintercept=0)+geom_vline(xintercept=0)+
#  scale_y_continuous(expand = c(0,0.00))+scale_x_continuous(expand=c(0,0))+
#  theme(axis.line = element_line(colour = "black"),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        panel.border = element_blank(),
#        panel.background = element_blank()) 
		
hw_wlines_combined<- hw_wlines_1 | hw_wlines_2

hw_wlines_plotname<- paste0("hw_wlines",plotname,".png")
ggsave(hw_wlines_plotname, hw_wlines_combined,
       width = 8, height = 5, dpi = 300, units = "in", device='png')
