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

  return(freqcomp)
  #return(trimmed)
  
}

hw_wlines_func<-function(freq_1, freq_2) {
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
  return(hw_wlines_combined)
}
freq_1<-freq_func(afreq_1, afreq_2, f1_genocounts)
freq_2<-freq_func(afreq_1_652, afreq_2_652, f1_genocounts_652)
freq_3<-freq_func(read.delim("../F1_just652_AF.afreq",header=TRUE), read.delim("../F1_just652_AF.afreq",header=TRUE), read.delim("../F1_just652.recode.vcf_genocounts.gcount",header=TRUE))
freq_4<-freq_func(read.delim("../F1_just652_AF.afreq",header=TRUE), read.delim("../F1_just652_AF.afreq",header=TRUE), read.delim("../F1_just652.recode.vcf_genocounts.gcount",header=TRUE))
freq_5<-freq_func(read.delim("../F1_just652_AF.afreq",header=TRUE), read.delim("../F1_just652_AF.afreq",header=TRUE), read.delim("../F1_just652.recode.vcf_genocounts.gcount",header=TRUE))
freq_6<-freq_func(read.delim("../F1_just652_AF.afreq",header=TRUE), read.delim("../F1_just652_AF.afreq",header=TRUE), read.delim("../F1_just652.recode.vcf_genocounts.gcount",header=TRUE))
freq_7<-freq_func(read.delim("../F1_just652_AF.afreq",header=TRUE), read.delim("../F1_just652_AF.afreq",header=TRUE), read.delim("../F1_just652.recode.vcf_genocounts.gcount",header=TRUE))
freq_8<-freq_func(read.delim("../F1_just652_AF.afreq",header=TRUE), read.delim("../F1_just652_AF.afreq",header=TRUE), read.delim("../F1_just652.recode.vcf_genocounts.gcount",header=TRUE))

plot_1<-hw_wlines_func(freq_1, freq_2)
plot_2<-hw_wlines_func(freq_3, freq_4)
plot_3<-hw_wlines_func(freq_5, freq_6)
plot_4<-hw_wlines_func(freq_7, freq_8)

hw_wlines_combined<- plot_1 / plot_2 / plot_3 / plot_4

hw_wlines_plotname<- paste0("hw_wlines",plotname,".png")
ggsave(hw_wlines_plotname, hw_wlines_combined,
       width = 8, height = 5, dpi = 300, units = "in", device='png')
