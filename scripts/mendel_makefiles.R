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
#when the "return" is "trimmed":
mendel_filename1<-paste0("mendel_all_",plotname,".txt")
write.table(freq_1,file=mendel_filename1,quote=FALSE, row.names=FALSE, sep = "\t")

mendel_filename2<-paste0("mendel_652_",plotname,".txt")
write.table(freq_2,file=mendel_filename2,quote=FALSE, row.names=FALSE, sep = "\t")

