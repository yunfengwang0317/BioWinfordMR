# install.packages("devtools")
#devtools::install_github("NightingaleHealth/ggforestplot")
# Load tidyverse and ggforestplot
# install.packages("tidyverse")
#library(tidyverse)
library(ggforestplot)
source('loaddat/Mendelian_Randomization.r')
drawforest<-function(df,namecol='id.exposure',outdir='./'){
  if (length(unique(df$exposure))==1) names(df)[1:4]=c('id.outcome','id.exposure','exposure','outcome')
  df$name=df[,namecol]
  df$name=shortenname(df$name,30,uni=F)
  poslist=df$id.exposure[df$pval<0.05 & df$method %in% c('IVW','Inverse variance weighted','MR IVW')]
  df=subset(df,id.exposure %in% poslist)
  # Forestplot
  p=forestplot(
    df = df,
    estimate = b,
    logodds = TRUE,
    colour = method,
    shape = method,
    title = unique(df$outcome)[1],
    xlab = "95% Confidential Interval"
  ) 
  # You may also want to add a manual shape scale to mark meta-analysis with a
  # diamond shape
  #p=p+ggplot2::scale_shape_manual(
  #  values = c(23L, 21L, 21L, 21L),
  #  labels = c("Inverse variance weighted", "Simple mode", "MR Egger", "Weighted mode")
  #)
  pdf(paste0(outdir,'/MRforest_',namecol,'.pdf'),wi=8.5,he=max(4,length(unique(df$name))/4))
  print (p)
  dev.off()
  print (p)
}
