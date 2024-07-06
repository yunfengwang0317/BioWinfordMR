library(meta)
source('loaddat/Mendelian_Randomization.r')
drawforest<-function(dt,outdir){
                library(forestploter)
                dt=dt[!is.infinite(dt$or),]
                dt$or_uci95[dt$or_uci95>100]=dt$or[dt$or_uci95>100]*1.5
                dt$` ` <- paste(rep(" ",15),collapse=" ")
                dt$`OR (95% CI)` <- ifelse(is.na(dt$or), "",
                           sprintf("%.3f (%.3f to %.3f)",
                                   dt$or, dt$or_lci95, dt$or_uci95))

                tm<-forest_theme(base_size=12,
                        refline_col='red',
                        footnote_col='#636363',
                        footnote_fontface='italic')

                dt$pval=round(dt$pval,3)
                dt[is.na(dt)]=1
                dt[,2]=as.integer(dt[,2])
		dt[grepl('Random|Fixed',dt[,1]),2]=''

		p<-forestploter::forest(
                        dt[,c(1,2,10,11,6)],
                        est=dt$or,lower=dt$or_lci95,upper=dt$or_uci95,ci_column=3,
                        xlim=c(round(.01+min(dt$or_lci95)-abs(min(dt$or_lci95))*.025,2),round(max(dt$or_uci95)+abs(max(dt$or_uci95))*.025,2)),x_trans='log',
                        ticks_at=c(.01+round(min(dt$or_lci95)-abs(min(dt$or_lci95))*.025,2),1,round(max(dt$or_uci95)+abs(max(dt$or_uci95))*.025,2)),
                        theme=tm
                )

                pdf(paste0(outdir,'/MR_forest.pdf'),wi=10,he=max(4,dim(dt)[1]*.6))
                plot(p)
                dev.off()
                return (p)
        }

MR_meta<-function(metafile=NULL,outdir='output/result'){
	system(paste0('rm ',outdir,'/*rds'))
    dat=read.delim(metafile,header=T)
    dat=dat[dat$method=='MR IVW'|dat$method=='Inverse variance weighted',]
    dat$study=paste(dat$id.exposure,dat$id.outcome,sep=' | ')
    names(dat)[grepl('^beta|^b$',names(dat))]='beta'
    names(dat)[grepl('^se',names(dat))]='se'
    names(dat)[grepl('^lo_ci',names(dat))]='lo_ci'
    names(dat)[grepl('^up_ci',names(dat))]='up_ci'
    names(dat)[grepl('^pval',names(dat))]='pval'
    dat$study=shortenname(dat$study,30)
    meta_analysis <- metagen(
	  TE = dat$beta,  # 效应量
	  seTE = dat$se,  # 标准误差
	  studlab = dat$study,  # 研究名称
	  fixed = T  # 是否进行随机效应模型的meta分析
	)
    saveRDS(meta_analysis,file=paste0(outdir,'/meta_analysis.rds'))
    dt=dat[,c('study','nsnp','beta','lo_ci','up_ci','pval')]
    dt=rbind(dt,c('Fixed effect model',0,meta_analysis$TE.fixed,meta_analysis$lower.fixed,meta_analysis$upper.fixed,meta_analysis$pval.fixed))
    dt=rbind(dt,c('Random effect model',0,meta_analysis$TE.random,meta_analysis$lower.random,meta_analysis$upper.random,meta_analysis$pval.random))
    dt=data.frame(dt)
    for (i in 2:6)dt[,i]=as.numeric(dt[,i])
    dt$or=exp(dt$b)
    dt$or_lci95=exp(dt$lo_ci)
    dt$or_uci95=exp(dt$up_ci)
    write.table(dt,paste0(outdir,'/MR_meta_result.txt'),row.names=F,col.names=T,quote=F,sep='\t')
    drawforest(dt,outdir)
}
