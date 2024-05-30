set.seed(1234)
source('loaddat/Mendelian_Randomization.r')
bubble_plot<-function(res,outdir){
        library(ggrepel)
        library(ggplot2)
        library(RColorBrewer)
        res$neglogP=-log(res$pval)
        labels=res$outcome
        labels[res$pval>0.05]=''
        p=ggplot(res,aes(x=b,
              y=neglogP,
              size=neglogP,
              fill=b))+
                  geom_point(alpha=0.75,shape=21,color='black')+
                  xlab(paste0('BETA'))+ylab(paste0('-log(pvalue)'))+
                  theme_bw()+ggtitle('')+
                  theme(plot.title = element_text(hjust = .5,size=16),
                       legend.position = "right",
                        panel.grid=element_blank(),
                        axis.title = element_text(size = 16),
                        axis.text = element_text(size = 14))+
                  scale_fill_gradient2(guide='none',
                       low = "#0571b0",
                       mid = "white",
                       high = "#ca0020",
                       midpoint = 0)+
                  scale_size_continuous(range = c(2, 10))+
                  geom_text_repel(aes(b, neglogP, label = labels),size=res$neglogP,max.overlaps=20)
        pdf(paste0(outdir,'/bubble_plot.pdf'),wi=8,he=6)
        print (p)
        dev.off()
        return(p)
}
MR_reverse<-function(exposure_id,expfile,pvalue,method,output,dir2load,outdir){
	tmpdir=gsub('result','tmp',dirname(dirname(outdir)))
        pvalue=as.numeric(pvalue)
        #load outcome data
        system(paste0('rm ',outdir,'/msg.log'))
	exposure=retrieve_exposure(exposure_id,expfile,pvalue,outdir=outdir)
	# R2 and F statistics
        if (all(is.na(exposure$samplesize)))exposure$samplesize=100000
        Fstat<-do.call(rbind,lapply(1:dim(exposure)[1],function(idx){
                N=exposure$samplesize[idx]
                sd=exposure$se[idx]*sqrt(N)
                #R2=2*(1-exposure$eaf.exposure[idx])*exposure$eaf.exposure[idx]*exposure$beta.exposure[idx]/sd
                R2=TwoSampleMR::get_r_from_bsen(exposure$beta[idx],exposure$se[idx],N)**2
                F=(N-2)*R2/(1-R2)
                c(exposure$SNP[idx],sd,R2,F)
        }))
        Fstat=data.frame(Fstat)
        names(Fstat)=c('SNP','sd','R2','F')
	exposure=merge(exposure,Fstat,by='SNP')
        write.table(exposure,paste0(outdir,'/exposure.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	snplist=c('SNP',exposure$SNP[!is.na(exposure$SNP)])
	write.table(snplist,paste0(tmpdir,'/snplist.txt'),row.names=F,col.names=F,quote=F)
	if (is.null(exposure)) return ()
	if (dim(exposure)[1]<3){
		msg=paste0('Not enough exposure SNP obtained')
                cat (msg,file=paste0(outdir,'/msg.log'))
                return ()
	}
	if (file.exists(paste0('/data1/MR_reverse/',dir2load))){
		outcomedir=paste0('/data1/MR_reverse/',dir2load)
	}else if (file.exists(paste0('/Volumes/TOSHIBA_EXT/resource/',dir2load))){
		outcomedir=paste0('/Volumes/TOSHIBA_EXT/resource/',dir2load)
	}else return ()
	final=c()
	presso_res=c()
	heter_res=c()
        direc_res=c()
	pleio_res=c()
	harm_res=c()
	if (file.exists(paste0(tmpdir,'/mr_result_rev.rda'))) load(paste0(tmpdir,'/mr_result_rev.rda'))
        if (file.exists(paste0(tmpdir,'/presso_res_rev.rda'))) load(paste0(tmpdir,'/presso_res_rev.rda'))
        if (file.exists(paste0(tmpdir,'/heter_res_rev.rda'))) load(paste0(tmpdir,'/heter_res_rev.rda'))
        if (file.exists(paste0(tmpdir,'/pleio_res_rev.rda'))) load(paste0(tmpdir,'/pleio_res_rev.rda'))
        if (file.exists(paste0(tmpdir,'/direc_res_rev.rda'))) load(paste0(tmpdir,'/direc_res_rev.rda'))

	library(data.table)
	library(parallel)
	myApply<-mclapply
	options("mc.cores"=8)
	options(java.parameters = "-Xmx4g" )
	options(java.parameters = c("-XX:+UseConcMarkSweepGC", paste0("-XX:ParallelGCThreads=",1), "-Xmx3072m"))
	outfiles=dir(outcomedir)
	outfiles=outfiles[!gsub('.gz','',outfiles) %in% gsub('outcome_','',dir(tmpdir))]
	tmp<-myApply(outfiles,function(outfile){
		system(paste0("loaddat/quickfilter ",tmpdir,'/snplist.txt ',outcomedir,'/',outfile,' ',tmpdir,'/outcome_',gsub('.gz','',outfile)))
	})


	for (outfile in list.files(outcomedir)){
		print (outfile)
		id=gsub('.tsv.gz','',gsub('.txt.gz','',outfile))
		if (file.exists(paste0(tmpdir,'/MR_rev_done'))){
			done=read.table(paste0(tmpdir,'/MR_rev_done'),header=F)
			if (outfile %in% done[,1]) next
		}
		print (paste0(tmpdir,'/outcome_',gsub('.gz','',outfile)))
		outcome=read.delim(paste0(tmpdir,'/outcome_',gsub('.gz','',outfile)),header=T,sep='\t')
		names(outcome)[1:10]=c('SNP','effect_allele','other_allele','beta','se','pval','samplesize','outcome','eaf','id')
		if (names(outcome)[2]!='effect_allele.outcome') names(outcome)[c(2,3,4,5,6,7,9,10)]=paste0(names(outcome)[c(2,3,4,5,6,7,9,10)],'.outcome')
		dat=merge(exposure,outcome,by='SNP')
                dat=TwoSampleMR::harmonise_data(exposure_dat=exposure,outcome_dat=outcome)
                dat=unique(dat)
                if (dim(dat)[1]<3) next
		harm_res=rbind(harm_res,dat)
                mr_res=TwoSampleMR::mr(dat,method_list=method)
                if (dim(mr_res)[1]==0) next
		if (is.na(mr_res$pval[mr_res$method=='Inverse variance weighted'])) next
                mr_result=TwoSampleMR::generate_odds_ratios(mr_res)
                final=rbind(final,mr_result)
		save(final,file=paste0(tmpdir,'/mr_result_rev.rda'))
		if (mr_res$pval[mr_res$method=='Inverse variance weighted']<0.05) {
			heter_res=rbind(heter_res,TwoSampleMR::mr_heterogeneity(dat,method_list=c('mr_ivw','mr_egger_regression')))
                        pleiotropy=TwoSampleMR::mr_pleiotropy_test(dat)
			pleio_res=rbind(pleio_res,pleiotropy)
                        direc_res=rbind(direc_res,TwoSampleMR::directionality_test(dat))
			suboutdir=paste0(outdir,'/',id)
			dir.create(suboutdir)
			write.table(dat,paste0(suboutdir,'/MR_harmonize.txt'),row.names=F,col.names=T,quote=F,sep='\t')
                        #final=rbind(final,mr_result)
                        if ('html' %in% output){
                                tryCatch({
                                        TwoSampleMR::mr_report(dat,output_path=suboutdir,output_type = "html")
                                },error=function(e){return (NULL)})
                        }
                        if ('scatter' %in% output){
                                pdf(paste0(suboutdir,'/',id,'_scatter_plot.pdf'),wi=8,he=6)
                                p=TwoSampleMR::mr_scatter_plot(mr_result,dat)
                                for (i in names(p)) print(p[[i]])
                                dev.off()
                        }
                        if ('forest' %in% output){
                                pdf(paste0(suboutdir,'/',id,'_forest_plot.pdf'),wi=6,he=8)
                                print (TwoSampleMR::mr_forest_plot(TwoSampleMR::mr_singlesnp(dat)))
                                dev.off()
                        }
                        if ('funnel' %in% output){
                                pdf(paste0(suboutdir,'/',id,'_funnel_plot.pdf'),wi=8,he=6)
                                print(TwoSampleMR::mr_funnel_plot(TwoSampleMR::mr_singlesnp(dat)))
                                dev.off()
                        }
                        if ('leaveoneout' %in% output){
                                pdf(paste0(suboutdir,'/',id,'_leaveOneOut_plot.pdf'),wi=6,he=8)
                                print(TwoSampleMR::mr_leaveoneout_plot(TwoSampleMR::mr_leaveoneout(dat)))
                                dev.off()
                        }
			if ('MRpresso' %in% output & dim(dat)[1]>3){
                                presso=MRPRESSO::mr_presso(dat,BetaOutcome='beta.outcome',BetaExposure='beta.exposure',SdOutcome='se.outcome',SdExposure='se.exposure',OUTLIERtest = TRUE, DISTORTIONtest = TRUE,NbDistribution=max(100,dim(dat)[1]+10))
				presso[[1]]$id=id
                                presso[[1]]$RSSobs=presso[[2]][[1]]$RSSobs
                                presso[[1]]$Global.Test.Pvalue=presso[[2]][[1]]$Pvalue
                                print (presso)
                                presso_res=rbind(presso_res,presso[[1]])
				if (!is.na(presso$`Main MR results`$`P-value`[2])) {
                                        outlier=data.frame(cbind(dat$SNP,presso$`MR-PRESSO results`$`Outlier Test`))
                                        outlier$outlier=FALSE
                                        if (presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`=='All SNPs considered as outliers'){ 
                                                outlier$outlier=TRUE
                                        }else outlier$outlier[presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`]=TRUE
                                        write.table(outlier,paste0(suboutdir,'/MR_presso_outlier.txt'),row.names=F,col.names=T,quote=F,sep='\t')
                                }
                        }
			save(presso_res,file=paste0(tmpdir,'/presso_res_rev.rda'))
                        save(heter_res,file=paste0(tmpdir,'/heter_res_rev.rda'))
                        save(pleio_res,file=paste0(tmpdir,'/pleio_res_rev.rda'))
                        save(direc_res,file=paste0(tmpdir,'/direc_res_rev.rda'))
                }
		cat(id,sep='\n',file=paste0(tmpdir,'/MR_rev_done'),append=T)
                #write.table(mr_result,paste0(outdir,'/MR_result.txt'),row.names=F,col.names=T,quote=F,sep='\t',append=T)
        }
	write.table(presso_res,paste0(outdir,'/MR_presso_result.txt'),row.names=F,col.names=T,quote=F,sep='\t')
        write.table(heter_res,paste0(outdir,'/heterogeneity_result.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	write.table(pleio_res,paste0(outdir,'/pleiotropy_result.txt'),row.names=F,col.names=T,quote=F,sep='\t')
        write.table(direc_res,paste0(outdir,'/steiger_result.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	write.table(harm_res,paste0(outdir,'/MR_harmonize.txt'),row.names=F,col.names=T,quote=F,sep='\t')
        mr_result=data.frame(unique(final))
	if (length(mr_result)==0) {
		msg=paste0('Not enough shared SNP obtained')
                cat (msg,file=paste0(outdir,'/msg.log'))
		return ()
	}
	
	mr_result$qvalue=p.adjust(mr_result$pval,method='fdr')
        write.table(mr_result[order(mr_result$pval),],paste0(outdir,'/MR_result.txt'),row.names=F,col.names=T,quote=F,sep='\t')

	if (sum(mr_result$pval<0.05 & mr_result$method=='Inverse variance weighted')==0){
                msg=paste0('No significant result obtained.')
                cat (msg,file=paste0(outdir,'/msg.log'))
		system(paste0("rm ",tmpdir,'/MR_rev_done ',tmpdir,'/*.rda ',tmpdir,'/outcome*'))
                return ()
        }

        drawforest<-function(dt){
                dt=dt[!is.infinite(dt$or),]
                dt=subset(dt,or_uci95<100 & or_lci95>0)
		if (dim(dt)[1]<2){
			msg=paste0('Your analysis has finished')
	                cat (msg,file=paste0(outdir,'/msg.log'))
			return ()
		}
                if (dim(dt)[1]>30)dt=dt[1:30,]
                library(forestploter)
                dt$or_uci95[dt$or_uci95>100]=dt$or[dt$or_uci95>100]*1.5
                dt$` ` <- paste(rep(" ",15),collapse=" ")
                dt$`OR (95% CI)` <- ifelse(is.na(dt$or), "",
                           sprintf("%.3f (%.3f to %.3f)",
                                   dt$or, dt$or_lci95, dt$or_uci95))

                tm<-forest_theme(base_size=12,
                        refline_col='red',
                        footnote_col='#636363',
                        footnote_fontface='italic')

                dt$pval=format(dt$pval, scientific = T,digits=2)
                dt$qvalue=format(dt$qvalue,scientific=T,digits=2)
                dt$method='mr_ivw'
                p<-forestploter::forest(
                        dt[,c(2,6,15,16,9)],
                        est=dt$or,lower=dt$or_lci95,upper=dt$or_uci95,ci_column=3,
                        xlim=c(round(min(dt$or_lci95)-abs(min(dt$or_lci95))*.025,2),round(max(dt$or_uci95)+abs(max(dt$or_uci95))*.025,2)),x_trans='log',
                        ticks_at=c(.01+round(min(dt$or_lci95)-abs(min(dt$or_lci95))*.025,2),1,round(max(dt$or_uci95)+abs(max(dt$or_uci95))*.025,2)),
                        theme=tm
                )
                pdf(paste0(outdir,'/MR_forest.pdf'),wi=10,he=max(4,dim(dt)[1]*.4))
                plot(p)
                dev.off()
                return (p)
        }
        mr_result$id.outcome=shortenname(mr_result$id.outcome,30,uni=F)
	plt=bubble_plot(mr_result[mr_result$method=='Inverse variance weighted',],outdir)
	plt=drawforest(mr_result[mr_result$pval<0.05 & mr_result$method=='Inverse variance weighted',names(mr_result)!='qvalue'])
        save(plt,file=paste0(outdir,'/plt.rda'))
	system(paste0("rm ",tmpdir,'/MR_rev_done ',tmpdir,'/*.rda ',tmpdir,'/outcome*'))
}
