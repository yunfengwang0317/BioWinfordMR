library(readxl)
library(gwasrapidd)
library(dplyr)
library(ggpubr)
set.seed(1234)
options(bitmapType='cairo')
if (file.exists('MR_available_outcomes.rda')){
        load('MR_available_outcomes.rda')
}else load('loaddat/MR_available_outcomes.rda')
data(gwas_catalog,package='MRInstruments')
source('loaddat/Mendelian_Randomization.r')
getGWAScatalog<-function(IDs,pvalue,type='exposure'){
        exposure_gwas <- subset(gwas_catalog, STUDY.ACCESSION %in% IDs)
        samplesize=do.call(c,lapply(unique(exposure_gwas$STUDY.ACCESSION),function(id) max(get_studies(study_id =id)@ancestries$number_of_individuals)))
        samplesize=data.frame(STUDY.ACCESSION=unique(exposure_gwas$STUDY.ACCESSION),samplesize=samplesize)
        exposure_gwas=merge(exposure_gwas,samplesize,by='STUDY.ACCESSION')
        exposure=exposure_gwas[,c('SNP','effect_allele','other_allele','beta','se','pval','samplesize','Phenotype','eaf','STUDY.ACCESSION')]
        names(exposure)[c(8,10)]=c(type,'id')
        if (type=='exposure')exposure<-exposure[exposure$pval<pvalue,]
        if (type=='outcome')exposure<-exposure[exposure$pval>pvalue,]
        #exposure<-format_data(exposure)
        return (exposure)
}

maxmin<-function(v)return ((v-min(v))/(max(v)-min(v)))
bubble_plot<-function(res,outdir){
	library(ggrepel)
	library(ggplot2)
	library(RColorBrewer)
	res$neglogP=-log(res$pval)
	labels=res$exposure
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
                  geom_text_repel(aes(b, neglogP, label = labels),size=maxmin(res$neglogP)*3,max.overlaps=20)
	pdf(paste0(outdir,'/bubble_plot.pdf'),wi=8,he=6)
	print (p)
	dev.off()
        return(p)
}

MR_gut<-function(outcome_id,outfile,pvalue,method,output,dir2load,outdir){
	method=head(method,5)
	pvalue=as.numeric(pvalue)
	#load outcome data
	system(paste0('rm ',outdir,'/msg.log'))
	load(paste0('loaddat/',dir2load,'_SNP.rda'))
	exp_tot=subset(exp_tot,pval.exposure<pvalue)
	exposure_SNP=unique(exp_tot$SNP)
	outcome=retrieve_outcome(exp_tot,outcome_id,outfile,outdir=outdir)
	#outcome=outcome[outcome$pval.outcome>pvalue,]#rm later
	if (is.null(outcome)) return ()

	if (max(table(exp_tot$exposure[exp_tot$SNP %in% outcome$SNP]))<3){
		msg=paste0('Few shared SNP were obtained\nTry to modify your Pval')
                cat (msg,file=paste0(outdir,'/msg.log'))
                return ()
	}

	outcome=data.frame(outcome %>% group_by(SNP) %>% filter(pval.outcome %in% min(sort(pval.outcome))) %>% ungroup)
	write.table(outcome,paste0(outdir,'/outcome.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	final=c()
	Fstat_sig=c()
	presso_res=c()
	heter_res=c()
	pleio_res=c()
	direc_res=c()
	tmpdir=gsub('result','tmp',dirname(outdir))
	if (file.exists(paste0(tmpdir,'/mr_result.rda'))) load(paste0(tmpdir,'/mr_result.rda'))
	if (file.exists(paste0(tmpdir,'/Fstat_sig.rda'))) load(paste0(tmpdir,'/Fstat_sig.rda'))
	if (file.exists(paste0(tmpdir,'/presso_res.rda'))) load(paste0(tmpdir,'/presso_res.rda'))
	if (file.exists(paste0(tmpdir,'/heter_res.rda'))) load(paste0(tmpdir,'/heter_res.rda'))
	if (file.exists(paste0(tmpdir,'/pleio_res.rda'))) load(paste0(tmpdir,'/pleio_res.rda'))
	if (file.exists(paste0(tmpdir,'/direc_res.rda'))) load(paste0(tmpdir,'/direc_res.rda'))
	for (id in unique(exp_tot$id.exposure)){
		if (file.exists(paste0(tmpdir,'/MR_',dir2load,'_done'))){
			done=read.table(paste0(tmpdir,'/MR_',dir2load,'_done'),header=F)
			if (id %in% done[,1])next
		}
		#if (as.integer(system(paste0("grep ",id,' ',tmpdir,'/MR_',dir2load,'_done|wc -l'),intern=T))>0) next
		exposure=exp_tot[exp_tot$pval.exposure<pvalue & exp_tot$id.exposure==id,]
		if (dim(exposure)[1]==0)next
		#merging and rm incompatable SNPs
		#exposure$id.exposure=exposure$exposure
	        dat=merge(exposure,outcome,by='SNP')
        	dat=dat[dat$effect_allele.exposure==dat$effect_allele.outcome & dat$other_allele.exposure==dat$other_allele.outcome,]
	        dat=TwoSampleMR::harmonise_data(exposure_dat=exposure,outcome_dat=outcome)
        	dat=unique(dat)
	        if (sum(dat$mr_keep)<3) next
	        dat[,'r.exposure']=TwoSampleMR::get_r_from_bsen(dat$beta.exposure,dat$se.exposure,max(exposure$samplesize.exposure))**2
        	dat[,'r.outcome']=TwoSampleMR::get_r_from_bsen(dat$beta.outcome,dat$se.outcome,max(outcome$samplesize.outcome))**2
		dat$r.exposure[is.na(dat$r.exposure)]=0
		dat$r.outcome[is.na(dat$r.outcome)]=0
	        #write.table(dat,paste0(outdir,'/MR_harmonise.txt'),row.names=F,col.names=T,quote=F,sep='\t')
        	#MR procedure
	        #proxy=paste(unlist(lapply(unique(dat$exposure),function(x) x)),collapse='; ')  
        	#dat$exposure=paste(unlist(lapply(unique(dat$exposure),function(x)strsplit(x,' \\|\\| ')[[1]][1])),collapse='; ')
	        #if (unique(dat$exposure)=='')dat$exposure=proxy
        	mr_res=TwoSampleMR::mr(dat,method_list=method)
		if (dim(mr_res)[1]==0) next
		if (all(is.na(mr_res$pval)))next
		if (is.na(mr_res$pval[mr_res$method=='Inverse variance weighted'])) next
	        mr_result=TwoSampleMR::generate_odds_ratios(mr_res)
		final=rbind(final,mr_result)
		save(final,file=paste0(tmpdir,'/mr_result.rda'))


        	if (mr_res$pval[mr_res$method=='Inverse variance weighted']<0.05) {
       			heter_res=rbind(heter_res,TwoSampleMR::mr_heterogeneity(dat,method_list=c('mr_ivw','mr_egger_regression')))
			pleiotropy=TwoSampleMR::mr_pleiotropy_test(dat)
        		pleio_res=rbind(pleio_res,pleiotropy)
       			direc_res=rbind(direc_res,TwoSampleMR::directionality_test(dat))
			suboutdir=paste0(outdir,'/',id)
			dir.create(paste0(suboutdir))
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
			save(presso_res,file=paste0(tmpdir,'/presso_res.rda'))
			save(heter_res,file=paste0(tmpdir,'/heter_res.rda'))
			save(pleio_res,file=paste0(tmpdir,'/pleio_res.rda'))
			save(direc_res,file=paste0(tmpdir,'/direc_res.rda'))
			save(Fstat_sig,file=paste0(tmpdir,'/Fstat_sig.rda'))
		}
		Fstat<-do.call(rbind,lapply(1:dim(dat)[1],function(idx){
                        N=dat$samplesize.exposure[idx]
                        sd=dat$se.exposure[idx]*sqrt(N)
                        #R2=2*(1-exposure$eaf.exposure[idx])*exposure$eaf.exposure[idx]*exposure$beta.exposure[idx]/sd
                        R2=TwoSampleMR::get_r_from_bsen(dat$beta.exposure[idx],dat$se.exposure[idx],N)**2
                        Fvalue=(N-2)*R2/(1-R2)
                        c(sd,R2,Fvalue)
                }))
		Fstat=cbind(dat,Fstat)
                Fstat_sig=rbind(Fstat_sig,Fstat)
		cat(id,sep='\n',file=paste0(tmpdir,'/MR_',dir2load,'_done'),append=T)
	}
	write.table(presso_res,paste0(outdir,'/MR_presso_result.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	write.table(heter_res,paste0(outdir,'/heterogeneity_result.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	write.table(pleio_res,paste0(outdir,'/pleiotropy_result.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	write.table(direc_res,paste0(outdir,'/steiger_result.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	
        Fstat=data.frame(Fstat_sig)
	NC=dim(Fstat)[2]
	colnames(Fstat)[(NC-2):NC]=c('sd','R2','F')
        write.table(Fstat,paste0(outdir,'/MR_harmonize.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	mr_result=data.frame(unique(final))
	mr_result$qvalue=p.adjust(mr_result$pval,method='fdr')
	write.table(mr_result[order(mr_result$pval),],paste0(outdir,'/MR_result.txt'),row.names=F,col.names=T,quote=F,sep='\t')


	if (sum(final[,9]<0.05 & final$method=='Inverse variance weighted')==0){
                msg=paste0('No significant result obtained.')
                cat (msg,file=paste0(outdir,'/msg.log'))
                return ()
        }



	drawforest<-function(dt){
		dt=dt[!is.infinite(dt$or),]
		dt=subset(dt,or_uci95<100 & or_lci95>0)
		if (dim(dt)[1]<2){
			msg=paste0('Your analysis has finished.')
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
                        dt[,c(1,6,15,16,9)],
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
	mr_result$id.exposure=shortenname(mr_result$id.exposure,30,uni=F)
	plt=bubble_plot(mr_result[mr_result$method=='Inverse variance weighted',],outdir)
	plt=drawforest(mr_result[mr_result$pval<0.05 & mr_result$method=='Inverse variance weighted',names(mr_result)!='qvalue'])
	save(plt,file=paste0(outdir,'/plt.rda'))
	system(paste0("rm ",tmpdir,'/MR_',dir2load,'_done ',tmpdir,'/*.rda'))
#	data<-reshape2::dcast(mr_result[,c('id.exposure','method','pval')],id.exposure~method)
#	CN=names(data)
#	data=data.frame(data[,-1],row.names=data[,1])
#	names(data)=CN[-1]
#	source('loaddat/MR_circle.r')
#	heatmap_circle(data,attr=NULL,showname=T,outdir=outdir)
}

#load('params_gut.rda')
#MRgut(params[['outcome_id']],params[['outcome_file']],as.numeric(params[['pvalue']]),params[['method']],params[['output']],params[['dir2load']],params[['outdir']])
