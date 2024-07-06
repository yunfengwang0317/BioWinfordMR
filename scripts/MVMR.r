source('loaddat/readinput.r')
source('loaddat/Mendelian_Randomization.r')
load('loaddat/MR_available_outcomes.rda')
MVMR<-function(exp_id,out_id,exp_file,out_file,pvalue=5e-8,clump_kb=10000,clump_r2=0.001,lasso=F,outdir){
	tmpdir=gsub('result','tmp',outdir)
	if (file.exists(outdir)) system(paste0('rm ',outdir,'/msg.log'))
	dir.create(tmpdir,recursive=T)
	library(rtracklayer)
	library(parallel)
	myApply<-mclapply
	options("mc.cores"=4)
	options(java.parameters = "-Xmx4g" )
	options(java.parameters = c("-XX:+UseConcMarkSweepGC", paste0("-XX:ParallelGCThreads=",1), "-Xmx3072m"))
	load('loaddat/aval_outs.rda')
	exp_id=strsplit(exp_id,';')[[1]]
	
	if (length(exp_id)>0 & all(do.call(c,lapply(exp_id,function(outcome_id) any(grepl(outcome_id,aval_outs)))))){
		exp_file=do.call(c,lapply(exp_id,function(outcome_id) aval_outs[grepl(outcome_id,aval_outs)][1]))
	}
	if (length(exp_id)>1 & is.null(exp_file)) {
		if (!all(exp_id %in% ao$id)){
			msg=paste0('Only ieugwas ID is supported for MVMR')
                        cat (msg,file=paste0(outdir,'/msg.log'))
                        return ()
		}
		tryCatch({
			exposure <- TwoSampleMR::mv_extract_exposures(exp_id,pval_threshold=pvalue,clump_kb=clump_kb,clump_r2=clump_r2)
			tmp=exposure
                },error=function(e){})
                if (!exists('tmp')){
                        msg=paste0('Server code 502; Server is possibly\nexperiencing traffic, trying again...')
                        cat (msg,file=paste0(outdir,'/msg.log'))
                        return ()
                }
	}else if (!is.null(exp_file)){
		N=length(exp_file)
		#tmp<-myApply(1:N,function(fi){
		for (fi in 1:N){				     
			readinput(exp_file[fi],paste0(tmpdir,'/',fi),paste0(tmpdir,'/exposure',fi,'.txt'))
		}
		snp=c()
		for (fi in 1:N){
			exposure <- data.frame(fread(paste0(tmpdir,'/exposure',fi,'.txt')))
			names(exposure)[1:10]=c('SNP','effect_allele.exposure','other_allele.exposure','beta.exposure','se.exposure','pval.exposure','samplesize.exposure','exposure','eaf.exposure','id.exposure')
			for (i in 4:6) exposure[,i]=as.numeric(gsub(',','.',exposure[,i]))
			d=exposure[exposure$pval.exposure<pvalue,1:10]
			for (i in 2:3) d[,i]=toupper(d[,i])
			d=local_clump(d,clump_p=1,clump_kb=clump_kb,clump_r2=clump_r2)
			snp=c(snp,d$SNP)
		}
		write.table(snp,paste0(tmpdir,'/snplst.txt'),row.names=F,col.names=F,quote=F)
		exposure=do.call(rbind,lapply(1:N,function(fi){
			d=fread(paste0(tmpdir,'/exposure',fi,'.txt'))
			names(d)[1:10]=c('SNP','effect_allele.exposure','other_allele.exposure','beta.exposure','se.exposure','pval.exposure','samplesize.exposure','exposure','eaf.exposure','id.exposure')
			d=d[d$SNP %in% snp,1:10]
			d
		}))	
		exposure=data.frame(exposure)
		exposure=exposure[!is.na(exposure$SNP),]
		write.table(exposure,paste0(outdir,'/exposure.txt'),row.names=F,col.names=T,quote=F,sep='\t')
		exposure=read.table(paste0(outdir,'/exposure.txt'),header=T,sep='\t')
	}else return ()
	
	if (sum(table(exposure$SNP)==length(unique(exposure$id.exposure)))<3){
		msg=paste0('No shared SNPs in exposures')
		cat (msg,file=paste0(outdir,'/msg.log'))	
		return ()
	}
	if (nchar(out_id)>0){
		outcome=retrieve_outcome(exposure,out_id,out_file,outdir=outdir)
		write.table(outcome,paste0(outdir,'/outcome.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	}else if (!is.null(out_file)){
		readinput(out_file,outdir,paste0(outdir,'/outcome.txt'))
                outcome <- read.delim(paste0(outdir,'/outcome.txt'),header=T,sep='\t')
		names(outcome)[1:10]=c('SNP','effect_allele.outcome','other_allele.outcome','beta.outcome','se.outcome','pval.outcome','samplesize.outcome','outcome','eaf.outcome','id.outcome')
		outcome=outcome[outcome$SNP %in% exposure$SNP,1:10]
		outcome=outcome[!duplicated(paste(outcome$SNP,outcome$outcome)),]
		for (i in 4:6) outcome[,i]=as.numeric(gsub(',','.',outcome[,i]))
		write.table(outcome,paste0(outdir,'/outcome.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	}else return ()

	dat=TwoSampleMR::mv_harmonise_data(exposure_dat=exposure,outcome_dat=outcome)
	if (dim(dat$exposure_pval)[1]<3){
		msg=paste0('No enough SNPs shared by exposure and outcome')
		cat (msg,file=paste0(outdir,'/msg.log'))	
		return ()
	}
	write.table(exposure,paste0(outdir,'/exposure.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	write.table(outcome,paste0(outdir,'/outcome.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	saveRDS(dat,paste0(outdir,'/harmonize.rds'))
	if (lasso==T) {
		fs=TwoSampleMR::mv_lasso_feature_selection(dat)
		if (dim(fs)[1]==0) fs=data.frame(exposure=dat$expname$exposure[1])
		write.table(fs,paste0(outdir,'/lasso_fs.txt'),row.names=F,col.names=T,quote=F,sep='\t')
		exposure_fs=exposure[exposure$id.exposure %in% fs$exposure,]
		if (length(unique(exposure_fs$id.exposure))>1) {
			dat=TwoSampleMR::mv_harmonise_data(exposure_dat=exposure_fs,outcome_dat=outcome)
		}else {
			dat=TwoSampleMR::harmonise_data(exposure_dat=exposure_fs,outcome_dat=outcome)
			mr_res=TwoSampleMR::mr(dat)
		        mr_result=TwoSampleMR::generate_odds_ratios(mr_res)
		        mr_result$method[mr_result$method=='Inverse variance weighted']='MR IVW'
		        write.table(mr_result,paste0(outdir,'/MR_result.txt'),row.names=F,col.names=T,quote=F,sep='\t')
			return (mr_result)
		}

	}
	tryCatch({
	mr_res=TwoSampleMR::mv_multiple(dat,intercept = T,pval_threshold=pvalue,plots=T)
	},error=function(e){})
	if (!exists('mr_res')) mr_res=TwoSampleMR::mv_multiple(dat,intercept = F,pval_threshold=pvalue,plots=T)
	write.table(mr_res$result,paste0(outdir,'/MV_multiple.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	pdf(paste0(outdir,'/MR_multiple_plot.pdf'),wi=6,he=4)
	print(mr_res$plots)
	dev.off()
	system(paste0('rm ',outdir,'/tmpfile.txt ',outdir,'/msg.log'))	
	mr_res$result$b=round(mr_res$result$b,3)
	mr_res$result$se=round(mr_res$result$se,3)
	mr_res$result$OR=exp(mr_res$result$b)
	mr_res$result$low95CI=exp(mr_res$result$b-1.96*mr_res$result$se)
	mr_res$result$up95CI=exp(mr_res$result$b+1.96*mr_res$result$se)
	return (mr_res$result)
}
