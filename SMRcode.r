source('loaddat/Mendelian_Randomization.r')
source("loaddat/plot_SMR.r")
SMR<-function(exp_id='',exp_file=NULL,smr_pval=5e-8,smr_tissues='',smr_mqtl='',smr_pqtl=F,probe='',outdir='./'){
	#system(paste0('rm ',outdir,'/*'))
	tmpdir=gsub('result','tmp',outdir)
	dir.create(tmpdir,recursive=T)
	sysname=Sys.info()["sysname"]
	load('loaddat/GTEx_V8_cis_eqtl_summary_lite/smr_snp.rda')
	if (!file.exists(paste0(tmpdir,'/smr_input.txt'))){
	d=retrieve_exposure(exposure_id=exp_id,expfile=exp_file,pvalue=smr_pval,clump_kb=10000,clump_r2=0.001,clump=F,outdir=tmpdir)
	#d=retrieve_outcome(smr_snp,outcome_id=exp_id,outfile=exp_file,outdir=outdir)
	if (is.null(d))return ()
	d=d[,c(1:3,9,4:7)]
	names(d)=c('SNP','A1','A2','Freq','b','se','p','n')
	d=subset(d,p<smr_pval)
	d=d[!grepl(',',d$SNP),]
	write.table(d,paste0(tmpdir,'/smr_input.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	}else d=read.table(paste0(tmpdir,'/smr_input.txt'),header=T,sep='\t')
	for (smr_tissue in smr_tissues){
		if (smr_tissue=='') next
		if (!any(d$SNP %in% system(paste0("cut -f 2 loaddat/GTEx_V8_cis_eqtl_summary_lite/eQTL_besd_lite/",smr_tissue,".lite.esi"),intern=T))) next
		job=paste0("loaddat/GTEx_V8_cis_eqtl_summary_lite/smr_",sysname," --bfile loaddat/1kg.v3/EUR --gwas-summary ",tmpdir,"/smr_input.txt --beqtl-summary loaddat/GTEx_V8_cis_eqtl_summary_lite/eQTL_besd_lite/",smr_tissue,".lite --out ",outdir,"/smr_",smr_tissue," --heidi-mtd 1 --thread-num 10 --diff-freq-prop 1")
		if (nchar(probe)>0) job=paste0(job,paste0(" --plot --probe-wind 500 --gene-list loaddat/GTEx_V8_cis_eqtl_summary_lite/glist-hg19 --probe ",probe))
		print (job)
		system(job)
	}
	for (mqtl in smr_mqtl){
		if (mqtl=='')next
                if (!any(d$SNP %in% system(paste0("less loaddat/GTEx_V8_cis_eqtl_summary_lite/LBC_BSGS_meta_lite/",mqtl,".esi|awk -F ' ' '{print $2}'"),intern=T))) next
                job=paste0("loaddat/GTEx_V8_cis_eqtl_summary_lite/smr_",sysname," --bfile loaddat/1kg.v3/EUR --gwas-summary ",tmpdir,"/smr_input.txt --beqtl-summary loaddat/GTEx_V8_cis_eqtl_summary_lite/LBC_BSGS_meta_lite/",mqtl," --out ",outdir,"/smr_",mqtl," --heidi-mtd 1 --thread-num 10 --diff-freq-prop 1")
                if (nchar(probe)>0) job=paste0(job,paste0(" --plot --probe-wind 500 --gene-list loaddat/GTEx_V8_cis_eqtl_summary_lite/glist-hg19 --probe ",probe))
                print (job)
                system(job)
        }
	if (smr_pqtl==T){
		job=paste0("loaddat/GTEx_V8_cis_eqtl_summary_lite/smr_",sysname," --bfile loaddat/1kg.v3/EUR --gwas-summary ",tmpdir,"/smr_input.txt --beqtl-summary loaddat/GTEx_V8_cis_eqtl_summary_lite/pQTL_decode/decode --out ",outdir,"/smr_pqtl --heidi-mtd 1 --thread-num 10 --diff-freq-prop 1")
                if (nchar(probe)>0) job=paste0(job,paste0(" --plot --probe-wind 500 --gene-list loaddat/GTEx_V8_cis_eqtl_summary_lite/glist-hg19 --probe ",probe))
                print (job)
                system(job)	
	}
	outputs=list.files(outdir)
	outputs=outputs[grepl('.smr$',outputs)]
	if (length(outputs)==0){
		msg=data.frame(msg='No shared SNP were found')
		write.table(msg,paste0(outdir,'/smr_result.txt'),row.names=F,quote=F)
		return (NULL)
	}
	res=c()
	for (fi in outputs){
		#if (!gsub('smr_|.smr','',fi) %in% smr_tissues) next
		d=read.table(paste0(outdir,'/',fi),header=T,sep='\t')
		if (dim(d)[1]==0) next	
		d$tissue=gsub('smr_','',gsub('.smr','',fi))
		res=rbind(res,d)
	}
	#write.table(res,paste0(outdir,'/smr_result.smr'),row.names=F,col.names=T,quote=F,sep='\t')
	res$FDR=p.adjust(res$p_SMR,method='fdr')
	write.table(res,paste0(outdir,'/smr_result.txt'),row.names=F,col.names=T,quote=F,sep='\t')
}


drawsmr<-function(probe,plotWindow=100,smrpval=5e-6,smr_tissue,outdir){
	probefile=paste0(outdir,'/plot/smr_',smr_tissue,'.',probe,'.txt')
	if (!file.exists(probefile)) {
		plotEmpty('no SNPs in common between \nreference data and GWAS data')
		return (NULL)
	}
	if (as.integer(system(paste0("grep GWAS ",probefile,"|awk -F ' ' '{print $2}'"),intern=T))==0) return ()
	SMRData = ReadSMRData(paste0(outdir,'/plot/smr_',smr_tissue,'.',probe,'.txt'))
	#SMRData[[2]]=SMRData[[2]][!is.na(SMRData[[2]][,4]),]
	SMRData[[2]][is.na(SMRData[[2]][,4]),4:6]=SMRData[[2]][which(is.na(SMRData[[2]][,4]))+1,4:6]	
	# Plot the SMR results in a genomic region centred around a probe:
	pdf(paste0(outdir,'/plot/',probe,'_plot.pdf'),wi=8,he=min(8+length(SMRData$probeID),15))
	SMRLocusPlot(data=SMRData, smr_thresh=smrpval, heidi_thresh=0.05, plotWindow=plotWindow, max_anno_probe=16)
	dev.off()
	SMRLocusPlot(data=SMRData, heidi_thresh=0.05, plotWindow=plotWindow, smr_thresh=smrpval,max_anno_probe=16)
}
