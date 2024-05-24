#https://github.com/lb664/MR2/
library(MR2)
library(reshape2)
library(ggplot2)
source('loaddat/readinput.r')
source('loaddat/Mendelian_Randomization.r')
load('loaddat/MR_available_outcomes.rda')
MVMR2<-function(exp_id='',out_id='',exp_file=NULL,out_file=NULL,pvalue=5e-8,clump_kb=10000,clump_r2=0.001,outdir){
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
        out_id=strsplit(out_id,';')[[1]]
	if (length(exp_id)+length(exp_file)<3){
		msg=paste0('At least 3 exposures are needed')
                cat (msg,file=paste0(outdir,'/msg.log'))
                return ()
	}
	if (length(out_id)+length(out_file)<3){
                msg=paste0('At least 3 outcomes are needed')
                cat (msg,file=paste0(outdir,'/msg.log'))
                return ()
        }


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
                if (N<3) {
			msg=paste0('At least 3 exposures are needed')
                        cat (msg,file=paste0(outdir,'/msg.log'))
			return ()
		}
                tmp<-myApply(1:N,function(fi){
                        readinput(exp_file[fi],paste0(tmpdir,'/',fi),paste0(tmpdir,'/exposure',fi,'.txt'))
                })
                snp=c()
                for (fi in 1:N){
                        exposure <- data.frame(fread(paste0(tmpdir,'/exposure',fi,'.txt')))
                        names(exposure)[1:10]=c('SNP','effect_allele.exposure','other_allele.exposure','beta.exposure','se.exposure','pval.exposure','samplesize.exposure','exposure','eaf.exposure'
,'id.exposure')
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
        
        if (sum(table(exposure$SNP)>1)<3){
                msg=paste0('No overlap SNPs in exposures')
                cat (msg,file=paste0(outdir,'/msg.log'))        
                return ()
        }

        if (length(out_id)>0 & all(do.call(c,lapply(out_id,function(outcome_id) any(grepl(outcome_id,aval_outs)))))){
                out_file=do.call(c,lapply(out_id,function(outcome_id) aval_outs[grepl(outcome_id,aval_outs)][1]))
        }
        if (is.null(out_file)){
                outcome=do.call(rbind,lapply(out_id,function(id) retrieve_outcome(exposure,id,out_file,outdir=outdir)))
                write.table(outcome,paste0(outdir,'/outcome.txt'),row.names=F,col.names=T,quote=F,sep='\t')
        }else if (!is.null(out_file)){
        	if (length(out_file)<3) {
			msg=paste0('At least 3 outcomes are needed')
                        cat (msg,file=paste0(outdir,'/msg.log'))
			return ()
		}
                outcome=do.call(rbind,lapply(out_file,function(fi) {
	               	out=retrieve_outcome(exposure,outcome_id='',outfile=fi,outdir=outdir)
        	       	out=out[,1:10]
               		names(out)=c('SNP','effect_allele.outcome','other_allele.outcome','beta.outcome','se.outcome','pval.outcome','samplesize.outcome','outcome','eaf.outcome','id.outcome')
	               	out
               	}))
                outcome=outcome[!duplicated(paste(outcome$SNP,outcome$outcome)),]
                for (i in 4:6) outcome[,i]=as.numeric(gsub(',','.',outcome[,i]))
                write.table(outcome,paste0(outdir,'/outcome.txt'),row.names=F,col.names=T,quote=F,sep='\t')
        }else return ()



    beta_Y <- reshape(outcome[,c('SNP','id.outcome','beta.outcome')], idvar = "SNP", timevar = "id.outcome", direction = "wide")
    beta_X <- reshape(exposure[,c('SNP','id.exposure','beta.exposure')], idvar = "SNP", timevar = "id.exposure", direction = "wide")
    rownames(beta_Y)=beta_Y$SNP
    rownames(beta_X)=beta_X$SNP
    beta_X=beta_X[,-1]
    beta_Y=beta_Y[,-1]
    MR2_output <- MR2(beta_Y, beta_X, EVgamma = c(1, 1),niter=200,burnin = 100, thin = 10, monitor = 100,nonDecomp = TRUE)
    PostProc_output <- PostProc(MR2_output, beta_Y, beta_X)
    BMM_EM_output <- BMM_EM(PostProc_output$gammaPost)
    FDR_BMM_output <- FDR_BMM(BMM_EM_output$PPI, BMM_EM_output$w)
    cutoff=FDR_BMM_output$cutoff
    PostProc_output$gammaPost > cutoff
    df=PostProc_output$gammaPost
    write.table(df,paste0(outdir,'/GammaPost.txt'),row.names=T,col.names=T,quote=F,sep='\t')
    res=melt(df)
    names(res)=c('exposure','outcome','gamma')
    res$exposure=gsub('beta.exposure.','',res$exposure)
    res$outcome=gsub('beta.outcome.','',res$outcome)
    res$sig <- ifelse(res$gamma >= cutoff,paste0("gamma≥",cutoff),paste0("gamma<",cutoff))

    p=ggplot(res, aes(exposure,outcome,shape = sig,color = gamma)) +
      geom_point(size = 5) +
      scale_shape_manual(values = c(15,7)) +
      scale_color_gradient2(limits = c(-1, 1),low = "#2b8cbe",mid = "white",high = "#e41a1c") +
      theme_bw() + labs(size= 3,x = "exposure",y = "outcome",title = "Multi-response result")+
      theme(plot.title = element_text(hjust = 0.5),
      	axis.text.x = element_text(angle = 45,
        hjust = 0,vjust = 0,family = "sans"))

    pdf(paste0(outdir,'/multi-response.pdf'),wi=nrow(df),he=ncol(df))
    print (p)
    dev.off()
    return (p)
}
