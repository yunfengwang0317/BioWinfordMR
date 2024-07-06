library(boot)
library(ggpubr)
source('loaddat/Mendelian_Randomization.r')
source('loaddat/readinput.r')
source('loaddat/MVMR.r')
calPval<-function(beta,n,se){
    # Calculate t-value
    t_value <- beta / se
    # Calculate degrees of freedom
    df <- n - 1
    # Calculate p-value
    p_value <- 2 * pt(-abs(t_value), df)
    # Print the p-value
    return(p_value)
}

mr_function <- function(data, indices) {  
    d <- data[indices,]   
    jieguo <- TwoSampleMR::mr(d,method_list=c('mr_ivw'))   
    return(jieguo  %>% pull(b)) 
}

drawforest_fix<-function(dt,figname,outdir){
                library(forestploter)
                dt=dt[!is.infinite(dt$or),]
		dt$id.exposure=shortenname(dt$id.exposure,20,uni=F)
		dt$id.outcome=shortenname(dt$id.outcome,20,uni=F)
                dt$or_uci95[dt$or_uci95>100]=dt$or[dt$or_uci95>100]*1.5
                dt$` ` <- paste(rep(" ",15),collapse=" ")
                dt$`OR (95% CI)` <- ifelse(is.na(dt$or), "",
                           sprintf("%.2f (%.2f to %.2f)",
                                   dt$or, dt$or_lci95, dt$or_uci95))

                tm<-forest_theme(base_size=12,
                        refline_col='red',
                        footnote_col='#636363',
                        footnote_fontface='italic')

                dt$pval=round(dt$pval,3)
                p<-forestploter::forest(
                        dt[,c(1,2,15,16,9)],
                        est=dt$or,lower=dt$or_lci95,upper=dt$or_uci95,ci_column=3,
                        xlim=c(round(.01+min(dt$or_lci95)-abs(min(dt$or_lci95))*.025,2),round(max(dt$or_uci95)+abs(max(dt$or_uci95))*.025,2)),x_trans='log',
                        #ticks_at=c(.01+round(min(dt$or_lci95)-abs(min(dt$or_lci95))*.025,2),1,round(max(dt$or_uci95)+abs(max(dt$or_uci95))*.025,2)),
                        theme=tm
                )
                pdf(paste0(outdir,'/',figname),wi=10,he=max(4,dim(dt)[1]*.6))
                plot(p)
                dev.off()
                return (p)
        }

plotEmpty<-function(txtlbl='Please wait until you see the list of samples done!'){
                 plot(c(0,1), c(0,1), type='n',axes=FALSE,ann=F)
                    text(x=0.5, y=0.5, labels = txtlbl,cex=2)
        }

twostep<-function(exposure_id='',exposure_file=NULL,mediator_id='',mediator_file=NULL,outcome_id='',outcome_file=NULL,pval=5e-8,tmpdir,outdir){
    if (file.exists(outdir)) system(paste0('rm ',outdir,'/*/msg.log'))
    exposure_id=gsub(' ','',exposure_id)
    mediator_id=gsub(' ','',mediator_id)
    outcome_id=gsub(' ','',outcome_id)
    #system(paste0('rm -r ',tmpdir))
    dir.create(outdir)
    dir.create(tmpdir)
    snplst=c()
    aval_exps=list.files('loaddat')
    aval_exps=aval_exps[grepl('_SNP.rda',aval_exps)]
    for (aval in aval_exps){
    	load(paste0('loaddat/',aval))
        if (gsub('ebi-a-','',exposure_id) %in% exp_tot$id.exposure & pval<max(exp_tot$pval.exposure)) {
        	print (paste0('retrieve ',exposure_id,' from ',aval))
                d=subset(exp_tot,id.exposure==gsub('ebi-a-','',exposure_id) & pval.exposure<pval)
		#write.table(d,exposure_file,row.names=F,col.names=T,quote=F,sep='\t')
	        #exposure_id=d[1,10]
        	snplst=c(snplst,d$SNP)
           	print ('exposure done')
	}
	if (gsub('ebi-a-','',mediator_id) %in% exp_tot$id.exposure & pval<max(exp_tot$pval.exposure)) {
                print (paste0('retrieve ',mediator_id,' from ',aval))
                d=subset(exp_tot,id.exposure==gsub('ebi-a-','',mediator_id) & pval.exposure<pval)
                #write.table(d,mediator_file,row.names=F,col.names=T,quote=F,sep='\t')
                #mediator_id=d[1,10]
                snplst=c(snplst,d$SNP)
                print ('mediator done')
        }
    }


    if (!is.null(exposure_file)){
	    print ('readinput exposure')
	    readinput(exposure_file,tmpdir,paste0(tmpdir,'/exposure.txt'))
            exposure_file=paste0(tmpdir,'/exposure.txt')
	    d=data.frame(fread(exposure_file,header=T,sep='\t'))
	    d=local_clump(d)
	    write.table(d,exposure_file,row.names=F,col.names=T,quote=F,sep='\t')
	    exposure_id=d[1,10]
	    snplst=c(snplst,d$SNP)
	    print ('exposure done')
    }
    if (!is.null(mediator_file)){
	    print ('readinput mediator')
	    readinput(mediator_file,tmpdir,paste0(tmpdir,'/mediator.txt'),keepallsnp=T)
            mediator_file=paste0(tmpdir,'/mediator.txt')
            d=data.frame(fread(mediator_file,header=T,sep='\t'))
	    d_clump=local_clump(d[d[,grepl('pval',names(d))]<pval,])
	    snplst=unique(c(snplst,d_clump$SNP))
	    d=d[d$SNP %in% snplst,]
	    write.table(d,mediator_file,row.names=F,col.names=T,quote=F,sep='\t')
            mediator_id=d[1,10]
	    print ('mediator done')
    }
    if (!is.null(outcome_file)){
	    print ('readinput outcome')
	    readinput(outcome_file,tmpdir,paste0(tmpdir,'/outcome.txt'),snplst=NULL,keepallsnp=T)
            outcome_file=paste0(tmpdir,'/outcome.txt')
            d=data.frame(fread(outcome_file,header=T,sep='\t'))
            outcome_id=d[1,10]
	    print ('outcome done')
    }
    if (file.exists(paste0(outdir,'/',exposure_id,'.txt'))){
	    MR(exposure_id='',outcome_id=outcome_id,expfile=paste0(outdir,'/',exposure_id,'.txt'),outfile=outcome_file,pvalue=pval,clump_kb=10000,clump_r2=0.001,methods='mr_ivw',runPheno=F,runMRpresso=F,outdir=paste0(outdir,'/',exposure_id,'_',outcome_id))
    }else{
	    MR(exposure_id=exposure_id,outcome_id=outcome_id,expfile=exposure_file,outfile=outcome_file,pvalue=pval,clump_kb=10000,clump_r2=0.001,methods='mr_ivw',runPheno=F,runMRpresso=F,outdir=paste0(outdir,'/',exposure_id,'_',outcome_id))
    }
    if (file.exists(paste0(outdir,'/',exposure_id,'_',outcome_id,'/msg.log'))){
	    msg=paste0(system(paste0('cat ',outdir,'/',exposure_id,'_',outcome_id,'//msg.log'),intern=T),collapse='\n')
	    return (plotEmpty(msg))
    }
    system(paste0('cp ',outdir,'/',exposure_id,'_',outcome_id,'/exposure.txt ',outdir,'/',exposure_id,'.txt'))
    med=mediator_id
    MR(exposure_id='',outcome_id=mediator_id,expfile=paste0(outdir,'/',exposure_id,'.txt'),outfile=mediator_file,pvalue=pval,clump_kb=10000,clump_r2=0.001,methods='mr_ivw',runPheno=F,runMRpresso=F,outdir=paste0(outdir,'/',exposure_id,'_',med))
    if (file.exists(paste0(outdir,'/',exposure_id,'_',med,'/msg.log'))){
            msg=paste0(system(paste0('cat ',outdir,'/',exposure_id,'_',med,'//msg.log'),intern=T),collapse='\n')
            return (plotEmpty(msg))
    }
    MR(exposure_id=med,outcome_id=outcome_id,expfile=mediator_file,outfile=outcome_file,pvalue=pval,clump_kb=10000,clump_r2=0.001,methods='mr_ivw',runPheno=F,runMRpresso=F,outdir=paste0(outdir,'/',med,'_',outcome_id))
	if (file.exists(paste0(outdir,'/',med,'_',outcome_id,'/msg.log'))){
		msg=paste0(system(paste0('cat ',outdir,'/',med,'_',outcome_id,'//msg.log'),intern=T),collapse='\n')
		return (plotEmpty(msg))
    	}


    mr_res=do.call(rbind,lapply(mediator_id,function(med){
        mr1=read.delim(paste0(outdir,'/',exposure_id,'_',med,'/MR_result.txt'),header=T,sep='\t')
	mr1$id.exposure=exposure_id
	mr1$id.outcome=med
        mr2=read.delim(paste0(outdir,'/',med,'_',outcome_id,'/MR_result.txt'),header=T,sep='\t')
	mr2$id.exposure=med
	mr2$id.outcome=outcome_id
        rbind(mr1,mr2)
    }))
    mr_tot=read.delim(paste0(outdir,'/',exposure_id,'_',outcome_id,'/MR_result.txt'),header=T,sep='\t')
    mr_tot$id.exposure=exposure_id
    mr_tot$id.outcome=outcome_id
    mr_res=rbind(mr_res,mr_tot)
    write.table(mr_res,paste0(outdir,'/before_removing_mediator_effect.txt'),row.names=F,col.names=T,quote=F,sep='\t')
    p1=drawforest_fix(mr_res,'before.pdf',outdir)
    for (med in mediator_id){
        beta_0 = mr_res %>% filter(id.exposure==exposure_id & id.outcome==outcome_id) %>% pull(b)
        beta_1 = mr_res %>% filter(id.exposure==exposure_id & id.outcome==med) %>% pull(b)
        beta_2 = mr_res %>% filter(id.exposure==med & id.outcome==outcome_id) %>% pull(b)
        #mr_res$b[mr_res$id.exposure==med & mr_res$id.outcome==outcome_id]=beta_2*beta_1
        mr_res$b[mr_res$id.exposure==exposure_id & mr_res$id.outcome==outcome_id]=beta_0-beta_2*beta_1
    }

    mr_res$or=exp(mr_res$b)
    set.seed(1234)
    sysname=Sys.info()["sysname"]
    for (i in 1:dim(mr_res)[1]){
	if (mr_res$id.exposure[i] != exposure_id | mr_res$id.outcome[i]!=outcome_id) next
        dat=read.delim(paste0(outdir,'/',mr_res$id.exposure[i],'_',mr_res$id.outcome[i],'/MR_harmonise.txt'),header=T,sep='\t')
        if (sysname!='Linux') {
		reps_mr <- summary(boot(data=dat, statistic=mr_function, R=100))
	}else{
		while (!exists('reps')){
	                tryCatch({
        	        reps=boot(data=dat, statistic=mr_function, R=100)
                	},error=function(e){})
                }
		reps_mr=data.frame(original=mr_res$b[i],bootSE=sd(reps$t))
	}
        mr_res$or_lci95[i]=exp(reps_mr$original-2*reps_mr$bootSE)
        mr_res$or_uci95[i]=exp(reps_mr$original+2*reps_mr$bootSE)
	mr_res$pval[i]=calPval(mr_res$b[i],mr_res$nsnp[i],mr_res$se[i])
    }
    write.table(mr_res,paste0(outdir,'/after_removing_mediator_effect.txt'),row.names=F,col.names=T,quote=F,sep='\t')
    p2=drawforest_fix(mr_res,'after.pdf',outdir)
    q=ggarrange(p1,p2,
                        labels = c("before removing mediator effect", "after removing mediator effect"),
                        ncol = 1, nrow = 2)
    pdf(paste0(outdir,'/Mediator_forest.pdf'),wi=8,he=6)
    print (q)
    dev.off()
    print (q)
}

MVMR_mediator<-function(exposure_id='',exposure_file=NULL,mediator_id='',mediator_file=NULL,outcome_id='',outcome_file=NULL,pval=5e-8,tmpdir,outdir){
	res=MVMR(exp_id=paste(exposure_id,mediator_id,collapse=';'),out_id=outcome_id,exp_file=c(exposure_file,mediator_file),out_file=outcome_file,pvalue=pval,clump_kb=10000,clump_r2=0.001,lasso=F,outdir)
	MR(exposure_id=exposure_id,outcome_id=mediator_id,expfile=exposure_file,outfile=mediator_file,pvalue=pval,methods='mr_ivw',runPheno=F,runMRpresso=F,outdir=paste0(outdir,'/exp_med'))
	exp_med=read.table(paste0(outdir,'/exp_med/MR_result.txt'),header=T,sep='\t')
	mv_res=read.table(paste0(outdir,'/MV_multiple.txt'),header=T,sep='\t')
	mr_res=data.frame(rbind(mv_res,exp_med[,names(mv_res)]))
	beta_0 = mr_res %>% filter(id.exposure==exposure_id & id.outcome==outcome_id) %>% pull(b)
	beta_1 = mr_res %>% filter(id.exposure==exposure_id & id.outcome==mediator_id) %>% pull(b)
	beta_2 = mr_res %>% filter(id.exposure==mediator_id & id.outcome==outcome_id) %>% pull(b)
	mr_res$b[mr_res$id.exposure==exposure_id & mr_res$id.outcome==outcome_id]=beta_0-beta_2*beta_1
	return (mr_res)
}
