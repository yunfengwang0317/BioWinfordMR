library(stringr)
library(gwasrapidd)
library(dplyr)
library(ggpubr)
library(MRPRESSO)
library(ieugwasr)
library(plinkbinr)
library(data.table)
set.seed(1234)
options(bitmapType='cairo')
if (file.exists('MR_available_outcomes.rda')){
	load('MR_available_outcomes.rda')
}else load('loaddat/MR_available_outcomes.rda')
data(gwas_catalog,package='MRInstruments')
FinnGen_R8=read.delim('loaddat/finngene_database_R8.txt',header=T,sep='\t')
FinnGen_R9=read.delim('loaddat/finngene_database_R9.txt',header=T,sep='\t')
FinnGen=data.frame(FinnGen_R9)
source('loaddat/readinput.r')
shortenname<-function(xx,L,uni=T){
    res=do.call(c,lapply(xx,function(x){
    x=strsplit(x,'__')[[1]][length(strsplit(x,'__')[[1]])]
    if (nchar(x)>L) {
            paste0(substr(x,1,L),'...')
    }else x
    }))
    if (uni==T) res[res %in% res[duplicated(res)]]=paste0(res[res %in% res[duplicated(res)]],which(res %in% res[duplicated(res)]))
    return (res)
}
local_clump<-function(exp_dat,clump_p=1,clump_kb=10000,clump_r2=0.001){
	names(exp_dat)[1:10]=c('SNP','effect_allele.exposure','other_allele.exposure','beta.exposure','se.exposure','pval.exposure','samplesize.exposure','exposure','eaf.exposure','id.exposure')
	tryCatch({
		exp_clp <-ld_clump_local(dplyr::tibble(rsid=exp_dat$SNP, pval=exp_dat$pval.exposure, 
		id=exp_dat$exposure),clump_p=clump_p,clump_kb=clump_kb,clump_r2=clump_r2,
		bfile="loaddat/1kg.v3/EUR", #https://mrcieu.github.io/ieugwasr/articles/local_ld.html
		#bfile='loaddat//g1000_eur/EUR',
		plink_bin=get_plink_exe())
	},error=function(e){})
	if (!exists('exp_clp')){
	    return (exp_dat)
	}else return (exp_dat[exp_dat$SNP %in% exp_clp$rsid,])
}

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

getFinnGen<-function(IDs,pvalue){
	exposure=c()
	for (ID in IDs){
		filedir=paste0('loaddat/FinnGen_db_R9/',ID,'.txt.gz')
		if (file.exists(filedir)){
			tmp=read.delim(filedir,header=T,sep='\t')
		}else if (file.exists(paste0('loaddat/FinnGen_db_R8/',ID,'.txt.gz'))){
			tmp=read.delim(paste0('loaddat/FinnGen_db_R8/',ID,'.txt.gz'),header=T,sep='\t')
		}else tmp=NULL
		exposure=data.frame(rbind(exposure,tmp[tmp$pval<pvalue,]))
	}
	return (exposure)
}

drawforest<-function(dt,outdir){
                library(forestploter)
                dt=dt[!is.infinite(dt$or),]
		if (dim(dt)[1]==0) return ()
                #dt$or_uci95[dt$or_uci95>100]=dt$or[dt$or_uci95>100]*1.5
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
                p<-forestploter::forest(
                        dt[,c(1,5,15,16,9)],
                        est=dt$or,lower=dt$or_lci95,upper=dt$or_uci95,ci_column=3,
                        #xlim=c(round(.01+min(dt$or_lci95)-abs(min(dt$or_lci95))*.025,2),round(max(dt$or_uci95)+abs(max(dt$or_uci95))*.025,2)),x_trans='log',
                        #ticks_at=c(.01+round(min(dt$or_lci95)-abs(min(dt$or_lci95))*.025,2),1,round(max(dt$or_uci95)+abs(max(dt$or_uci95))*.025,2)),
                        theme=tm
                )
                pdf(paste0(outdir,'/MR_forest.pdf'),wi=10,he=max(4,dim(dt)[1]*.6))
                plot(p)
                dev.off()
                return (p)
        }

check_avail_exposure<-function(exposure_id,pval,outdir){
	exposure_id=strsplit(sub(' ','',exposure_id),';')[[1]]
	aval_exps=list.files('loaddat')
        aval_exps=aval_exps[grepl('_SNP.rda',aval_exps)]
        exposure=NULL
	for (aval in aval_exps){
		if (grepl('original',aval)) next
                load(paste0('loaddat/',aval))
		exp_i=do.call(rbind,lapply(exposure_id,function(id){
                	if (gsub('ebi-a-|ebi-b-','',id) %in% exp_tot$id.exposure) {
                        print (paste0('retrieve ',id,' from ',aval))
                        exp_=subset(exp_tot,id.exposure==gsub('ebi-a-|ebi-b-','',id))
                        exp_=exp_[,c('SNP','effect_allele.exposure','other_allele.exposure','beta.exposure','se.exposure','pval.exposure','samplesize.exposure','exposure','eaf.exposure','id.exposure')]
			exp_[exp_$pval.exposure<pval,]
			}
		}))
		exposure=rbind(exposure,exp_i)
	}
	exposure=unique(exposure)
        write.table(exposure,paste0(outdir,'/exposure.txt'),row.names=F,col.names=T,quote=F,sep='\t')
        return (exposure)
}


check_avail_outcome<-function(outcome_id){
	#load outcome data
        load('loaddat/aval_outs.rda')
        outcome_id=gsub(' |ebi-a-|ebi-b-','',outcome_id)
        if (any(grepl(outcome_id,aval_outs))){
                print ('extract outcome from availabel files')
                outfile=aval_outs[grepl(outcome_id,aval_outs)][1]
                return (outfile)
        }
	return (NULL)
}

retrieve_exposure<-function(exposure_id='',expfile=NULL,pvalue=5e-08,clump_kb=10000,clump_r2=0.001,clump=T,outdir='./'){
	pvalue=as.numeric(pvalue)
	tmpdir=gsub('result','tmp',outdir)
	dir.create(tmpdir,recursive=T)
	if (exposure_id!='' & clump==T){
		exposure=check_avail_exposure(exposure_id,pvalue,outdir)
		if (!is.null(exposure)) return (exposure)
	}
	if (exposure_id!='' & clump==F){
		rawfi=check_avail_outcome(exposure_id)
		if (!is.null(rawfi)) expfile=rawfi
	}


	#load exposure data
	if (is.null(expfile) & exposure_id!='') {
		exposure_id=strsplit(gsub(' ','',exposure_id),';')[[1]]
		if (length(exposure_id)>5) {
			msg=paste0('exposure should be no more than 5')
			cat (msg,file=paste0(outdir,'/msg.log'))
			return ()
		}
		ieu_id=unique(ao$id[tolower(ao$id) %in% tolower(exposure_id)])
		catalog_id=unique(gwas_catalog$STUDY.ACCESSION[tolower(gwas_catalog$STUDY.ACCESSION) %in% tolower(exposure_id)])
		FG_id=FinnGen$phenocode[tolower(FinnGen$phenocode) %in% tolower(exposure_id)]
		exposure=c()
		tries=0
		#load from open GWAS
		if (length(ieu_id)>0){
			print (paste0(length(ieu_id),' IDs from opengwas'))
			tryCatch({
				#exposure=TwoSampleMR::extract_instruments(outcomes=ieu_id,p1=min(5e-4,pvalue),clump=clump)
				obj=list()
				obj$ieu_id=ieu_id
				obj$snplst=NULL
				obj$pvalue=pvalue
				obj$r2=clump_r2
				obj$kb=clump_kb
				obj$clump=clump
				obj$outdir=outdir
				objpath=paste0(tmpdir,'/params_exposure.rds')
				saveRDS(obj,file=objpath)
				job=read.table('MRjob_lst',header=T,sep='\t')
				job=rbind(job[job$obj!=objpath,],c(objpath,'wait',as.character(Sys.time())))
				names(job)=c('obj','status','time')
				write.table(job,'MRjob_lst',row.names=F,col.names=T,quote=F,sep='\t')
				system(paste0('rm ',outdir,'/exposure.txt'))
				t0=Sys.time()
				while (as.numeric(difftime(Sys.time(), t0, units = "secs"))<300){
					if (file.exists(paste0(outdir,'/exposure.txt'))){
						exposure=read.table(paste0(outdir,'/exposure.txt'),header=T,sep='\t')
						break
					}else Sys.sleep(10)
				}
				tmp=exposure
			},error=function(e){})
			if (!exists('tmp')){
				msg=paste0('Server code 502; Server is possibly\nexperiencing traffic, trying again...')
                                cat (msg,file=paste0(outdir,'/msg.log'))
                                return ()
			}
			if (is.null(exposure)){
                                msg=paste0('Not enough exposure SNP obtained.\nPlease modify the Pvalue.')
                                cat (msg,file=paste0(outdir,'/msg.log'))
                                return ()
                        }
			if (is.null(exposure$samplesize.exposure)){
				exposure$samplesize.exposure=100000
			}
			
			if (all(is.na(exposure$samplesize.exposure))){
				exposure$samplesize.exposure=100000
			}
			
			exposure$samplesize.exposure=max(exposure$samplesize.exposure,na.rm=T)
			exposure=exposure[,c('SNP','effect_allele.exposure','other_allele.exposure','beta.exposure','se.exposure','pval.exposure','samplesize.exposure','exposure','eaf.exposure','id.exposure')]
			names(exposure)=c('SNP','effect_allele','other_allele','beta','se','pval','samplesize','exposure','eaf','id')
		}
		#load from GWAS catalog
		if (length(catalog_id)>0){
			print (paste0(length(catalog_id),' IDs from gwas_catalog'))
			exposure=data.frame(rbind(exposure,getGWAScatalog(catalog_id,pvalue,'exposure')))
		}

		#load from FinnGen
		if (length(FG_id)>0){
			print (paste0(length(FG_id),' IDs from FinnGen'))
			exposure=data.frame(rbind(exposure,getFinnGen(FG_id,pvalue)))
		}
		exposure=unique(exposure)	

		if (length(ieu_id)+length(catalog_id)+length(FG_id)>0){
			if (dim(exposure)[1]<2){
				msg=paste0('Not enough exposure SNP obtained.\nPlease modify the Pvalue.')
				cat (msg,file=paste0(outdir,'/msg.log'))
				return ()
			}
		}else if (length(ieu_id)+length(catalog_id)+length(FG_id)==0 & length(exposure_id)>0){
			msg=paste0('Your exposure_id is not available.\nPlease check it.')
			cat (msg,file=paste0(outdir,'/msg.log'))
			return ()
		}
	}else{
		tmpdir=gsub('result','tmp',dirname(outdir))
		readinput(expfile,tmpdir,paste0(outdir,'/exposure.txt'))	
		d=fread(paste0(outdir,'/exposure.txt'),header=T,sep='\t')
		d=data.frame(d)
		names(d)[1:10]=c('SNP','effect_allele','other_allele','beta','se','pval','samplesize','exposure','eaf','ID')
		if (length(unique(d$ID))!=length(unique(d$exposure))) d$ID=d$exposure
		d=d[d$SNP!='',]
		if (any(is.na(d$samplesize))) d$samplesize[is.na(d$samplesize)]=100000
		d=subset(d,pval<pvalue)
		if (length(d$SNP)<2){
        	        msg=paste0('Not enough exposure SNP obtained.\nPlease modify the Pvalue')
                	cat (msg,file=paste0(outdir,'/msg.log'))
	                return ()
	        }

		if (any(grepl('eaf',tolower(names(d))))) names(d)[grepl('eaf',tolower(names(d)))]='eaf'
		write.table(d[d$pval<pvalue,],paste0(outdir,'/exposure.txt'),row.names=F,col.names=T,quote=F,sep='\t')
		#read exposure data
		exposure=TwoSampleMR::read_exposure_data(filename=paste0(outdir,'/exposure.txt'),
        	   clump = FALSE,
                   sep = "\t",
                   snp_col = "SNP",
                   beta_col = "beta",
                   se_col = "se",
                   eaf_col = "eaf",
                   effect_allele_col = "effect_allele",
                   other_allele_col = "other_allele",phenotype_col='exposure',
                   pval_col = "pval",samplesize_col = "samplesize",id_col='ID')
		if (!'id.exposure' %in% names(exposure))exposure$id.exposure=exposure$exposure
		exposure$samplesize.exposure=max(exposure$samplesize.exposure,na.rm=T)
	}
	#rm na lines and modify the header
	#exposure=na.omit(exposure)
	exposure=exposure[!is.na(exposure$effect_allele) & !is.na(exposure$other_allele),]
	if (names(exposure)[2]!='effect_allele.exposure') names(exposure)[c(2,3,4,5,6,7,9,10)]=paste0(names(exposure)[c(2,3,4,5,6,7,9,10)],'.exposure')
	#exposure=data.frame(exposure %>% group_by(SNP) %>% filter(pval.exposure %in% min(sort(pval.exposure))) %>% ungroup)
	exposure=subset(exposure,pval.exposure<pvalue)
	if (clump==T) {
		exposure_clump=c()
		if (length(unique(exposure$id.exposure))>30){
			msg=paste0('Your exposure contains too many gwas')
                        cat (msg,file=paste0(outdir,'/msg.log'))
                        return ()
		}
		for (exp_id in unique(exposure$id.exposure)){
			exposure_clump=rbind(exposure_clump,local_clump(exposure[exposure$id.exposure==exp_id,],1,clump_kb,clump_r2))
		}
		exposure=data.frame(exposure_clump)
	}
	write.table(exposure,paste0(outdir,'/exposure.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	return (exposure)
}


retrieve_outcome<-function(exposure,outcome_id='',outfile=NULL,outdir='./'){
	tmpdir=gsub('result','tmp',outdir)
	dir.create(tmpdir,recursive=T)
	if (length(exposure$SNP)<2){
		msg=paste0('Not enough outcome SNP obtained.\nPlease modify the Pvalue')
		cat (msg,file=paste0(outdir,'/msg.log'))
		return ()
	}
	#load outcome data
	if (is.null(outfile) & outcome_id!=''){
        	outfile=check_avail_outcome(outcome_id)
	}
        

	if (is.null(outfile) & outcome_id!=''){
		outcome_id=strsplit(gsub(' ','',outcome_id),';')[[1]]
		ieu_id=unique(ao$id[tolower(ao$id) %in% tolower(outcome_id)])
		outcome=c()
		tries=0
		#load from open GWAS
		if (length(ieu_id)>0){
			print (paste0(length(ieu_id),' IDs from opengwas as outcome'))
			tryCatch({
				#outcome=TwoSampleMR::extract_outcome_data(snps=unique(exposure$SNP),outcomes=ieu_id)
				obj=list()
                                obj$ieu_id=ieu_id
                                obj$snplst=exposure$SNP
                                obj$outdir=outdir
                                objpath=paste0(tmpdir,'/params_outcome.rds')
                                saveRDS(obj,file=objpath)
                                job=read.table('MRjob_lst',header=T,sep='\t')
                                job=rbind(job[job$obj!=objpath,],c(objpath,'wait',as.character(Sys.time())))
                                names(job)=c('obj','status','time')
				print (job)
                                write.table(job,'MRjob_lst',row.names=F,col.names=T,quote=F,sep='\t')
                                system(paste0('rm ',outdir,'/outcome.txt'))
				t0=Sys.time()
                                while (as.numeric(difftime(Sys.time(), t0, units = "secs"))<300){
                                        if (file.exists(paste0(outdir,'/outcome.txt'))){
                                                outcome=read.table(paste0(outdir,'/outcome.txt'),header=T,sep='\t')
						break
                                        }else Sys.sleep(10)
                                }
				tmp1=outcome
			},error=function(e){})
			if (!exists('tmp1')){
                                msg=paste0('Server code 502; Server is possibly\nexperiencing traffic, trying again...')
                                cat (msg,file=paste0(outdir,'/msg.log'))
                                return ()
                        }

			if (is.null(outcome)){
                                msg=paste0('Not enough outcome SNP obtained.\nPlease modify the Pvalue.')
                                cat (msg,file=paste0(outdir,'/msg.log'))
                                return ()
                        }
        	        if (is.null(outcome$samplesize.outcome)){
                                outcome$samplesize.outcome=100000
                        }
                        
                        if (all(is.na(outcome$samplesize.outcome))){
                                outcome$samplesize.outcome=100000
                        }
                        
			outcome=outcome[,c('SNP','effect_allele.outcome','other_allele.outcome','beta.outcome','se.outcome','pval.outcome','samplesize.outcome','outcome','eaf.outcome','id.outcome')]
                	names(outcome)=c('SNP','effect_allele','other_allele','beta','se','pval','samplesize','outcome','eaf','id')
		}

		if (length(ieu_id)>0){
                	if (dim(outcome)[1]<2){
                        	msg=paste0('Not enough outcome SNP obtained.\nPlease modify the Pvalue.')
	                        cat (msg,file=paste0(outdir,'/msg.log'))
        	                return ()
                	}	
		}else if (length(ieu_id)==0 & length(outcome_id)>0){
			msg=paste0('Only ieugwas ID \ncan be used as outcome')
			cat (msg,file=paste0(outdir,'/msg.log'))
			return ()
		}
	}else if (!is.null(outfile)){
		tmpdir=gsub('result','tmp',dirname(outdir))
		readinput(outfile,tmpdir,paste0(outdir,'/outcome.txt'),snplst=exposure$SNP)
		if (!file.exists(paste0(outdir,'/outcome.txt'))){
			msg=paste0('Not enough shared SNP obtained.\nPlease modify the Pvalue.')
                        cat (msg,file=paste0(outdir,'/msg.log'))
                        return ()
		}
		d=fread(paste0(outdir,'/outcome.txt'),header=T,sep='\t')
		d=data.frame(d)
		names(d)[1:10]=c('SNP','effect_allele','other_allele','beta','se','pval','samplesize','outcome','eaf','id.outcome')
		if (length(unique(d$outcome))!=length(unique(d$id.outcome))) d$id.outcome=d$outcome
		write.table(d[d$SNP %in% exposure$SNP,],paste0(outdir,'/outcome.txt'),row.names=F,col.names=T,quote=F,sep='\t')
		if (length(intersect(exposure$SNP,d$SNP))==0){
			msg=paste0('No overlap SNPs were found.\nPlease check it.')
                        cat (msg,file=paste0(outdir,'/msg.log'))
                        return ()
		}
		#read outcome data
		outcome=TwoSampleMR::read_outcome_data(snps=exposure$SNP,filename=paste0(outdir,'/outcome.txt'),
                   sep = "\t",
                   snp_col = "SNP",
                   beta_col = "beta",
                   se_col = "se",
                   eaf_col = "eaf",
                   effect_allele_col = "effect_allele",
                   other_allele_col = "other_allele",phenotype_col='outcome',
                   pval_col = "pval",id_col='id.outcome')
		if (!'id.outcome' %in% names(outcome)) outcome$id.outcome=outcome$outcome
		if (is.null(outcome$samplesize.outcome)){
			outcome$samplesize.outcome=100000
		}

		if (all(is.na(outcome$samplesize.outcome))){
			outcome$samplesize.outcome=100000
		}
	}else return ()
	outcome=outcome[!is.na(outcome$effect_allele) & !is.na(outcome$other_allele),]
        if (names(outcome)[2]!='effect_allele.outcome') names(outcome)[c(2,3,4,5,6,7,9,10)]=paste0(names(outcome)[c(2,3,4,5,6,7,9,10)],'.outcome')
	outcome=data.frame(outcome %>% group_by(SNP) %>% filter(pval.outcome %in% min(sort(pval.outcome))) %>% ungroup)
        write.table(outcome,paste0(outdir,'/outcome.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	return (outcome)
}

MR<-function(exposure_id='',outcome_id='',expfile=NULL,outfile=NULL,pvalue=5e-08,clump_kb=10000,clump_r2=0.001,methods,runPheno=F,runMRpresso=F,outdir){
        cat (exposure_id,outcome_id,expfile,outfile,pvalue,clump_kb,clump_r2,runPheno,runMRpresso,outdir,'\n')
        system(paste0('rm -r ',outdir))
        if (!file.exists(outdir)) dir.create(outdir)
        #ao <- available_outcomes()
        outcomes=unique(ao$id[tolower(ao$id) %in% tolower(strsplit(gsub(' ','',outcome_id),';')[1])])

	exposure=retrieve_exposure(exposure_id,expfile,pvalue,clump_kb,clump_r2,outdir=outdir)
	if (is.null(exposure)) return ()
	outcome=retrieve_outcome(exposure,outcome_id,outfile,outdir=outdir)
	if (is.null(outcome)) {
		msg=paste0('Not enough shared SNP obtained')
		cat (msg,file=paste0(outdir,'/msg.log'))
		return ()
	}
	if (dim(outcome)[1]<3) {
		msg=paste0('Not enough shared SNP obtained.')
                cat (msg,file=paste0(outdir,'/msg.log'))
		return ()
	}
	if (length(unique(exposure$id.exposure))>30|length(unique(outcome$id.outcome))>30){
		msg=paste0('Your data contain too many gwas')
		cat (msg,file=paste0(outdir,'/msg.log'))
		return ()
	}
	# R2 and F statistics
	if (all(is.na(exposure$samplesize)))exposure$samplesize=100000
	if (all(is.na(outcome$samplesize)))outcome$samplesize=100000
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
	write.table(Fstat,paste0(outdir,'/SNP_F_statistic.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	#if (length(exposure_id)>1 & sum(table(exposure$SNP)>1)>1){
        #        print ('Running MR multiple')
	#	tryCatch({
        #        dat=TwoSampleMR::mv_harmonise_data(exposure_dat=exposure,outcome_dat=outcome)
        #        mr_res=TwoSampleMR::mv_multiple(dat)
        #        write.table(mr_res$result,paste0(outdir,'/MV_multiple.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	#	},error=function(e){cat('MR multiple oncur error')})
        #}

	if (dim(exposure)[1]>20) exposure=exposure[exposure$pval.exposure<pvalue,]
	#merging and rm incompatable SNPs
	dat=merge(exposure,outcome,by='SNP')
	dat=dat[dat$effect_allele.exposure==dat$effect_allele.outcome & dat$other_allele.exposure==dat$other_allele.outcome,]
	dat=TwoSampleMR::harmonise_data(exposure_dat=exposure,outcome_dat=outcome)
	dat=unique(dat)
	if (dim(dat)[1]<3){
                msg=paste0('Not enough exposure SNP obtained.\nPlease modify the Pvalue.')
                cat (msg,file=paste0(outdir,'/msg.log'))
                return ()
        }

	dat[,'r.exposure']=TwoSampleMR::get_r_from_bsen(dat$beta.exposure,dat$se.exposure,max(exposure$samplesize.exposure))**2
	dat[,'r.outcome']=TwoSampleMR::get_r_from_bsen(dat$beta.outcome,dat$se.outcome,max(outcome$samplesize.outcome))**2
	write.table(dat,paste0(outdir,'/MR_harmonise.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	if (sum(!dat$remove)<3){
		msg=paste0('Few SNP pass Harmonize\nYou can download and check details')
                cat (msg,file=paste0(outdir,'/msg.log'))
		return ()
	}
	tryCatch({
		TwoSampleMR::mr_report(dat,output_path=outdir,output_type = "html")	
		},error=function(e){cat('error ocurred')})
	if (runPheno==T){
		#FastDownloader::install_pkg("FastTraitR")
		library(FastTraitR) #引用R包
		mr_scanner =look_trait(file_name=paste0(outdir,'/MR_harmonise.txt'), pval=pvalue)
		#mr_scanner=MendelianRandomization::phenoscanner(snpquery=dat$SNP[i:min(dim(dat)[1],i+90)])
		write.table(mr_scanner,paste0(outdir,'/PhenoScanner_result.txt'),row.names=F,col.names=T,quote=F,sep='\t',append=T)
	}
	#MR procedure
	#proxy=paste(unlist(lapply(unique(dat$exposure),function(x) x)),collapse='; ')	
	#dat$exposure=paste(unlist(lapply(unique(dat$exposure),function(x)strsplit(x,' \\|\\| ')[[1]][1])),collapse='; ')
	#if (unique(dat$exposure)=='')dat$exposure=proxy
	mr_res=TwoSampleMR::mr(dat,method_list=methods)
	mr_result=TwoSampleMR::generate_odds_ratios(mr_res)
	mr_result$method[mr_result$method=='Inverse variance weighted']='MR IVW'
	write.table(mr_result,paste0(outdir,'/MR_result.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	p2=drawforest(mr_result,outdir)
	di_test=do.call(rbind,lapply(1:dim(dat)[1],function(i) TwoSampleMR::directionality_test(dat[i,])))
	newdat=data.frame(cbind(dat,di_test[,5:8]))
	write.table(newdat,paste0(outdir,'/MR_harmonise.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	if (length(unique(mr_result$id.exposure))>1 | length(unique(mr_result$id.outcome))>1) {
		msg=paste0('graph failed due to multiple exposure\nYou can download the other results')
		print (msg)
		plt=p2
		save(plt,file=paste0(outdir,'/plt.rda'))
		return ()	
	}
	pdf(paste0(outdir,'/MR_scatter_plot.pdf'),wi=8,he=6)
	p=TwoSampleMR::mr_scatter_plot(mr_result,dat)
	for (i in names(p)) print(p[[i]])
	dev.off()
	
	#sensitivity test
	tryCatch({
	heter_res=TwoSampleMR::mr_heterogeneity(dat,method_list=c('mr_ivw','mr_egger_regression'))
	print (heter_res)
	write.table(heter_res,paste0(outdir,'/MR_heterogeneity.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	},error=function(e){cat(paste("method causes Error\n"))})

	#run_mr_presso(dat,NbDistribution=3000)
	pdf(paste0(outdir,'/MR_funnel_plot.pdf'),wi=8,he=6)
	print(TwoSampleMR::mr_funnel_plot(TwoSampleMR::mr_singlesnp(dat)))
	dev.off()

	pleiotropy=TwoSampleMR::mr_pleiotropy_test(dat)
	print (pleiotropy)
	write.table(pleiotropy,paste0(outdir,'/MR_pleiotropy_test.txt'),row.names=F,col.names=T,quote=F,sep='\t')

	if (runMRpresso==T & dim(dat)[1]>3){
		presso=mr_presso(dat,BetaOutcome='beta.outcome',BetaExposure='beta.exposure',SdOutcome='se.outcome',SdExposure='se.exposure',NbDistribution=max(100,dim(dat)[1]+10),OUTLIERtest=T,DISTORTIONtest=T)
		presso_res=presso$`Main MR results`
		presso_res$Global_Test_Pvalue=presso$`MR-PRESSO results`$`Global Test`$Pvalue
		write.table(presso_res,paste0(outdir,'/MR_presso_test.txt'),row.names=F,col.names=T,quote=F,sep='\t')
		if (!is.na(presso$`Main MR results`$`P-value`[2])) {
			outlier=data.frame(cbind(dat$SNP,presso$`MR-PRESSO results`$`Outlier Test`))
			outlier$outlier=FALSE
			if (presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`=='All SNPs considered as outliers'){ 
				outlier$outlier=TRUE
			}else outlier$outlier[presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`]=TRUE
			write.table(outlier,paste0(outdir,'/MR_presso_outlier.txt'),row.names=F,col.names=T,quote=F,sep='\t')
		}
	}
	
	pdf(paste0(outdir,'/MR_forest_plot.pdf'),wi=6,he=min(20,max(3,log(dim(dat)[1],2))))
	print (TwoSampleMR::mr_forest_plot(TwoSampleMR::mr_singlesnp(dat)))
	dev.off()

	pdf(paste0(outdir,'/MR_leaveOneOut_plot.pdf'),wi=6,he=min(20,max(3,log(dim(dat)[1],2))))
	print(TwoSampleMR::mr_leaveoneout_plot(TwoSampleMR::mr_leaveoneout(dat)))
	dev.off()
	
	p1=TwoSampleMR::mr_scatter_plot(mr_result,dat)
	plt=ggarrange(p1[[1]],p2, 
                                labels = c("A", "B"),
                                ncol = 1, nrow = 2)
	save(plt,file=paste0(outdir,'/plt.rda'))
}
#load('params.rda')
#MR(params[['exposure_id']],params[['outcome_id']],params[['exposure_file']],params[['outcome_file']],as.numeric(params[['pvalue']]),as.numeric(params[['clump_kb']]),as.numeric(params[['clump_r2']]),params[['methods']],params[['runPheno']],params[['runMRpresso']],params[['outdir']])
