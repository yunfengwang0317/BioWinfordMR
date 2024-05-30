library(data.table)
filleaf<-function(dat){
	frq=fread('loaddat/1kg.v3/EUR.frq',header=T)
	dat=merge(dat,frq,by='SNP',all.x=T)
	dat$eaf[is.na(dat$eaf)]=dat$MAF[is.na(dat$eaf)]
	dat$samplesize[is.na(dat$samplesize)]=dat$NCHROBS[is.na(dat$samplesize)]
	return (dat[,1:10])
}
clean_format<-function(inputfile='dirtyfile.txt',IVs=NULL,gwas_type='exposure',snpcol='SNP',EAcol='effect_allele',OAcol='other_allele',betacol='beta',secol='se',pvalcol='pval',Ncol='samplesize',phenocol='exposure',eafcol='eaf',idcol='ID',fill_eaf=F,outdir='output/Mendelian_Randomization'){
	dat=data.frame(fread(inputfile))
	if (idcol!='' & !idcol %in% names(dat)) outdir=paste0(dirname(outdir),'/',idcol)
	dir.create(outdir)

	if (!Ncol %in% names(dat)) {
		Ncol='samplesize'
		dat$samplesize=''
	}
	if (!eafcol %in% names(dat)){
		eafcol='eaf'
		dat$eaf=''
	}
	if (!phenocol %in% names(dat)){
		dat$trait=phenocol
		phenocol='trait'
	}
	if (!idcol %in% names(dat)){
		dat$id=idcol
		idcol='id'
	}
	if (!all(c(snpcol,EAcol,OAcol,betacol,secol,pvalcol,Ncol,phenocol,eafcol,idcol) %in% c('',names(dat)))){
		miscol=paste(c(snpcol,EAcol,OAcol,betacol,secol,pvalcol,Ncol,phenocol,eafcol,idcol)[!c(snpcol,EAcol,OAcol,betacol,secol,pvalcol,Ncol,phenocol,eafcol,idcol)  %in% names(dat)],',')
		dat=data.frame(msg=paste0(miscol,' may be misspelled'))
		return (dat)
	}
	if (any(grepl(',',dat[,pvalcol]))|is.character(dat[1,pvalcol])) dat[,pvalcol]=as.numeric(gsub(',','.',dat[,pvalcol]))
	if (min(dat[,betacol],na.rm=T)>0) dat[,betacol]=log(dat[,betacol])
        if (max(dat[,pvalcol],na.rm=T)>1) dat[,pvalcol]=exp(-dat[,pvalcol])
	if (max(dat[,pvalcol],na.rm=T)<=0) dat[,pvalcol]=exp(dat[,pvalcol])
	if (secol=='') {
		secol='se'
		dat$se=abs(dat[,betacol]/qnorm(dat[,pvalcol]/2))
	}
	dat=dat[,c(snpcol,EAcol,OAcol,betacol,secol,pvalcol,Ncol,phenocol,eafcol,idcol)]
	if (!is.null(IVs)) dat=dat[dat[,snpcol] %in% IVs[IVs!=''] | dat[,pvalcol]<5e-4,]
	if (gwas_type=='exposure' & is.null(IVs)) dat=dat[dat[,pvalcol]<5e-4,]
	if (gwas_type=='exposure') names(dat)[1:10]=c('SNP','effect_allele','other_allele','beta','se','pval','samplesize','exposure','eaf','ID')
	if (gwas_type=='outcome') names(dat)[1:10]=c('SNP','effect_allele','other_allele','beta','se','pval','samplesize','outcome','eaf','ID')
	dat$eaf=as.numeric(dat$eaf)
	dat$samplesize=as.integer(dat$samplesize)
	if (fill_eaf==T) dat=filleaf(dat)
	dat[dat[,8]=='',8]='textfile'
	dat[dat[,10]=='',10]='ID'
	write.table(dat,paste0(outdir,'/',gwas_type,'_clean_format.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	return (head(dat,10000))
}
