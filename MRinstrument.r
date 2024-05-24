library(MRInstruments)
load('loaddat/MR_available_outcomes.rda')
load('loaddat/pQTL.rda')
load('loaddat/metab_qtls.rda')
MRinstru<-function(sourcename){
	if (sourcename=='ieu open gwas') {
		d=ao
	}else if (sourcename=='finngen_R8'){
		d=read.delim('loaddat/finngene_database_R8.txt',header=T,sep='\t')
	}else if (sourcename=='finngen_R9'){
		d=read.delim('loaddat/finngene_database_R9.txt',header=T,sep='\t')
	}else if (sourcename=='finngen_R10'){
                d=read.delim('loaddat/finngene_database_R10.txt',header=T,sep='\t')
	}else if (sourcename=='proteomic_qtls'){
		d=pQTL
	}else if (sourcename=='metab_qtls'){
		d=metab_qtls
	}else {
		d=eval(parse(text=sourcename))
	}
	return (d)
}
