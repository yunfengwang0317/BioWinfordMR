alignseq<-function(fa,GFFmap,tmpdir,outdir){
	if (as.integer(system(paste0("head -1 ",fa,"|grep '>'|wc -l"),intern=T))==1) {
		system(paste0("cp ",fa,' ',tmpdir,'/input.fa'))#fa format
        }else{
                print ('TXT 2 FA')
                seqs=read.delim(fa,header=F,sep='\t')
                for (i in 1:dim(seqs)[1]){
                	s=seqs[i,1]
                	ss=seqs[i,2]
                        cat (paste0('>',s),file=paste0(tmpdir,'/input.fa'),sep='\n',append=T)
                        cat (ss,file=paste0(tmpdir,'/input.fa'),sep='\n',append=T)
                }
        }
        system(paste0("sh loaddat/singleline.sh -m many -i ",tmpdir,'/input.fa'," -o ",tmpdir,"/input_oneline.fa"))
        system(paste0("gsnap -t 4 -A sam -N 1 -D ./loaddat/reference/ -d ",GFFmap," -w 50000 ",tmpdir,"/input_oneline.fa 2>/dev/null |grep -v @ |cut -f 1,3,4 > ",tmpdir,"/mapping_result.txt"))
	d=read.delim(paste0(tmpdir,'/mapping_result.txt'),sep='\t',header=F)
        d=d[d[,2]!='*',]
        if (GFFmap=='gencode'){
		d$ENSEMBL_ID=unlist(lapply(d$V2,function(x) strsplit(x,'\\|')[[1]][2]))
                d$SYMBOL=unlist(lapply(d$V2,function(x) strsplit(x,'\\|')[[1]][6]))
                d$region=unlist(lapply(d$V2,function(x) strsplit(x,'\\|')[[1]][8]))
                write.table(d[,c('V1','ENSEMBL_ID','SYMBOL','region')],paste0(outdir,'/sequence/anno_filter.tsv'),row.names=F,col.names=T,quote=F,sep='\t')
                annotres=read.delim(paste0(outdir,'/sequence/anno_filter.tsv'),header=T,sep='\t')
	}else if (GFFmap=='refseq'){
		refseq_hash=readRDS('loaddat/reference/refseq/refseq_gene.rds')
		d$genbank=unlist(lapply(d[,2],function(x) paste(strsplit(x,'_')[[1]][4:5],collapse='_')))	
		d=merge(d,refseq_hash,by='genbank')
		d=d[,c('V1','genbank','symbol')]	
		write.table(d[,c('V1','genbank','symbol')],paste0(outdir,'/sequence/anno_filter.tsv'),row.names=F,col.names=T,quote=F,sep='\t')
                annotres=read.delim(paste0(outdir,'/sequence/anno_filter.tsv'),header=T,sep='\t')
	}
	return (annotres)
}
