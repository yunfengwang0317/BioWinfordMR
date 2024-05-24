library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(dendsort)
source('loaddat/Mendelian_Randomization.r')
heatmap_circle<-function(mr_result,namecol='id.exposure',value='pval', attr=NULL,showname=T,col_low,col_high,clust='YES',outdir,savefig=F){
	if (length(unique(mr_result$exposure))==1) names(mr_result)[1:4]=c('id.outcome','id.exposure','exposure','outcome')
	rawidx=unique(mr_result$id.exposure)
	if (namecol=='id.exposure') data<-reshape2::dcast(mr_result[,c(namecol,'method','pval')],id.exposure~method,mean)
	if (namecol=='exposure') {
		data<-reshape2::dcast(mr_result[,c(namecol,'method','pval')],exposure~method,mean)
		data$exposure=gsub('class |genus |family |order |phylum ','',data$exposure)
	}
	CN=names(data)
	data=data.frame(data[,-1],row.names=data[,1])
	names(data)=CN[-1]
        names(data)=gsub('Inverse variance weighted','IVW',names(data))
        rownames(data)=shortenname(rownames(data),30)
        data[is.na(data)]=1
	
	if (namecol=='id.exposure')data<-reshape2::dcast(mr_result[,c(namecol,'method',value)],id.exposure~method,mean)
	if (namecol=='exposure')data<-reshape2::dcast(mr_result[,c(namecol,'method',value)],exposure~method,mean)

	df=unique(mr_result[,c('id.exposure','exposure')])
	hash=df$id.exposure
	names(hash)=df$exposure
	if (namecol=='exposure'){
		rownames(data)=hash[data$exposure]
		data=data[unique(mr_result$id.exposure),]
		data$exposure=gsub('class |genus |family |order |phylum ','',data$exposure)
	}else{
		rownames(data)=data[,1]
		data=data[unique(mr_result$id.exposure),]
	}
	CN=names(data)
	data=data.frame(data[,-1],row.names=data[,1])
	names(data)=CN[-1]
	names(data)[names(data)=='Inverse variance weighted']='IVW'
	rownames(data)=shortenname(rownames(data),25)
	data[is.na(data)]=1
	sig=as.matrix(data)
	madt2<-as.matrix(data)
	#重新定义热图颜色梯度：
	mycol = colorRamp2(c(quantile(unlist(c(data)),.05), quantile(unlist(c(data)),.5), quantile(unlist(c(data)),.95)), c(col_low, "#FFFFFF", col_high))
	if (value=='pval') mycol = colorRamp2(c(0, 0.05,1), c(col_low, "#FFFFFF", col_high))
	if (savefig==T) pdf(paste0(outdir,'/circos_',namecol,'.pdf'),wi=8,he=8,onefile=F)
	split=NULL
	if (!is.null(attr)){
		attr=read.table(attr,header=T,sep='\t')
		names(attr)=c('exposure','group')
		attr$id=hash[attr$exposure]
		if (namecol=='id.exposure') {
			rownames(attr)=attr$id
			attr=attr[rownames(madt2),'group']
		}else {
			rownames(attr)=shortenname(attr$exposure)
			attr=attr[rownames(madt2),'group']
		}
		split=factor(attr)
	}
	plot.new()
	circos.clear()
	range(madt2)	
	circos.par(gap.after=c(50))
	
	cols=rep('black',dim(data)[1])
	fonts=rep(1,dim(data)[1])
	cexs=rep(.5,dim(data)[1])
	cols[sig[,1]<0.05]='red'
	fonts[sig[,1]<0.05]=2
	cexs[sig[,1]<0.05]=.6
	if (length(cexs)<60 & max(nchar(rownames(madt2)))<15) cexs=cexs*2
	showcluster=T
	if (clust=='NO') showcluster=F
	print (split)
	circos.heatmap(madt2,split=split,
		show.sector.labels = ifelse(!is.null(split),T,F),
		col=mycol,dend.side="inside",
		rownames.side="outside",
		track.height = min(0.2,1.5/dim(madt2)[2]),
		rownames.col=cols,
		rownames.cex=cexs,
		rownames.font=fonts,
		cluster=showcluster,
		dend.track.height=0.18,
		dend.callback=function(dend,m,si) {
			color_branches(dend,k=4,col=2:5)
		}
	)
	lg=Legend(title=value,col_fun=mycol,direction = c("vertical"))
	grid.draw(lg)

	#添加列名：
	if (showname==T){
		circos.track(track.index=2,
			panel.fun=function(x,y){
				if(CELL_META$sector.numeric.index==1){
					cn=colnames(madt2)
					n=length(cn)
					circos.text(rep(CELL_META$cell.xlim[2],n)+convert_x(0.8,"mm"),#x坐标
					1:n - convert_y(0.5, "mm"), ##y轴坐标
					rev(cn),cex=1,adj=c(0,1),facing="inside")
				}
			},bg.border=NA
		)
	}
	circos.clear()
	print ('basegraph finished')
	#circos.par()调整圆环首尾间的距离，数值越大，距离越宽
	#dend.side：控制行聚类树的方向，inside为显示在圆环内圈，outside为显示在圆环外圈
	#rownames.side：控制矩阵行名的方向,与dend.side相同；但注意二者不能在同一侧，必须一内一外
	#cluster=TRUE为对行聚类，cluster=FALSE则不显示聚类
	#track.height：轨道的高度，数值越大圆环越粗
	#dend.track.height：调整行聚类树的高度
	#dend.callback：用于聚类树的回调，当需要对聚类树进行重新排序，或者添加颜色时使用
	#包含的三个参数：dend：当前扇区的树状图；m：当前扇区对应的子矩阵；si：当前扇区的名称
	#color_branches():修改聚类树颜色

	if (savefig==T)dev.off()

}
