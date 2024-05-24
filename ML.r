library(MASS)
library(foreign)
library(rms)
library(rsample)
library(glmnet)
#library(tidyverse)
library(tidyr)
library(e1071)
library(caret)
library(randomForest)
library(plotROC)
library(pROC)
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
smooth<-function(v){
	loess(v~c(1:length(v)),span=.1)$fitted
}
ML<-function(dat,datapath,tmpdir,outdir){
		RN=unlist(lapply(names(dat),function(x)gsub('\\.','-',x)))
		dat=as.data.frame(scale(apply(t(dat),2,as.numeric)))
		names(dat)=unlist(lapply(names(dat),function(x)gsub('-','.',x)))
		rownames(dat)=RN
		
		con=read.delim(datapath,header=T,sep='\t')
                con=unique(con[con[,1] %in% rownames(dat),])
		dat=dat[rownames(dat) %in% con[,1],]
                g1=unique(con[,2])[1]
                g2=unique(con[,2])[2]
                condition=c(sum(con[,2]==g1),sum(con[,2]==g2))

		label=rep(g1,dim(dat)[1])
		label[rownames(dat) %in% con[con[,2]==g2,1]]=g2
                Y=as.factor(label)
                perm=1
		roc_nb=list()
		roc_svm=list()
		roc_rf=list()
		dat[is.na(dat)]=0

	while (perm<11){
		N=sample(1:dim(dat)[1],as.integer(dim(dat)[1]*0.75))
		train_dat=dat[N,]
                test_dat=dat[-N,]
	      	env <- parent.frame()
                #formula=as.formula(paste0('Y[N]~.'))
                #boruta_output <- Boruta('Y[N]~.', data=na.omit(train_dat), doTrace=0)
                #plot(boruta_output, cex.axis=.7, las=2, xlab="", main="Variable Importance")

                env <- parent.frame()
                formula=as.formula(paste0('Y[N]~.'))
                nb_model = naiveBayes(formula,data=train_dat,type='raw')
                predValid_nb <- predict(nb_model, test_dat,type='raw')[,1]
                
                env <- parent.frame()
                formula=as.formula(paste('Y[N]~.'))
                rf_model <- randomForest(formula, data = train_dat, ntree = 500, mtry = 6,importance = TRUE)
                predValid_rf <- predict(rf_model, test_dat, type = "prob",probability=T)[,2]

                env <- parent.frame()
                formula=as.formula(paste0('Y[N]~.'))
                svm_model <- svm(formula, data=train_dat,probability = TRUE)
                predValid_svm <- attr(predict(svm_model, test_dat, probability = TRUE,type='raw'),'probabilities')[,2]
                
                roc_nb[[perm]]=roc(as.integer(Y[-N]),as.numeric(predValid_nb))
	       	roc_svm[[perm]]=roc(as.integer(Y[-N]), as.numeric(predValid_svm))
                roc_rf[[perm]]=roc(as.integer(Y[-N]), as.numeric(predValid_rf))
                
		roc_nb[[perm]]$sensitivities=roc_nb[[perm]]$sensitivities[as.integer(seq(1,length(roc_nb[[perm]]$sensitivities),length.out=10))]
		roc_svm[[perm]]$sensitivities=roc_svm[[perm]]$sensitivities[as.integer(seq(1,length(roc_svm[[perm]]$sensitivities),length.out=10))]
		roc_rf[[perm]]$sensitivities=roc_rf[[perm]]$sensitivities[as.integer(seq(1,length(roc_rf[[perm]]$sensitivities),length.out=10))]	
		roc_nb[[perm]]$specificities=roc_nb[[perm]]$specificities[as.integer(seq(1,length(roc_nb[[perm]]$specificities),length.out=10))]
                roc_svm[[perm]]$specificities=roc_svm[[perm]]$specificities[as.integer(seq(1,length(roc_svm[[perm]]$specificities),length.out=10))]
                roc_rf[[perm]]$specificities=roc_rf[[perm]]$specificities[as.integer(seq(1,length(roc_rf[[perm]]$specificities),length.out=10))]        
		perm=perm+1
	}

		sen_nb=apply(matrix(unlist(lapply(roc_nb,function(x) x$sensitivities)),ncol=10),1,mean)
		spe_nb=apply(matrix(unlist(lapply(roc_nb,function(x) x$specificities)),ncol=10),1,mean)
		auc_nb=mean(unlist(lapply(roc_nb,function(x) x$auc)))
		sen_svm=apply(matrix(unlist(lapply(roc_svm,function(x) x$sensitivities)),ncol=10),1,mean)
                spe_svm=apply(matrix(unlist(lapply(roc_svm,function(x) x$specificities)),ncol=10),1,mean)
                auc_svm=mean(unlist(lapply(roc_svm,function(x) x$auc)))
		sen_rf=apply(matrix(unlist(lapply(roc_rf,function(x) x$sensitivities)),ncol=10),1,mean)
                spe_rf=apply(matrix(unlist(lapply(roc_rf,function(x) x$specificities)),ncol=10),1,mean)
                auc_rf=mean(unlist(lapply(roc_rf,function(x) x$auc)))
		
		plot(1-spe_nb, smooth(sen_nb),type='l',lty=1,xlim=c(0,1),ylim=c(0,1),main='ROC of 10-fold CV',xlab='specificity',ylab='sensitivity')
                lines(1-spe_svm,smooth(sen_svm), col="red", lty=2)
                lines(sort(1-spe_rf),smooth(sen_rf[order(1-spe_rf)]), col='blue', lty=3)
                text(0.45,0.5,paste0('AUC=',round(auc_nb,2)),col='black')
		text(0.45,0.55,paste0('AUC=',round(auc_svm,2)),col='red')
                text(0.45,0.6, paste0('AUC=',round(auc_rf,2)),col='blue')
                legend('bottomright',legend=c('naive bayes','SVM','RF'),col=c('black','red','blue'),lwd=4,lty=c(1,2,3))
		
		pdf(paste0(outdir,'/machine_learning/ROC.pdf'),wi=6,he=6)
                par(pty='s')
		plot(1-spe_nb, smooth(sen_nb),type='l',lty=1,xlim=c(0,1),ylim=c(0,1),main='ROC of 10-fold CV',xlab='specificity',ylab='sensitivity')
                lines(1-spe_svm,smooth(sen_svm), col="red", lty=2)
                lines(sort(1-spe_rf),smooth(sen_rf[order(1-spe_rf)]), col='blue', lty=3)
                text(0.45,0.5,paste0('AUC=',round(auc_nb,2)),col='black')
                text(0.45,0.55,paste0('AUC=',round(auc_svm,2)),col='red')
                text(0.45,0.6, paste0('AUC=',round(auc_rf,2)),col='blue')
                legend('bottomright',legend=c('naive bayes','SVM','RF'),col=c('black','red','blue'),lwd=4,lty=c(1,2,3))
		dev.off()
}



ML_py<-function(dat,datapath,models,X_val=NULL,y_val=NULL,metrics='Acc',tmpdir,outdir){
	source('loaddat/readinput.r')
	dat=t(dat)
        if (dim(dat)[2]>50){
		print ('too many genes are used')
		vars=apply(dat,2,var)
		vars_=rev(sort(vars))
		dat=dat[,vars>vars_[50]]
		print (paste0(dim(dat)[2],' genes were selected as HVG'))
	}
	con=read.delim(datapath,header=T,sep='\t',row.names=1)
	rownames(con)=unlist(lapply(rownames(con),function(x)gsub('\\.','-',x)))
	if (length(unique(con[,1]))>10){#continuous values exist instead of labels
		MED=median(con[,1])
		tmp=con[,1]
		con[tmp<MED,1]='low'
		con[tmp>=MED,1]='high'
	}
	dat=merge(dat,con,by='row.names')
	rownames(dat)=dat[,1]
	dat=dat[,-1]
	feature=colnames(dat)
        g1=unique(dat[,dim(dat)[2]])
        g2=unique(dat[,dim(dat)[2]])
	write.table(dat,paste0(tmpdir,'/input.txt'),row.names=T,col.names=NA,quote=F,sep='\t')
	source('loaddat/nomogram_clf.r')
	tryCatch({
		nomoclf(paste0(tmpdir,'/input.txt'),paste0(outdir,'/machine_learning'))
	},error=function(cond) {print('failed')})
	py='/Users/yunfengwang/miniconda3/envs/ML/bin/python'
	if (!file.exists(py)) py='/usr/bin/anaconda3/bin/python'
	write.table(models,paste0(outdir,'/machine_learning/models.txt'),row.names=F,col.names=F,quote=F)
	
	job=paste0(py,' loaddat/Machine_Learning_Models.py ',tmpdir,'/input.txt ',outdir,'/machine_learning')
	print (job)
	system(job)
	res=read.table(paste0(outdir,'/machine_learning/accuracy.txt'),header=T,sep='\t')
	
	if ('ggradar' %in% installed.packages()){
		library(ggradar)
		dat=res
		input=data.frame(do.call(rbind,lapply(unique(dat$type),function(x) c(x,dat$TestSet.Acc[dat$type==x]))))
		rownames(input)=unique(dat$type)
		colnames(input)=c('type',dat$model[dat$type=='Acc'])
		for (i in 2:dim(input)[2]) input[,i]=as.numeric(input[,i])
		pdf(paste0(outdir,'/machine_learning/radar_graph.pdf'),wi=8,he=6)
		print (ggradar(plot.data = input,axis.label.size=3))
		dev.off()
	}
	
	res=res[res$type==metrics,]
	#return (round(res,3))
	res$Avg.Acc=(res$TrainSet.Acc+res$TestSet.Acc)/2
	res=res[res$model %in% models,]
	rownames(res)=res$model
	res=res[,-1]
	print (res)
	p=ggplot(res,aes(x=TestSet.Acc,
              y=TrainSet.Acc,
              size=Avg.Acc,
              fill=Avg.Acc))+
		  geom_point(alpha=0.75,shape=21,color='black')+
		  coord_flip()+xlab(paste0('TestSet ',metrics))+ylab(paste0('TrainSet ',metrics))+
		  theme_bw()+ggtitle(metrics)+
		  theme(plot.title = element_text(hjust = .5,size=16),
                       legend.position = "right",
                        panel.grid=element_blank(),
                        axis.title = element_text(size = 16),
                        axis.text = element_text(size = 14))+
		  scale_fill_gradient2(guide='none',
                       low = "#0571b0",
                       mid = "white",
                       high = "#ca0020",
                       midpoint = median(res$Avg.Acc))+
		  scale_size_continuous(range = c(2, 10))+
		  #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .1)) + 
		  #geom_hline(yintercept = median(res$Avg.Acc), linetype = 2, color = 'gray55')+
		  geom_text_repel(aes(TestSet.Acc, TrainSet.Acc, label = rownames(res)),size=res$Avg.Acc*5,max.overlaps=20)
	print (p)
	save(p,file=paste0(tmpdir,'/MLfig.rda'))
	
	if (!is.null(X_val) & !is.null(y_val)){
	
	dir.create(paste0(outdir,'/machine_learning/validate'))
	write.table(models,paste0(outdir,'/machine_learning/validate/models.txt'),row.names=F,col.names=F,quote=F)
	readinput(X_val,tmpdir,paste0(tmpdir,'/X_val.txt'))
	readinput(y_val,tmpdir,paste0(tmpdir,'/y_val.txt'))
	X_val=read.delim(paste0(tmpdir,'/X_val.txt'),header=T,row.names=1,sep='\t',check.names=F)
	dat=t(X_val)
        con=read.delim(paste0(tmpdir,'/y_val.txt'),header=T,sep='\t',row.names=1)
	names(con)[1]='condition'
        dat=merge(dat,con,by='row.names')
	if (mode(dat$condition)=='numeric')dat$condition=paste0('group',dat$condition)
	RN=dat[,1]
	dat=dat[,-1]

	dat=dat[,colnames(dat) %in% c(feature,'condition')]
	#dat=dat[,match(feature,colnames(dat))]
        g1=unique(dat[,dim(dat)[2]])[1]
        g2=unique(dat[,dim(dat)[2]])[2]
	rownames(dat)=RN
        write.table(dat,paste0(tmpdir,'/validate.txt'),row.names=T,col.names=NA,quote=F,sep='\t')
	job=paste0(py,' loaddat/Machine_Learning_Models.py ',tmpdir,'/input.txt ',outdir,'/machine_learning/validate ',tmpdir,'/validate.txt')
        print (job)
        system(job)
	cols=colorRampPalette(brewer.pal(12,'Paired'))(20)
        dat=read.table(paste0(outdir,'/machine_learning/validate/probares.txt'),header=T,sep='\t',row.names=1)
        truelabel=dat[,1]
        clf=names(dat)[2]
        roc1=roc(as.factor(truelabel), as.numeric(dat[,clf]))
        #roc1$sensitivities=smooth(roc1$sensitivities)
        #roc1$specificities=smooth(roc1$specificities)
        plot(roc1, print.auc=F, col=cols[1], legacy.axes=T, main="ROC of validation dataset")
        text(0.45, 0.03, paste0('AUC of ',clf,":",round(roc1$auc,3)), col=cols[1],pos=4,cex=1)
        if (dim(dat)[2]>2){
            for (i in 3:dim(dat)[2]){
                clf=names(dat)[i]
                roc1=roc(as.factor(truelabel), as.numeric(dat[,clf]))
                #roc1$sensitivities=smooth(roc1$sensitivities)
                #roc1$specificities=smooth(roc1$specificities)
                lines(roc1, print.auc=F, col=cols[i], legacy.axes=T)
                text(0.45, 0.03+(i-2)/20, paste0("AUC of ",clf,':',round(roc1$auc,3)),pos=4, col=cols[i],cex=1)
            }
        }	
	
	}else print(p)
}




Reg_py<-function(inputfile,models,tmpdir,outdir){
	dir.create(outdir,recursive=T)
	py='/Users/yunfengwang/miniconda3/envs/ML/bin/python'
    if (!file.exists(py)) py='/usr/bin/anaconda3/bin/python'
    write.table(models,paste0(outdir,'/machine_learning/models.txt'),row.names=F,col.names=F,quote=F)
    job=paste0(py,' loaddat/Regression_Models.py ',inputfile,' ',outdir,'/machine_learning')
    print (job)
    system(job)
		
    library(ggplot2)
    d=read.table(paste0(outdir,'/machine_learning/predict_Y_models.txt'),header=T,sep='\t')
    models=names(d)[-1]
    d=d[order(d$Y_test),]
    df=do.call(rbind,lapply(names(d)[2:dim(d)[2]],function(x) {
	 tmp=d[,c('Y_test',x)]
	 names(tmp)=c('Y_test','Y_pred')
	 tmp$model=x
	 tmp
	 }))
	#Code
	library(ggpubr)
	p=ggscatter(df, x = 'Y_test', y = 'Y_pred',
	          add = "reg.line",size=1,
	          add.params = list(color = "blue", fill = "lightgray"),
	          color = "black", palette = "jco", fill = "lightgray",
	          fullrange = TRUE,
	          rug = TRUE, facet.by = "model", cor.coef = T,
	          title = paste0('Correlation between Y_test and Y_pred'),
	          conf.int = TRUE,
	          cor.method = "spearman",
	          cor.coef.coord = c(NULL, NULL),
	          cor.coef.size = 3,
	          cor.coeff.args = list(method = "spearman",col='red', label.x = 3, label.sep = " ")
	)

	pdf(paste0(outdir,'/machine_learning/correlation_Y_test_Y_pred.pdf'),wi=10,he=10)
	print (p)
	dev.off()
	return (p)
}
