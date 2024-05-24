library(data.table)
library(tidyverse)
library(dplyr)
library(ldscr)
source('loaddat/Mendelian_Randomization.r')
ldsc_rg_fix<-function (munged_sumstats, ancestry, sample_prev = NA, population_prev = NA, 
    ld, wld, n_blocks = 200, chisq_max = NA, chr_filter = seq(1, 
        22, 1)) 
{
    if (missing(ancestry)) {
        cli::cli_progress_step("No ancestry specified, checking for user-specified `ld` and `wld`")
        checkmate::assert_directory_exists(ld)
        checkmate::assert_directory_exists(wld)
    }else {
        checkmate::assert_choice(ancestry, c("AFR", "AMR", "CSA", 
            "EAS", "EUR", "MID"), null.ok = FALSE)
        cli::cli_progress_step("Using {ancestry} reference from Pan-UKB")
    }
    checkmate::assert_list(munged_sumstats)
    if (missing(sample_prev)) {
        cli::cli_alert_info("No sample prevalence data provided. Estimating heritabilities on the observed scale.")
    }
    checkmate::assert_number(n_blocks)
    n.blocks <- n_blocks
    n.traits <- length(munged_sumstats)
    n.V <- n.traits * (n.traits + 1)/2
    if (n.traits > 18) {
        n.blocks <- (((n.traits + 1) * (n.traits + 2))/2) + 1
        cli::cli_alert_info("Setting the number of blocks used to perform the block jacknife used to estimate the sampling covariance matrix (V) to {n.blocks}")
        if (n.blocks > 1000) {
            cli_alert_warning("The number of blocks needed to estimate V is > 1000, which may result in sampling dependencies across the blocks used to estimate standard errors and can bias results.")
        }
    }
    cov <- matrix(NA, nrow = n.traits, ncol = n.traits)
    V.hold <- matrix(NA, nrow = n.blocks, ncol = n.V)
    N.vec <- matrix(NA, nrow = 1, ncol = n.V)
    Liab.S <- rep(1, n.traits)
    I <- matrix(NA, nrow = n.traits, ncol = n.traits)
    h2_res <- tibble()
    cli::cli_progress_step("Reading LD Scores")
    x <- ldscr:::read_ld(ancestry, ld)
    x$CM <- x$MAF <- NULL
    cli::cli_progress_step("Reading weights")
    w <- ldscr:::read_wld(ancestry, wld)
    w$CM <- w$MAF <- NULL
    colnames(w)[ncol(w)] <- "wLD"
    cli::cli_progress_step("Reading M")
    m <- ldscr:::read_m(ancestry, ld)
    M.tot <- sum(m)
    m <- M.tot
    cli::cli_progress_step("Reading summary statistics")
    all_y <- purrr::imap(munged_sumstats, ~{
        if (is.character(.x)) {
            cli::cli_progress_step("Reading summary statistics for '{.y}' from {.x}")
            sumstats_df <- vroom::vroom(.x, col_types = vroom::cols())
        }else {
            cli::cli_progress_step("Reading summary statistics for '{.y}' from dataframe")
            sumstats_df <- .x
        }
        cli::cli_progress_step("Merging '{.y}' with LD-score files")
        merged <- ldscr:::merge_sumstats(sumstats_df, w, x, chr_filter)
        cli::cli_alert_info(glue::glue("{nrow(merged)}/{nrow(sumstats_df)} SNPs remain after merging '{.y}' with LD-score files"))
        if (is.na(chisq_max)) {
            chisq_max <- max(0.001 * max(merged$N), 80)
        }
        rm <- (merged$Z^2 > chisq_max)
        merged <- merged[!rm, ]
        cli::cli_alert_info(glue::glue("Removed {sum(rm)} SNPs with Chi^2 > {chisq_max} from '{.y}'; {nrow(merged)} SNPs remain"))
        return(merged)
    })
    s <- 1
    for (j in 1:n.traits) {
        y1 <- all_y[[j]]
        y1$chi1 <- y1$Z^2
        for (k in j:n.traits) {
            if (j == k) {
                trait <- names(munged_sumstats[j])
                cli::cli_progress_step("Estimating heritability for '{trait}'")
                samp.prev <- sample_prev[j]
                pop.prev <- population_prev[j]
                merged <- y1
                n.snps <- nrow(merged)
                merged$intercept <- 1
                merged$x.tot <- merged$L2
                merged$x.tot.intercept <- 1
                initial.w <- ldscr:::make_weights(chi1 = merged$chi1, 
                  L2 = merged$L2, wLD = merged$wLD, N = merged$N, 
                  M.tot)
                merged$weights <- initial.w/sum(initial.w)
                N.bar <- mean(merged$N)
                weighted.LD <- as.matrix(cbind(merged$L2, merged$intercept) * 
                  merged$weights)
                weighted.chi <- as.matrix(merged$chi1 * merged$weights)
                analysis_res <- ldscr:::perform_analysis(n.blocks, n.snps, 
                  weighted.LD, weighted.chi, N.bar, m)
                V.hold[, s] <- analysis_res$pseudo.values
                N.vec[1, s] <- analysis_res$N.bar
                lambda.gc <- median(merged$chi1)/qchisq(0.5, 
                  df = 1)
                mean.Chi <- mean(merged$chi1)
                ratio <- (analysis_res$intercept - 1)/(mean.Chi - 
                  1)
                ratio.se <- analysis_res$intercept.se/(mean.Chi - 
                  1)
                if (is.na(population_prev) == F & is.na(sample_prev) == 
                  F) {
                  h2_lia <- h2_liability(h2 = analysis_res$reg.tot, 
                    sample_prev, population_prev)
                  h2_res <- h2_res %>% bind_rows(tibble(trait = trait, 
                    mean_chisq = mean.Chi, lambda_gc = lambda.gc, 
                    intercept = analysis_res$intercept, intercept_se = analysis_res$intercept.se, 
                    ratio = ratio, ratio_se = ratio.se, h2_observed = analysis_res$reg.tot, 
                    h2_observed_se = analysis_res$tot.se, h2_Z = analysis_res$reg.tot/analysis_res$tot.se, 
                    h2_p = 2 * pnorm(abs(h2_Z), lower.tail = FALSE), 
                    h2_liability = h2_lia, h2_liability_se = h2_lia/h2_Z))
                }else {
                  h2_res <- h2_res %>% bind_rows(tibble(trait = trait, 
                    mean_chisq = mean.Chi, lambda_gc = lambda.gc, 
                    intercept = analysis_res$intercept, intercept_se = analysis_res$intercept.se, 
                    ratio = ratio, ratio_se = ratio.se, h2_observed = analysis_res$reg.tot, 
                    h2_observed_se = analysis_res$tot.se, h2_Z = analysis_res$reg.tot/analysis_res$tot.se, 
                    h2_p = 2 * pnorm(abs(h2_Z), lower.tail = FALSE)))
                }
                cov[j, j] <- analysis_res$reg.tot
                I[j, j] <- analysis_res$intercept
            }
            if (j != k) {
                trait1 <- names(munged_sumstats[j])
                trait2 <- names(munged_sumstats[k])
                cli::cli_progress_step("Estimating genetic covariance for for '{trait1}' and '{trait2}'")
                y2 <- all_y[[k]]
                y <- merge(y1, y2[, c("SNP", "N", "Z", "A1")], 
                  by = "SNP", sort = FALSE)
                y$Z.x <- ifelse(y$A1.y == y$A1.x, y$Z.x, -y$Z.x)
                y$ZZ <- y$Z.y * y$Z.x
                y$chi2 <- y$Z.y^2
                merged <- na.omit(y)
                n.snps <- nrow(merged)
                merged$intercept <- 1
                merged$x.tot <- merged$L2
                merged$x.tot.intercept <- 1
                initial.w <- ldscr:::make_weights(chi1 = merged$chi1, 
                  L2 = merged$L2, wLD = merged$wLD, N = merged$N.x, 
                  M.tot)
                initial.w2 <- ldscr:::make_weights(chi1 = merged$chi2, 
                  L2 = merged$L2, wLD = merged$wLD, N = merged$N.y, 
                  M.tot)
                merged$weights_cov <- (initial.w + initial.w2)/sum(initial.w + 
                  initial.w2)
                N.bar <- sqrt(mean(merged$N.x) * mean(merged$N.y))
                weighted.LD <- as.matrix(cbind(merged$L2, merged$intercept) * 
                  merged$weights)
                weighted.chi <- as.matrix(merged$ZZ * merged$weights_cov)
                covariance_res <- ldscr:::perform_analysis(n.blocks, 
                  n.snps, weighted.LD, weighted.chi, N.bar, m)
                V.hold[, s] <- covariance_res$pseudo.values
                N.vec[1, s] <- covariance_res$N.bar
                cov[k, j] <- cov[j, k] <- covariance_res$reg.tot
                I[k, j] <- I[j, k] <- covariance_res$intercept
            }
            s <- s + 1
        }
    }
    v.out <- cov(V.hold)/crossprod(N.vec * (sqrt(n.blocks)/m))
    ratio <- tcrossprod(sqrt(Liab.S))
    S <- cov * ratio
    scaleO <- gdata::lowerTriangle(ratio, diag = TRUE)
    V <- v.out * tcrossprod(scaleO)
    colnames(S) <- if (is.null(names(munged_sumstats))) { 
        paste0("V", 1:ncol(S))
    }else names(munged_sumstats)
    rownames(S) <- if (is.null(names(munged_sumstats))) {
        paste0("V", 1:ncol(S))
    }else names(munged_sumstats)
    if (mean(Liab.S) != 1) {
        r <- nrow(S)
        SE <- matrix(0, r, r)
        SE[lower.tri(SE, diag = TRUE)] <- sqrt(diag(V))
        colnames(SE) <- colnames(S)
        rownames(SE) <- rownames(S)
    }
    if (all(diag(S) > 0)) {
        ratio <- tcrossprod(1/sqrt(diag(S)))
        S_Stand <- S * ratio
        scaleO <- gdata::lowerTriangle(ratio, diag = TRUE)
        V_Stand <- V * tcrossprod(scaleO)
        r <- nrow(S)
        SE_Stand <- matrix(0, r, r)
        SE_Stand[lower.tri(SE_Stand, diag = TRUE)] <- sqrt(diag(V_Stand))
        colnames(SE_Stand) <- colnames(S)
        rownames(SE_Stand) <- rownames(S)
    }else {
        cli::cli_alert_warning("Your genetic covariance matrix includes traits estimated to have a negative heritability.")
        return ()
    }
    ind <- which(lower.tri(S, diag = F), arr.ind = TRUE)
    rg_res <- tibble(trait1 = dimnames(S_Stand)[[2]][ind[, 2]], 
        trait2 = dimnames(S_Stand)[[1]][ind[, 1]], rg = S_Stand[ind], 
        rg_se = SE_Stand[ind], rg_p = 2 * pnorm(abs(rg/rg_se), 
            lower.tail = FALSE))
    output <- list(h2 = h2_res, rg = rg_res, raw = list(V = V, 
        S = S, I = I, N = N.vec, m = m, V_Stand = V_Stand, S_Stand = S_Stand, 
        SE_Stand = SE_Stand))
    class(output) <- c("ldscr_list", "list")
    return(output)
}

LDSC_rg<-function(expo,outcome,an='EUR',sample_prev=NA, population_prev=NA,ld,wld,chr_filter=c(1:22),n_blocks=10){
    id.o<-outcome$id.outcome[1]
    id.e<-expo$id.exposure[1]
    
    expo<-expo%>%mutate(Z=beta.exposure/se.exposure)
    expo<-expo%>%select(SNP=SNP,N=samplesize.exposure,Z=Z
                        ,A1=effect_allele.exposure
                        ,A2=other_allele.exposure)
    expo<-as_tibble(expo)
    
    outcome<-outcome%>%mutate(Z=beta.outcome/se.outcome)
    outcome<-outcome%>%select(SNP=SNP,N=samplesize.outcome,Z=Z
                              ,A1=effect_allele.outcome
                              ,A2=other_allele.outcome)
    outcome<-as_tibble(outcome)
    dat<-list(expo,outcome)
    names(dat)<-c(id.e,id.o)  
    rm(expo,outcome)
    res<-try(ldsc_rg_fix(dat,ancestry = an,sample_prev=sample_prev,
                            population_prev=population_prev,ld=ld,wld=wld,
                            n_blocks=n_blocks,chr_filter=chr_filter))
    return(res)
  }
  


runLDSC_old<-function(exp_id=NULL,exp_file=NULL,out_id=NULL,out_file=NULL,pval=5e-6,show='h2',outdir){
	pval=as.numeric(pval)
	dir.create(outdir)
	ld='loaddat/ldscr/inst/extdata/EUR/'
	if (!is.null(exp_file))exp_id='file'
	if (!file.exists(paste0(outdir,'/',exp_id,'_h2.txt'))){
		exposure=retrieve_exposure(exposure_id=exp_id,expfile=exp_file,pvalue=pval,clump_kb=10000,clump_r2=0.001,clump=F,outdir=outdir)
		if (is.null(exposure)) return (data.frame(msg="Not enough exposure SNP obtained"))
		outcome=retrieve_outcome(exposure,outcome_id=out_id,outfile=out_file,outdir=outdir)
		if (is.null(outcome)) return (data.frame(msg="Not enough shared SNP obtained"))
		LDSC_res<-LDSC_rg(exposure,outcome,n_blocks=2,ld=ld,wld=ld)
		if (is.null(LDSC_res)) return (data.frame(msg='Your genetic covariance matrix includes traits estimated to have a negative heritability.'))
	  	h2 <- LDSC_res[["h2"]]
		rg <- LDSC_res[["rg"]]
		write.table(h2,file = paste0(outdir,"/",exp_id,"_h2.txt"), row.names = FALSE,quote=F,sep='\t')
		write.table(rg,file = paste0(outdir,"/",exp_id,"_rg.txt"), row.names = FALSE,quote=F,sep='\t')
	}else{
		h2=read.delim(paste0(outdir,"/",exp_id,"_h2.txt"),header=T,sep='\t')
		rg=read.delim(paste0(outdir,"/",exp_id,"_rg.txt"),header=T,sep='\t')
	}
	if (show=='h2') return (h2)
	if (show=='rg') return (rg)
}

source('loaddat/Mendelian_Randomization.r')
runLDSC<-function(exp_id='',exp_file=NULL,out_id='',out_file=NULL,outdir){
	tmpdir=gsub('result','tmp',outdir)
	dir.create(tmpdir,recursive=T)
	if (exp_id!='' | out_id!='') outdir=paste0(outdir,'/',exp_id,'_',out_id)
	dir.create(outdir,recursive=T)
	if (is.null(exp_file) & exp_id!=''){
                exp_file=check_avail_outcome(exp_id)
        }
	if (is.null(exp_file))return ()
	readinput(exp_file,tmpdir,paste0(tmpdir,'/gwas1.txt'),snplst=NULL,keepallsnp=T)
	trait1=fread(paste0(tmpdir,'/gwas1.txt'))
	names(trait1)[1:10]=c('SNP','ALT','REF','BETA','SE','P','N','exposure','Frq','ID')
	trait1$N[is.na(trait1$N)]=10000
	trait1$N[trait1$N=='']=10000
	name1=trait1$ID[1]
	if (!file.exists(paste0(tmpdir,'/',name1,'.txt.gz'))){
		write.table(trait1[,c(1,2,3,4,5,6,7,9)],paste0(tmpdir,'/',name1,'.txt'),sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
		system(paste0('gzip ',tmpdir,'/',name1,'.txt'))
	}
	if (is.null(out_file) & out_id!=''){
                out_file=check_avail_outcome(out_id)
        }
	if (is.null(out_file)) return ()
	readinput(out_file,tmpdir,paste0(tmpdir,'/gwas2.txt'),snplst=NULL,keepallsnp=T)
        trait2=fread(paste0(tmpdir,'/gwas2.txt'))
        names(trait2)[1:10]=c('SNP','ALT','REF','BETA','SE','P','N','exposure','Frq','ID')
	trait2$N[is.na(trait2$N)]=10000
	trait2$N[trait2$N=='']=10000
	name2=trait2$ID[1]
	if (!file.exists(paste0(tmpdir,'/',name2,'.txt.gz'))){
		write.table(trait2[,c(1,2,3,4,5,6,7,9)],paste0(tmpdir,'/',name2,'.txt'),sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
		system(paste0('gzip ',tmpdir,'/',name2,'.txt'))
	}

	system(paste0('/usr/bin/anaconda3/envs/ldsc/bin/python /home/shiny/ldsc-master/munge_sumstats.py --sumstats ',tmpdir,'/',name1,'.txt.gz --merge-alleles /home/shiny/ldsc-master/resource/w_hm3.snplist.bz2 --a1 ALT --a2 REF --chunksize 500000 --out ',tmpdir,'/',name1))
	system(paste0('/usr/bin/anaconda3/envs/ldsc/bin/python /home/shiny/ldsc-master/munge_sumstats.py --sumstats ',tmpdir,'/',name2,'.txt.gz --merge-alleles /home/shiny/ldsc-master/resource/w_hm3.snplist.bz2 --a1 ALT --a2 REF --chunksize 500000 --out ',tmpdir,'/',name2))
	system(paste0('/usr/bin/anaconda3/envs/ldsc/bin/python /home/shiny/ldsc-master/ldsc.py --h2 ',tmpdir,'/',name1,'.sumstats.gz --ref-ld-chr /home/shiny/ldsc-master/resource/LDscores/ --w-ld-chr /home/shiny/ldsc-master/resource/LDscores/ --out ',outdir,'/',name1))
	system(paste0('/usr/bin/anaconda3/envs/ldsc/bin/python /home/shiny/ldsc-master/ldsc.py --h2 ',tmpdir,'/',name2,'.sumstats.gz --ref-ld-chr /home/shiny/ldsc-master/resource/LDscores/ --w-ld-chr /home/shiny/ldsc-master/resource/LDscores/ --out ',outdir,'/',name2))
	system(paste0('/usr/bin/anaconda3/envs/ldsc/bin/python /home/shiny/ldsc-master/ldsc.py --rg ',tmpdir,'/',name1,'.sumstats.gz,',tmpdir,'/',name2,'.sumstats.gz --ref-ld-chr /home/shiny/ldsc-master/resource/LDscores/ --w-ld-chr /home/shiny/ldsc-master/resource/LDscores/ --out ',outdir,'/',name1,'_',name2))
	res=NULL
	logfi=paste0(outdir,'/',name1,'_',name2,'.log')
	if (file.exists(logfi)) res=system(paste0('cat ',logfi),intern=T)
	return (res)
}
