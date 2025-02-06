#!/usr/bin/env Rscript

library(e1071)
library(optparse)

args <- commandArgs(trailingOnly = TRUE)
N_s <-as.numeric(args[2])

final.sumstats <- c()
#Initialize Fst
for (s in 1:(N_s-1)){
  tmp <- sprintf('Fst.s%s.s%s', as.character(s), as.character((s+1):N_s))
  final.sumstats[tmp] <- NA
}
tmp <- sprintf('Fst.s%s.adm', as.character(1:N_s))
final.sumstats[tmp] <- NA

#Initialize het
tmp <- sprintf('mean.het.%s', c(paste0('s', as.character(1:N_s)),'adm'))
final.sumstats[tmp] <- NA
tmp <- sprintf('var.het.%s', c(paste0('s', as.character(1:N_s)),'adm'))
final.sumstats[tmp] <- NA

#Initialize F
tmp <- sprintf('mean.F.%s', c(paste0('s', as.character(1:N_s)),'adm'))
final.sumstats[tmp] <- NA
tmp <- sprintf('var.F.%s', c(paste0('s', as.character(1:N_s)),'adm'))
final.sumstats[tmp] <- NA

#Initialize f3
for (s in 1:(N_s-1)){
  tmp <- sprintf('f3.s%s.s%s', as.character(s), as.character((s+1):N_s))
  final.sumstats[tmp] <- NA
}

#Initialize ASD
for (s in 1:(N_s-1)){
  tmp <- sprintf('mean.ASD.s%s.s%s', as.character(s), as.character((s+1):N_s))
  final.sumstats[tmp] <- NA
}
tmp <- sprintf('mean.ASD.s%s.adm', as.character(1:N_s))
final.sumstats[tmp] <- NA

for (s in 1:(N_s-1)){
  tmp <- sprintf('var.ASD.s%s.s%s', as.character(s), as.character((s+1):N_s))
  final.sumstats[tmp] <- NA
}
tmp <- sprintf('var.ASD.s%s.adm', as.character(1:N_s))
final.sumstats[tmp] <- NA

tmp <- sprintf('mean.ASD.%s', c(paste0('s', as.character(1:N_s)),'adm'))
final.sumstats[tmp] <- NA
tmp <- sprintf('var.ASD.%s', c(paste0('s', as.character(1:N_s)),'adm'))
final.sumstats[tmp] <- NA


#Initialize pairwise admixture proportions
for(s1 in 1:(N_s-1)){
  for(s2 in (s1+1):N_s){
    tmp <- sprintf('mean.adm.prop.s%s.s%s', as.character(s1), as.character(s2))
    final.sumstats[tmp] <- NA
    tmp <- sprintf('var.adm.prop.s%s.s%s', as.character(s1), as.character(s2))
    final.sumstats[tmp] <- NA
    tmp <- sprintf('skew.adm.prop.s%s.s%s', as.character(s1), as.character(s2))
    final.sumstats[tmp] <- NA
    tmp <- sprintf('kurt.adm.prop.s%s.s%s', as.character(s1), as.character(s2))
    final.sumstats[tmp] <- NA
    tmp <- sprintf('mode.adm.prop.s%s.s%s', as.character(s1), as.character(s2))
    final.sumstats[tmp] <- NA
    tmp <- sprintf('perc%d.adm.prop.s%s.s%s', seq(0,100,10), as.character(s1), as.character(s2))
    final.sumstats[tmp] <- NA
    
  }
}


## pairwise Fst computed by vcftools
get.Fsts <- function(){
    d <- read.table('result.fst')
    r <- 1
    for (s1 in 1:(N_s-1)){
      final.sumstats[paste0('Fst.s', as.character(s1), '.adm')] <- d[r,]
      r <- r+1
      for (s2 in (s1+1):N_s){
        final.sumstats[paste0('Fst.s', as.character(s1),'.s', as.character(s2))] <- d[r,]
        r <- r+1
      }
    }
    final.sumstats[paste0('Fst.s', as.character(N_s), '.adm')] <- d[r,]
    final.sumstats
}


## we use 'tryCatch' tricks to prevent problems if files do not exist
tryCatch(final.sumstats <- get.Fsts())#, error=function(e){invisible()}, warning=function(w){invisible()})



## F3 statistic, following Patterson 2012
compute.F3 <- function(){
    p.adm <- read.table('adm.frq', skip=1)[,6]
    d <- data.frame(p.adm)
    ## get allelic frequencies computed by vcftools
    for (s in 1:N_s){
      tmp <- sprintf('p.s%s', as.character(s))
      d[tmp] <- read.table(paste0('s', as.character(s), '.frq'), skip=1)[,6]
    }
    for (s1 in 1:(N_s-1)){
      for (s2 in (s1+1):N_s){
        ## compute numerators and denominators of statistic
        num <- sum((d$p.adm - d[paste0('p.s', as.character(s1))]) * (d$p.adm - d[paste0('p.s', as.character(s2))]))
        denom <- sum(2 * d$p.adm * (1-d$p.adm))
        ## compute statistic
        final.sumstats[paste0('f3.s', as.character(s1), '.s', as.character(s2))] <- num/denom
      }
    }
    final.sumstats
}

tryCatch(final.sumstats <- compute.F3())#, error=function(e){invisible()}, warning=function(w){invisible()})

compute.inbreeding.het <- function(){
    for (s in 1:N_s){
        ## inbreeding coefficients from vcftools
        d <- read.table(sprintf('s%s.het', as.character(s)), header=TRUE)
        final.sumstats[sprintf('mean.F.s%s', as.character(s))] <- mean(d$F)
        final.sumstats[sprintf('var.F.s%s', as.character(s))] <- var(d$F)
        ## we compute heterozigosity using allelic frequencies computed by vcftools
        d <- read.table(sprintf('s%s.frq', as.character(s)), skip=1)
        freqs <- d[,6]
        nb.sample <- d[1,4]
        tmp <- freqs*freqs + (1-freqs)*(1-freqs)
        tmp <- 1 - tmp
        tmp <- tmp * nb.sample / (nb.sample - 1)
        final.sumstats[sprintf('mean.het.s%s', as.character(s))] <- mean(tmp)
        final.sumstats[sprintf('var.het.s%s', as.character(s))] <- var(tmp)
    }
    ## inbreeding coefficients from vcftools
    d <- read.table('adm.het', header=TRUE)
    final.sumstats['mean.F.adm'] <- mean(d$F)
    final.sumstats['var.F.adm'] <- var(d$F)
    ## we compute heterozigosity using allelic frequencies computed by vcftools
    d <- read.table('adm.frq', skip=1)
    freqs <- d[,6]
    nb.sample <- d[1,4]
    tmp <- freqs*freqs + (1-freqs)*(1-freqs)
    tmp <- 1 - tmp
    tmp <- tmp * nb.sample / (nb.sample - 1)
    final.sumstats['mean.het.adm'] <- mean(tmp)
    final.sumstats['var.het.adm'] <- var(tmp)
    final.sumstats
}

tryCatch(final.sumstats <- compute.inbreeding.het())#, error=function(e){invisible()}, warning=function(w){invisible()})

## compute statistics based on ASD (which was computed by asd)
compute.ASD.stats <- function(){
    d <- read.table('data.asd.dist', header=TRUE, row.names=1)
    colnames(d) <- rownames(d)          #prevents problems if indiv names starts with [0-9]
    all.pops <- c(paste0('s', as.character(1:N_s)), 'adm')
    indivs <- NULL
    for (pop in all.pops){
        ids <- read.table(sprintf('%s.het', pop), header=TRUE, as.is=TRUE)[,1]
        indivs <- rbind(indivs, data.frame(ids=ids, pop=pop, stringsAsFactors=FALSE))
    }
    for (i in 1:(N_s+1)){
        pop1 <- all.pops[i]
        wanted.ids.1 <- indivs$ids[indivs$pop == pop1]
        for (j in i:(N_s+1)){
            pop2 <- all.pops[j]
            wanted.ids.2 <- indivs$ids[indivs$pop == pop2]
            tmp <- as.matrix(d[wanted.ids.1, wanted.ids.2])
            
            if (i == j){
                tmp <- tmp[lower.tri(tmp)]
                pop <- pop1
                final.sumstats[sprintf('mean.ASD.%s', pop)] <- mean(tmp)
                final.sumstats[sprintf('var.ASD.%s', pop)] <- var(tmp)
            } else {
                tmp <- as.vector(tmp)
                final.sumstats[sprintf('mean.ASD.%s.%s', pop1, pop2)] <- mean(tmp)
                final.sumstats[sprintf('var.ASD.%s.%s', pop1, pop2)] <- var(tmp)
            }
        }
    }
    final.sumstats
}


tryCatch(final.sumstats <- compute.ASD.stats())#, error=function(e){invisible()}, warning=function(e){invisible()})

## computes admixture proportions stats using projection in MDS
compute.adm.props <- function(){
    d <- read.table('data.asd.dist', header=TRUE, row.names=1)
    colnames(d) <- rownames(d)          #prevents problems if indiv names starts with [0-9]
    all.pops <- c(paste0('s', as.character(1:N_s)),'adm')
    indivs <- NULL
    for (pop in all.pops){
        ids <- read.table(sprintf('%s.het', pop), header=TRUE, as.is=TRUE)[,1]
        indivs <- rbind(indivs, data.frame(ids=ids, pop=pop, stringsAsFactors=FALSE))
    }
    mds <- suppressWarnings(cmdscale(as.matrix(d), k=2))
    for (s1 in 1:(N_s-1)){
      for(s2 in (s1+1):N_s){
        ids1 <- indivs$ids[indivs$pop == paste0('s', as.character(s1))]
        ids2 <- indivs$ids[indivs$pop == paste0('s', as.character(s2))]
        centroid1 <- colMeans(mds[ids1,])
        centroid2 <- colMeans(mds[ids2,])
        centroids.dist <- dist(rbind(centroid1,centroid2))[1]
        ids.adm <- indivs$ids[indivs$pop == 'adm']
        adm.props <- rep(0, length(ids.adm))
        vec.b <- centroid2 - centroid1
        for (i in 1:length(adm.props)){
            idx.adm <- which(rownames(mds) == ids.adm[i])
            vec.a <- mds[idx.adm,] - centroid1
            proj <- (sum(vec.a*vec.b) / sum(vec.b*vec.b)) * vec.b
            new.point <- centroid1 + proj
            adm.props[i] <- 1 - dist(rbind(centroid1, new.point))[1] / centroids.dist
        }
        print(adm.props)
        final.sumstats[paste0('mean.adm.prop.s', as.character(s1),'.s', as.character(s2))] <- mean(adm.props)
        final.sumstats[paste0('var.adm.prop.s', as.character(s1),'.s', as.character(s2))] <- var(adm.props)
        final.sumstats[paste0('skew.adm.prop.s', as.character(s1),'.s', as.character(s2))] <- skewness(adm.props)
        final.sumstats[paste0('kurt.adm.prop.s', as.character(s1),'.s', as.character(s2))] <- kurtosis(adm.props)
        tmp <- density(adm.props)
        final.sumstats[paste0('mode.adm.prop.s', as.character(s1),'.s', as.character(s2))] <- tmp$x[which.max(tmp$y)][1]
        tmp <- sprintf('perc%d.adm.prop.s%s.s%s', seq(0,100,10), as.character(s1), as.character(s2))
        final.sumstats[tmp] <- quantile(adm.props, seq(0,1,.1))
      }
    }
    final.sumstats
}

tryCatch(final.sumstats <- compute.adm.props())#, error=function(e){invisible()}, warning=function(e){invisible()})

compute.adm.angles <- function(){
    d <- read.table('data.asd.dist', header=TRUE, row.names=1)
    colnames(d) <- rownames(d)          #prevents problems if indiv names starts with [0-9]
    all.pops <- c(paste0('s', as.character(1:N_s)),'adm')
    indivs <- NULL
    for (pop in all.pops){
        ids <- read.table(sprintf('%s.het', pop), header=TRUE, as.is=TRUE)[,1]
        indivs <- rbind(indivs, data.frame(ids=ids, pop=pop, stringsAsFactors=FALSE))
    }
    mds <- suppressWarnings(cmdscale(as.matrix(d), k=2))
    for (s1 in 1:(N_s-1)){
      for (s2 in (s1+1):N_s){
        ids1 <- indivs$ids[indivs$pop == paste0('s', as.character(s1))]
        ids2 <- indivs$ids[indivs$pop == paste0('s', as.character(s2))]
        centroid1 <- colMeans(mds[ids1,])
        centroid2 <- colMeans(mds[ids2,])
        ids.adm <- indivs$ids[indivs$pop == 'adm']
        adm.angles <- rep(0, length(ids.adm))
        for (i in 1:length(adm.angles)){
            idx.adm <- which(rownames(mds) == ids.adm[i])
            vec.a <- mds[idx.adm,] - centroid1
            vec.b <- mds[idx.adm,] - centroid2
            cosine <- sum(vec.a * vec.b) / (sqrt(sum(vec.a*vec.a)) * sqrt(sum(vec.b*vec.b)))
            adm.angles[i] <- acos(cosine)
        }
        final.sumstats[paste0('mean.adm.angles.s', as.character(s1), '.s', as.character(s2))] <- mean(adm.angles)
        final.sumstats[paste0('var.adm.angles.s', as.character(s1), '.s', as.character(s2))] <- var(adm.angles)
        final.sumstats[paste0('skew.adm.angles.s', as.character(s1), '.s', as.character(s2))] <- skewness(adm.angles)
        final.sumstats[paste0('kurt.adm.angles.s', as.character(s1), '.s', as.character(s2))] <- kurtosis(adm.angles)
        tmp <- density(adm.angles)
        final.sumstats[paste0('mode.adm.angles.s', as.character(s1), '.s', as.character(s2))] <- tmp$x[which.max(tmp$y)][1]
        tmp <- sprintf('perc%d.adm.angles.s%s.s%s', seq(0,100,10), as.character(s1), as.character(s2))
        final.sumstats[tmp] <- quantile(adm.angles, seq(0,1,.1))
      }
    }
    final.sumstats
}


tryCatch(final.sumstats <- compute.adm.angles())#, error=function(e){invisible()}, warning=function(e){invisible()})

write.table(t(final.sumstats),args[1], quote=FALSE, row.names=FALSE)
