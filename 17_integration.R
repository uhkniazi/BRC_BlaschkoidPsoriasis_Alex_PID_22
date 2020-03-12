# File: 17_integration.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: merging the data from current and previous studies
# Date: 11/03/2020

source('header.R')

## load the data
library(RMySQL)

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# get the query
g_did
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 43) AND (MetaFile.comment like "%count%")')
dfSample = dbGetQuery(db, q)
dfSample
n = paste0(dfSample$location, dfSample$name)
load(n[2])

## load the metadata i.e. covariates
q = paste0('select Sample.* from Sample where Sample.idData = 43')
dfSample = dbGetQuery(db, q)
dim(dfSample)
dfSample
# close connection after getting data
dbDisconnect(db)

## make count matrix
names(lCounts)
mCounts = do.call(cbind, lCounts)
colnames(mCounts) = names(lCounts)

# sanity check
identical(dfSample$id, as.integer(colnames(mCounts)))

mData = mCounts
dim(mData)

## some EDA on raw data before merging replicates
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

## use combination of batch and biological source as identifier for technical replicates
fReplicates = factor(dfSample$group1):factor(dfSample$group2)
levels(fReplicates)
dfSample$fReplicates = factor(fReplicates)
# combine the technical replicates
i = seq_along(1:ncol(mData))
m = tapply(i, dfSample$fReplicates, function(x) {
  return(x)
})

mData = sapply(m, function(x){
  return(rowSums(mCounts[,x]))
})

# get a shorter version of dfSample after adding technical replicates
dfSample.2 = dfSample[sapply(m, function(x) return(x[1])), ]
identical(colnames(mData), as.character(dfSample.2$fReplicates))
dim(dfSample.2)
dfSample.2 = droplevels.data.frame(dfSample.2)

# drop the rows where average across rows is less than 3
i = rowMeans(mData)
table( i < 3)
mData = mData[!(i< 3),]
dim(mData)

ivProb = apply(mData, 1, function(inData) {
  inData[is.na(inData) | !is.finite(inData)] = 0
  inData = as.logical(inData)
  lData = list('success'=sum(inData), fail=sum(!inData))
  return(mean(rbeta(1000, lData$success + 0.5, lData$fail + 0.5)))
})

hist(ivProb)

## convert the gene names to symbols
library(org.Hs.eg.db)
df = AnnotationDbi::select(org.Hs.eg.db, keys = as.character(rownames(mData)), columns = 'SYMBOL', keytype = 'ENTREZID')
i = match(rownames(mData), df$ENTREZID)
df = df[i,]
identical(rownames(mData), df$ENTREZID)
rownames(mData) = df$SYMBOL

############### load the second data set
dfStudy.meta = read.csv('dataExternal/E-GEOD-41745/E-GEOD-41745-experiment-design.tsv', header=T, sep='\t')
dfStudy.counts = read.csv('dataExternal/E-GEOD-41745/E-GEOD-41745-raw-counts.tsv', header=T, sep='\t')
mData.study = as.matrix(dfStudy.counts[,-c(1,2)])
rownames(mData.study) = as.character(dfStudy.counts$Gene.Name)

# remove low count genes
i = rowMeans(mData.study)
table(i < 3)
mData.study = mData.study[!(i < 3),]
dim(mData.study)

## put samples in same order in both metadata and count tables
dfSample.study = dfStudy.meta#[,c(1:4)]
head(dfSample.study)
dfSample.study$Run = as.character(dfSample.study$Run)
str(dfSample.study)
table(dfSample.study$Run %in% colnames(mData.study))
mData.study = mData.study[,dfSample.study$Run]
identical(dfSample.study$Run, colnames(mData.study))

## match the genes in the two data sets
table(rownames(mData) %in% rownames(mData.study))
mData = mData[rownames(mData) %in% rownames(mData.study), ]
table(rownames(mData) %in% rownames(mData.study))
length(unique(rownames(mData.study))); dim(mData.study)
mData.study = mData.study[rownames(mData), ]
dim(mData.study); dim(mData)
identical(rownames(mData), rownames(mData.study))

## normalise the 2 data sets separately
library(DESeq2)
sf = estimateSizeFactorsForMatrix(mData)
sf2 = estimateSizeFactorsForMatrix(mData.study)

mData.norm.1 = sweep(mData, 2, sf, '/')
mData.norm.2 = sweep(mData.study, 2, sf2, '/')

################## merge the two samples 
mData.merged = cbind(mData.norm.1, mData.norm.2)
fGroups = c(as.character(dfSample.2$group1), as.character(dfSample.study$Factor.Value.phenotype.))
fGroups = factor(fGroups)
levels(fGroups)
## check matrix structure
oDiag.study = CDiagnosticPlots(log(mData.norm.2+1), 'normalised study')
l = CDiagnosticPlotsGetParameters(oDiag.study)
l$PCA.jitter = F
l$HC.jitter = F

oDiag.study = CDiagnosticPlotsSetParameters(oDiag.study, l)
fBatch = dfSample.study$Factor.Value.phenotype.

par(mfrow=c(1,1))
boxplot.median.summary(oDiag.study, fBatch, legend.pos = 'topright', axis.label.cex = 0.5)
plot.mean.summary(oDiag.study, fBatch, axis.label.cex = 0.5)
plot.sigma.summary(oDiag.study, fBatch, axis.label.cex = 0.5)
plot.missing.summary(oDiag.study, fBatch, axis.label.cex = 0.5, cex.main=1)
plot.PCA(oDiag.study, fBatch, csLabels = '')
plot.dendogram(oDiag.study, fBatch, labels_cex = 0.5)

## choose some house keeping genes and plot
cvGeneList = c('UBC', 'SDHA', 'CYC1', 'CANX', 'TBP',
               'YWHAZ', 'B2M', 'RPLP2', 'GAPDH')

cvGeneList = scan(what=character())
mData.sub = mData.merged[rownames(mData.merged) %in% cvGeneList, ]
dim(mData.sub)
cvGeneList = rownames(mData.sub)
## plot these genes
library(lattice)
df = data.frame(t(log(mData.sub+1)))
df = stack(df)
levels(fGroups)
f = as.character(fGroups)
f = gsub('Lesional', 'L', f)
f = gsub('Non-lesional', 'NL', f)
f = gsub('normal', 'nor', f)
f = gsub('psoriasis', 'ps', f)
df$fGroups = factor(f, levels = c('L', 'NL', 'ps', 'nor'))
dotplot(values ~ fGroups | ind, data=df, panel=function(x, y, ...) panel.bwplot(x, y, pch='|',...), type='b',
        par.strip.text=list(cex=0.7), scales=list(relation='free', x=list(cex=0.7), y=list(cex=0.7)))

detach('package:org.Hs.eg.db', unload=T)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='nbResponsePartialPooling.stan')

############### create data for input
# dfData = data.frame(y = mData.norm.1[cvGeneList[1],])
# str(dfData)
# str(dfSample.2)
# dfData$fTreatment = factor(dfSample.2$group1)
# dfData$fPatient = factor(dfSample.2$group2)
# 
# str(dfData)
# m1 = model.matrix(y ~ fTreatment - 1, data=dfData)
# m2 = model.matrix(y ~ fPatient - 1, data=dfData)
# m = cbind(m1, m2)
# 
# lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
#                  NscaleBatches=2, NBatchMap=c(rep(1, times=nlevels(dfData$fTreatment)),
#                                               rep(2, times=nlevels(dfData$fPatient))),
#                  y=as.integer(dfData$y))
# 
# fit.stan.1 = sampling(stanDso, data=lStanData, iter=1000, chains=2, pars=c('betas', 'populationMean', 'sigmaRan',
#                                                                             'phi', 'mu'),
#                       cores=2, control=list(adapt_delta=0.99, max_treedepth = 12))
# print(fit.stan.1, c('betas', 'populationMean', 'sigmaRan', 'phi'), digits=3)
# 
# traceplot(fit.stan.1, 'populationMean')
# traceplot(fit.stan.1, 'betas')
# traceplot(fit.stan.1, 'sigmaRan')
# 
# ######## second data set
# dfData = data.frame(y = mData.norm.2[cvGeneList[1],])
# str(dfData)
# str(dfSample.study)
# dfData$fTreatment = factor(dfSample.study$Sample.Characteristic.disease.)
# 
# str(dfData)
# m = model.matrix(y ~ fTreatment - 1, data=dfData)
# 
# lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
#                  NscaleBatches=1, NBatchMap=c(rep(1, times=nlevels(dfData$fTreatment))),
#                  y=as.integer(dfData$y))
# 
# fit.stan.2 = sampling(stanDso, data=lStanData, iter=1000, chains=2, pars=c('betas', 'populationMean', 'sigmaRan',
#                                                                            'phi', 'mu'),
#                       cores=2, control=list(adapt_delta=0.99, max_treedepth = 12))
# print(fit.stan.2, c('betas', 'populationMean', 'sigmaRan', 'phi'), digits=3)

## function to calculate statistics for differences between coefficients
getDifference = function(ivData, ivBaseline){
  stopifnot(length(ivData) == length(ivBaseline))
  # get the difference vector
  d = ivData - ivBaseline
  # get the z value
  z = mean(d)/sd(d)
  # get 2 sided p-value
  p = pnorm(-abs(mean(d)/sd(d)))*2
  return(list(z=z, p=p))
}

getDifferenceVector = function(ivData, ivBaseline){
  stopifnot(length(ivData) == length(ivBaseline))
  # get the difference vector
  d = ivData - ivBaseline
  return(d)
}

# mCoef1 = extract(fit.stan.1)$betas[,1:2]
# mCoef2 = extract(fit.stan.2)$betas
# 
# colnames(mCoef1) = levels(factor(dfSample.2$group1))
# colnames(mCoef2) = levels(factor(dfSample.study$Sample.Characteristic.disease.))
# 
# dif1 = getDifferenceVector(ivData = mCoef1[,'Lesional'], ivBaseline = mCoef1[,'Non-lesional'])
# dif2 = getDifferenceVector(ivData = mCoef2[,'psoriasis'], ivBaseline = mCoef2[,'normal'])
# dif = getDifference(ivData = dif1, ivBaseline = dif2)
# r = data.frame(ind= cvGeneList[1], lVSnl=mean(dif1), 
#                  psoVSnor=mean(dif2), difference=mean(dif1-dif2),
#                  zscore=dif$z, pvalue=dif$p)
# dfData = data.frame(blas=dif1, psor=dif2)
# dfData = stack(dfData)
# 
# bwplot(values ~ ind, data=dfData, panel=panel.violin, type='b',
#        par.strip.text=list(cex=0.7), varwidth=F, main='Log Fold difference - disease vs healthy skin')
# 
# bwplot(values ~ ind, data=dfData, panel=function(x, y, ...) panel.bwplot(x, y, pch='|',...), type='b',
#        par.strip.text=list(cex=0.7), varwidth=T, main='Log Fold difference - disease vs healthy skin')


dfResults = data.frame()
pdf('temp/figs.pdf')
############################################ repeat this analysis in a loop for various genes
for (i in seq_along(cvGeneList)){
  ############### create data for input
  dfData = data.frame(y = mData.norm.1[cvGeneList[i],])
  dfData$fTreatment = factor(dfSample.2$group1)
  dfData$fPatient = factor(dfSample.2$group2)
  
  m1 = model.matrix(y ~ fTreatment - 1, data=dfData)
  m2 = model.matrix(y ~ fPatient - 1, data=dfData)
  m = cbind(m1, m2)
  
  lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
                   NscaleBatches=2, NBatchMap=c(rep(1, times=nlevels(dfData$fTreatment)),
                                                rep(2, times=nlevels(dfData$fPatient))),
                   y=as.integer(dfData$y))
  
  fit.stan.1 = sampling(stanDso, data=lStanData, iter=1000, chains=2, pars=c('betas', 'populationMean', 'sigmaRan',
                                                                             'phi', 'mu'),
                        cores=2, control=list(adapt_delta=0.99, max_treedepth = 12))
  print(fit.stan.1, c('betas', 'populationMean', 'sigmaRan', 'phi'), digits=3)
  
  ######## second data set
  dfData = data.frame(y = mData.norm.2[cvGeneList[i],])
  dfData$fTreatment = factor(dfSample.study$Factor.Value.phenotype.)
  dfData$fPatient = factor(dfSample.study$Sample.Characteristic.individual.)
  
  m1 = model.matrix(y ~ fTreatment - 1, data=dfData)
  m2 = model.matrix(y ~ fPatient - 1, data=dfData)
  m = cbind(m1, m2)
  
  lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
                   NscaleBatches=2, NBatchMap=c(rep(1, times=nlevels(dfData$fTreatment)),
                                                rep(2, times=nlevels(dfData$fPatient))),
                   y=as.integer(dfData$y))
  
  fit.stan.2 = sampling(stanDso, data=lStanData, iter=1000, chains=2, pars=c('betas', 'populationMean', 'sigmaRan',
                                                                             'phi', 'mu'),
                        cores=2, control=list(adapt_delta=0.99, max_treedepth = 12))
  print(fit.stan.2, c('betas', 'populationMean', 'sigmaRan', 'phi'), digits=3)
  
  mCoef1 = extract(fit.stan.1)$betas[,1:2]
  mCoef2 = extract(fit.stan.2)$betas[,1:2]
  
  colnames(mCoef1) = levels(factor(dfSample.2$group1))
  colnames(mCoef2) = levels(factor(dfSample.study$Factor.Value.phenotype.))
  
  dif1 = getDifferenceVector(ivData = mCoef1[,'Lesional'], ivBaseline = mCoef1[,'Non-lesional'])
  dif2 = getDifferenceVector(ivData = mCoef2[,'lesional'], ivBaseline = mCoef2[,'non-lesional'])
  dif = getDifference(ivData = dif1, ivBaseline = dif2)
  r = data.frame(ind= cvGeneList[i], lVSnl=mean(dif1), 
                 psoVSnor=mean(dif2), difference=mean(dif1-dif2),
                 zscore=dif$z, pvalue=dif$p)
  dfData = data.frame(blas=dif1, psor=dif2)
  dfData = stack(dfData)
  
  dfResults = rbind(dfResults, r)
  yl = unlist(tapply(dfData$values, INDEX = dfData$ind, quantile, prob=c(0.05, 0.95)))
  print(bwplot(values ~ ind, data=dfData, panel=function(x, y, ...) panel.bwplot(x, y, pch='|',...), type='b',
         par.strip.text=list(cex=0.7), varwidth=F, do.out=F,
         ylim=c(min(yl)-0.5, max(yl)+0.5),
         main=paste('Log Fold difference - disease vs healthy skin', cvGeneList[i])))
}
dev.off(dev.cur())
###################

