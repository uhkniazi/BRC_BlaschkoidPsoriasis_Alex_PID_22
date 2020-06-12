# File: 18_linePlots.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: plotting results from the current and previous studies
# Date: 29/5/2020

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

## convert the gene names to symbols
library(org.Hs.eg.db)
df = AnnotationDbi::select(org.Hs.eg.db, keys = as.character(rownames(mData)), columns = 'SYMBOL', keytype = 'ENTREZID')
i = match(rownames(mData), df$ENTREZID)
df = df[i,]
identical(rownames(mData), df$ENTREZID)
rownames(mData) = df$SYMBOL

############### load the second data set
dfStudy.meta = read.csv('dataExternal/E-GEOD-41745-experiment-design.tsv', header=T, sep='\t')
dfStudy.counts = read.csv('dataExternal/E-GEOD-41745-raw-counts.tsv', header=T, sep='\t')
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
fStudy = factor(c(rep('Bl', times=4), rep('s417', times=6)))
# sanity check
table(fGroups, fStudy)

cvGeneList = scan(what=character())
mData.sub = mData.merged[rownames(mData.merged) %in% cvGeneList, ]
dim(mData.sub)
cvGeneList = rownames(mData.sub)
library(lattice)
# df = data.frame(t(log(mData.sub+1)))
# df = stack(df)
# levels(fGroups)
# f = as.character(fGroups)
# f = gsub('^Lesional', 'L', f, ignore.case = T)
# f = gsub('Non-lesional', 'NL', f, ignore.case = T)
# 
# df$fGroups = factor(f, levels = c('L', 'NL'))
# df$fStudy = fStudy
# 
# dotplot(values ~ fGroups | ind, data=df, groups=df$fStudy, 
#         panel=function(x, y, ...) panel.bwplot(x, y, pch='|',...), type='b',
#         par.strip.text=list(cex=0.7), scales=list(relation='free', x=list(cex=0.7), y=list(cex=0.7)))

detach('package:org.Hs.eg.db', unload=T)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='nbResponsePartialPooling.stan')

dfResults = data.frame(row.names = 1:4)
dfResults.raw = data.frame(row.names=1:4)

############################################ repeat this analysis in a loop for various genes
for (i in seq_along(cvGeneList)){
  ############### create data for input
  
  dfData = data.frame(y = mData.sub[cvGeneList[i],])
  f1 = c(dfSample.2$group1, as.character(dfSample.study$Factor.Value.phenotype.))
  f1[grep('^lesional', f1, ignore.case = T)] = 'lesional'
  f1[grep('^non', f1, ignore.case = T)] = 'non-lesional'
  dfData$fTreatment = factor(f1)
  f2 = c(dfSample.2$group2, as.character(dfSample.study$Sample.Characteristic.individual.))
  dfData$fPatient = factor(f2)
  f3 = c(rep(1, times=4), rep(2, times=6))
  dfData$batch = factor(f3)
  
  m1 = model.matrix(y ~ fTreatment:batch - 1, data=dfData)
  m2 = model.matrix(y ~ fPatient - 1, data=dfData)
  m = cbind(m1, m2)
  
  lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
                   NscaleBatches=4, NBatchMap=c(rep(1, times=2),
                                                rep(2, times=2),
                                                rep(3, times=2),
                                                rep(4, times=3)),
                   y=as.integer(dfData$y))
  
  fit.stan.1 = sampling(stanDso, data=lStanData, iter=1000, chains=2, pars=c('betas', 'populationMean', 'sigmaRan',
                                                                             'phi', 'mu'),
                        cores=2, control=list(adapt_delta=0.99, max_treedepth = 12))
  print(fit.stan.1, c('betas', 'populationMean', 'sigmaRan', 'phi'), digits=3)
  
   
  mCoef1 = extract(fit.stan.1)$betas[,1:4]
  mCoef1 = sweep(mCoef1, 1, extract(fit.stan.1)$populationMean, '+')
  # mCoef2 = extract(fit.stan.2)$betas[,1:2]
  
  colnames(mCoef1) = gsub('fTreatment', '', colnames(m1))
  r = data.frame(colMeans(mCoef1))
  colnames(r) = cvGeneList[i]
  r2 = data.frame(tapply(log(lStanData$y+1), dfData$fTreatment:dfData$batch, mean))
  colnames(r2) = cvGeneList[i]
   
  dfResults = cbind(dfResults, r)
  rownames(dfResults) = rownames(r)
  dfResults.raw = cbind(dfResults.raw, r2)
  rownames(dfResults.raw) = rownames(r2)
  
}

pdf('temp/figs.pdf')
df = stack(dfResults)
df$f1 = factor(c('L', 'NL', 'L', 'NL'))
df$f2 = factor(c('bl', 'bl', 's417', 's417'))

xyplot(values ~ f1 | ind, data=df, groups=f2, type=c('l', 'p'), scales=list(relation='free', x=list(cex=0.7), y=list(cex=0.7)), 
       ylab='Model Estimated log Average', main=list(label='profile of Selected genes', cex=0.8),
       xlab='Condition', auto.key = list(columns=2))

df = stack(dfResults.raw)
df$f1 = factor(c('L', 'L', 'NL', 'NL'))
df$f2 = factor(c('bl', 's417', 'bl', 's417'))

xyplot(values ~ f1 | ind, data=df, groups=f2, type=c('l', 'p'), scales=list(relation='free', x=list(cex=0.7), y=list(cex=0.7)), 
       ylab='log Average', main=list(label='profile of Selected genes', cex=0.8),
       xlab='Condition', auto.key = list(columns=2))


dev.off(dev.cur())
###################
