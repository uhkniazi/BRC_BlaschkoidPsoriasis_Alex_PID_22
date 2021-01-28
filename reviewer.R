# File: reviewer.R
# Auth: umar.niazi@kcl.ac.uk
# Desc: answers to various comments and questions by reviewers
# Date: 25/01/2021

source('header.R')

############################################################################
######### Rev 1 - Q2 - showing heatmaps of selected genes
############################################################################
load('reviewer_data/bp_counts_norm.rds')
cvSel = c("CARD14", "CARD6", "FUT2", "GJB2", "LCE3D", "PRSS53", "RPS6KA4", "STAT3")

mData = mData.norm.bp[cvSel,]
library(NMF)
library(RColorBrewer)
mData = log(mData+1)

aheatmap(mData, annRow = NA, scale = 'row', Rowv = T, 
         Colv=T, cexRow=1, cexCol = 1, col='-RdBu:5')

######### make line plots
library(lattice)

dfResults.raw = data.frame(t(mData)) 
df = stack(dfResults.raw)
df$f1 = factor(c('L', 'L', 'NL', 'NL'), levels = c('NL', 'L'))
df$f2 = factor(c('P2', 'P1', 'P2', 'P1'))

xyplot(values ~ f1 | ind, data=df, groups=f2, type=c('l', 'p'), scales=list(relation='free', x=list(cex=0.7), y=list(cex=0.7)), 
       ylab='log Expression', main=list(label='profile of Selected genes', cex=0.8),
       xlab='Condition', auto.key = list(columns=2))


######### perform a PCA of the data matrix with a plot
m = log(mData.norm.bp+1)
m = scale(m)
dim(m)
mPC = prcomp(t(m), scale = T)
plot(mPC$x[,1:2], pch=20, cex=2, main='First 2 principal components',
     xlim=c(-110, 110), ylim=c(-100, 100))
text(mPC$x[,1:2], labels=colnames(mData.norm.bp), pos=1)

############################################################################

############################################################################
######### Rev 1 - Q4 - Reporting DEGs at various cutoffs
############################################################################
df = read.csv('reviewer_data/LF_pvals.csv', header=T)
sapply(c('10%'=0.1, '5%'=0.05, '1%'=0.01), function(x) table(df$adj.P.Val < x))

############################################################################

############################################################################
######### Rev 1 - Reporting adjusted p-values for divergent genes
############################################################################
df = read.csv(file.choose(), header=T)
head(df)
df$adj.P.Val = p.adjust(df$pvalue, 'BH')
table(df$adj.P.Val < 0.1)
table(df$pvalue < 0.05)
write.csv(df, file='temp/divergent-geo41745.csv')
############################################################################

############################################################################
######### Rev 1 - Testing results with a new dataset
############################################################################
source('header.R')

# ## load the data
# library(RMySQL)
# 
# db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
# dbListTables(db)
# 
# # get the query
# g_did
# q = paste0('select MetaFile.* from MetaFile
#            where (MetaFile.idData = 43) AND (MetaFile.comment like "%count%")')
# dfSample = dbGetQuery(db, q)
# dfSample
# n = paste0(dfSample$location, dfSample$name)
# load(n[2])
# 
# ## load the metadata i.e. covariates
# q = paste0('select Sample.* from Sample where Sample.idData = 43')
# dfSample = dbGetQuery(db, q)
# dim(dfSample)
# dfSample
# # close connection after getting data
# dbDisconnect(db)
# 
# ## make count matrix
# names(lCounts)
# mCounts = do.call(cbind, lCounts)
# colnames(mCounts) = names(lCounts)
# 
# # sanity check
# identical(dfSample$id, as.integer(colnames(mCounts)))
# 
# mData = mCounts
# dim(mData)

## some EDA on raw data before merging replicates
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

# ## use combination of batch and biological source as identifier for technical replicates
# fReplicates = factor(dfSample$group1):factor(dfSample$group2)
# levels(fReplicates)
# dfSample$fReplicates = factor(fReplicates)
# # combine the technical replicates
# i = seq_along(1:ncol(mData))
# m = tapply(i, dfSample$fReplicates, function(x) {
#   return(x)
# })
# 
# mData = sapply(m, function(x){
#   return(rowSums(mCounts[,x]))
# })
# 
# # get a shorter version of dfSample after adding technical replicates
# dfSample.2 = dfSample[sapply(m, function(x) return(x[1])), ]
# identical(colnames(mData), as.character(dfSample.2$fReplicates))
# dim(dfSample.2)
# dfSample.2 = droplevels.data.frame(dfSample.2)
# 
# # drop the rows where average across rows is less than 3
# i = rowMeans(mData)
# table( i < 3)
# mData = mData[!(i< 3),]
# dim(mData)
# 
# ivProb = apply(mData, 1, function(inData) {
#   inData[is.na(inData) | !is.finite(inData)] = 0
#   inData = as.logical(inData)
#   lData = list('success'=sum(inData), fail=sum(!inData))
#   return(mean(rbeta(1000, lData$success + 0.5, lData$fail + 0.5)))
# })
# 
# hist(ivProb)
# 
# ## convert the gene names to symbols
# library(org.Hs.eg.db)
# df = AnnotationDbi::select(org.Hs.eg.db, keys = as.character(rownames(mData)), columns = 'SYMBOL', keytype = 'ENTREZID')
# i = match(rownames(mData), df$ENTREZID)
# df = df[i,]
# identical(rownames(mData), df$ENTREZID)
# rownames(mData) = df$SYMBOL
# 
# ############### load the second data set
# dfStudy.counts = read.csv('reviewer_data/GSE121212_readcount.txt', header=T, sep='\t',
#                           stringsAsFactors = F)
# f = duplicated(dfStudy.counts$X)
# table(f)
# # drop duplicate genes
# dfStudy.counts = dfStudy.counts[!f,]
# rownames(dfStudy.counts) = dfStudy.counts$X
# # keep only Psoriasis samples
# mData.study = as.matrix(dfStudy.counts[,-c(1)])
# i = grep('PSO', colnames(mData.study))
# mData.study = mData.study[,i]
# # create factors
# fPID = factor(gsub('PSO_(\\d+)_.+', '\\1', colnames(mData.study)))
# nlevels(fPID)
# fLesional = rep('L', times=length(fPID))
# fLesional[grep('non', x = colnames(mData.study))] = 'NL'
# fLesional = factor(fLesional, levels = c('NL', 'L'))
# dfSample.study = data.frame(fPID, fLesional)
# 
# # remove low count genes
# i = rowMeans(mData.study)
# table(i < 3)
# mData.study = mData.study[!(i < 3),]
# dim(mData.study)


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
fGroups = c(as.character(dfSample.2$group1), 
            as.character(dfSample.study$fLesional))
fGroups = factor(fGroups)
levels(fGroups)
## check matrix structure
oDiag.study = CDiagnosticPlots(log(mData.norm.2+1), 'normalised study')
l = CDiagnosticPlotsGetParameters(oDiag.study)
l$PCA.jitter = F
l$HC.jitter = F

oDiag.study = CDiagnosticPlotsSetParameters(oDiag.study, l)
fBatch = dfSample.study$fLesional

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
f = gsub('Lesional', 'L-BP', f)
f = gsub('Non-lesional', 'NL-BP', f)
df$fGroups = factor(f, levels = c('L', 'NL', 'L-BP', 'NL-BP'))
dotplot(values ~ fGroups | ind, data=df, panel=function(x, y, ...) panel.bwplot(x, y, pch='|',...), type='b',
        par.strip.text=list(cex=0.7), scales=list(relation='free', x=list(cex=0.7), y=list(cex=0.7)))

detach('package:org.Hs.eg.db', unload=T)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='nbResponsePartialPooling.stan')


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

dfResults = data.frame()
#pdf('temp/figs.pdf')
############################################ repeat this analysis in a loop for various genes
for (i in seq_along(cvGeneList)){
  ############### create data for input
  
  dfData = data.frame(y = mData.sub[cvGeneList[i],])
  f1 = c(dfSample.2$group1, as.character(dfSample.study$fLesional))
  f1[grep('^lesional', f1, ignore.case = T)] = 'L'
  f1[grep('^non', f1, ignore.case = T)] = 'NL'
  dfData$fTreatment = factor(f1)
  f2 = c(dfSample.2$group2, 
         as.character(dfSample.study$fPID))
  dfData$fPatient = factor(f2)
  f3 = c(rep(1, times=4), rep(2, times=55))
  dfData$batch = factor(f3)
  
  m1 = model.matrix(y ~ fTreatment:batch - 1, data=dfData)
  m2 = model.matrix(y ~ fPatient - 1, data=dfData)
  m = cbind(m1, m2)
  
  lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
                   NscaleBatches=4, NBatchMap=c(rep(1, times=2),
                                                rep(2, times=2),
                                                rep(3, times=28),
                                                rep(4, times=2)),
                   y=as.integer(dfData$y))
  
  fit.stan.1 = sampling(stanDso, data=lStanData, iter=1000, chains=2, pars=c('betas', 'populationMean', 'sigmaRan',
                                                                             'phi', 'mu'),
                        cores=2, control=list(adapt_delta=0.99, max_treedepth = 12))
  print(fit.stan.1, c('betas', 'populationMean', 'sigmaRan', 'phi'), digits=3)
  
  # ######## second data set
  # dfData = data.frame(y = mData.norm.2[cvGeneList[i],])
  # dfData$fTreatment = factor(dfSample.study$Factor.Value.phenotype.)
  # dfData$fPatient = factor(dfSample.study$Sample.Characteristic.individual.)
  # 
  # m1 = model.matrix(y ~ fTreatment - 1, data=dfData)
  # m2 = model.matrix(y ~ fPatient - 1, data=dfData)
  # m = cbind(m1, m2)
  # 
  # lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m,
  #                  NscaleBatches=2, NBatchMap=c(rep(1, times=nlevels(dfData$fTreatment)),
  #                                               rep(2, times=nlevels(dfData$fPatient))),
  #                  y=as.integer(dfData$y))
  # 
  # fit.stan.2 = sampling(stanDso, data=lStanData, iter=1000, chains=2, pars=c('betas', 'populationMean', 'sigmaRan',
  #                                                                            'phi', 'mu'),
  #                       cores=2, control=list(adapt_delta=0.99, max_treedepth = 12))
  # #print(fit.stan.2, c('betas', 'populationMean', 'sigmaRan', 'phi'), digits=3)
  # 
  mCoef1 = extract(fit.stan.1)$betas[,1:4]
  # mCoef2 = extract(fit.stan.2)$betas[,1:2]
  
  colnames(mCoef1) = gsub('fTreatment', '', colnames(m1))
  # colnames(mCoef2) = levels(factor(dfSample.study$Factor.Value.phenotype.))
  # 
  dif1 = getDifferenceVector(ivData = mCoef1[,'lesional:batch1'], ivBaseline = mCoef1[,'non-lesional:batch1'])
  dif2 = getDifferenceVector(ivData = mCoef1[,'lesional:batch2'], ivBaseline = mCoef1[,'non-lesional:batch2'])
  dif = getDifference(ivData = dif1, ivBaseline = dif2)
  r = data.frame(ind= cvGeneList[i], lVSnl=mean(dif1), 
                 psoVSnor=mean(dif2), difference=mean(dif1-dif2),
                 zscore=dif$z, pvalue=dif$p)
  # dfData = data.frame(blas=dif1, psor=dif2)
  # dfData = stack(dfData)
  # 
  dfResults = rbind(dfResults, r)
  # yl = unlist(tapply(dfData$values, INDEX = dfData$ind, quantile, prob=c(0.05, 0.95)))
  # print(bwplot(values ~ ind, data=dfData, panel=function(x, y, ...) panel.bwplot(x, y, pch='|',...), type='b',
  #        par.strip.text=list(cex=0.7), varwidth=F, do.out=F,
  #        ylim=c(min(yl)-0.5, max(yl)+0.5),
  #        main=paste('Log Fold difference - disease vs healthy skin', cvGeneList[i])))
}
#dev.off(dev.cur())
write.csv(dfResults, file='results/Alex_DEgenesAt10per_list_differences_in_dataset_E-GEOD-41745_3.csv')
library(org.Hs.eg.db)
df = AnnotationDbi::select(org.Hs.eg.db, keys = as.character(dfResults$ind), 
                           keytype = 'SYMBOL', columns = 'ENTREZID')
i = match(as.character(dfResults$ind), df$SYMBOL)
df = df[i,]
identical(as.character(dfResults$ind), df$SYMBOL)
dfResults$ENTREZID = df$ENTREZID
write.csv(dfResults, file='results/Alex_DEgenesAt10per_list_differences_in_dataset_E-GEOD-41745_3_.csv')

###################



############################################################################