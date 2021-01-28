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