# File: 17_integration.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: merging the data from current and previous studies, look at other branches
#       where integration was performed, we merge the result tables from those
#       2 branches called dataset_**
# Date: 02/04/2020

source('header.R')

# load the datasets
dfCombined = read.csv('results/divergent_GSE121212_GSE41745_combinedModel.xls', 
                 header=T, stringsAsFactors = F)
dfCombined = dfCombined[,-1]
rownames(dfCombined) = dfCombined$ENTREZID

df544 = read.csv('results/Alex_DEgenesAt10per_list_differences_in_dataset_E-GEOD-54456_2.csv',
                 header=T, stringsAsFactors = F)

df544 = df544[,-1]
rownames(df544) = df544$ENTREZID

dfBP = read.csv('results/Transcriptomics of Blaschkoid psoriasis - REVERSE - Alex - Project ID 22.csv',
                header=T, stringsAsFactors = F)

rownames(dfBP) = dfBP$ENTREZ.ID

df544 = df544[rownames(dfCombined),]
dfBP = dfBP[rownames(dfCombined), ]
identical(rownames(dfBP), rownames(dfCombined))

### correlations between results
plot(dfCombined$BP_logFC_lVSnl, dfBP$logFC, pch=20)
cor(dfCombined$BP_logFC_lVSnl, dfBP$logFC)
d = dfCombined$BP_logFC_lVSnl - dfBP$logFC
summary(round(d, 2))
hist(d)

plot(dfCombined$GSE41745_logFC, dfCombined$GSE121212_logFC, pch=20)
cor(dfCombined$GSE41745_logFC, dfCombined$GSE121212_logFC)

dim(dfCombined)
dim(df544)

## genes at various cutoffs
sapply(c('10%'=0.1, '5%'=0.05, '1%'=0.01), function(x) table(dfCombined$GSE41745_adj.P.Val < x))
sapply(c('10%'=0.1, '5%'=0.05, '1%'=0.01), function(x) table(dfCombined$GSE121212_adj.P.Val < x))

## create one matrix
mData = cbind(BP=dfCombined$BP_logFC_lVSnl,
              GSE41745=dfCombined$GSE41745_logFC,
              GSE54456=df544$psoVSnor,
              GSE121212=dfCombined$GSE121212_logFC)
rownames(mData) = dfCombined$ind

## make appropriate figures
range(mData)
quantile(as.vector(mData), 0:20/20)
mData.p = mData

mData.p[mData < -2.5] = -2.5
mData.p[mData > 3.5] = 3.5
pairs(mData.p, pch=20, las=0.5)

library(NMF)
library(RColorBrewer)
library(amap)

d = Dist(t(mData), method='correlation')
hc = hclust(d)
plot(hc)

aheatmap(mData.p, annRow = NA, scale = 'none', Rowv = T, breaks=0,
         Colv=hc, cexRow=5, cexCol = 1, distfun = 'correlation', 
         col=rev(brewer.pal(9, 'RdBu')))

cvSig = dfCombined$ind[dfCombined$GSE121212_adj.P.Val < 0.01]
length(cvSig)

hc = hclust(Dist(t(mData[cvSig,]), method = 'correlation'))
plot(hc)
aheatmap(mData[cvSig,], annRow = NA, scale = 'none', Rowv = T, breaks=0,
         Colv=hc, cexRow=5, cexCol = 1, distfun = 'correlation', 
         col=rev(brewer.pal(9, 'RdBu')))

m = cor(mData)
aheatmap(m, annRow = NA, scale = 'none', Rowv = NA, breaks=0.7,
         Colv=NA, cexRow=1, cexCol = 1,  
         col=(brewer.pal(9, 'RdYlBu')))


## select a list of genes here
i = order(dfBP$adj.P.Val, decreasing = F)
cvSig = dfBP$SYMBOL[i[1:14]]
#cvSig = scan(what=character())
length(cvSig)
cvSig = cvSig[cvSig %in% rownames(mData)]

hc = hclust(Dist(t(mData[cvSig,]), method = 'correlation'))
plot(hc)

aheatmap(mData[cvSig,], annRow = NA, scale = 'none', Rowv = T, #breaks=0,
         Colv=hc, cexRow=1, cexCol = 1, distfun = 'correlation', 
         col=rev(brewer.pal(9, 'RdBu')))

m = cor(mData[cvSig, ])
aheatmap(m, annRow = NA, scale = 'none', Rowv = NA, breaks=0.7,
         Colv=NA, cexRow=1, cexCol = 1,  
         col=(brewer.pal(9, 'RdYlBu')))

# 
# ## genes different in both comparisons
# cvSig.121 = dfCombined$ind[dfCombined$GSE121212_adj.P.Val < 0.1]
# cvSig.417 = dfCombined$ind[dfCombined$GSE41745_adj.P.Val < 0.1]
# cvSig = Reduce(intersect, list(cvSig.121, cvSig.417))
# 
# round(cor(mData[cvSig,]), 2)
# aheatmap(mData[cvSig,], annRow = NA, scale = 'none', Rowv = T, 
#          Colv=T, cexRow=5, cexCol = 1, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
#          col=c('white', brewer.pal(9, 'YlOrRd')))


### cluster using a subset of genes
library(amap)
d = Dist(t(mData[cvSig,]), method='correlation')
hc = hclust(d)
plot(hc, main='clustering of bl', sub='')
# hcr = hclust(Dist(mData[cvSig,], method='correlation'))
# aheatmap(mData[cvSig,], annRow = NA, scale = 'none', Rowv = hcr, 
#          Colv=hc, cexRow=5, cexCol = 1, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
#          col=c('white', brewer.pal(9, 'YlOrRd')))

round(d, 3)

## returns TRUE if BP shares a cluster with any other
## member otherwise FALSE
getCluster = function(m){
   d = Dist(t(m), method='correlation')
   c = cutree(hclust(d), k = 2)
   i = which(names(c) == 'BP')
   return(any(c[i] == c[-i]))
}

# getDistance = function(m){
#   d = as.matrix(Dist(t(m), method='correlation'))
#   d = max(d['BP',])
#   return(d)
# }


simulateOne = function(){
  c = sample(rownames(mData), size = length(cvSig), replace = F)
  return(getCluster(mData[c,]))# < 0)
}

set.seed(123)
fCluster = replicate(1000, simulateOne())
table(fCluster)
sum(!fCluster)/length(fCluster)
sum(fCluster)/length(fCluster)

## See bayesian computations with r page 24
b = rbeta(1000, 2, 8)
hist(b)
1 - pbeta(0.5, 2, 8)
sum(b >= 0.5)/length(b)




# ivDist = replicate(1000, simulateOne())
# iObs = getDistance(mData[cvSig,])
# ## calculate bayesian p-value for this test statistic
# getPValue = function(Trep, Tobs){
#   left = sum(Trep <= Tobs)/length(Trep)
#   right = sum(Trep >= Tobs)/length(Trep)
#   return(min(left, right))
# }
# 
# hist(ivDist, main='Simulated distances based on resampling', 
#      xlab='', xlim=c(0, 0.6), xaxt='n')
# axis(1, at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6))
# points(iObs, 0, pch=20, col=2, cex=2)
# getPValue(ivDist, iObs)

################################################
## check this with a selected  list of genes
################################################

cvSel = scan(what=character())
cvSel = unique(cvSel)
length(cvSel)
cvSel = cvSel[cvSel %in% rownames(mData)]
length(cvSel)
aheatmap(mData[cvSel,], annRow = NA, scale = 'none', Rowv = T, 
         Colv=T, cexRow=5, cexCol = 1, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
         col=c('white', brewer.pal(9, 'YlOrRd')))

dfSel = cbind(df417[,c(2, 3, 4, 8)], df544[, c(4, 9)])
dfSel = dfSel[, -c(4)]
colnames(dfSel) = c('symbol', 'bl', 's417', 's544', 'entrezid')
dfSel = dfSel[dfSel$symbol %in% cvSel,]
write.csv(dfSel, file='results/alex_gwas_common.xls')

