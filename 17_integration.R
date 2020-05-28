# File: 17_integration.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: merging the data from current and previous studies, look at other branches
#       where integration was performed, we merge the result tables from those
#       2 branches called dataset_**
# Date: 02/04/2020

source('header.R')

# load the first dataset
df417 = read.csv('results/Alex_DEgenesAt10per_list_differences_in_dataset_E-GEOD-41745_3_.csv', 
                 header=T, stringsAsFactors = F)
df544 = read.csv('results/Alex_DEgenesAt10per_list_differences_in_dataset_E-GEOD-54456_2.csv',
                 header=T, stringsAsFactors = F)

dim(df417)
dim(df544)

table(df417$ind %in% df544$ind)
i = match(df544$ind, df417$ind)
df417 = df417[i,]
identical(df417$ind, df544$ind)

## create one matrix
mData = cbind(bl=df417$lVSnl, s417=df417$psoVSnor, s544=df544$psoVSnor)
rownames(mData) = df417$ind

## make appropriate figures
range(mData)
quantile(as.vector(mData), 0:20/20)

# mData[mData < -2] = -2
# mData[mData > 3] = 3

library(NMF)
library(RColorBrewer)

aheatmap(mData, annRow = NA, scale = 'none', Rowv = T, 
         Colv=T, cexRow=5, cexCol = 1, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
         col=c('white', brewer.pal(9, 'YlOrRd')))

cvSig = df417$ind[df417$pvalue < 0.01]
length(cvSig)
aheatmap(mData[cvSig,], annRow = NA, scale = 'none', Rowv = T, 
         Colv=T, cexRow=5, cexCol = 1, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
         col=c('white', brewer.pal(9, 'YlOrRd')))

round(cor(mData),2)
pdf('temp/heatmaps.pdf')

aheatmap(mData, annRow = NA, scale = 'none', Rowv = T, 
         Colv=T, cexRow=5, cexCol = 1, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
         col=c('white', brewer.pal(9, 'YlOrRd')))

cvSig = df417$ind[df417$pvalue < 0.1]
length(cvSig)
aheatmap(mData[cvSig,], annRow = NA, scale = 'none', Rowv = T, 
         Colv=T, cexRow=2, cexCol = 2, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
         col=c('white', brewer.pal(9, 'YlOrRd')))
dev.off(dev.cur())

cvSig = df417$ind[df417$pvalue < 0.001]
length(cvSig)
aheatmap(mData[cvSig,], annRow = NA, scale = 'none', Rowv = T, 
         Colv=T, cexRow=5, cexCol = 1, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
         col=c('white', brewer.pal(9, 'YlOrRd')))

## select a list of genes here
cvSig = scan(what=character())
length(cvSig)
cvSig = cvSig[cvSig %in% rownames(mData)]
round(cor(mData[cvSig,]), 2)
aheatmap(mData[cvSig,], annRow = NA, scale = 'none', Rowv = T, 
         Colv=T, cexRow=5, cexCol = 1, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
         col=c('white', brewer.pal(9, 'YlOrRd')))

### cluster using a subset of genes
hc = hclust(dist(t(mData[cvSig,])))
plot(hc, main='clustering of bl', sub='')
round(as.matrix(dist(t(mData[cvSig,]))),3)

getDistance = function(m){
   d = as.matrix(dist(t(m)))
   d1 = (d['bl', 's544'])
   d2 = (d['s417', 's544'])
   return(d1 - d2)
}

simulateOne = function(){
  c = sample(rownames(mData), size = length(cvSig), replace = F)
  return(getDistance(mData[c,]))# < 0)
}

ivDist = replicate(1000, simulateOne())
iObs = getDistance(mData[cvSig,])
## calculate bayesian p-value for this test statistic
getPValue = function(Trep, Tobs){
  left = sum(Trep <= Tobs)/length(Trep)
  right = sum(Trep >= Tobs)/length(Trep)
  return(min(left, right))
}

hist(ivDist, main='simulated difference of distance', 
     xlab='(bl:s544) - (s417:544)', xlim=c(-3, 5))
points(iObs, 0, pch=20, col=2, cex=2)
getPValue(ivDist, iObs)

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

### make some line plots
library(lattice)

df = 