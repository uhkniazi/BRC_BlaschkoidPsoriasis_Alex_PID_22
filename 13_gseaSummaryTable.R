# File: 13_gseaSummaryTable.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: merge the gsea results for all contrasts in one table
# Date: 18/12/2019


lFiles = list.files('results/', pattern='*pathways_mSigDb_c2_*', full.names = T, ignore.case = F)

# load the files
ldfData = lapply(lFiles, function(x) read.csv(x, row.names=1))
names(ldfData) = lFiles
sapply(ldfData, nrow)
# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

# get the upregulated/downregulated in order
ldfData.up = ldfData[grepl(lFiles, pattern = 'upregulated')]
ldfData.down = ldfData[grepl(lFiles, pattern = 'downregulated')]

## set the names for each contrast
sn = gsub('results//(.+Vs.+)_up.+', '\\1', names(ldfData.up))
names(ldfData.up) = sn

sn = gsub('results//(.+Vs.+)_down.+', '\\1', names(ldfData.down))
names(ldfData.down) = sn

## create a table/matrix of p-values
mMerged.up = sapply(ldfData.up, function(x) x$p.val)
rownames(mMerged.up) = rownames(ldfData.up[[1]])
## create similar table for downregulated
mMerged.down = sapply(ldfData.down, function(x) x$p.val)
rownames(mMerged.down) = rownames(ldfData.down[[1]])

# sanity check
identical(rownames(mMerged.up), rownames(mMerged.down))
colnames(mMerged.up) = paste(colnames(mMerged.up), 'up', sep='-')
colnames(mMerged.down) = paste(colnames(mMerged.down), 'down', sep='-')

mMerged.c2 = cbind(mMerged.up, mMerged.down)
# reorder the columns
colnames(mMerged.c2)
o = c(1, 2)
## sanity check
matrix(colnames(mMerged.c2)[o], ncol = 2, byrow = T)
colnames(mMerged.c2)[o]
mMerged.c2 = mMerged.c2[,o]

# remove na sections
dim(mMerged.c2)
mMerged.c2 = na.omit(mMerged.c2)
dim(mMerged.c2)
head(mMerged.c2)

### create a binary matrix based on cutoffs
getBinaryMatrix = function(mat, cutoff=0.01){
  mat2 = apply(mat, 2, function(x) round(x, 3) <= cutoff)
}

mMerged.c2.bin = getBinaryMatrix(mMerged.c2)

## group this matrix into combinations
mMerged.c2.bin.grp = mMerged.c2.bin
set.seed(123)
dm = dist(mMerged.c2.bin.grp, method='binary')
hc = hclust(dm)
plot(hc, labels=F)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mMerged.c2.bin.grp = cbind(mMerged.c2.bin.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mMerged.c2.bin.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## map these names to the cp
groups = cp
sig.pvals = rowSums(mMerged.c2.bin)
dfMerged.c2 = data.frame(round(mMerged.c2, 3), sig.pvals, groups, DB='C2')
str(dfMerged.c2)
head(dfMerged.c2)
tail(dfMerged.c2)

### repeat for the C5 database
lFiles = list.files('results/', pattern='*pathways_mSigDb_c5.xls', full.names = T, ignore.case = F)

# load the files
ldfData = lapply(lFiles, function(x) read.csv(x, row.names=1))
names(ldfData) = lFiles
sapply(ldfData, nrow)
# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

# get the upregulated/downregulated in order
ldfData.up = ldfData[grepl(lFiles, pattern = 'upregulated')]
ldfData.down = ldfData[grepl(lFiles, pattern = 'downregulated')]

## set the names for each contrast
sn = gsub('results//(.+Vs.+)_up.+', '\\1', names(ldfData.up))
names(ldfData.up) = sn

sn = gsub('results//(.+Vs.+)_down.+', '\\1', names(ldfData.down))
names(ldfData.down) = sn

## create a table/matrix of p-values
mMerged.up = sapply(ldfData.up, function(x) x$p.val)
rownames(mMerged.up) = rownames(ldfData.up[[1]])
## create similar table for downregulated
mMerged.down = sapply(ldfData.down, function(x) x$p.val)
rownames(mMerged.down) = rownames(ldfData.down[[1]])

# sanity check
identical(rownames(mMerged.up), rownames(mMerged.down))
colnames(mMerged.up) = paste(colnames(mMerged.up), 'up', sep='-')
colnames(mMerged.down) = paste(colnames(mMerged.down), 'down', sep='-')

mMerged.c5 = cbind(mMerged.up, mMerged.down)
# reorder the columns
colnames(mMerged.c5)
o = c(1, 2)
## sanity check
matrix(colnames(mMerged.c5)[o], ncol = 2, byrow = T)
colnames(mMerged.c5)[o]
mMerged.c5 = mMerged.c5[,o]
identical(colnames(mMerged.c2), colnames(mMerged.c5))

# remove na sections
dim(mMerged.c5)
mMerged.c5 = na.omit(mMerged.c5)
dim(mMerged.c5)
head(mMerged.c5)

### create a binary matrix based on cutoffs
getBinaryMatrix = function(mat, cutoff=0.01){
  mat2 = apply(mat, 2, function(x) round(x, 3) <= cutoff)
}

mMerged.c5.bin = getBinaryMatrix(mMerged.c5)

## group this matrix into combinations
mMerged.c5.bin.grp = mMerged.c5.bin
set.seed(123)
dm = dist(mMerged.c5.bin.grp, method='binary')
hc = hclust(dm)
plot(hc, labels=F)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mMerged.c5.bin.grp = cbind(mMerged.c5.bin.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mMerged.c5.bin.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## map these names to the cp
groups = cp
sig.pvals = rowSums(mMerged.c5.bin)
table(sig.pvals)
dfMerged.c5 = data.frame(round(mMerged.c5, 3), sig.pvals, groups, DB='C5')
str(dfMerged.c5)
head(dfMerged.c5)
tail(dfMerged.c5)

identical(colnames(dfMerged.c2), colnames(dfMerged.c5))

write.csv(dfMerged.c2, file='results/gsea_msigdb_c2_merged.xls')
write.csv(dfMerged.c5, file='results/gsea_msigdb_c5_merged.xls')

## merge together into one dataframe
# drop the group with most zeros
table(dfMerged.c2$groups)
t = rowSums(mMerged.c2.bin)
table(t, dfMerged.c2$groups)
dfMerged.c2.sub = dfMerged.c2[dfMerged.c2$groups != 2,]

table(dfMerged.c5$groups)
dfMerged.c5.sub = dfMerged.c5[dfMerged.c5$groups != 2,]

dfMerged = rbind(dfMerged.c2.sub, dfMerged.c5.sub)
dfMerged = droplevels.data.frame(dfMerged)
dim(dfMerged)

write.csv(dfMerged, file='results/gsea_msigdb_significant_c2_c5_merged.xls')

### heatmaps
### just for a quick visual check, do not use for results
df = dfMerged
head(df)
dim(df)
mMat = as.matrix(df[,c(1:2)])
head(mMat)
mMat = -1*log(mMat+1e-16)
g1 = df[,'groups']
g1 = factor(as.character(g1))
levels(g1)
g2 = df[,'DB']
g2 = factor(as.character(g2))
levels(g2)

ann = data.frame(DB=g2, Group=g1 )
#ann = data.frame(Group=g1 )
range(mMat)
quantile(as.vector(mMat), 0:20/20)
mMat[mMat > 15] = 15 
#mMat[mMat > 8] = 8

library(NMF)
library(RColorBrewer)

aheatmap(mMat, annRow = NA, scale = 'none', Rowv = order(g2:g1), Colv=NA, cexRow=5, cexCol = 0.6, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
         col=c('white', brewer.pal(9, 'YlOrRd')))

pdf('results/gsea_msigdb_significant_merged.pdf')
aheatmap(mMat, annRow = NA, scale = 'none', Rowv = order(g2:g1), Colv=NA, cexRow=5, cexCol = 0.6, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
         col=c('white', brewer.pal(9, 'YlOrRd')))
aheatmap(mMat, annRow = ann, scale = 'none', Rowv = order(g2:g1), Colv=NA, cexRow=5, cexCol = 0.6, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
         col=c('white', brewer.pal(9, 'YlOrRd')))

dev.off(dev.cur())

