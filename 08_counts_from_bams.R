# File: 08_counts_from_bams.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: generate count tables for transcripts from bam files
# Date: 11/12/2019


## set variables and source libraries
source('header.R')

## load the transcript db objects
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicAlignments)

# get the exons into GRangesList object
oGRLgenes = exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = 'gene')

## create the bamfile list from database
## connect to mysql database to get sample information
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'File')
# get the query
g_did
q = paste0('select Sample.id as sid, Sample.group1, Sample.group2, Sample.title, File.* from Sample, File
           where (Sample.idData = 43) AND (File.idSample = Sample.id AND File.type like "%duplicates removed%")')
dfSample = dbGetQuery(db, q)
nrow(dfSample)
dfSample
# close connection after getting data
dbDisconnect(db)
# remove any whitespace from the names
dfSample$name = gsub(" ", "", dfSample$name, fixed = T)
dfSample$title = gsub(" ", "", dfSample$title, fixed = T)
dfSample$group1 = gsub(" ", "", dfSample$group1, fixed = T)
# create a new file path variable 
dfSample$fp = paste0(dfSample$name)

#### set working directory to appropriate location with bam files
setwd('dataExternal/bams/')
csFiles = list.files('.', pattern = '*.bam$', recursive = F)
# check if these files match the file names in database
table(dfSample$fp %in% csFiles)
dim(dfSample)
csFiles = dfSample$fp

## order the data frame by sample
dfSample = dfSample[order(dfSample$id),]
# split the files by ids to reduce memory usage
lFiles = split(dfSample$fp, dfSample$sid)

# for each of these bam file lists do the counting
lCounts = lapply(lFiles, function(bfl){
  ## create a bamfiles list object
  oBamFiles = BamFileList(bfl, index=paste0(bfl, '.bai'))
  return(assays(summarizeOverlaps(oGRLgenes, oBamFiles, ignore.strand = F, singleEnd=F, preprocess.reads=invertStrand))$counts)
})

# sanity checks
identical(as.character(sapply(lCounts, colnames)), dfSample$name)
identical(names(lCounts), as.character(dfSample$sid))
names(lCounts) = dfSample$sid

## save the summarized experiment object
setwd(gcswd)
n = make.names(paste('lCounts rd did 43 alex reverse rds'))
n2 = paste0('~/Data/MetaData/', n)
save(lCounts, file=n2)

## comment out after first time execution
library('RMySQL')
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/',
                comment='list of Count matrix from alex human blaschkoid psoriasis project with quality 10 and duplicates removed and strands reversed')
dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
dbDisconnect(db)
