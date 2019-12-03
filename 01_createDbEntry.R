# File: 01_createDbEntry.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: list the samples and create appropriate db entries
# Date: 3/12/2019


## set variables and source libraries
source('header.R')

## connect to mysql database 
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# sample and file table
dbGetQuery(db, paste('describe Sample;'))
cSampleCol = dbGetQuery(db, paste('describe Sample;'))$Field[-1]

dbGetQuery(db, paste('describe File;'))
cFileCol = dbGetQuery(db, paste('describe File;'))$Field[-1]

# setwd(gcRemoteDir)
setwd('dataExternal/')
setwd('fastq/')

# list the files
cvFiles = list.files(pattern = 'fastq.gz', recursive = T)

# each sample has 2 files 
fSplit = gsub('_R[1|2]', '', cvFiles)

lFiles = split(cvFiles, fSplit)

setwd(gcswd)
## load the metadata file
dfMeta = read.csv('dataExternal/metadata.csv', header=T, stringsAsFactors = F)
str(dfMeta)

cn = colnames(dfMeta)
for (i in seq_along(1:ncol(dfMeta))) dfMeta[,cn[i]] = gsub(' ', '', dfMeta[,cn[i]])

# sanity check
table(as.character(dfMeta$Name) %in% cvFiles)

## order the table in the same sequence as file names
i = match(cvFiles, dfMeta$Name)
dfMeta = dfMeta[i,]
identical(as.character(dfMeta$Name), cvFiles)
dfMeta$fSplit = fSplit

## extract reduced sample table by removing one pair of the fastq file
cvSample = unique(dfMeta$fSplit)
i = match(cvSample, dfMeta$fSplit)
dfMeta.sam = dfMeta[i,]
str(dfMeta.sam)
xtabs( ~ Sample_ID + Lesional + Patient_ID, data=dfMeta.sam)
xtabs( ~ Sequencing_Run + Lesional + Patient_ID, data=dfMeta.sam)
xtabs( ~ Sequencing_Run + Sequencing_Lane + Patient_ID, data=dfMeta.sam)

## create the entry for samples
cSampleCol

dfSamples = data.frame(idProject=g_pid, idData=g_did, title=dfMeta.sam$Sample_ID, 
                       description= paste('sequencing lane', as.character(dfMeta.sam$Sequencing_Lane),
                                          'group1 is Treatment',
                                          'group2 is patient batch',
                                          'group3 is sequencing run',
                                          'use replicate pairs as technical replicates identifier', sep=';'),
                       group1 = dfMeta.sam$Lesional, group2= dfMeta.sam$Patient_ID, group3=dfMeta.sam$Sequencing_Run)
# write this data to the database
rownames(dfSamples) = NULL

### NOTE: Do not execute this anymore as entry created
# write this table to database
#dbWriteTable(db, name='Sample', value=dfSamples, append=T, row.names=F)
# get this table again from database with ids added
g_did
dfSamples = dbGetQuery(db, paste('select * from Sample where Sample.idData = 43;'))

# create entries for these files in the database
dbListTables(db)
cn = dbListFields(db, 'File')[-1]
cn
identical(names(lFiles), dfMeta.sam$fSplit)
names(lFiles) = dfSamples$id

# get the names of the samples
temp = lapply(as.character(dfSamples$id), function(x){
  # get the file names
  df = data.frame(name=lFiles[[x]], type='fastq', idSample=dfSamples[as.character(dfSamples$id) == x, 'id'])
  return(df)
})

dfFiles = do.call(rbind, temp)
rownames(dfFiles) = NULL

# write this table to database
## note: do not execute as it is already done
#dbWriteTable(db, name='File', value=dfFiles, append=T, row.names=F)

dbDisconnect(db)
