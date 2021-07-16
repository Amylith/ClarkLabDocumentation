#clean the exon coordinates downloaded from UCSC for use in "ExtractGenesFrommaf.R"

#what you need:
#exon coordinates: UCSC Table Browser: hg38 --> Genes and Gene Predictions --> GENCODE V36 --> knownGene

#what you need to change in this script:
#filenames


coords=read.table("hg38genes.txt", header = T, stringsAsFactors = F)
#keep only canonical genes:
coords=coords[grepl("canonical", coords$tier),]
#keep only protein coding:
coords=coords[coords$transcriptType=="protein_coding",]
nrow(coords)

notunique=names(sort(table(coords$geneName),decreasing = T))[sort(table(coords$geneName),decreasing = T)>1]

#fixing not unique - genes with multiple entries manually corrected
#####################################################################################

#ACTL10 - do not use retrogene
coords=coords[!(coords$geneName=="ACTL10" & grepl("retrogene", coords$tag)),]

#AKAP17A - use X chr, discard Y
coords=coords[!(coords$geneName=="AKAP17A" & coords$chrom=="chrY"),]

#ASMT - use X chr, discard Y
coords=coords[!(coords$geneName=="ASMT" & coords$chrom=="chrY"),]

#ASMTL - use X chr, discard Y
coords=coords[!(coords$geneName=="ASMTL" & coords$chrom=="chrY"),]

#CD99 - use X chr, discard Y
coords=coords[!(coords$geneName=="CD99" & coords$chrom=="chrY"),]

#CRLF2 - use X chr, discard Y
coords=coords[!(coords$geneName=="CRLF2" & coords$chrom=="chrY"),]

#CSF2RA - use X chr, discard Y
coords=coords[!(coords$geneName=="CSF2RA" & coords$chrom=="chrY"),]

#DHRSX - use X chr, discard Y
coords=coords[!(coords$geneName=="DHRSX" & coords$chrom=="chrY"),]

#GTPBP6 - use X chr, discard Y
coords=coords[!(coords$geneName=="GTPBP6" & coords$chrom=="chrY"),]

#IL3RA - use X chr, discard Y
coords=coords[!(coords$geneName=="IL3RA" & coords$chrom=="chrY"),]

#IL9R - use X chr, discard Y
coords=coords[!(coords$geneName=="IL9R" & coords$chrom=="chrY"),]

#MATR3 - longer overall length transcript (not MANE_Select)
# starts=coords$chromStart[coords$geneName=="MATR3"]
# ends=coords$chromEnd[coords$geneName=="MATR3"]
coords=coords[!(coords$geneName=="MATR3" & grepl("MANE_Select", coords$tag)),]

#P2RY8 - use X chr, discard Y
coords=coords[!(coords$geneName=="P2RY8" & coords$chrom=="chrY"),]

#PDE11A - no readthrough transcript
coords=coords[!(coords$geneName=="PDE11A" & grepl("readthrough_transcript", coords$tag)),]

#PLCXD1 - use X chr, discard Y
coords=coords[!(coords$geneName=="PLCXD1" & coords$chrom=="chrY"),]

#POLR2J3 - no readthrough transcript
coords=coords[!(coords$geneName=="POLR2J3" & grepl("readthrough_transcript", coords$tag)),]

#PPP2R3B - use X chr, discard Y
coords=coords[!(coords$geneName=="PPP2R3B" & coords$chrom=="chrY"),]

#SHOX - use X chr, discard Y
coords=coords[!(coords$geneName=="SHOX" & coords$chrom=="chrY"),]

#SLC25A6 - use X chr, discard Y
coords=coords[!(coords$geneName=="SLC25A6" & coords$chrom=="chrY"),]

#TMSB15B - longer transcript (yes CCDS) 
# starts=coords$chromStart[coords$geneName=="TMSB15B"]
# ends=coords$chromEnd[coords$geneName=="TMSB15B"]
coords=coords[!(coords$geneName=="TMSB15B" & !grepl("CCDS", coords$tag)),]

#VAMP7 - use X chr, discard Y
coords=coords[!(coords$geneName=="VAMP7" & coords$chrom=="chrY"),]

#WASH6P - use X chr, discard Y
coords=coords[!(coords$geneName=="WASH6P" & coords$chrom=="chrY"),]

#ZBED1 - use X chr, discard Y
coords=coords[!(coords$geneName=="ZBED1" & coords$chrom=="chrY"),]


#####################################################################################

#save table of genes used
write.table(coords, file="hg38genes_clean.txt", row.names = F, col.names = F, quote = F, sep = "\t")



#reformat coordinates so each exon is one line
count=1
startcoords=c()
endcoords=c()
chrs=c()
name=c()
score=c()
strands=c()

while(count<=nrow(coords)){
  start=coords[count,]$chromStart
  end=coords[count,]$chromEnd
  chr=coords[count,]$chrom
  n=coords[count,]$geneName
  strand=coords[count,]$strand
  
  exstart=as.numeric(strsplit(coords[count,]$chromStarts, split=",")[[1]])
  exsize=as.numeric(strsplit(coords[count,]$blockSizes, split=",")[[1]])
  
  excount=1
  while(excount<=length(exstart)){
    startcoords=c(startcoords, start+exstart[excount])
    endcoords=c(endcoords, start+exstart[excount]+exsize[excount])
    excount=excount+1
  }
  
  
  chrs=c(chrs, rep(chr, excount-1))
  name=c(name, rep(n, excount-1))
  score=c(score, rep(0, excount-1))
  strands=c(strands, rep(strand, excount-1))
  
  count=count+1
  print(count)
}

exondf=data.frame(chrs, startcoords, endcoords, name, score, strands)
colnames(exondf)=c("chr", "start", "end", "name", "value", "strand")
write.table(exondf, file="hg38exoncoords.txt", row.names = F, col.names = F, quote = F, sep = "\t")



