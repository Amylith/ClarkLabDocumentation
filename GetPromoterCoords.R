#this script shows how to get promoter coordinates from gene coordinates
#using gene coordinates from hg19

#what you need:
#gene coordinates
#(if necessary) file mapping gene codes to gene names

#what you need to change in this script:
#filenames
#(optional) distance-based definition of promoter



#gene coordinates and map of gene names to gene codes
genecoords=read.table("C:/Users/apacz/Desktop/PromoterAnalysis/GeneCoords.txt", header = T, sep = "\t", stringsAsFactors = F)
namemap=read.table("C:/Users/apacz/Desktop/PromoterAnalysis/ucidCorrect_commonName.map", stringsAsFactors = F, header=T)

sum(namemap$ucid %in% genecoords$name)
nrow(namemap)

#get gene names
subcoords=genecoords[genecoords$name %in% namemap$ucid,]
subcoords$genename=namemap$gene[match(subcoords$name,namemap$ucid)]

#functions - get start and end coordinate for promoter
#####################################################################
#plus strand - txstart-2000 (or 0) and txstart
#minus strand - txend (or 0) and txend+2000
getpromoterstart=function(start, end, strand, distance=2000){
  startprom=0
  endprom=0
  if(strand=="+"){
    startprom=max(c(start-distance, 0))
    endprom=start
  }
  if(strand=="-"){
    startprom=max(c(end,0))
    endprom=start+distance
  }
  return(startprom)
}

getpromoterend=function(start, end, strand, distance=2000){
  startprom=0
  endprom=0
  if(strand=="+"){
    startprom=max(c(start-distance, 0))
    endprom=start
  }
  if(strand=="-"){
    startprom=max(c(end,0))
    endprom=end+distance
  }
  return(endprom)
}
#####################################################################


#get promoter using distance=100bases upstream of transcription start 
subcoords$promoterstart=mapply(getpromoterstart, subcoords$txStart, subcoords$txEnd, subcoords$strand, 100)
subcoords$promoterend=mapply(getpromoterend, subcoords$txStart, subcoords$txEnd, subcoords$strand, 100)
#make bed
#chr start end name value
promoterbed=data.frame(subcoords$chrom, subcoords$promoterstart, subcoords$promoterend, subcoords$genename, rep(0, nrow(subcoords)))
colnames(promoterbed)=c("chr", "start", "end", "name", "value")
write.table(promoterbed, file="C:/Users/apacz/Desktop/PromoterAnalysis/PromoterCoordsMultipleDistances/PromoterCoords100.txt", row.names = F, col.names = F, quote = F, sep = "\t")


