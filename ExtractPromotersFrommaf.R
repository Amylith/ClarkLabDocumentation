#this script is an example of how to extract promoters from a maf file
#assuming you have human (hg19) coordinates for exons
#note that this is simplified from "ExtractGenesFrommaf.R"
#because promoters are like simpler (single exon) genes

#what you need:
#promoter coordinates
#maf files

#what you need to change in this script:
#filenames
#copy/paste code and change chrom to run other chromosomes


library(rphast)
promoterbed=read.table("/home/kowaae22/100way/PromoterCoords.txt", stringsAsFactors = F)


#function - extract alignment from maf
#####################################################################################
makealnsfrommaf=function(maf, start, end, fnames, message="done"){
  count=1
  while(count<=length(fnames)){
    newmsa=sub.msa(maf, start.col=start[count]+1, end.col=end[count]+1, refseq = "hg19")
    write.msa(newmsa, file=fnames[count], format="FASTA")
    newmsa=NULL
    gc()
    count=count+1
    # print(count)
  }
  return(message)
}
#####################################################################################


#chr22 - creates .fasta alignments for genes on chromosome 22
#####################################################################################
# sum(promoterbed$V1=="chr22") #421
maf=read.msa("/home/kowaae22/100way/maf/chr22.maf")
#get alns
start=promoterbed$V2[promoterbed$V1=="chr22"]
end=promoterbed$V3[promoterbed$V1=="chr22"]
names=promoterbed$V4[promoterbed$V1=="chr22"]
fnames=unlist(lapply(names, function(x){paste0("/home/kowaae22/100way/promoteralns/", x, ".fasta", collapse="")}))
makealnsfrommaf(maf, start, end, fnames, message="done promoterchr22")
maf=NULL
gc()
#####################################################################################




