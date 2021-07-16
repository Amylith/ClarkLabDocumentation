#this script is an example of how to extract genes from a maf file
#assuming you have human (hg38) coordinates for exons
#coords file should have separate entries for each exon
#where all exons for one gene have the same gene name
#exons are extracted from the maf file and concatenated into one .fasta file
#for genes with many exons, a temporary file is created to conserve memory

#what you need:
#exon coordinates
#maf files

#what you need to change in this script:
#filenames
#copy/paste code and change chrom to run other chromosomes

library(rphast)
coords=read.table("/home/kowaae22/120way/hg38exoncoords.txt", stringsAsFactors = F)


#function - split coords per gene to create batches of exons
#each batch runs together and is saved in a temporary file before running the next batch
#this uses less memory than storing all exons in working memory
#####################################################################################

splitcoords=function(x, n=2){
  indices=seq(from=1, to=length(x), by=n)
  indices[length(indices)]=length(x)
  vallist=vector(mode = "list", length=length(indices)-1)
  count=1
  while(count<=length(vallist)){
    if(count==1){
      start=indices[count]
    }else{
      start=indices[count]+1
    }
    end=indices[count+1]
    vallist[[count]]=x[start:end]
    count=count+1
  }
  return(vallist)
}

#####################################################################################


#chrY - creates .fasta alignments for genes on chromosome Y
#####################################################################################
chrom="chrY"
maf=read.msa("/home/kowaae22/120way/maffiles/chrY.maf")
#get alns
names=coords$V4[coords$V1==chrom]
uniquenames=unique(names)
count=1
while(count<=length(uniquenames)){
  fn=paste0("/home/kowaae22/120way/NewCodingAlns/", uniquenames[count], ".fasta", collapse="")
  start=coords$V2[coords$V1==chrom&coords$V4==uniquenames[count]]
  end=coords$V3[coords$V1==chrom&coords$V4==uniquenames[count]]
  if(sum(uniquenames[count]==names)==1){
    newmsa=sub.msa(maf, start.col=start+1, end.col=end+1, refseq = "hg38")
    write.msa(newmsa, file=fn, format="FASTA")
    newmsa=NULL
    gc()
  }else{
    if(length(start)<=5){
      excount=1
      exlist=vector(mode="list", length=length(start))
      while(excount<=length(start)){
        exlist[[excount]]=sub.msa(maf, start.col=start[excount]+1, end.col=end[excount]+1, refseq = "hg38")
        excount=excount+1
        # print(excount)
      }
      newmsa=concat.msa(exlist)
      write.msa(newmsa, file=fn, format="FASTA")
    }else{
      startsplits=splitcoords(start)
      endsplits=splitcoords(end)
      
      splitcount=1
      while(splitcount<=length(startsplits)){
        substarts=startsplits[[splitcount]]
        subends=endsplits[[splitcount]]
        
        excount=1
        exlist=vector(mode="list", length=length(substarts))
        while(excount<=length(substarts)){
          exlist[[excount]]=sub.msa(maf, start.col=substarts[excount]+1, end.col=subends[excount]+1, refseq = "hg38")
          excount=excount+1
          # print(excount)
        }
        # newmsa=concat.msa(exlist)
        # exlist=NULL
        # gc()
        if(splitcount==1){
          write.msa(concat.msa(exlist), file=paste0("/home/kowaae22/temp",chrom,".fasta"), format="FASTA")
        }else if(splitcount<length(startsplits)){
          # newmsa=concat.msa(list(read.msa(file=paste0("/home/kowaae22/temp",chrom,".fasta"), format="FASTA"),newmsa))
          write.msa(concat.msa(list(read.msa(file=paste0("/home/kowaae22/temp",chrom,".fasta"), format="FASTA"),concat.msa(exlist))), file=paste0("/home/kowaae22/temp",chrom,".fasta"), format="FASTA")
        }else{
          # newmsa=concat.msa(list(read.msa(file=paste0("/home/kowaae22/temp",chrom,".fasta"), format="FASTA"),newmsa))
          write.msa(concat.msa(list(read.msa(file=paste0("/home/kowaae22/temp",chrom,".fasta"), format="FASTA"),concat.msa(exlist))), file=fn, format="FASTA")
        }
        newmsa=NULL
        gc()
        splitcount=splitcount+1
        print(splitcount)
      }
      
    }
    newmsa=NULL
    exlist=NULL
    gc()
  }
  
  count=count+1
  print(paste0(chrom, ": ", count, " out of ", length(uniquenames)))
}
#####################################################################################


#helper code - run this to see if any genes didn't write to files and gives their names
#####################################################################################
#ones that got skipped
chrom="chrY"
names=coords$V4[coords$V1==chrom]
uniquenames=unique(names)
fns=unlist(strsplit(list.files("/home/kowaae22/120way/NewCodingAlns/"), ".fasta"))

inds=which(!(uniquenames %in% fns))
uniquenames[which(!(uniquenames %in% fns))]
#####################################################################################


