#this script contains code to calculate TFBS statistics across 350k noncoding regions and the UCSC mammals
#coordinates for 350k regions in hg19 must first be lifted over to other speices
#those coordinates must then be combined and removed if their lengths are very incorrect
#note that this code has already run (to summarize TFs over CNEs), so you can use those and start at the PGLS step
#all of this is stored on gobie/discus

#what you need:
#noncoding region trees
#noncoding region hg19 coordinates
#TFBS calls across species for all TFBS motifs
#liftover files
#phenotype

#what you need to change in this script:
#filenames





#creating TF coords
###########################################################################
#make directories for summaries by TF
specs=list.dirs("/home/Genomes/HOCOMOCO/Organism/", full.names = F)
specs=specs[-1]
dir.create("/home/kowaae22/TFcalls/TFsummaries")

for(s in specs){
  dname=paste0("/home/kowaae22/TFcalls/TFsummaries/", s)
  dir.create(dname)
}

#make correct coordinates per species
dir.create("/home/kowaae22/TFcalls/coords")

#check for missing chain files/incorrect names
# allchain=list.files("/home/kowaae22/TFcalls/chainfiles/", full.names = F)
# for(s in specs){
#   temp=strsplit(s, split="")[[1]]
#   temp[1]=toupper(temp[1])
#   temp=paste0(temp, collapse = "")
#   chain=paste0("hg19To", temp, ".over.chain.gz")
#   if(!(chain %in% allchain)){
#     print(s)
#   }
# }

#liftover coordinates to all other species
for(s in specs){
  if(s!="hg19"){
    newfile=paste0("/home/kowaae22/TFcalls/coords/", s, "coords")
    temp=strsplit(s, split="")[[1]]
    temp[1]=toupper(temp[1])
    temp=paste0(temp, collapse = "")
    chain=paste0("/home/kowaae22/TFcalls/chainfiles/hg19To", temp, ".over.chain.gz")
    # system(paste0("CrossMap.py bed ",chain," /home/kowaae22/TFcalls/phastconsbedsimple.bed ",newfile))
    print(paste0("CrossMap.py bed ",chain," /home/kowaae22/TFcalls/phastconsbedsimple.bed ",newfile))
  }
}
#can't allocate memory - put in TFcallsshell and ran directly in terminal


#reorder coordinates 
for(s in specs){
  fn=paste0("/home/kowaae22/TFcalls/coords/", s, "coords")
  newfn=paste0("/home/kowaae22/TFcalls/orderedcoords/", s, "coords")
  # system(paste0("sort -k 1,1 -k2,2n ",fn," > ",newfn))
  print(paste0("sort -k 1,1 -k2,2n ",fn," > ",newfn))
}
#error -put in sortshell and ran in terminal

#skip this - it combines too many things
# #merge bookended reads
# 
# for(s in specs){
#   fn=paste0("/home/kowaae22/TFcalls/orderedcoords/", s, "coords")
#   newfn=paste0("/home/kowaae22/TFcalls/finalcoords/", s, "coords")
#   # system(paste0("sort -k 1,1 -k2,2n ",fn," > ",newfn))
#   print(paste0("bedtools merge -i ",fn," > ",newfn))
# }
# #in mergeshell


#doesn't work - need to do with TF coords
# for(s in specs){
#   coords=read.table(paste0("/home/kowaae22/TFcalls/orderedcoords/",s,"coords"), stringsAsFactors=F)
#   chr=unique(coords[,1])
#   len=rep(1000000, length(chr))
#   df=data.frame(chr, len)
#   write.table(df, file=paste0("/home/kowaae22/TFcalls/genomefiles/",s), quote=F, col.names=F, row.names=F, sep="\t")
# }


#get genome order based on TF calls AND genome coordinates
#this defines the order of the files to use for mapping
for(s in specs){
  coords=read.table(paste0("/home/kowaae22/TFcalls/orderedcoords/",s,"coords"), stringsAsFactors=F)
  chr=unique(coords[,1])
  tffiles=list.files(paste0("/home/Genomes/HOCOMOCO/Organism/", s), full.names=T)
  count=1
  for(t in tffiles){
    calls=read.table(t, stringsAsFactors=F)
    cvals=unique(calls[,1])
    chr=c(chr, cvals)
    chr=unique(chr)
    print(count)
    count=count+1
  }
  chr=sort(chr)
  df=data.frame(chr, rep(1e6, length(chr)))
  write.table(df, file=paste0("/home/kowaae22/TFcalls/genomefiles/",s), quote=F, col.names=F, row.names=F, sep="\t")
  print(s)
}


#sort TF calls
for(s in specs){
  dname=paste0("/home/kowaae22/TFcalls/rawTFcalls/", s)
  dir.create(dname)
}

sortlines=c()
for(s in specs){
  tffiles=list.files(paste0("/home/Genomes/HOCOMOCO/Organism/", s))
  for(f in tffiles){
    test=paste0("/home/Genomes/HOCOMOCO/Organism/", s, "/", f)
    outfiles=paste0("/home/kowaae22/TFcalls/rawTFcalls/", s, "/", f)
    l=paste0("sort -k 1,1 -k2,2n ",test," > ",outfiles)
    sortlines=c(sortlines, l)
  }
  print(s)
}
conn=file("/home/kowaae22/TFcalls/sortrawTF")
writeLines(sortlines, conn)
close(conn)



#STOP HERE and go to "custom merge" - this code (using unfiltered/unmerged coords) is only included for completeness

#summarize over each species coordinates over each TF
# bedtools map -a test.merged.out -b /home/Genomes/HOCOMOCO/Organism/mm10/TYY1_HUMAN.H11MO.0.A.bed -c 5 -o mean,median,sum,max,min,count > testsummary.out

#output: TFsummaries, input: orderedcoords

alllines=c()
for(s in specs){
  tffiles=list.files(paste0("/home/kowaae22/TFcalls/rawTFcalls/", s))
  coords=paste0("/home/kowaae22/TFcalls/orderedcoords/", s, "coords")
  gfile=paste0("/home/kowaae22/TFcalls/genomefiles/",s)
  for(f in tffiles){
    test=paste0("/home/kowaae22/TFcalls/rawTFcalls/", s, "/", f)
    outfiles=paste0("/home/kowaae22/TFcalls/TFsummaries/", s, "/", f)
    l=paste0("bedtools map -a ", coords, " -b ", test, " -c 5 -o mean,median,sum,max,min,count -g ", gfile," > ", outfiles)
    alllines=c(alllines, l)
  }
}

conn=file("/home/kowaae22/TFcalls/runmap")
writeLines(alllines, conn)
close(conn)

###########################################################################



#custom merge

#custom merging - done!
###########################
specs=list.dirs("/home/Genomes/HOCOMOCO/Organism/", full.names = F)
specs=specs[-1]


for(s in specs){
  f=read.table(paste0("/home/kowaae22/TFcalls/orderedcoords/",s, "coords"), stringsAsFactors=F)
  c=table(f$V4)
  c=data.frame(c)
  c=setNames(c$Freq, c$Var1)
  c=sort(c, decreasing=T)
  dups=names(c)[c>1]
  count=1
  while(count<=length(dups)){
    curdups=f[f$V4==dups[count],]
    locs=which(f$V4==dups[count])
    if(sum((locs-locs[1])==0:(length(locs)-1))==length(locs)){
      f$V3[locs[1]]=curdups$V3[nrow(curdups)]
      f=f[-locs[-1],]
    }else{
      f=f[-locs,]
    }
    count=count+1
    if(count %% 1000 ==0){
      print(count/length(dups))
    }
  }
  write.table(f, file=paste0("/home/kowaae22/TFcalls/custommergedcoords/",s, "coords"), quote=F, col.names=F, row.names=F, sep="\t")
  print(s) 
}
###########################

#check custom merge length - done!
#removed if more than 100 bases different in length than hg19
#removed if less than 30 bases total length
###########################
specs=list.dirs("/home/Genomes/HOCOMOCO/Organism/", full.names = F)
specs=specs[-1]
specs=specs[-23]
f=read.table("/home/kowaae22/TFcalls/custommergedcoords/hg19coords", stringsAsFactors=F)
lengths=data.frame(f$V4, f$V3-f$V2)
n=c("name", "hg19", specs)

for(s in specs){
  f=read.table(paste0("/home/kowaae22/TFcalls/custommergedcoords/",s,"coords"), stringsAsFactors=F)
  f=data.frame(f$V4, f$V3-f$V2)
  colnames(f)=c("name", "length")
  lengths=merge(lengths, f, by.x="f.V4", by.y="name", all=T)
  print(s)
}
colnames(lengths)=n

# count=3
# while(count<=ncol(lengths)){
#   print(sum(abs(lengths$hg19-lengths[,count])>100, na.rm=T))
#   count=count+1
# }

# plot(lengths$hg19, lengths$ailMel1)
# hist(lengths$hg19-lengths$ailMel1)
# sum(abs(lengths$hg19-lengths$ailMel1)>100, na.rm=T)


# filteredcustommergedcoords

# count=46
# while(count<=63){
#   print(sum(abs(lengths$hg19-lengths[,count])>100, na.rm=T))
#   count=count+1
# }


f=read.table("/home/kowaae22/TFcalls/custommergedcoords/hg19coords", stringsAsFactors=F)
write.table(f, file=paste0("/home/kowaae22/TFcalls/filteredcustommergedcoords/","hg19", "coords"), quote=F, col.names=F, row.names=F, sep="\t")

count=3
while(count<=ncol(lengths)){
  s=colnames(lengths)[count]
  f=read.table(paste0("/home/kowaae22/TFcalls/custommergedcoords/",s,"coords"), stringsAsFactors=F)
  
  bad=as.character(na.omit(lengths$name[abs(lengths$hg19-lengths[,count])>100]))
  short=as.character(na.omit(lengths$name[lengths[,count]<30]))
  bad=unique(c(bad, short))
  
  f=f[!(f$V4 %in% bad),]
  
  write.table(f, file=paste0("/home/kowaae22/TFcalls/filteredcustommergedcoords/",s, "coords"), quote=F, col.names=F, row.names=F, sep="\t")
  print(s)
  count=count+1
}

###########################

#get custom merge values (makes a bash file you need to run on command line)
###########################
#summarize over each species coordinates over each TF
#FOR CUSTOM MERGED PHASTCONS COORDS
specs=list.dirs("/home/Genomes/HOCOMOCO/Organism/", full.names = F)
specs=specs[-1]

for(s in specs){
  dname=paste0("/home/kowaae22/TFcalls/TFsummariescustommerge/", s)
  dir.create(dname)
}

alllines=c()
for(s in specs){
  tffiles=list.files(paste0("/home/kowaae22/TFcalls/rawTFcalls/", s))
  coords=paste0("/home/kowaae22/TFcalls/filteredcustommergedcoords/", s, "coords")
  gfile=paste0("/home/kowaae22/TFcalls/genomefiles/",s)
  for(f in tffiles){
    test=paste0("/home/kowaae22/TFcalls/rawTFcalls/", s, "/", f)
    outfiles=paste0("/home/kowaae22/TFcalls/TFsummariescustommerge/", s, "/", f)
    l=paste0("bedtools map -a ", coords, " -b ", test, " -c 5 -o mean,median,sum,max,min,count -g ", gfile," > ", outfiles)
    alllines=c(alllines, l)
  }
}

conn=file("/home/kowaae22/TFcalls/runmapcustommerge")
writeLines(alllines, conn)
close(conn)
###########################


#summarize with custom merge
###########################

#one df for each: mean,median,sum,max,min,count
#rows as CNEs and cols as species, one file per TF

specs=list.dirs("/home/Genomes/HOCOMOCO/Organism/", full.names = F)
specs=specs[-1]

allhoco=unique(list.files("/home/kowaae22/TFcalls/TFsummariescustommerge/hg19"))
coords=read.table("/home/kowaae22/TFcalls/filteredcustommergedcoords/hg19coords", stringsAsFactors=F)

n=c("name", "chr", "start", "end", specs)

count=1
for(h in allhoco){
  newfn=paste0("/home/kowaae22/TFcalls/groupspecscustommergedups/", h)
  
  dfmean=coords[,1:4]
  dfmedian=coords[,1:4]
  dfsum=coords[,1:4]
  dfmax=coords[,1:4]
  dfmin=coords[,1:4]
  dfcount=coords[,1:4]
  for(s in specs){
    f=read.table(paste0("/home/kowaae22/TFcalls/TFsummariescustommerge/", s, "/", h), stringsAsFactors=F)
    f[f=="."]=0
    
    #remove duplicates - this shouldn't do anything
    c=table(f$V4)
    c=data.frame(c)
    c=setNames(c$Freq, c$Var1)
    c=sort(c, decreasing=T)
    dups=names(c)[c>1]
    f=f[!(f$V4 %in% dups),]
    
    #summarize
    f1=data.frame(f$V4, f$V7)
    dfmean=merge(dfmean, f1, by.x="V4", by.y="f.V4", all=T)
    
    f2=data.frame(f$V4, f$V8)
    dfmedian=merge(dfmedian, f2, by.x="V4", by.y="f.V4", all=T)
    
    f3=data.frame(f$V4, f$V9)
    dfsum=merge(dfsum, f3, by.x="V4", by.y="f.V4", all=T)
    
    f4=data.frame(f$V4, f$V10)
    dfmax=merge(dfmax, f4, by.x="V4", by.y="f.V4", all=T)
    
    f5=data.frame(f$V4, f$V11)
    dfmin=merge(dfmin, f5, by.x="V4", by.y="f.V4", all=T)
    
    f6=data.frame(f$V4, f$V12)
    dfcount=merge(dfcount, f6, by.x="V4", by.y="f.V4", all=T)
    
    print(s)
  }
  colnames(dfmean)=n
  colnames(dfmedian)=n
  colnames(dfsum)=n
  colnames(dfmax)=n
  colnames(dfmin)=n
  colnames(dfcount)=n
  
  write.table(dfmean, file=paste0(newfn, "mean"), quote=F, row.names=F, sep="\t")
  write.table(dfmedian, file=paste0(newfn, "median"), quote=F, row.names=F, sep="\t")
  write.table(dfsum, file=paste0(newfn, "sum"), quote=F, row.names=F, sep="\t")
  write.table(dfmax, file=paste0(newfn, "max"), quote=F, row.names=F, sep="\t")
  write.table(dfmin, file=paste0(newfn, "min"), quote=F, row.names=F, sep="\t")
  write.table(dfcount, file=paste0(newfn, "count"), quote=F, row.names=F, sep="\t")
  
  count=count+1
  print(paste0("allhoco: ", count))
}

###########################

#Analysis with custom merge and region-specific trees (binary phenotype)
###########################

specs=list.dirs("/home/Genomes/HOCOMOCO/Organism/", full.names = F)
specs=specs[-1]

foreground=readRDS("/home/kowaae22/AnalysisWithThreeTrees/hairlessSpecs.rds")

# alltrees=readRDS("/home/kowaae22/AnalysisWithThreeTrees/finaltreesv4.rds")
# class(alltrees)[2]="treesObj"
library(nlme)
library(ape)
library(phytools)

allhoco=unique(list.files("/home/kowaae22/TFcalls/TFsummariescustommerge/hg19"))
# folder="groupspecsremovedups" #duplicates removed
folder="groupspecscustommergedups" #custom merge
stat="median"
# stat="mean"
# stat="min"
# stat="max"
# stat="count"
# stat="sum"

trees1=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part1.trees.rds")
trees2=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part2.trees.rds")
trees3=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part3.trees.rds")
trees4=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part4.trees.rds")
trees5=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part5.trees.rds")
trees6=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part6.trees.rds")
trees7=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part7.trees.rds")

getTree=function(regionname){
  if(regionname %in% names(trees1$trees)){
    ind=which(names(trees1$trees)==regionname)
    return(trees1$trees[[ind]])
  }
  if(regionname %in% names(trees2$trees)){
    ind=which(names(trees2$trees)==regionname)
    return(trees2$trees[[ind]])
  }
  if(regionname %in% names(trees3$trees)){
    ind=which(names(trees3$trees)==regionname)
    return(trees3$trees[[ind]])
  }
  if(regionname %in% names(trees4$trees)){
    ind=which(names(trees4$trees)==regionname)
    return(trees4$trees[[ind]])
  }
  if(regionname %in% names(trees5$trees)){
    ind=which(names(trees5$trees)==regionname)
    return(trees5$trees[[ind]])
  }
  if(regionname %in% names(trees6$trees)){
    ind=which(names(trees6$trees)==regionname)
    return(trees6$trees[[ind]])
  }
  if(regionname %in% names(trees7$trees)){
    ind=which(names(trees7$trees)==regionname)
    return(trees7$trees[[ind]])
  }
}
# mt=trees$masterTree

pheno=readRDS("/home/kowaae22/AnalysisWithThreeTrees/hairlessSpecs.rds")
phenvec=rep(0, length(mt$tip.label))
names(phenvec)=mt$tip.label
phenvec[names(phenvec) %in% pheno]=1

full=read.table("/home/kowaae22/TFcalls/filteredcustommergedcoords/hg19coords", stringsAsFactors=F)

resultsdf=data.frame(matrix(nrow=nrow(full), ncol=length(allhoco)))
# snames=unlist(lapply(strsplit(allhoco, split="_"), function(x){return(x[1])}))
colnames(resultsdf)=allhoco
rownames(resultsdf)=full$V4

allresults=list(resultsdf,
                resultsdf,
                resultsdf,
                resultsdf,
                resultsdf,
                resultsdf)
names(allresults)=c("PGLSp", "PGLSstat", "PGLSbinp", "PGLSbinstat", "PGLSPAp", "PGLSPAstat")

#PGLS
# allresults=readRDS("allresultsfirst5.rds")
# allhoco=allhoco[6:26]
hcount=1
for(h in allhoco){
  fn=paste0("/home/kowaae22/TFcalls/", folder, "/", h, stat)
  data=read.table(fn, stringsAsFactors = F)
  colnames(data)=data[1,]
  data=data[-1,]
  colnames(data)[colnames(data)=="odoRosDiv1"]="odoRosDi"
  
  rownames(allresults$PGLSp)=data$name
  rownames(allresults$PGLSstat)=data$name
  rownames(allresults$PGLSbinp)=data$name
  rownames(allresults$PGLSbinstat)=data$name
  rownames(allresults$PGLSPAp)=data$name
  rownames(allresults$PGLSPAstat)=data$name
  
  count=1
  while(count<=nrow(data)){
    curcne=data$name[count]
    TF=setNames(data[count,], colnames(data))
    TF=TF[-c(1:4)]
    TF=TF[match(names(phenvec), names(TF))]
    TF=as.numeric(TF)
    df=data.frame(TF, phenvec)
    df2=na.omit(df)
    mt2=getTree(curcne)
    
    keep=intersect(rownames(df2), mt2$tip.label)
    
    if(!is.null(mt2)){
      df2=df2[rownames(df2) %in% keep,]
      mt2=keep.tip(mt2, keep)
    }
    
    df$PA=ifelse(is.na(df$TF), 0, 1)
    df2bin=df2
    df2bin$TF=ifelse(df2$TF==0,0,1)
    df=df[,2:3]
    #continuous PGLS
    if(length(unique(df2$TF))!=1 & length(unique(df2$phenvec))!=1 & nrow(unique(df2))>2 & !is.null(mt2)){
      
      pgls=gls(TF~phenvec, correlation = corBrownian(phy=mt2), data=df2)
      pvalPGLS=summary(pgls)$tTable[2,4]
      statPGLS=summary(pgls)$tTable[2,3]
      allresults$PGLSp[count, hcount]=pvalPGLS
      allresults$PGLSstat[count, hcount]=statPGLS 
    }else{
      allresults$PGLSp[count, hcount]=NA
      allresults$PGLSstat[count, hcount]=NA
    }
    
    #binary PGLS
    if(length(unique(df2bin$TF))!=1 & length(unique(df2bin$phenvec))!=1 & nrow(unique(df2bin))>2 & !is.null(mt2)){
      pglsbin=gls(TF~phenvec, correlation = corBrownian(phy=mt2), data=df2bin)
      pvalPGLSbin=summary(pgls)$tTable[2,4]
      statPGLSbin=summary(pgls)$tTable[2,3]
      allresults$PGLSbinp[count, hcount]=pvalPGLSbin
      allresults$PGLSbinstat[count, hcount]=statPGLSbin
    }else{
      allresults$PGLSbinp[count, hcount]=NA
      allresults$PGLSbinstat[count, hcount]=NA
    }
    
    #presence/absence of CNE
    if(length(unique(df$PA))!=1 & length(unique(df$phenvec))!=1 & nrow(unique(df))>2){
      pglsPA=gls(PA~phenvec, correlation = corBrownian(phy=mt), data=df)
      pvalPGLSPA=summary(pglsPA)$tTable[2,4]
      statPGLSPA=summary(pglsPA)$tTable[2,3]
      allresults$PGLSPAp[count, hcount]=pvalPGLSPA
      allresults$PGLSPAstat[count, hcount]=statPGLSPA
    }else{
      allresults$PGLSPAp[count, hcount]=NA
      allresults$PGLSPAstat[count, hcount]=NA
    }
    
    if(count %% 10000==0){
      print(paste0("CNE count: ", count)) #345786
    }
    count=count+1
  }
  print(paste0("TF count: ", hcount)) #771
  saveRDS(allresults, paste0("/home/kowaae22/TFcalls/allresultsmergedcustomtrees/allresultsmergedCTfirst",hcount,".rds"))
  hcount=hcount+1
}

###########################








#BONUS - another PGLS example

#PC1 longevity analysis with custom merge and region-specific trees
###########################

specs=list.dirs("/home/Genomes/HOCOMOCO/Organism/", full.names = F)
specs=specs[-1]

# foreground=readRDS("/home/kowaae22/AnalysisWithThreeTrees/hairlessSpecs.rds")
PC1=readRDS("/home/kowaae22/AnalysisWithThreeTrees/PC1.rds")

# alltrees=readRDS("/home/kowaae22/AnalysisWithThreeTrees/finaltreesv4.rds")
# class(alltrees)[2]="treesObj"
library(nlme)
library(ape)
library(phytools)

allhoco=unique(list.files("/home/kowaae22/TFcalls/TFsummariescustommerge/hg19"))
# folder="groupspecsremovedups" #duplicates removed
folder="groupspecscustommergedups" #custom merge
stat="median"
# stat="mean"
# stat="min"
# stat="max"
# stat="count"
# stat="sum"

trees1=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part1.trees.rds")
trees2=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part2.trees.rds")
trees3=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part3.trees.rds")
trees4=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part4.trees.rds")
trees5=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part5.trees.rds")
trees6=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part6.trees.rds")
trees7=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part7.trees.rds")

getTree=function(regionname){
  if(regionname %in% names(trees1$trees)){
    ind=which(names(trees1$trees)==regionname)
    return(trees1$trees[[ind]])
  }
  if(regionname %in% names(trees2$trees)){
    ind=which(names(trees2$trees)==regionname)
    return(trees2$trees[[ind]])
  }
  if(regionname %in% names(trees3$trees)){
    ind=which(names(trees3$trees)==regionname)
    return(trees3$trees[[ind]])
  }
  if(regionname %in% names(trees4$trees)){
    ind=which(names(trees4$trees)==regionname)
    return(trees4$trees[[ind]])
  }
  if(regionname %in% names(trees5$trees)){
    ind=which(names(trees5$trees)==regionname)
    return(trees5$trees[[ind]])
  }
  if(regionname %in% names(trees6$trees)){
    ind=which(names(trees6$trees)==regionname)
    return(trees6$trees[[ind]])
  }
  if(regionname %in% names(trees7$trees)){
    ind=which(names(trees7$trees)==regionname)
    return(trees7$trees[[ind]])
  }
}
# mt=trees$masterTree

# pheno=readRDS("/home/kowaae22/AnalysisWithThreeTrees/hairlessSpecs.rds")
PC1=readRDS("/home/kowaae22/AnalysisWithThreeTrees/PC1.rds")
phenvec=PC1
# phenvec=rep(0, length(mt$tip.label))
# names(phenvec)=mt$tip.label
# phenvec[names(phenvec) %in% pheno]=1

full=read.table("/home/kowaae22/TFcalls/filteredcustommergedcoords/hg19coords", stringsAsFactors=F)

resultsdf=data.frame(matrix(nrow=nrow(full), ncol=length(allhoco)))
# snames=unlist(lapply(strsplit(allhoco, split="_"), function(x){return(x[1])}))
colnames(resultsdf)=allhoco
rownames(resultsdf)=full$V4

allresults=list(resultsdf,
                resultsdf,
                resultsdf,
                resultsdf,
                resultsdf,
                resultsdf)
names(allresults)=c("PGLSp", "PGLSstat", "PGLSbinp", "PGLSbinstat", "PGLSPAp", "PGLSPAstat")

#PGLS
# allresults=readRDS("allresultsfirst5.rds")
# allhoco=allhoco[6:26]
hcount=1
for(h in allhoco){
  fn=paste0("/home/kowaae22/TFcalls/", folder, "/", h, stat)
  data=read.table(fn, stringsAsFactors =F)
  colnames(data)=data[1,]
  data=data[-1,]
  colnames(data)[colnames(data)=="odoRosDiv1"]="odoRosDi"
  
  rownames(allresults$PGLSp)=data$name
  rownames(allresults$PGLSstat)=data$name
  rownames(allresults$PGLSbinp)=data$name
  rownames(allresults$PGLSbinstat)=data$name
  rownames(allresults$PGLSPAp)=data$name
  rownames(allresults$PGLSPAstat)=data$name
  
  count=1
  while(count<=nrow(data)){
    curcne=data$name[count]
    TF=setNames(data[count,], colnames(data))
    TF=TF[-c(1:4)]
    TF=TF[match(names(phenvec), names(TF))]
    TF=as.numeric(TF)
    df=data.frame(TF, phenvec)
    df2=na.omit(df)
    mt2=getTree(curcne)
    
    keep=intersect(rownames(df2), mt2$tip.label)
    
    if(!is.null(mt2)){
      df2=df2[rownames(df2) %in% keep,]
      mt2=keep.tip(mt2, keep)
    }
    
    df$PA=ifelse(is.na(df$TF), 0, 1)
    df2bin=df2
    df2bin$TF=ifelse(df2$TF==0,0,1)
    df=df[,2:3]
    #continuous PGLS
    if(length(unique(df2$TF))!=1 & length(unique(df2$phenvec))!=1 & nrow(unique(df2))>2 & !is.null(mt2)){
      
      pgls=gls(TF~phenvec, correlation = corBrownian(phy=mt2), data=df2)
      pvalPGLS=summary(pgls)$tTable[2,4]
      statPGLS=summary(pgls)$tTable[2,3]
      allresults$PGLSp[count, hcount]=pvalPGLS
      allresults$PGLSstat[count, hcount]=statPGLS 
    }else{
      allresults$PGLSp[count, hcount]=NA
      allresults$PGLSstat[count, hcount]=NA
    }
    
    #binary PGLS
    if(length(unique(df2bin$TF))!=1 & length(unique(df2bin$phenvec))!=1 & nrow(unique(df2bin))>2 & !is.null(mt2)){
      pglsbin=gls(TF~phenvec, correlation = corBrownian(phy=mt2), data=df2bin)
      pvalPGLSbin=summary(pgls)$tTable[2,4]
      statPGLSbin=summary(pgls)$tTable[2,3]
      allresults$PGLSbinp[count, hcount]=pvalPGLSbin
      allresults$PGLSbinstat[count, hcount]=statPGLSbin
    }else{
      allresults$PGLSbinp[count, hcount]=NA
      allresults$PGLSbinstat[count, hcount]=NA
    }
    
    if(count %% 10000==0){
      print(paste0("CNE count: ", count)) #345786
    }
    count=count+1
  }
  print(paste0("TF count: ", hcount)) #771
  saveRDS(allresults, paste0("/home/kowaae22/TFcalls/allresultsmergedcustomtreesLongevityPC1/allresultsmergedCTfirst",hcount,".rds"))
  hcount=hcount+1
}

###########################

#PC2 longevity analysis with custom merge and region-specific trees
###########################

specs=list.dirs("/home/Genomes/HOCOMOCO/Organism/", full.names = F)
specs=specs[-1]

# foreground=readRDS("/home/kowaae22/AnalysisWithThreeTrees/hairlessSpecs.rds")
PC2=readRDS("/home/kowaae22/AnalysisWithThreeTrees/PC2.rds")

# alltrees=readRDS("/home/kowaae22/AnalysisWithThreeTrees/finaltreesv4.rds")
# class(alltrees)[2]="treesObj"
library(nlme)
library(ape)
library(phytools)

allhoco=unique(list.files("/home/kowaae22/TFcalls/TFsummariescustommerge/hg19"))
# folder="groupspecsremovedups" #duplicates removed
folder="groupspecscustommergedups" #custom merge
stat="median"
# stat="mean"
# stat="min"
# stat="max"
# stat="count"
# stat="sum"

trees1=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part1.trees.rds")
trees2=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part2.trees.rds")
trees3=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part3.trees.rds")
trees4=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part4.trees.rds")
trees5=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part5.trees.rds")
trees6=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part6.trees.rds")
trees7=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part7.trees.rds")

getTree=function(regionname){
  if(regionname %in% names(trees1$trees)){
    ind=which(names(trees1$trees)==regionname)
    return(trees1$trees[[ind]])
  }
  if(regionname %in% names(trees2$trees)){
    ind=which(names(trees2$trees)==regionname)
    return(trees2$trees[[ind]])
  }
  if(regionname %in% names(trees3$trees)){
    ind=which(names(trees3$trees)==regionname)
    return(trees3$trees[[ind]])
  }
  if(regionname %in% names(trees4$trees)){
    ind=which(names(trees4$trees)==regionname)
    return(trees4$trees[[ind]])
  }
  if(regionname %in% names(trees5$trees)){
    ind=which(names(trees5$trees)==regionname)
    return(trees5$trees[[ind]])
  }
  if(regionname %in% names(trees6$trees)){
    ind=which(names(trees6$trees)==regionname)
    return(trees6$trees[[ind]])
  }
  if(regionname %in% names(trees7$trees)){
    ind=which(names(trees7$trees)==regionname)
    return(trees7$trees[[ind]])
  }
}
# mt=trees$masterTree

# pheno=readRDS("/home/kowaae22/AnalysisWithThreeTrees/hairlessSpecs.rds")
PC2=readRDS("/home/kowaae22/AnalysisWithThreeTrees/PC2.rds")
phenvec=PC2
# phenvec=rep(0, length(mt$tip.label))
# names(phenvec)=mt$tip.label
# phenvec[names(phenvec) %in% pheno]=1

full=read.table("/home/kowaae22/TFcalls/filteredcustommergedcoords/hg19coords", stringsAsFactors=F)

resultsdf=data.frame(matrix(nrow=nrow(full), ncol=length(allhoco)))
# snames=unlist(lapply(strsplit(allhoco, split="_"), function(x){return(x[1])}))
colnames(resultsdf)=allhoco
rownames(resultsdf)=full$V4

allresults=list(resultsdf,
                resultsdf,
                resultsdf,
                resultsdf,
                resultsdf,
                resultsdf)
names(allresults)=c("PGLSp", "PGLSstat", "PGLSbinp", "PGLSbinstat", "PGLSPAp", "PGLSPAstat")

#PGLS
# allresults=readRDS("allresultsfirst5.rds")
# allhoco=allhoco[6:26]
hcount=1
for(h in allhoco){
  fn=paste0("/home/kowaae22/TFcalls/", folder, "/", h, stat)
  data=read.table(fn, stringsAsFactors = F)
  colnames(data)=data[1,]
  data=data[-1,]
  colnames(data)[colnames(data)=="odoRosDiv1"]="odoRosDi"
  
  rownames(allresults$PGLSp)=data$name
  rownames(allresults$PGLSstat)=data$name
  rownames(allresults$PGLSbinp)=data$name
  rownames(allresults$PGLSbinstat)=data$name
  rownames(allresults$PGLSPAp)=data$name
  rownames(allresults$PGLSPAstat)=data$name
  
  count=1
  while(count<=nrow(data)){
    curcne=data$name[count]
    TF=setNames(data[count,], colnames(data))
    TF=TF[-c(1:4)]
    TF=TF[match(names(phenvec), names(TF))]
    TF=as.numeric(TF)
    df=data.frame(TF, phenvec)
    df2=na.omit(df)
    mt2=getTree(curcne)
    
    keep=intersect(rownames(df2), mt2$tip.label)
    
    if(!is.null(mt2)){
      df2=df2[rownames(df2) %in% keep,]
      mt2=keep.tip(mt2, keep)
    }
    
    df$PA=ifelse(is.na(df$TF), 0, 1)
    df2bin=df2
    df2bin$TF=ifelse(df2$TF==0,0,1)
    df=df[,2:3]
    #continuous PGLS
    if(length(unique(df2$TF))!=1 & length(unique(df2$phenvec))!=1 & nrow(unique(df2))>2 & !is.null(mt2)){
      
      pgls=gls(TF~phenvec, correlation = corBrownian(phy=mt2), data=df2)
      pvalPGLS=summary(pgls)$tTable[2,4]
      statPGLS=summary(pgls)$tTable[2,3]
      allresults$PGLSp[count, hcount]=pvalPGLS
      allresults$PGLSstat[count, hcount]=statPGLS 
    }else{
      allresults$PGLSp[count, hcount]=NA
      allresults$PGLSstat[count, hcount]=NA
    }
    
    #binary PGLS
    if(length(unique(df2bin$TF))!=1 & length(unique(df2bin$phenvec))!=1 & nrow(unique(df2bin))>2 & !is.null(mt2)){
      pglsbin=gls(TF~phenvec, correlation = corBrownian(phy=mt2), data=df2bin)
      pvalPGLSbin=summary(pgls)$tTable[2,4]
      statPGLSbin=summary(pgls)$tTable[2,3]
      allresults$PGLSbinp[count, hcount]=pvalPGLSbin
      allresults$PGLSbinstat[count, hcount]=statPGLSbin
    }else{
      allresults$PGLSbinp[count, hcount]=NA
      allresults$PGLSbinstat[count, hcount]=NA
    }
    
    if(count %% 10000==0){
      print(paste0("CNE count: ", count)) #345786
    }
    count=count+1
  }
  print(paste0("TF count: ", hcount)) #771
  saveRDS(allresults, paste0("/home/kowaae22/TFcalls/allresultsmergedcustomtreesLongevityPC2/allresultsmergedCTfirst",hcount,".rds"))
  hcount=hcount+1
}

###########################










#code you probably don't want to use, but I'll include it for the sake of completeness
#duplicates removed - included for completeness, but suboptimal strategy

#summarize by removing duplicates (rather than custom merging)
###########################

#one df for each: mean,median,sum,max,min,count
#rows as CNEs and cols as species, one file per TF

specs=list.dirs("/home/Genomes/HOCOMOCO/Organism/", full.names = F)
specs=specs[-1]

allhoco=unique(list.files("/home/kowaae22/TFcalls/TFsummaries/hg19"))
coords=read.table("/home/kowaae22/TFcalls/orderedcoords/hg19coords", stringsAsFactors=F)

n=c("name", "chr", "start", "end", specs)

count=1
for(h in allhoco){
  newfn=paste0("/home/kowaae22/TFcalls/groupspecsremovedups/", h)
  
  dfmean=coords[,1:4]
  dfmedian=coords[,1:4]
  dfsum=coords[,1:4]
  dfmax=coords[,1:4]
  dfmin=coords[,1:4]
  dfcount=coords[,1:4]
  for(s in specs){
    f=read.table(paste0("/home/kowaae22/TFcalls/TFsummaries/", s, "/", h), stringsAsFactors=F)
    f[f=="."]=0
    
    #remove duplicates
    c=table(f$V4)
    c=data.frame(c)
    c=setNames(c$Freq, c$Var1)
    c=sort(c, decreasing=T)
    dups=names(c)[c>1]
    f=f[!(f$V4 %in% dups),]
    
    #summarize
    f1=data.frame(f$V4, f$V7)
    dfmean=merge(dfmean, f1, by.x="V4", by.y="f.V4", all=T)
    
    f2=data.frame(f$V4, f$V8)
    dfmedian=merge(dfmedian, f2, by.x="V4", by.y="f.V4", all=T)
    
    f3=data.frame(f$V4, f$V9)
    dfsum=merge(dfsum, f3, by.x="V4", by.y="f.V4", all=T)
    
    f4=data.frame(f$V4, f$V10)
    dfmax=merge(dfmax, f4, by.x="V4", by.y="f.V4", all=T)
    
    f5=data.frame(f$V4, f$V11)
    dfmin=merge(dfmin, f5, by.x="V4", by.y="f.V4", all=T)
    
    f6=data.frame(f$V4, f$V12)
    dfcount=merge(dfcount, f6, by.x="V4", by.y="f.V4", all=T)
    
    print(s)
  }
  colnames(dfmean)=n
  colnames(dfmedian)=n
  colnames(dfsum)=n
  colnames(dfmax)=n
  colnames(dfmin)=n
  colnames(dfcount)=n
  
  write.table(dfmean, file=paste0(newfn, "mean"), quote=F, row.names=F, sep="\t")
  write.table(dfmedian, file=paste0(newfn, "median"), quote=F, row.names=F, sep="\t")
  write.table(dfsum, file=paste0(newfn, "sum"), quote=F, row.names=F, sep="\t")
  write.table(dfmax, file=paste0(newfn, "max"), quote=F, row.names=F, sep="\t")
  write.table(dfmin, file=paste0(newfn, "min"), quote=F, row.names=F, sep="\t")
  write.table(dfcount, file=paste0(newfn, "count"), quote=F, row.names=F, sep="\t")
  
  count=count+1
  print(paste0("allhoco: ", count))
}



# f=read.table(paste0("/home/kowaae22/TFcalls/TFsummaries/", s, "/", h), stringsAsFactors=F)
# c=table(f$V4)
# c=data.frame(c)
# c=setNames(c$Freq, c$Var1)
# c=sort(c, decreasing=T)
# dups=names(c)[c>1]
# 
# count=1
# s=0
# while(count<=length(dups)){
#   temp=f[f$V4==dups[count],]
#   s=s+sum(temp$V12)
#   print(count)
#   count=count+1
# }

###########################

#Analysis with removed duplicates
###########################

specs=list.dirs("/home/Genomes/HOCOMOCO/Organism/", full.names = F)
specs=specs[-1]

foreground=readRDS("/home/kowaae22/AnalysisWithThreeTrees/hairlessSpecs.rds")

alltrees=readRDS("/home/kowaae22/AnalysisWithThreeTrees/finaltreesv4.rds")
class(alltrees)[2]="treesObj"
library(nlme)
library(ape)
library(phytools)

allhoco=unique(list.files("/home/kowaae22/TFcalls/TFsummaries/hg19"))
folder="groupspecsremovedups" #duplicates removed
# folder="groupspecscustommergedups" #custom merge
stat="median"
# stat="mean"
# stat="min"
# stat="max"
# stat="count"
# stat="sum"

trees=readRDS("/home/kowaae22/AllEnhancers/finaltreesv4.rds")
mt=trees$masterTree

pheno=readRDS("/home/kowaae22/AnalysisWithThreeTrees/hairlessSpecs.rds")
phenvec=rep(0, length(mt$tip.label))
names(phenvec)=mt$tip.label
phenvec[names(phenvec) %in% pheno]=1

full=read.table("/home/kowaae22/TFcalls/orderedcoords/hg19coords", stringsAsFactors=F)

resultsdf=data.frame(matrix(nrow=nrow(full), ncol=length(allhoco)))
# snames=unlist(lapply(strsplit(allhoco, split="_"), function(x){return(x[1])}))
colnames(resultsdf)=allhoco
rownames(resultsdf)=full$V4

allresults=list(resultsdf,
                resultsdf,
                resultsdf,
                resultsdf,
                resultsdf,
                resultsdf)
names(allresults)=c("PGLSp", "PGLSstat", "PGLSbinp", "PGLSbinstat", "PGLSPAp", "PGLSPAstat")

#PGLS
# allresults=readRDS("allresultsfirst5.rds")
# allhoco=allhoco[6:26]
hcount=1
for(h in allhoco){
  fn=paste0("/home/kowaae22/TFcalls/", folder, "/", h, stat)
  data=read.table(fn, stringsAsFactors = F)
  colnames(data)=data[1,]
  data=data[-1,]
  colnames(data)[colnames(data)=="odoRosDiv1"]="odoRosDi"
  
  rownames(allresults$PGLSp)=data$name
  rownames(allresults$PGLSstat)=data$name
  rownames(allresults$PGLSbinp)=data$name
  rownames(allresults$PGLSbinstat)=data$name
  rownames(allresults$PGLSPAp)=data$name
  rownames(allresults$PGLSPAstat)=data$name
  
  count=1
  while(count<=nrow(data)){
    TF=setNames(data[count,], colnames(data))
    TF=TF[-c(1:4)]
    TF=TF[match(names(phenvec), names(TF))]
    TF=as.numeric(TF)
    df=data.frame(TF, phenvec)
    df2=na.omit(df)
    mt2=keep.tip(mt, rownames(df2))
    df$PA=ifelse(is.na(df$TF), 0, 1)
    df2bin=df2
    df2bin$TF=ifelse(df2$TF==0,0,1)
    df=df[,2:3]
    #continuous PGLS
    if(length(unique(df2$TF))!=1 & length(unique(df2$phenvec))!=1 & nrow(unique(df2))>2){
      
      pgls=gls(TF~phenvec, correlation = corBrownian(phy=mt2), data=df2)
      pvalPGLS=summary(pgls)$tTable[2,4]
      statPGLS=summary(pgls)$tTable[2,3]
      allresults$PGLSp[count, hcount]=pvalPGLS
      allresults$PGLSstat[count, hcount]=statPGLS 
    }else{
      allresults$PGLSp[count, hcount]=NA
      allresults$PGLSstat[count, hcount]=NA
    }
    
    #binary PGLS
    if(length(unique(df2bin$TF))!=1 & length(unique(df2bin$phenvec))!=1 & nrow(unique(df2bin))>2){
      pglsbin=gls(TF~phenvec, correlation = corBrownian(phy=mt2), data=df2bin)
      pvalPGLSbin=summary(pgls)$tTable[2,4]
      statPGLSbin=summary(pgls)$tTable[2,3]
      allresults$PGLSbinp[count, hcount]=pvalPGLSbin
      allresults$PGLSbinstat[count, hcount]=statPGLSbin
    }else{
      allresults$PGLSbinp[count, hcount]=NA
      allresults$PGLSbinstat[count, hcount]=NA
    }
    
    #presence/absence of CNE
    if(length(unique(df$PA))!=1 & length(unique(df$phenvec))!=1 & nrow(unique(df))>2){
      pglsPA=gls(PA~phenvec, correlation = corBrownian(phy=mt), data=df)
      pvalPGLSPA=summary(pglsPA)$tTable[2,4]
      statPGLSPA=summary(pglsPA)$tTable[2,3]
      allresults$PGLSPAp[count, hcount]=pvalPGLSPA
      allresults$PGLSPAstat[count, hcount]=statPGLSPA
    }else{
      allresults$PGLSPAp[count, hcount]=NA
      allresults$PGLSPAstat[count, hcount]=NA
    }
    
    if(count %% 10000==0){
      print(paste0("CNE count: ", count)) #345786
    }
    count=count+1
  }
  print(paste0("TF count: ", hcount)) #771
  saveRDS(allresults, paste0("/home/kowaae22/TFcalls/allresultsNEWfirst",hcount,".rds"))
  hcount=hcount+1
}

###########################

#Analysis with custom merge (but NOT region-specific trees)
###########################


specs=list.dirs("/home/Genomes/HOCOMOCO/Organism/", full.names = F)
specs=specs[-1]

foreground=readRDS("/home/kowaae22/AnalysisWithThreeTrees/hairlessSpecs.rds")

alltrees=readRDS("/home/kowaae22/AnalysisWithThreeTrees/finaltreesv4.rds")
class(alltrees)[2]="treesObj"
library(nlme)
library(ape)
library(phytools)

allhoco=unique(list.files("/home/kowaae22/TFcalls/TFsummariescustommerge/hg19"))
# folder="groupspecsremovedups" #duplicates removed
folder="groupspecscustommergedups" #custom merge
stat="median"
# stat="mean"
# stat="min"
# stat="max"
# stat="count"
# stat="sum"

trees=readRDS("/home/kowaae22/AllEnhancers/finaltreesv4.rds")
mt=trees$masterTree

pheno=readRDS("/home/kowaae22/AnalysisWithThreeTrees/hairlessSpecs.rds")
phenvec=rep(0, length(mt$tip.label))
names(phenvec)=mt$tip.label
phenvec[names(phenvec) %in% pheno]=1

full=read.table("/home/kowaae22/TFcalls/filteredcustommergedcoords/hg19coords", stringsAsFactors=F)

resultsdf=data.frame(matrix(nrow=nrow(full), ncol=length(allhoco)))
# snames=unlist(lapply(strsplit(allhoco, split="_"), function(x){return(x[1])}))
colnames(resultsdf)=allhoco
rownames(resultsdf)=full$V4

allresults=list(resultsdf,
                resultsdf,
                resultsdf,
                resultsdf,
                resultsdf,
                resultsdf)
names(allresults)=c("PGLSp", "PGLSstat", "PGLSbinp", "PGLSbinstat", "PGLSPAp", "PGLSPAstat")

#PGLS
# allresults=readRDS("allresultsfirst5.rds")
# allhoco=allhoco[6:26]
hcount=1
for(h in allhoco){
  fn=paste0("/home/kowaae22/TFcalls/", folder, "/", h, stat)
  data=read.table(fn, stringsAsFactors = F)
  colnames(data)=data[1,]
  data=data[-1,]
  colnames(data)[colnames(data)=="odoRosDiv1"]="odoRosDi"
  
  rownames(allresults$PGLSp)=data$name
  rownames(allresults$PGLSstat)=data$name
  rownames(allresults$PGLSbinp)=data$name
  rownames(allresults$PGLSbinstat)=data$name
  rownames(allresults$PGLSPAp)=data$name
  rownames(allresults$PGLSPAstat)=data$name
  
  count=1
  while(count<=nrow(data)){
    TF=setNames(data[count,], colnames(data))
    TF=TF[-c(1:4)]
    TF=TF[match(names(phenvec), names(TF))]
    TF=as.numeric(TF)
    df=data.frame(TF, phenvec)
    df2=na.omit(df)
    mt2=keep.tip(mt, rownames(df2))
    df$PA=ifelse(is.na(df$TF), 0, 1)
    df2bin=df2
    df2bin$TF=ifelse(df2$TF==0,0,1)
    df=df[,2:3]
    #continuous PGLS
    if(length(unique(df2$TF))!=1 & length(unique(df2$phenvec))!=1 & nrow(unique(df2))>2){
      
      pgls=gls(TF~phenvec, correlation = corBrownian(phy=mt2), data=df2)
      pvalPGLS=summary(pgls)$tTable[2,4]
      statPGLS=summary(pgls)$tTable[2,3]
      allresults$PGLSp[count, hcount]=pvalPGLS
      allresults$PGLSstat[count, hcount]=statPGLS 
    }else{
      allresults$PGLSp[count, hcount]=NA
      allresults$PGLSstat[count, hcount]=NA
    }
    
    #binary PGLS
    if(length(unique(df2bin$TF))!=1 & length(unique(df2bin$phenvec))!=1 & nrow(unique(df2bin))>2){
      pglsbin=gls(TF~phenvec, correlation = corBrownian(phy=mt2), data=df2bin)
      pvalPGLSbin=summary(pgls)$tTable[2,4]
      statPGLSbin=summary(pgls)$tTable[2,3]
      allresults$PGLSbinp[count, hcount]=pvalPGLSbin
      allresults$PGLSbinstat[count, hcount]=statPGLSbin
    }else{
      allresults$PGLSbinp[count, hcount]=NA
      allresults$PGLSbinstat[count, hcount]=NA
    }
    
    #presence/absence of CNE
    if(length(unique(df$PA))!=1 & length(unique(df$phenvec))!=1 & nrow(unique(df))>2){
      pglsPA=gls(PA~phenvec, correlation = corBrownian(phy=mt), data=df)
      pvalPGLSPA=summary(pglsPA)$tTable[2,4]
      statPGLSPA=summary(pglsPA)$tTable[2,3]
      allresults$PGLSPAp[count, hcount]=pvalPGLSPA
      allresults$PGLSPAstat[count, hcount]=statPGLSPA
    }else{
      allresults$PGLSPAp[count, hcount]=NA
      allresults$PGLSPAstat[count, hcount]=NA
    }
    
    if(count %% 10000==0){
      print(paste0("CNE count: ", count)) #345786
    }
    count=count+1
  }
  print(paste0("TF count: ", hcount)) #771
  saveRDS(allresults, paste0("/home/kowaae22/TFcalls/allresultsmerged/allresultsmergedfirst",hcount,".rds"))
  hcount=hcount+1
}

###########################


















