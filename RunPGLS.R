#this script demonstrates how to run PGLS to scan for coevolution between TFBS score and a phenotype
#across all ~350k conserved noncoding regions in the UCSC 100way alignment mammals
#and a chosen TFBS motif

#what you need:
#trees from RERconverge readTrees (for 350k regions, they are split into several batches)
#phenotype
#output from TFcalls calculated over regions

#what you need to change in this script:
#filenames
#alter loop (if desired) to perform calculation over all TFs



#function - get the set of trees that contains a given CNE
###############################
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

###############################


#load packages
library(nlme)
library(phytools)
library(RERconverge)

#get species to use based on those with TF calls
specs=list.dirs("/home/Genomes/HOCOMOCO/Organism/", full.names = F)
specs=specs[-1]

#read in phenotype
PC1=readRDS("/home/kowaae22/AnalysisWithThreeTrees/PC1.rds")
phenvec=PC1

#name of noncoding region to investigate
# statfns=c("STAT1_HUMAN.H11MO.0.A.bed",
#           "STAT1_HUMAN.H11MO.1.A.bed",
#           "STAT2_HUMAN.H11MO.0.A.bed",
#           "STAT3_HUMAN.H11MO.0.A.bed",
#           "STAT4_HUMAN.H11MO.0.A.bed",
#           "STAT6_HUMAN.H11MO.0.B.bed")
statfns=c("STAT2_HUMAN.H11MO.0.A.bed")

#folder identity that contains TFBS scores + which statistic to use
folder="groupspecscustommergedups" #custom merge
# stat="median"
# stat="mean"
# stat="min"
# stat="max"
stat="count"
# stat="sum"

#read in tree batches for noncoding regions
trees1=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part1.trees.rds")
trees2=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part2.trees.rds")
trees3=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part3.trees.rds")
trees4=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part4.trees.rds")
trees5=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part5.trees.rds")
trees6=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part6.trees.rds")
trees7=readRDS("/home/kowaae22/AllEnhancers/trees/phastcons46way.final.wSpalax.part7.trees.rds")



#run PGLS
########################################

#create a dataframe of the correct size to contain results
full=read.table("/home/kowaae22/TFcalls/filteredcustommergedcoords/hg19coords", stringsAsFactors=F)
full=full[full$V1=="chr1",]
resultsdf=data.frame(matrix(nrow=nrow(full), ncol=length(statfns)))
colnames(resultsdf)=statfns
rownames(resultsdf)=full$V4

#create list to store results
allresults=list(resultsdf,
                resultsdf)
names(allresults)=c("PGLSp", "PGLSstat")




start=Sys.time()
hcount=1
for(h in statfns){
  #read in TF call information
  fn=paste0("/home/kowaae22/TFcalls/", folder, "/", h, stat)
  data=read.table(fn, stringsAsFactors =F)
  colnames(data)=data[1,]
  data=data[-1,]
  colnames(data)[colnames(data)=="odoRosDiv1"]="odoRosDi"
  data=data[data$chr=="chr1",]
  
  
  rownames(allresults$PGLSp)=data$name
  rownames(allresults$PGLSstat)=data$name
  
  #loop over all CNEs:
  count=1
  while(count<=nrow(data)){
    #match species content in CNE species tree and phenvec
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
    
    #run PGLS if there are enough species + enough diversity in phenotype and TF across those species
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
    
    if(count %% 10000==0){
      print(paste0("CNE count: ", count)) #345786
    }
    count=count+1
  }
  print(paste0("TF count: ", hcount)) #771
  # saveRDS(allresults, paste0("/home/kowaae22/TFcalls/allresultsmergedcustomtreesLongevityPC1/allresultsmergedCTfirst",hcount,".rds"))
  hcount=hcount+1
}
end=Sys.time()
end-start

saveRDS(allresults, "/home/kowaae22/permPGLSrealSTAT2PC1count.rds")

########################################