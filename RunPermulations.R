#this script includes examples of running permulations for RERconverge and PGLS
#notes and descriptions in this script are limited
#because full vignettes are available on GitHub

#functions
#####################################
#function to permulate phenotype (continuous)
simpermvec=function(namedvec, treewithbranchlengths){
  #returns sim/perm vec
  #tree must be rooted and fully dichotomous
  #species in tree must match species in vec
  #simulate vector
  vec=simulatevec(namedvec, treewithbranchlengths)
  
  #assign real values to vec
  simsorted=sort(vec)
  realsorted=sort(namedvec)
  l=length(simsorted)
  c=1
  while(c<=l){
    simsorted[c]=realsorted[c]
    c=c+1
  }
  simsorted
}

#function to simulate phenotype (continuous)
simulatevec=function(namedvec, treewithbranchlengths){
  #returns simulated vec
  #tree must be rooted and fully dichotomous
  #species in tree must match species in vec
  library("geiger")
  rm=ratematrix(treewithbranchlengths, namedvec)
  sims=sim.char(treewithbranchlengths, rm, nsim = 1)
  nam=rownames(sims)
  s=as.data.frame(sims)
  simulatedvec=s[,1]
  names(simulatedvec)=nam
  vec=simulatedvec
  vec
}

#function to permute phenotype (continuous or binary)
permutevec=function(namedvec){
  #returns permuted vec
  n=names(namedvec)
  vec=sample(namedvec)
  names(vec)=n
  vec
}

#function to permulate phenotype (binary)
simBinPheno=function(trees, root, phenvec, fgnum=NULL, internal=0, drop=NULL){
  blsum=0
  if(is.null(fgnum)){
    fgnum=sum(phenvec)
  }
  tips=fgnum-internal
  while(blsum!=fgnum){
    t=root.phylo(trees$masterTree, root, resolve.root = T)
    t=drop.tip(t, drop)
    rm=ratematrix(t, phenvec)
    sims=sim.char(t, rm, nsim = 1)
    nam=rownames(sims)
    s=as.data.frame(sims)
    simulatedvec=s[,1]
    names(simulatedvec)=nam
    top=names(sort(simulatedvec, decreasing = TRUE))[1:tips]
    t=foreground2Tree(top, trees, clade="all", plotTree = F)
    blsum=sum(t$edge.length)
  }
  # plot(t)
  return(t)
}

#####################################


#permulations with RERconverge using a continuous phenotype
#####################################
trees=readRDS("/home/kowaae22/100way/promotertrees/promotertrees.rds")
RERs=readRDS("/home/kowaae22/100way/RERs/promoterRERslongevitySpecs.rds")
annots=readRDS("/home/kowaae22/Annotations/fullcodingannots.rds")
res=readRDS("/home/kowaae22/100way/RERanalysisresults/PC1cors.rds")
enrichment=readRDS("/home/kowaae22/100way/RERanalysisresults/PC1enrich.rds")
PC1=readRDS("/home/kowaae22/AnalysisWithThreeTrees/PC1.rds")

mt=trees$masterTree
mt=drop.tip(mt, "chrAsi1")
mt=root.phylo(mt, outgroup="ornAna1", resolve.root=T)

perms=RERconverge::getPermsContinuous(1000, PC1, RERs, annots, trees, mt)
saveRDS(perms, "/home/kowaae22/100way/RERanalysisresults/PC1perms.rds")

corpermpvals=RERconverge::permpvalcor(res, perms)
saveRDS(corpermpvals, "/home/kowaae22/100way/RERanalysisresults/PC1correlationpermp.rds")
enrichpermpvals=RERconverge::permpvalenrich(enrichment, perms)
saveRDS(enrichpermpvals, "/home/kowaae22/100way/RERanalysisresults/PC1enrichpermp.rds")

res$permpval=corpermpvals[match(rownames(res), names(corpermpvals))]
res$permpvaladj=p.adjust(res$permpval, method="BH")
saveRDS(res, "/home/kowaae22/100way/RERanalysisresults/PC1corwithpermp.rds")
count=1
while(count<=length(enrichment)){
  enrichment[[count]]$permpval=enrichpermpvals[[count]][match(rownames(enrichment[[count]]),
                                                              names(enrichpermpvals[[count]]))]
  enrichment[[count]]$permpvaladj=p.adjust(enrichment[[count]]$permpval, method="BH")
  count=count+1
}
saveRDS(enrichment, "/home/kowaae22/100way/RERanalysisresults/PC1enrichwithpermp.rds")
#####################################


#permulations with RERconverge using a binary phenotype
#####################################
trees=readRDS("/home/kowaae22/100way/promotertrees/promotertrees.rds")
RERs=readRDS("/home/kowaae22/100way/RERs/promoterRERsallSpecsweightresid.rds")
annots=readRDS("/home/kowaae22/Annotations/fullcodingannots.rds")
res=readRDS("/home/kowaae22/100way/RERanalysisresults/hairlesscors.rds")
enrichment=readRDS("/home/kowaae22/100way/RERanalysisresults/hairlessenrich.rds")

fg=readRDS("/home/kowaae22/AnalysisWithThreeTrees/hairlessSpecs.rds")
s=list(clade1=c("orcOrc1", "turTru2"))

perms=getPermsBinary(1000, fg, s, "ornAna1", RERs, trees, trees$masterTree, permmode="cc",calculateenrich=T,annotlist=annots)
saveRDS(perms, "/home/kowaae22/100way/RERanalysisresults/hairlessperms.rds")
permpcor = permpvalcor(res,perms)
saveRDS(permpcor, "/home/kowaae22/100way/RERanalysisresults/hairlesscorrelationpermp.rds")
enrichpermpvals=permpvalenrich(enrichment, perms)
saveRDS(enrichpermpvals, "/home/kowaae22/100way/RERanalysisresults/hairlessenrichpermp.rds")

res$permpval=permpcor[match(rownames(res), names(permpcor))]
res$permpvaladj=p.adjust(res$permpval, method="BH")
saveRDS(res, "/home/kowaae22/100way/RERanalysisresults/hairlesscorwithpermp.rds")
count=1
while(count<=length(enrichment)){
  enrichment[[count]]$permpval=enrichpermpvals[[count]][match(rownames(enrichment[[count]]),
                                                              names(enrichpermpvals[[count]]))]
  enrichment[[count]]$permpvaladj=p.adjust(enrichment[[count]]$permpval, method="BH")
  count=count+1
}
saveRDS(enrichment, "/home/kowaae22/100way/RERanalysisresults/hairlessenrichwithpermp.rds")

#####################################


#permulations with PGLS and a continuous phenotype - this will take a very long time
########################################

numperms=500

full=read.table("/home/kowaae22/TFcalls/filteredcustommergedcoords/hg19coords", stringsAsFactors=F)
full=full[full$V1=="chr1",]
resultsdf=data.frame(matrix(nrow=nrow(full), ncol=numperms))
rownames(resultsdf)=full$V4


allresults=list(resultsdf,
                resultsdf)
names(allresults)=c("PGLSp", "PGLSstat")




start=Sys.time()
hcount=1
while(hcount<=numperms){
  h=statfns[1]
  fn=paste0("/home/kowaae22/TFcalls/", folder, "/", h, stat)
  data=read.table(fn, stringsAsFactors =F)
  colnames(data)=data[1,]
  data=data[-1,]
  colnames(data)[colnames(data)=="odoRosDiv1"]="odoRosDi"
  data=data[data$chr=="chr1",]
  
  
  rownames(allresults$PGLSp)=data$name
  rownames(allresults$PGLSstat)=data$name
  
  
  #permulate phenotype
  rtmt=root.phylo(trees1$masterTree, outgroup = "ornAna1", resolve.root = T)
  rtmt=drop.tip(rtmt, c("sgal", "chrAsi1"))
  phenvec=simpermvec(PC1, rtmt)
  
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
    
    if(count %% 10000==0){
      print(paste0("CNE count: ", count)) #345786
    }
    count=count+1
  }
  print(paste0("perm count: ", hcount)) #771
  # saveRDS(allresults, paste0("/home/kowaae22/TFcalls/allresultsmergedcustomtreesLongevityPC1/allresultsmergedCTfirst",hcount,".rds"))
  hcount=hcount+1
}
end=Sys.time()
end-start

saveRDS(allresults, "/home/kowaae22/permPGLSpermsSTAT2PC1count.rds")

########################################





