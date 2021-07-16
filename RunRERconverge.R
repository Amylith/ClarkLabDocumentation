#this script includes several examples of running RERconverge
#notes and descriptions in this script are limited
#because full vignettes are available on GitHub

#topics included in this script:
#   *full continuous trait analysis
#   *creating RERs regressed by a continuous phenotype



#full continuous trait analysis
############################################################
#readTrees - this takes a long time
library(RERconverge)
treespromoters100way=readTrees("/home/kowaae22/100way/promotertrees/promotertrees.trees")
saveRDS(treespromoters100way, "/home/kowaae22/100way/promotertrees/promotertrees.rds")

#make RERs - this also takes a long time
trees=readRDS("/home/kowaae22/100way/promotertrees/promotertrees.rds")
RERallSpecs=getAllResiduals(trees, plot=F)
saveRDS(RERallSpecs, "/home/kowaae22/100way/RERs/promoterRERsallSpecs.rds")
pheno=readRDS("/home/kowaae22/AnalysisWithThreeTrees/PC1.rds")
longevityspecs=names(pheno)
RERlongevitySpecs=getAllResiduals(trees, plot=F, useSpecies = longevityspecs)
saveRDS(RERlongevitySpecs, "/home/kowaae22/100way/RERs/promoterRERslongevitySpecs.rds")

#run RERconverge - this is relatively faster
library(RERconverge)
trees=readRDS("/home/kowaae22/100way/promotertrees/promotertrees.rds")
RERallSpecs=readRDS("/home/kowaae22/100way/RERs/promoterRERsallSpecs.rds")
RERlongevitySpecs=readRDS("/home/kowaae22/100way/RERs/promoterRERslongevitySpecs.rds")
annots=readRDS("/home/kowaae22/Annotations/annotswithGO.RDS")
gtex=readRDS("/home/kowaae22/Annotations/collapsedTissueAnnots.RDS")
follicle=readRDS("/home/kowaae22/Annotations/expgenes.rds")
folliclecorrected=readRDS("/home/kowaae22/Annotations/expgenesnew.rds")
allexp=readRDS("/home/kowaae22/Annotations/allexpgenes.rds")
haircompartment=readRDS("/home/kowaae22/Annotations/haircompartmentgenes.rds")
annots$gtex=gtex$tissueannots
annots$expression=list(genesets=list(follicle=follicle,folliclecorrected=folliclecorrected, allexp=allexp), geneset.names=c("follicle", "folliclecorrected", "allexp"))
annots$haircompartments=haircompartment$haircompartment
saveRDS(annots, "/home/kowaae22/Annotations/fullcodingannots.rds")

PC1=readRDS("/home/kowaae22/AnalysisWithThreeTrees/PC1.rds")
PC1path=char2Paths(PC1, trees)
PC1cors=correlateWithContinuousPhenotype(RERlongevitySpecs, PC1path)
PC1enrich=fastwilcoxGMTall(getStat(PC1cors), annots, outputGeneVals=T)
saveRDS(PC1cors, "/home/kowaae22/100way/RERanalysisresults/PC1cors.rds")
saveRDS(PC1enrich, "/home/kowaae22/100way/RERanalysisresults/PC1enrich.rds")

PC2=readRDS("/home/kowaae22/AnalysisWithThreeTrees/PC2.rds")
PC2path=char2Paths(PC2, trees)
PC2cors=correlateWithContinuousPhenotype(RERlongevitySpecs, PC2path)
PC2enrich=fastwilcoxGMTall(getStat(PC2cors), annots, outputGeneVals=T)
saveRDS(PC2cors, "/home/kowaae22/100way/RERanalysisresults/PC2cors.rds")
saveRDS(PC2enrich, "/home/kowaae22/100way/RERanalysisresults/PC2enrich.rds")
############################################################



#binary analysis - RERs regressed based on a continuous phenotype
############################################################

#read in phenotype values and keep species contained in binary phenotype of interest
library(xlsx)
phen=read.xlsx("/home/kowaae22/AnalysisWithThreeTrees/SpeciesKeyWInfo.xlsx", sheetIndex = 1)
mass=phen$Adult.weight..g.
names(mass)=phen[,1]
mass=mass[!is.na(mass)]
statetree=readRDS("/home/kowaae22/AnalysisWithThreeTrees/statetree")
speciestouse=statetree$tip.label
mass=mass[speciestouse]

#read in RERs and trees
RERmat=readRDS("/home/kowaae22/100way/RERs/promoterRERsallSpecs.rds")
trees=readRDS("/home/kowaae22/100way/promotertrees/promotertrees.rds")

#create path based on continuous phenotype to regress out
library(RERconverge)
weightpath=char2Paths(mass, trees)

#regress out based on continuous phenotype
newRERmat=RERmat
count=1
while(count<=nrow(newRERmat)){
  if(sum(!is.na(newRERmat[count,]))>0){
    mod=lm(newRERmat[count,]~weightpath)
    inds=which(!is.na(newRERmat[count,]))
    inds=inds[!is.na(inds)]
    newvec=rep(NA, ncol(newRERmat))
    newvec[inds]=mod$residuals
    newRERmat[count,]=newvec
    print(count)
  }
  count=count+1
}
saveRDS(newRERmat,"/home/kowaae22/100way/RERs/promoterRERsallSpecsweightresid.rds")

#run RERconverge using regressed RERs:
RERwr=readRDS("/home/kowaae22/100way/RERs/promoterRERsallSpecsweightresid.rds")
hairlesstree=readRDS("/home/kowaae22/AnalysisWithThreeTrees/statetree")
hairlesspath=tree2Paths(hairlesstree, trees)
hairlesscors=correlateWithBinaryPhenotype(RERwr, hairlesspath)
hairlessenrich=fastwilcoxGMTall(getStat(hairlesscors), annots, outputGeneVals=T)
saveRDS(hairlesscors, "/home/kowaae22/100way/RERanalysisresults/hairlesscors.rds")
saveRDS(hairlessenrich, "/home/kowaae22/100way/RERanalysisresults/hairlessenrich.rds")

############################################################








