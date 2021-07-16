#this script shows how to calculate Bayes factors from RERs and a phenotype
#Bayes factors indicate the level of support for a particular model (between 0 and 1 indicates more support for the null model)
#ratios of Bayes factors indicate relative support for one model over another
#bigger Bayes factor = better; different sources give different thresholds for a "good" Bayes factor

#what you need:
#RERs
#phenotype tree
#trees

#what you need to change in this script:
#filenames
#phenotype names/labels


#read in phenotype trees
marinetree=readRDS("marinetree.rds")
hairlesstree=readRDS("statetree.rds")
hairlessnomarinetree=readRDS("hairlessnomarinetree.rds")
marinenohairlesstree=readRDS("marinenohairlesstree.rds")

#read in trees and RERs
trees=readRDS("finaltreesv4.rds")
class(trees)[2]="treesObj"
newRERmat=readRDS("weightresidRER.rds")
#create phenotype paths using trees and phenotype trees
marinepath=tree2Paths(marinetree, trees)
hairlesspath=tree2Paths(hairlesstree, trees)
hairlessnomarinepath=tree2Paths(hairlessnomarinetree, trees)
marinenohairlesspath=tree2Paths(marinenohairlesstree, trees)

#make dataframe to store Bayes factors
library(BayesFactor)
allBayesFactors=data.frame(matrix(nrow=nrow(RERmat), ncol=4))
colnames(allBayesFactors)=c("marine", "hairless", "hairlessnomarine", "marinenohairless")
rownames(allBayesFactors)=rownames(newRERmat)

#calculate Bayes factors per gene
count=1
while(count<=nrow(newRERmat)){
  #marine
  tempdf=data.frame(trait=marinepath, RER=newRERmat[count,])
  tempdf=na.omit(tempdf)
  bfmod=lmBF(trait~RER, data=tempdf)
  bf=as.numeric(extractBF(bfmod)[1])
  allBayesFactors$marine[count]=bf
  
  #hairless
  tempdf=data.frame(trait=hairlesspath, RER=newRERmat[count,])
  tempdf=na.omit(tempdf)
  bfmod=lmBF(trait~RER, data=tempdf)
  bf=as.numeric(extractBF(bfmod)[1])
  allBayesFactors$hairless[count]=bf 
  
  #hairlessnomarine
  tempdf=data.frame(trait=hairlessnomarinepath, RER=newRERmat[count,])
  tempdf=na.omit(tempdf)
  bfmod=lmBF(trait~RER, data=tempdf)
  bf=as.numeric(extractBF(bfmod)[1])
  allBayesFactors$hairlessnomarine[count]=bf 
  
  #marinenohairless
  tempdf=data.frame(trait=marinenohairlesspath, RER=newRERmat[count,])
  tempdf=na.omit(tempdf)
  bfmod=lmBF(trait~RER, data=tempdf)
  bf=as.numeric(extractBF(bfmod)[1])
  allBayesFactors$marinenohairless[count]=bf
  
  count=count+1
  print(count)
}
saveRDS(allBayesFactors, "allBayesFactors.rds")


