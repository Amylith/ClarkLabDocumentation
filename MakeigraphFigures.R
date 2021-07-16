#this script shows how to get make igraph clusters from RERconverge enrichment results
#this step is important because pathway annotations are extremely redundant
#clustering allows us to group together pathways that contain similar genes
#and identify overarching functional themes in a less biased way
#This process involves two steps:
#1) properly format node and edge weights/labels
#2) make the igraph figure
#IMPORTANT: there are points in this script where you need to stop and check something (they say CHECK)
#make sure to SAVE YOUR SCRIPTS AND WORKSPACE (save.image()) before running the code after the checks
#optimization can cause R to crash, and if you don't save before that, you will lose your work

#what you need:
#RERconverge enrichment results run with outputGeneVals=T
#(optional but strongly recommended) enrichment permulations
#OR
#any definition of nodes and edges you want to create a network from (like ERC, for example) - this will require more modification of this code

#what you need to change in this script:
#filenames
#cutoffs for curenrich to create nodes and edges (remove permulation check if you didn't use permulations)
#places where it says CHANGE THIS: edge weight scaling factor, node scaling factor


#functions to extract genes from RERconverge enrich + count matching genes between two lists
##################################

getGeneList=function(genestring){
  genes=unlist(strsplit(genestring, split=", "))
  ng=c()
  for(g in genes){
    nocol=strsplit(g, split=":")[[1]][[1]]
    ng=c(ng, nocol)
  }
  genes=list(ng)
  genes[[1]]
}

numberMatchingGenes=function(genelist1, genelist2){
  m=match(genelist1, genelist2)
  nummatch=length(m[!is.na(m)])
  nummatch
}
##################################


#function to make files with node vals and edge weights
##################################
writenodesandedges=function(curenrich, edgesfn, nodesfn){
  #make df without group numbers
  goodnum=nrow(curenrich)*nrow(curenrich)*2 #max number of edges if everything is connected to everything else
  #columns: stat, pval, simpermPval, p.adj, # of genes, top ten genes and rank
  groupedPaths=data.frame(matrix(ncol = 3, nrow = goodnum))
  colnames(groupedPaths)=c("Pathway1", "edgeweight", "Pathway2")
  enrichnames=rownames(curenrich)
  
  pathstats=data.frame(rownames(curenrich), curenrich$stat)
  
  notgrouped=c() #track unique pathways
  #for each pathway, look at every other pathway after the current pathway (so we don't get duplicates)
  dfcount=1
  count=1
  while(count<=nrow(curenrich)){
    curpathway=curenrich[count,]$gene.vals
    curpathway=getGeneList(curpathway)
    curstat=curenrich[count,]$stat
    curname=enrichnames[count]
    if(count==1){
      firstnm=curname
    }
    count2=count
    flag=F
    while(count2<=nrow(curenrich)){ 
      if(count2!=count){ #don't compare a pathway to itself
        comppathway=curenrich[count2,]$gene.vals
        comppathway=getGeneList(comppathway)
        compname=enrichnames[count2]
        compstat=curenrich[count2,]$stat
        denom=min(length(comppathway), length(curpathway))
        simscore=numberMatchingGenes(curpathway, comppathway)/denom
        if(simscore>0){ #only make edge if the pathways share genes
          flag=T
          groupedPaths[dfcount,1]=curname
          groupedPaths[dfcount,2]=simscore
          groupedPaths[dfcount,3]=compname
          dfcount=dfcount+1
        }
      }
      count2=count2+1
    }
    if(!flag){
      notgrouped=c(notgrouped, curname)
      groupedPaths[dfcount,1]=curname
      groupedPaths[dfcount,2]=0
      groupedPaths[dfcount,3]=firstnm
      dfcount=dfcount+1
    }
    message(paste0("finished pathway ", count," of ", nrow(curenrich)))
    count=count+1
  }
  
  groupedPaths=groupedPaths[!is.na(groupedPaths$Pathway1),]
  
  
  # edgesfn=paste0("./igraphfiles/edgeinfo",fn)
  # nodesfn=paste0("./igraphfiles/nodeinfo",fn)
  write.csv(groupedPaths, edgesfn, row.names = F)
  write.csv(pathstats, nodesfn, row.names = F)
  return(list(edges=groupedPaths, nodes=pathstats))
}
##################################


#example of getting nodes and edges from RERconverge enrichment:
##################################
curenrich=readRDS("RERconvergeResults/PC1enrich120.rds")
curenrich=curenrich$mgi

#MGI "mammalian phenotype" contains all genes annotated through MGI - remove it because it's not informative
if("mammalian phenotype" %in% rownames(curenrich)){
  r=which(rownames(curenrich)=="mammalian phenotype")
  curenrich=curenrich[-r,]
}

#subset enrichment to a smaller set, the actual points that you want to plot (for example, only significant results)
#if your network has too many points, clustering won't run and the plot will be a hairball
curenrich=curenrich[curenrich$stat<0 & curenrich$permpval<=0.02,]

nrow(curenrich) 
#CHECK - does your enrichment contain too many points? I usually use no more than 70
#you will have to play around with cutoffs for curenrich to see what works for your data
#also, if you have too many entries in curenrich, writenodesandedges will take too long to run

ne=writenodesandedges(curenrich, "igraphfiles/MGIedgesPC1neg.csv", "igraphfiles/MGInodesPC1neg.csv")
##################################



#example of making igraph figure with clustering:
##################################
#packages used to make networks:
library(igraph)
library(network)
library(sna)
library(stringr)
library(RColorBrewer)
#blog with more information/more examples:
#http://kateto.net/network-visualization
# for help do ?plot.igraph

par(mfrow=c(1,2))
#function to make simple igraph plot
plotNetSimple=function(net, ...){V(net)$label <- NA;plot(net, edge.arrow.size=0, edge.color="grey30",  vertex.label.color="#000000", ...)}

#read table of edges created previously
netTable=read.csv("igraphfiles/CANedgesPC1pos.csv",header=T)

#reformat table of edges to make into network
netTable=netTable[,c(1,3,2)]
net=graph_from_data_frame(netTable[,])

#CHANGE THIS: scaling factors for edge values (pick one so you can see all edges without them being too bulky)
E(net)$width=netTable$edgeweight^3

#read in node information created previously
statTable=read.csv("igraphfiles/CANnodesPC1pos.csv")

#create color maping from node values:
iim=match(names(V(net)), statTable[,1])
mycol=rev(brewer.pal(8,"Blues"))
V(net)$color=mycol[cut(statTable[iim,2],8)]

#create size mapping from node values
#CHANGE THIS: scaling factors for node sizes
V(net)$size=abs(statTable[iim,2])*20

#create cutoff to determine which edges will be contained in the figure
#CHANGE THIS: if you have too many edges, clustering will not work + you'll end up with a hairball
cut.off <- quantile(netTable$edgeweight, .95)

#remove edges below cutoff
net.sp <- delete_edges(net, E(net)[netTable$edgeweight<cut.off])

#make simple plot
plotNetSimple(net.sp)

#THIS IS VERY IMPORTANT:
#check the simple plot
#is it too complicated?
#if so, reduce number of nodes (from previous section) and/or edges (from this section) and rerun everything
#regardless, SAVE YOUR WORK before running the next part - save scripts and run save.image() to save workspace
#R may crash/get stuck on the next line

#runs optimal clustering:
clp <- cluster_optimal(net.sp)
#if optimal clustering doesn't work (for example if you want to include many nodes/edges), 
#you can try other clustering methods:
# clp=cluster_fast_greedy(net.sp) #fast and greedy clustering
# clp=cluster_louvain(net.sp) #louvain clustering
# clp=cluster_walktrap(net.sp) #random walk

#plots igraph network with clusters:
plot(clp, net.sp, vertex.label=NA)


#at this point, your network plot is made
#subsequent code is different visualization options
#and extracting cluster membership


#which nodes are in clusters?
tt=table(clp$membership)
iimodules=which(tt>2)
inmodules=clp$membership %in% iimodules
inmodules

#make new network that only contains clusters (no singletons)
net.sp.clust=delete_vertices(net.sp, V(net.sp)[which(!inmodules)])

#make cluster list structure to label with cluster membership (tell which color is which cluster number)
vl=list()
mem=clp$membership[inmodules]
mem=as.numeric(as.factor(mem))
for(i in 1:max(mem)){
  vl[[i]]=V(net.sp.clust)[mem==i]
}

#make singleton graph
net.sp.single=delete_vertices(net.sp, V(net.sp)[which(inmodules)])
#define a grid layout
n=length(V(net.sp.single))
l=matrix(nrow=n, ncol=2)
l[,1]=1
l[,2]=1:n
dev.off()
par(mfrow=c(1,2), mai=rep(0,4))
plotNetSimple(net.sp.clust, mark.groups=vl)
#with padding
plotNetSimple(net.sp.single, layout=l,edge.arrow.size=0, vertex.label=str_pad(names(V(net.sp.single)), width=180, side="left"), vertex.label.color="#000000")
dev.off()
#no padding
plotNetSimple(net.sp.single, layout=l,edge.arrow.size=0, vertex.label=str_pad(names(V(net.sp.single)), width=0, side="left"), vertex.label.color="#000000", xlim=c(-2,0))

names(mem)=names(V(net)[inmodules])
mem
plotNetSimple(net.sp.clust, mark.groups=vl, vertex.label=mem, vertex.label.dist=1)


#my preferred strategy for plotting
#I plot the clusters and singletons together and save as a PDF
#then I move to inkscape and adjust the organization
# pdf("igraphPlots/PC2mgipos.pdf", width = 15)
par(mfrow=c(1,2), mai=rep(0,4))
plotNetSimple(net.sp.clust, mark.groups=vl)
plotNetSimple(net.sp.single, layout=l,edge.arrow.size=0, vertex.label=str_pad(names(V(net.sp.single)), width=0, side="left"), vertex.label.color="#000000", xlim=c(-2,0))
# dev.off()


#mem is membership - fullmem lists all nodes (in clusters + singletons) with cluster membership
#use fullmem to make labels for your clusters (for example cell cycle, apoptosis, etc.)
nomemname=V(net.sp.single)$name
nomem=c(rep(NA, length(nomemname)))
names(nomem)=nomemname
fullmem=c(mem, nomem)
fullmem=data.frame(fullmem)
colnames(fullmem)=c("group membership")
# write.csv(fullmem, file="igraphPlots/PC2mgiposmembership.csv")


##################################










