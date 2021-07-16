#this (very short) script demonstrates how to use RERconverge tree-building functions
#note that these functions were directly built on top of phangorn tree-building functions

#what you need:
#master tree + preferred species to root on
#directory contained filtered + cleaned alignments

#what you need to change in this script:
#filenames
#(if necessary) changing substitution model depending on if you have nucleotides or amino acids



#make trees
mt=readRDS("/home/kowaae22/AnalysisWithThreeTrees/finaltreesv4.rds")
mt=mt$masterTree
mt$edge.length[]=1
mt=root.phylo(mt, "hg19", resolve.root=T)
write.tree(mt, file="/home/kowaae22/100way/mastertree.tree")

library(RERconverge)
#parameters for nucleotide sequences (uses GTR model)
estimatePhangornTreeAll(alndir="/home/kowaae22/100way/promoteralnsfiltered/",treefile="/home/kowaae22/100way/mastertree.tree", output.file = "/home/kowaae22/100way/promotertrees/promotertrees.trees", format = "fasta", type = "DNA", submodel = "GTR")

#parameters for amino acid sequences (uses LG model with k=4)
estimatePhangornTreeAll(alndir="/home/kowaae22/100way/promoteralnsfiltered/",treefile="/home/kowaae22/100way/mastertree.tree", output.file = "/home/kowaae22/100way/promotertrees/promotertrees.trees", format = "fasta")

