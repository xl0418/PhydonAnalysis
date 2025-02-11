library(ape)
library(phytools)
library(ggtree)
library(ggplot2)
library(DDD)
library(ggsci)
library(viridis)
library(wesanderson)
library(tidyverse)
library(castor)
library(ggtree)
get_short_accession <- function(accid){
  temp_acc <- strsplit(accid, "_")
  right_accid <- paste0(temp_acc[[1]][2], '_',temp_acc[[1]][3])
  return(right_accid)
}

## read gtdb subtree
fulltree <- read.tree("data_final/tree.madin.traits.tree")
print(paste("The number of species in tree.madin.traits.tree is", length(unique(fulltree$tip.label))))


## read trait data
species_data <- read_csv("data_final/GTDB_tax_trait_repGenome_in_tree_expanded.csv")
genomes_representative <- species_data[,c("genomes", "Representative.genome")]

species_data <- species_data[,c("species", "Representative.genome", "superkingdom", "doubling_h")]


colnames(species_data) <- c("species", "full_accession", "kingdoms", "d")
print(paste("The number of species in trait data is", length(unique(species_data$species))))


## ultrametic tree
fulltree_ddd <- phylo2L(fulltree)
fulltree_ddd[,4] <- -1
ss_tree <- L2phylo(fulltree_ddd)

##
fulltree_label <- fulltree$tip.label

ss_tree$tip.label <- fulltree_label

# ggtree(ss_tree) +theme_tree2()
# ggsave("figures/newtree.png")

# time slicing tree
obj<-ltt(ss_tree,plot=FALSE)

no.clades <- seq(10, 410, 10) # 410 for bac

no_clades <- c()
which_clade <- c()
no_tips <- c()
training_tips <- c()

####### cut the tree into n clades; This is used to determine the training data sets and the testing data sets ##########
for(k in no.clades){
  ## find the height when there are k lineages

  h<-mean(obj$times[c(which(obj$ltt==k),which(obj$ltt==(k+1)))])

  subtrees<-treeSlice(ss_tree,h,trivial=TRUE)
  tt_trees <- length(subtrees)
  for(kk in 1:tt_trees){
    nottraining_acc <- subtrees[[kk]]$tip.label

    training_acc_rep <- fulltree_label[-which(fulltree_label %in% nottraining_acc)]
    training_acc <- c()
    for(rep_genome in training_acc_rep){
      training_acc <- c(training_acc, genomes_representative[which(genomes_representative$Representative.genome == rep_genome),1]$genomes)
    }


    filename <- paste0('data_final/GR3subtrees/Accession_', k, '_clades_train_',kk,'.csv')
    write.table(training_acc, filename, row.names = FALSE, quote=FALSE,col.names=FALSE, sep = ',')

    no_clades <- c(no_clades, k)
    which_clade <- c(which_clade, kk)
    no_tips <- c(no_tips, length(nottraining_acc))
    training_tips <- c(training_tips, length(training_acc_rep))

  }

}

tree_cutting_info <- data.frame(no_clades, which_clade, no_tips, training_tips)
colnames(tree_cutting_info) <- c("no_clades", "which_clade", "nottraintips", "traintips")


hk <- c()
crown_t <- max(obj$times)
cuttime_clades <- seq(10, 410, 10)

for(k in cuttime_clades){
  ## find the height when there are k lineages

  h<-mean(obj$times[c(which(obj$ltt==k),which(obj$ltt==(k+1)))])
  hk <- c(hk, crown_t - h)
}
hk_df <- data.frame(hk, k = cuttime_clades)


hk_filename <- paste0('data_final/cuttime.csv')
write.table(hk_df, hk_filename, row.names = FALSE, quote=FALSE,col.names=TRUE, sep = ',')



