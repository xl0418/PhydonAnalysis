library(tidyr)
library(dplyr)
library(data.table)
library(readr)
library(ape)
library(castor)
library(foreach)
library(doParallel)


# read data: final_data_rawmadin.csv for bac data; final_data_rawmadin_archaea.csv for archaea data
final_data <- read.csv("GTDBr220/final_data_rawmadin_archaea.csv", header = TRUE, stringsAsFactors = FALSE)

# function expands rep genomes by taking out the clustered species
expand_rep_genomes <- function(df) {
  # set seeds
  set.seed(123)
  # get the species of rep genomes
  rep_genomes <- df$Representative.genome
  # get the species of clustered genomes
  clustered_genomes <- df$Clustered.genomes
  # split by ,
  clustered_genomes <- strsplit(clustered_genomes, ",")[[1]]
  if (length(clustered_genomes) < 5) {
    expand.df <- data.frame(
      species = df$GTDB.species,
      genomes = c(clustered_genomes, rep_genomes),
      superkingdom = df$superkingdom,
      doubling_h = df$doubling_h,
      growth_tmp = df$growth_tmp,
      Representative.genome = rep_genomes,
      phylum = df$phylum,
      stringsAsFactors = FALSE
    )
  } else {
    # sample 5 clustered genomes
    clustered_genomes <- sample(clustered_genomes, 5)
    expand.df <- data.frame(
      species = df$GTDB.species,
      genomes = unique(c(clustered_genomes, rep_genomes)),
      superkingdom = df$superkingdom,
      doubling_h = df$doubling_h,
      growth_tmp = df$growth_tmp,
      Representative.genome = rep_genomes,
      phylum = df$phylum,
      stringsAsFactors = FALSE
    )
  }

  return(expand.df)
}

# parallelize through all rows of final_data
#setup parallel backend to use many processors
cores = detectCores()
cl <- makeCluster(cores[1] - 1) #not to overload your computer
registerDoParallel(cl)

expanded_df <-
  foreach(
    i = 1:nrow(final_data),
    .combine = rbind
  ) %dopar% {
    outdf = expand_rep_genomes(final_data[i,]) #calling a function

    outdf
  }
#stop cluster
stopCluster(cl)

# drop duplicate for expanded_df$genomes
expanded_df <- expanded_df[!duplicated(expanded_df$genomes),]


# save   # GTDB_tax_trait_repGenome_in_tree_expanded.csv for bacteria
write.table(
  expanded_df,
  file = "data_final/GTDB_tax_trait_repGenome_in_tree_expanded_archaea.csv",
  row.names = FALSE,
  quote = FALSE,
  sep = ","
)


