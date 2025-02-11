library(tidyr)
library(dplyr)
library(data.table)
library(readr)
library(ape)
library(castor)
library(foreach)
library(doParallel)

# a function removes the first part of a string sep by '__'
rm_first_part <- function(ss) {
  temp <- strsplit(ss, "__")
  return(temp[[1]][2])
}


get_superkindom <- function(ss) {
  temp1 <- strsplit(ss, ";")[[1]][1]
  temp2 <- strsplit(temp1, "__")[[1]][2]
  return(temp2)
}


##### using raw madin data #####
raw_data_madin <-
  read.csv(
    "GTDBr220/condensed_traits_GTDB.csv",
    sep = ",",
    header = TRUE,
    stringsAsFactors = FALSE
  )
# remove doubling_h == NA
raw_data_madin <- raw_data_madin[!is.na(raw_data_madin$doubling_h), ]
raw_data_madin <- raw_data_madin[, c("species", "doubling_h", "growth_tmp", "phylum")]
# average doubling_h and growth_tmp for the same species
raw_data_madin <- aggregate(raw_data_madin,
                            by = list(raw_data_madin$species),
                            FUN = min)
raw_data_madin <- raw_data_madin[, -2]
colnames(raw_data_madin) <- c("Species", "d", "GrowthTemp", "phylum")
stat_data_madin <- raw_data_madin


######
num_species_in_madin <- length(unique(stat_data_madin$Species))
print(paste(
  "The number of unique species in stat_data_madin is:",
  num_species_in_madin
))

cleared_stat_data_madin <- stat_data_madin[, c("Species", "d", "GrowthTemp", "phylum")]
colnames(cleared_stat_data_madin) <- c("species", "doubling_h", "growth_tmp", "phylum")

num_species_in_madin_no_growth_temp <- length(cleared_stat_data_madin[is.na(cleared_stat_data_madin$growth_tmp), ]$species)
print(
  paste(
    "The number of unique species in stat_data_madin without growth_temp is:",
    num_species_in_madin_no_growth_temp
  )
)

############# read sp_clusters_r220.tsv get from GTDB ############
sp_clusters_r220 <-
  read.csv(
    "GTDBr220/sp_clusters_r220_arc.tsv",
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
  )
sp_clusters_r220$GTDB.species <-
  sapply(sp_clusters_r220$GTDB.species, rm_first_part)
colnames(sp_clusters_r220)[2] <- c("species")

# merge sp_clusters_r220 with cleared_stat_data_madin on species
GTDB_taxonomy_trait_repGenome <-
  merge(sp_clusters_r220,
        cleared_stat_data_madin,
        by = "species",
        all.x = FALSE)
colnames(GTDB_taxonomy_trait_repGenome)
cleared_GTDB_tax_trait_repGenome <-
  GTDB_taxonomy_trait_repGenome[, c(
    "species",
    "Representative.genome",
    "doubling_h",
    "growth_tmp",
    "Clustered.genomes",
    "GTDB.taxonomy",
    "phylum"
  )]
cleared_GTDB_tax_trait_repGenome$superkingdom <-
  sapply(cleared_GTDB_tax_trait_repGenome$GTDB.taxonomy,
         get_superkindom)

# remove the column GTDB.taxonomy
cleared_GTDB_tax_trait_repGenome <-
  cleared_GTDB_tax_trait_repGenome[, -6]


# species in stat_data_madin but not in cleared_GTDB_tax_trait_repGenome
species_not_in_cleared <-
  setdiff(stat_data_madin$Species,
          cleared_GTDB_tax_trait_repGenome$species)

print(
  paste(
    "The number of unique species in stat_data_madin after matching with GTDB is:",
    length(unique(
      cleared_GTDB_tax_trait_repGenome$species
    ))
  )
)

cleared_GTDB_tax_trait_repGenome$GTDB.species <- cleared_GTDB_tax_trait_repGenome$species




#### parallelize the above loop ####

cl <- makeCluster(4)
registerDoParallel(cl)

species_not_in_cleared <- as.character(species_not_in_cleared)

complementary_rows <- foreach(species = species_not_in_cleared, .combine =
                                rbind) %dopar% {
                                  species_p1 <- strsplit(species, " ")[[1]][1]
                                  species_p2 <- strsplit(species, " ")[[1]][2]

                                  complementary_row <- NULL

                                  for (species_gtdb in sp_clusters_r220$species) {
                                    species_gtdb_p1 <- strsplit(species_gtdb, " ")[[1]][1]
                                    species_gtdb_p2 <- strsplit(species_gtdb, " ")[[1]][2]

                                    if (grepl(tolower(species_p1), tolower(species_gtdb_p1), fixed = TRUE) &
                                        grepl(tolower(species_p2), tolower(species_gtdb_p2), fixed = TRUE)) {
                                      stat_madin <- stat_data_madin[stat_data_madin$Species == species, c("Species", "d", "GrowthTemp", "phylum")]
                                      colnames(stat_madin) <- c("species", "doubling_h", "growth_tmp", "phylum")
                                      stat_gtdb <- sp_clusters_r220[sp_clusters_r220$species == species_gtdb, c("Representative.genome",
                                                                                                                "species",
                                                                                                                "Clustered.genomes")]
                                      colnames(stat_gtdb) <- c("Representative.genome",
                                                               "GTDB.species",
                                                               "Clustered.genomes")
                                      superkingdom <- get_superkindom(sp_clusters_r220[sp_clusters_r220$species == species_gtdb, "GTDB.taxonomy"])

                                      complementary_row <- data.frame(
                                        species = species,
                                        Representative.genome = stat_gtdb$Representative.genome,
                                        doubling_h = stat_madin$doubling_h,
                                        growth_tmp = stat_madin$growth_tmp,
                                        Clustered.genomes = stat_gtdb$Clustered.genomes,
                                        superkingdom = superkingdom,
                                        GTDB.species = stat_gtdb$GTDB.species,
                                        phylum = stat_madin$phylum
                                      )

                                    }
                                  }

                                  return(complementary_row)
                                }

stopCluster(cl)



cleared_GTDB_tax_trait_repGenome <- rbind(cleared_GTDB_tax_trait_repGenome, complementary_rows)



### gtdb full tree: bac120_r220.tree for bac; ar53_r220.tree for archaea
gtdb_tree <- read.tree("GTDBr220/ar53_r220.tree")

representative_genomes_stat_madin <-
  unique(cleared_GTDB_tax_trait_repGenome$Representative.genome)

rep_genomes_intree <- gtdb_tree$tip.label






# get the intersection of representative_genomes_stat_madin and rep_genomes_intree
matched_genomes <- intersect(representative_genomes_stat_madin, rep_genomes_intree)

print(paste(
  "The number of matched genomes in gtdb_tree in the end is:",
  length(matched_genomes)
))

sub.gtdb.tree <-
  get_subtree_with_tips(gtdb_tree, only_tips = matched_genomes)
gtdb.subtree.madin <- sub.gtdb.tree$subtree

write.tree(gtdb.subtree.madin, file = "data_final/tree.madin.traits_archaea.tree") # tree.madin.traits.tree for bac




# the genomes in representative_genomes_stat_madin but not in matched_genomes
unmatched_genomes <- setdiff(representative_genomes_stat_madin, matched_genomes)

# the species names
species_unmatched_genomes <- cleared_GTDB_tax_trait_repGenome[cleared_GTDB_tax_trait_repGenome$Representative.genome %in% unmatched_genomes, "species"]


#### take the rows with representative.genomes in matched_genomes ####
final_data <- cleared_GTDB_tax_trait_repGenome[cleared_GTDB_tax_trait_repGenome$Representative.genome %in% matched_genomes, ]

# write the final_data to a csv file
write.csv(final_data,
          "GTDBr220/final_data_rawmadin_archaea.csv",
          row.names = FALSE) # final_data_rawmadin.csv for bac data;   final_data_rawmadin_archaea.csv for archaea data
