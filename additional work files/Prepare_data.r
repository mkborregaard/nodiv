#TODO replace with import
library(picante)

Definition of nodiv_data-class











Definition of nodiv-class



commat <- read.delim("../OCCUR_ALL.xls")
hcom <- matrix2sample(commat)

htree <- read.tree("../phylogeny.txt")

commat = commat[, match(htree$tip.label, colnames(commat))]

source("Node_based_analysis.R")

node_species = Create_node_by_species_matrix(htree)

dat.LL <- NA

save(hcom, dat.LL, htree, commat, node_species, file = "GlobMammals_data.RData")

