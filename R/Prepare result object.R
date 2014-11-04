######################################################
# Functions to interpret the representation analysis #
######################################################


## FUNCTIONS TO CREATE THE DATA OBJECT

Summarize.nodes <- function(dispersion, comm = hcom, tree = htree)
  #Summarizes the information in the representation matrix in a by-node manner
{
  # dispersion   	: a representation matrix with nodes by sites
  
  nodenums <- as.numeric(colnames(dispersion))
  overdisp <- apply(dispersion,2, function(x) mean(abs(x),na.rm = T))
  overdisp[1] <- 0
  overdisp[overdisp == Inf] <- 0
  nodesize <- sapply(nodenums, Node.size)
  nodeoccupancy <- sapply(nodenums, Node.occupancy, comm = comm, tree = tree)
  #node_overrep_overlap <- sapply(nodenums, Node.overrep.overlap, dats = dispersion)
  #nodeoverlap <- sapply(nodenums, Node.overlap)
  nodestats <- data.frame(nodenums, overdisp,  nodeoccupancy, nodesize )#,node_overrep_overlap , nodeoverlap, ) 
  return(nodestats)
}


Summarize.sites <- function(dispersion)
  #Summarizes the information in the representation matrix in a by-sites manner
{
  # dispersion 		: a representation matrix resulting from running measure_dispersion()
  # dat.LL			: the geographical coordinates of the sites
  
  cell <- rownames(dispersion)  
  nodeeff <- apply(dispersion, 1, function(x) mean(abs(x),na.rm = T))
  mapdata <- data.frame(cell, nodeeff)
  return(mapdata)
}

parent_representation = function(node_number, rep_matrix, tree = htree)
  # takes the representation matrix, and summarizes at the parent node (because sister species in the representation matrix are mirror images)
{
  desc = Descendants(node_number)
  if(desc[1] < basal.node(tree) | desc[2] < basal.node(tree))
    return(rep(NA, nrow(rep_matrix)))
  desc1row = which(colnames(rep_matrix) == as.character(desc[1]))
  desc2row = which(colnames(rep_matrix) == as.character(desc[2]))
  return(rowMeans(cbind(rep_matrix[,desc1row], -rep_matrix[, desc2row])))
}

parent_pval_representation = function(node_number, rep_matrix, tree = htree)
  # takes the representation matrix, and summarizes at the parent node (because sister species in the representation matrix are mirror images)
{
  desc = Descendants(node_number)
  if(desc[1] < basal.node(tree) | desc[2] < basal.node(tree))
    return(rep(NA, nrow(rep_matrix)))
  desc1row = which(colnames(rep_matrix) == as.character(desc[1]))
  desc2row = which(colnames(rep_matrix) == as.character(desc[2]))
  # return(rowMeans(cbind(rep_matrix[,desc1row], -rep_matrix[, desc2row])))
  return(rowMeans(cbind(rep_matrix[,desc1row], 1-rep_matrix[, desc2row])))
}

parent_pval_representation_matrix <- function(rep_matrix, tree = htree)
{
  retmat <- sapply(nodenumbers(tree), parent_pval_representation, rep_matrix = rep_matrix, tree = tree)
  colnames(retmat) <- colnames(rep_matrix)
  rownames(retmat) <- rownames(rep_matrix)
  return(retmat)
}

parent_representation_matrix <- function(rep_matrix, tree = htree)
{
  retmat <- sapply(nodenumbers(tree), parent_representation, rep_matrix = rep_matrix, tree = tree)
  colnames(retmat) <- colnames(rep_matrix)
  rownames(retmat) <- rownames(rep_matrix)
  return(retmat)
}

create.datalist <- function(dispersion, sitestatistics, coords)
{
  # a function to create the final result object after the representation analysis has completed
  library(ape)
  datalist <- list()
  dispersion[abs(dispersion) == Inf] <- NA
  datalist$siteresults <- merge(sitestatistics, Summarize.sites(dispersion), by = "cell")
  datalist$siteresults <- merge(datalist$siteresults, coords, by = "cell")
  datalist$siteresults <- datalist$siteresults[match(as.character(coords$cell), datalist$siteresults$cell),]
  datalist$noderesults <- Summarize.nodes(dispersion)
  datalist$rep_matrix <- dispersion[match(as.character(datalist$siteresults$cell), rownames(dispersion)),]
  datalist$parent_rep_matrix <- parent_representation_matrix(datalist$rep_matrix)
  datalist$parent_rep_matrix[abs(datalist$parent_rep_matrix) == Inf] <- NA
  pod <- apply(datalist$parent_rep_matrix,2, function(x) mean(abs(x),na.rm = T))  #I guess this really belongs in the summarize.nodes function
  pod[pod == Inf] <- NA
  datalist$noderesults$parent_overdisp <- pod
  return(datalist)
}


