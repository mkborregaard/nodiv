
##################################################################
# Functions for use in node-based analysis of clade distribution #
##################################################################

# most functions here require the 'ape' package
require(ape)

## Functions relating trees and sites
	 
# returns the internal node number of the basal node on the phylogeny
basal.node <- function(tree) return(Ntip(tree) + 1)

Descendants <- function(node, tree = htree) 
  return(tree$edge[ tree$edge[,1] == node , 2])

Parent <- function(node, tree = htree)
{
  if (node == Ntip(tree) +1 )   # If the node is the basal node it does not have a parent node
    return (NA)
  return(tree$edge[ tree$edge[,2] == node , 1])
}

Sister <- function(node, tree = htree) 
{
  if (node == Ntip(tree) +1 )   # If the node is the basal node it does not have a sister node
    return (NA)
  sisters = Descendants(Parent(node, tree) , tree)
  return(sisters[! sisters == node])
}	



Site.species <- function(site, comm = hcom)
{
	# returns a list of the species that occur in a site
	sitecom <- comm[comm$Plot == unique(comm$Plot)[site],]  # find all species records of the grid cell (called Plot)
	return( as.character(sitecom$Species))
}

Sitestats <- function(comm = hcom, tree = htree)
#summarizes statistics, such as richness, for each site
{
	# comm  	: the occurrence matrix
	# tree		: the tree
	
	require(ape)
	# some important statistics for each site
	cell <- sort(unique(comm$Plot))
	richness = numeric()
	mean_range_in_cell = numeric()
	
	# the range sizes of all species
	ranges = table(comm$Species)
	ranges = ranges[match(tree$tip.label,names(ranges))]

	for ( plotnum in seq_along(cell))
	{
		print(plotnum)
		
		testcom <- comm[comm$Plot == cell[plotnum],]  # find all species records of the grid cell (called Plot)
		emp <- as.character(testcom$Species)                       # a list of species present at the site

		# calculate mean range and richness
		mean_range_in_cell[plotnum] = mean(ranges[tree$tip.label %in% emp])
		richness[plotnum] = length(emp)
	}
	
	return(data.frame(cell, richness, mean_range_in_cell))
}



Node.size <- function(node)
# calculates the number of species that descend from a certain node
{
	# node : the internal (ape) number of the node 
	return(length(Node.species2(node)))
}

Node.species <- function(node, tree = htree)
# returns a character vector with names of species that descend from a node
{
	# node : the internal (ape) number of the node
	
	require(ape)
	nodetree <- extract.clade(tree, node)
	return(nodetree$tip.label)
}

Node.species2 <- function(node, tree = htree, nodespecmatrix = node_species)
# returns a character vector with names of species that descend from a node
{
	# node : the internal (ape) number of the node

	return(colnames(nodespecmatrix)[nodespecmatrix[node-Ntip(tree),] > 0])
}

Node.comm <- function(node, comm = hcom, tree = htree)
# returns a samplelist of sites occupied by at least one member of the node
{
	# node : the internal (ape) number of the node
	if (node < Ntip(tree)) #if it is in fact a tip
		nodespecs = tree$tip.label[node] else nodespecs <- Node.species2(node, tree) #how to pass the nodespecmatrix?
	
	nodecom <- subset(comm, Species %in% nodespecs)
	return(nodecom)
}

Node.sites <- function(node, comm = hcom, tree = htree)
# calculates which sites are occupied by at least one member of the node
{
	# node : the internal (ape) number of the node
	
	nodecom <- Node.comm(node, comm, tree)
	return(unique(nodecom$Plot))
}
	

Node.occupancy <- function(node, comm = hcom, tree = htree)
# calculates the number of sites occupied by at least one member of the node
{
	# node : the internal (ape) number of the node 
	return(length(Node.sites(node, comm, tree)))
}

Node.richness <- function(node_number, comm = hcom, tree = htree, coords = dat.LL)
# A summary function that calculates the species richness pattern of all species that descend from a certain node
{
	# node_number  :  the internal ape node number (the number showed when calling nodelabels() on a plot of the phylogenetic tree)

	require(ape)
	
	nodecom = Node.comm(node_number, comm, tree)
	richs <- data.frame(table(as.character(nodecom$Plot)), stringsAsFactors = FALSE)
	names(richs) <- c("cell", "richness")
	richs <- merge(coords, richs, all.x = T, by="cell")
	richs <- richs[match(as.character(coords$cell), as.character(richs$cell)),]

	return(richs)
}

# returns a vector with the internal numbers of all nodes on the tree
nodenumbers <- function(tree) {require(ape); return(1:Nnode(tree) + Ntip(tree))}


Create_node_by_species_matrix = function(tree = htree)
{
  # create a matrix with 0s and 1s indicating which species descend from each node
	require(ape)
	
	nodespecies <- matrix(0, nrow = Nnode(tree), ncol = Ntip(tree))
	colnames(nodespecies) <- htree$tip.label
	rownames(nodespecies) <- 1:Nnode(tree) + Ntip(tree)
	
	for ( i in 1:Nnode(tree))
	{
		nodespecies[i,which(colnames(nodespecies) %in% Node.species(nodenumbers(tree)[i]))] <- 1
	}
	
	return(nodespecies)
}

Node_representation = function(species_index, node_species_matrix = node_species)
{
	# returns a vector of how many times each node is represented in a community (species list)
	submatrix = node_species_matrix[,as.logical(species_index)]
	return(rowSums(submatrix))
}

	
######################################################
# Functions to interpret the representation analysis #
######################################################


Summarize.nodes <- function(dispersion, comm = hcom, tree = htree)
#Summarizes the information in the representation matrix in a by-node manner
{
	# dispersion 		: a representation matrix with nodes by sites
	
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

print("Succesfully read representation functions")
