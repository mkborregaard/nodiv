
##################################################################
# Functions for use in node-based analysis of clade distribution #
##################################################################


##############INTERNAL FUNCTIONS (NOT TO BE EXPORTED)###############


## Functions relating trees and sites
	 
# returns the internal node number of the basal node on the phylogeny
basal.node <- function(tree) return(Ntip(tree) + 1)

Descendants <- function(node, tree) 
  return(tree$edge[ tree$edge[,1] == node , 2])

Parent <- function(node, tree)
{
  if (node == Ntip(tree) +1 )   # If the node is the basal node it does not have a parent node
    return (NA)
  return(tree$edge[ tree$edge[,2] == node , 1])
}

Sister <- function(node, tree) 
{
  if (node == Ntip(tree) +1 )   # If the node is the basal node it does not have a sister node
    return (NA)
  sisters = Descendants(Parent(node, tree), tree)
  return(sisters[! sisters == node])
}	



Site_species <- function(site, comm)
{
	# returns a list of the species that occur in a site
	sitecom <- comm[comm$Plot == unique(comm$Plot)[site],]  # find all species records of the grid cell (called Plot)
	return( as.character(sitecom$Species))
}


Sitestats <- function(comm, tree )
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

Node_species <- function(node, tree = htree)
  # returns a character vector with names of species that descend from a node
{
  # node : the internal (ape) number of the node
  
  require(ape)
  nodetree <- extract.clade(tree, node)
  return(nodetree$tip.label)
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
    nodespecies[i,which(colnames(nodespecies) %in% Node_species(nodenumbers(tree)[i]))] <- 1
  }
  
  return(nodespecies)
}




#######EXPORTED FUNCTIONS################


## TODO these should all be available from the data object, S3
Node_size <- function(node)
# calculates the number of species that descend from a certain node
{
	# node : the internal (ape) number of the node 
	return(length(Node_species2(node)))
}

#TODO rename here and everywhere and make an S3 function
Node_species2 <- function(node, tree = htree, nodespecmatrix = node_species)
  # returns a character vector with names of species that descend from a node
{
  # node : the internal (ape) number of the node
  
  return(colnames(nodespecmatrix)[nodespecmatrix[node-Ntip(tree),] > 0])
}


Node_comm <- function(node, comm = hcom, tree = htree)
# returns a samplelist of sites occupied by at least one member of the node
{
	# node : the internal (ape) number of the node
	if (node < Ntip(tree)) #if it is in fact a tip
		nodespecs = tree$tip.label[node] else nodespecs <- Node_species2(node, tree)
	
	nodecom <- subset(comm, Species %in% nodespecs)
	return(nodecom)
}

Node_sites <- function(node, comm = hcom, tree = htree)
# calculates which sites are occupied by at least one member of the node
{
	# node : the internal (ape) number of the node
	
	nodecom <- Node_comm(node, comm, tree)
	return(unique(nodecom$Plot))
}
	

Node_occupancy <- function(node, comm = hcom, tree = htree)
# calculates the number of sites occupied by at least one member of the node
{
	# node : the internal (ape) number of the node 
	return(length(Node_sites(node, comm, tree)))
}

Node_richness <- function(node_number, comm = hcom, tree = htree, coords = dat.LL)
# A summary function that calculates the species richness pattern of all species that descend from a certain node
{
	# node_number  :  the internal ape node number (the number showed when calling nodelabels() on a plot of the phylogenetic tree)

	require(ape)
	
	nodecom = Node_comm(node_number, comm, tree)
	richs <- data.frame(table(as.character(nodecom$Plot)), stringsAsFactors = FALSE)
	names(richs) <- c("cell", "richness")
	richs <- merge(coords, richs, all.x = T, by="cell")
	richs <- richs[match(as.character(coords$cell), as.character(richs$cell)),]

	return(richs)
}

## TODO should just be deleted
# Node_representation = function(species_index, node_species_matrix = node_species)
# {
# 	# returns a vector of how many times each node is represented in a community (species list)
# 	submatrix = node_species_matrix[,as.logical(species_index)]
# 	return(rowSums(submatrix))
# }

	

