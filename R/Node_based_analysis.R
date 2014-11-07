
##################################################################
# Functions for use in node-based analysis of clade distribution #
##################################################################


##############INTERNAL FUNCTIONS (NOT TO BE EXPORTED)###############


## Functions relating trees and sites
	 
# returns the internal node number of the basal node on the phylogeny
basal_node <- function(tree) 
{
  if(inherits(tree, "nodiv_data"))
    tree <- tree$phylo
  if(!inherits(tree, "phylo"))
    stop("tree must be an object of type phylo or nodiv_data")
  return(Ntip(tree) + 1)
}

Descendants <- function(node, tree) 
{
  if(inherits(tree, "nodiv_data"))
    tree <- tree$phylo
  if(!inherits(tree, "phylo"))
    stop("tree must be an object of type phylo or nodiv_data")
  return(tree$edge[ tree$edge[,1] == node , 2])
}
  

Parent <- function(node, tree)
{
  if(inherits(tree, "nodiv_data"))
    tree <- tree$phylo
  if(!inherits(tree, "phylo"))
    stop("tree must be an object of type phylo or nodiv_data")  
  if (node == Ntip(tree) +1 )   # If the node is the basal node it does not have a parent node
    return (NA)
  return(tree$edge[ tree$edge[,2] == node , 1])
}

Sister <- function(node, tree) 
{
  if(inherits(tree, "nodiv_data"))
    tree <- tree$phylo
  if(!inherits(tree, "phylo"))
    stop("tree must be an object of type phylo or nodiv_data")
  if (node == Ntip(tree) +1 )   # If the node is the basal node it does not have a sister node
    return (NA)
  sisters = Descendants(Parent(node, tree), tree)
  return(sisters[! sisters == node])
}	

########

Site_species <- function(site, nodiv_data)
{
	# returns a list of the species that occur in a site
	spec_index <- nodiv_data$comm[site,] == 1
	return(nodiv_data$species[spec_index])
}

## TODO change this - richness is already calculated, but perhaps mean_range_in_cell is relevant?
Sitestats <- function(comm, tree )
#summarizes statistics, such as richness, for each site
{
	# comm  	: the occurrence matrix
	# tree		: the tree
	
	require(ape)
	# some important statistics for each site
	cell <- sort(unique(comm$plot))
	richness = numeric()
	mean_range_in_cell = numeric()
	
	# the range sizes of all species
	ranges = table(comm$id)
	ranges = ranges[match(tree$tip.label,names(ranges))]

	for ( plotnum in seq_along(cell))
	{
		print(plotnum)
		
		testcom <- comm[comm$plot == cell[plotnum],]  # find all species records of the grid cell (called plot)
		emp <- as.character(testcom$id)                       # a list of species present at the site

		# calculate mean range and richness
		mean_range_in_cell[plotnum] = mean(ranges[tree$tip.label %in% emp])
		richness[plotnum] = length(emp)
	}
	
	return(data.frame(cell, richness, mean_range_in_cell))
}

Node_spec <- function(node, tree)
  # returns a character vector with names of species that descend from a node
{
  # node : the internal (ape) number of the node
  
  require(ape)
  nodetree <- extract.clade(tree, node)
  return(nodetree$tip.label)
}



# returns a vector with the internal numbers of all nodes on the tree
nodenumbers <- function(x) 
{
  if(inherits(x, "nodiv_data"))
    x <- x$phylo
  if(!inherits(x, "phylo"))
    stop("x must be of type phylo or nodiv_data")
  return(1:Nnode(x) + Ntip(x))
}


Create_node_by_species_matrix = function(tree)
{
  # create a matrix with 0s and 1s indicating which species descend from each node
  nodespecies <- matrix(0, nrow = Nnode(tree), ncol = Ntip(tree))
  colnames(nodespecies) <- tree$tip.label
  rownames(nodespecies) <- 1:Nnode(tree) + Ntip(tree)
  
  for ( i in 1:Nnode(tree))
  {
    nodespecies[i,which(colnames(nodespecies) %in% Node_spec(nodenumbers(tree)[i], tree))] <- 1
  }
  
  return(nodespecies)
}




#######EXPORTED FUNCTIONS################


## TODO these should all be available from the data object, S3
Node_size <- function(node, nodiv_data)
# calculates the number of species that descend from a certain node
{
  if(node == 0)
    return(sapply(nodenumbers(nodiv_data), function(nod) Node_size(nodiv_data, nod))) else
      if(! node %in% nodenumbers(nodiv_data)) stop("node must be one of the internal nodes in nodiv_data$phylo")
	# node : the internal (ape) number of the node 
	return(sum(nodiv_data$node_species[node - Nspecies(nodiv_data),]))
}


Node_species <- function(node, nodiv_data)
  # returns a character vector with names of species that descend from a node
{
  # node : the internal (ape) number of the node
  if(!inherits(nodiv_data, "nodiv_data"))
    stop("nodiv_data must be an object of type nodiv_data")
  
  return(colnames(nodiv_data$node_species)[nodiv_data$node_species[node-Nspecies(nodiv_data),] > 0])
}


Node_comm <- function(node, nodiv_data)
# returns a samplelist of sites occupied by at least one member of the node
{
	# node : the internal (ape) number of the node
	if (node < Ntip(nodiv_data$tree)) #if it is in fact a tip
		nodespecs = nodiv_data$species[node] else nodespecs <- Node_species(node, nodiv_data)
	
	nodecom <- subset(nodiv_data$hcom, id %in% nodespecs)
	return(nodecom)
}

Node_sites <- function(node, nodiv_data)
# calculates which sites are occupied by at least one member of the node
{
	# node : the internal (ape) number of the node
	
	nodecom <- Node_comm(node, nodiv_data)
	return(unique(nodecom$plot))
}
	

Node_occupancy <- function(nodiv_data, node = 0)
# calculates the number of sites occupied by at least one member of the node
{
  if(node == 0)
    return(sapply(nodenumbers(nodiv_data), function(nod) Node_occupancy(nodiv_data, nod))) else
    if(! node %in% nodenumbers(nodiv_data)) stop("node must be one of the internal nodes in nodiv_data$phylo")
	
  return(length(Node_sites(node, nodiv_data)))
}

# Node_richness <- function(node_number, comm = hcom, tree = htree, coords = dat.LL)
# # A summary function that calculates the species richness pattern of all species that descend from a certain node
# {
# 	# node_number  :  the internal ape node number (the number showed when calling nodelabels() on a plot of the phylogenetic tree)
# 
# 	require(ape)
# 	
# 	nodecom = Node_comm(node_number, comm, tree)
# 	richs <- data.frame(table(as.character(nodecom$plot)), stringsAsFactors = FALSE)
# 	names(richs) <- c("cell", "richness")
# 	richs <- merge(coords, richs, all.x = T, by="cell")
# 	richs <- richs[match(as.character(coords$cell), as.character(richs$cell)),]
# 
# 	return(richs)
# }

## TODO almost all of the above functions can be changed to simple access functions 
## of a nodiv_data object that has been subsetted


	

