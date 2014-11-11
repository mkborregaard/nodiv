
##################################################################
# Functions for use in node-based analysis of clade distribution #
##################################################################


##############INTERNAL FUNCTIONS (NOT TO BE EXPORTED)###############


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

Node_spec <- function(node, tree)
  # returns a character vector with names of species that descend from a node
{
  # node : the internal (ape) number of the node
  if(!inherits(tree, "phylo"))
    stop("tree must be an object of type phylo")
  
  if(!node %in% nodenumbers(tree))
    stop("not a valid node number")
  
  nodetree <- extract.clade(tree, node)
  return(nodetree$tip.label)
}

Node_comm <- function(nodiv_data, node)
# returns a samplelist of sites occupied by at least one member of the node
{
  # node : the internal (ape) number of the node
  if (node < Ntip(nodiv_data$phylo)) #if it is in fact a tip
    nodespecs = nodiv_data$species[node] else nodespecs <- Node_species(nodiv_data, node)
  
  nodecom <- subset(nodiv_data$hcom, nodiv_data$hcom$id %in% nodespecs)
  return(nodecom)
}

identify_node <- function(node, tree)
{
  if(inherits(tree, "nodiv_data"))
    tree <- tree$phylo
  if(!inherits(tree, "phylo"))
    stop("tree must be an object of type phylo or nodiv_data")

  if(!is.vector(node)) stop("node must be either numeric or character")
  if(length(node)>1) warning("node was had length > 1 - only the first element was used")
  
  node <- node[1]
  
  if(is.character(node))
  {
    if(is.null(tree$node.label))
      stop("node could not be matched, as the phylogeny does not have node labels")
    node <- match(node, tree$node.label)
    if(is.na(node))
      stop("the node could not be matched to the node labels")
  }
  
  if(!node %in% nodenumbers(tree))
    stop("Undefined node")

  node
}

############ EXPORTED FUNCTIONS #############################

## Functions relating trees and nodes
	 
# returns the internal node number of the basal node on the phylogeny
basal_node <- function(tree) 
{
  if(inherits(tree, "nodiv_data"))
    tree <- tree$phylo
  if(!inherits(tree, "phylo"))
    stop("tree must be an object of type phylo or nodiv_data")
  return(Ntip(tree) + 1)
}

# returns a vector with the internal numbers of all nodes on the tree
nodenumbers <- function(tree) 
{
  if(inherits(tree, "nodiv_data"))
    tree <- tree$phylo
  if(!inherits(tree, "phylo"))
    stop("tree must be of type phylo or nodiv_data")
  return(1:Nnode(tree) + Ntip(tree))
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



####### Functions summarizing nodiv_data on nodes

Node_size <- function(nodiv_data, node = NULL)
# calculates the number of species that descend from a certain node
{
  if(!inherits(nodiv_data, "nodiv_data"))
    stop("nodiv_data must be an object of type nodiv_data or nodiv_result")

  if(is.null(node))
    return(sapply(nodenumbers(nodiv_data), function(nod) Node_size(nodiv_data, nod))) else
     
  node <- identify_node(node, nodiv_data)
	return(sum(nodiv_data$node_species[node - Nspecies(nodiv_data),]))
}


Node_species <- function(nodiv_data, node)
  # returns a character vector with names of species that descend from a node
{
  if(!inherits(nodiv_data, "nodiv_data"))
    stop("nodiv_data must be an object of type nodiv_data")
  
  node <- identify_node(node, nodiv_data)
  return(colnames(nodiv_data$node_species)[nodiv_data$node_species[node-Nspecies(nodiv_data),] > 0])
}


Node_sites <- function(nodiv_data, node)
# calculates which sites are occupied by at least one member of the node
{
	if(!inherits(nodiv_data, "nodiv_data"))
    stop("nodiv_data must be an object of type nodiv_data or nodiv_result")
  node <- identify_node(node, nodiv_data)
	nodecom <- Node_comm(nodiv_data, node)
	return(unique(nodecom$plot))
}
	

Node_occupancy <- function(nodiv_data, node = NULL)
# calculates the number of sites occupied by at least one member of the node
{
  if(!inherits(nodiv_data, "nodiv_data"))
    stop("nodiv_data must be an object of type nodiv_data or nodiv_result")
  if(is.null(node))
    return(sapply(nodenumbers(nodiv_data), function(nod) Node_occupancy(nodiv_data, nod))) 

  node <- identify_node(node, nodiv_data)
  return(length(Node_sites(nodiv_data, node)))
}


	

