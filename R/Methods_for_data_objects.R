

Nspecies <- function(x)
{
  if (!inherits(x, "distrib_data")) 
    stop("object is not of class \"distrib_data\"")
  length(x$species)
}

Nsites<- function(x)
{
  if (!inherits(x, "distrib_data")) 
    stop("object is not of class \"distrib_data\"")
  length(x$coords)
}

print.distrib_data <- function(x, printlen = 4, ...)
{
  cat(paste("Data object with", x$type," distributions of", Nspecies(x),"species in", Nsites(x),"sites\n\n"))
  cat("Species names:\n")
  cat(paste("\t", paste(x$species[1:printlen], collapse = ", "),", ...\n", sep = ""))
  cat("Site names:\n")
  cat(paste("\t", paste(x$coords$sites[1:printlen], collapse = ", "),", ...\n", sep = ""))
}

print.nodiv_data <- function(x, printlen = 4, ...)
{
  cat(paste("Data object with", x$type,"distributions and phylogenetic relationships of", Nspecies(x),"species in", Nsites(x),"sites\n\n"))
  cat("Species names:\n")
  cat(paste("\t", paste(x$species[1:printlen], collapse = ", "),", ...\n", sep = ""))
  cat("Site names:\n")
  cat(paste("\t", paste(x$coords$sites[1:printlen], collapse = ", "),", ...\n", sep = ""))
}

summary.distrib_data <- function(object, ...)
{
  richness <- if(object$type == "grid")
    SpatialPixelsDataFrame(SpatialPoints(object$coords), data.frame(richness = rowSums(object$comm))) else
    SpatialPointsDataFrame(SpatialPoints(object$coords), data.frame(richness = rowSums(object$comm))) 
      
  occupancy <- colSums(object$comm)
  ret <- list(species = object$species, coords = object$coords, richness = richness, occupancy = occupancy, type = object$type)
  class(ret) <- "summary_distrib_data"
  ret
}

print.summary_distrib_data <- function(x, printlen = 4, ...)
{
  cat(paste("Data object with", x$type,"distributions of", length(x$species),"species in", nrow(x$coords),"sites\n\n"))
  cat("Species names:\n")
  cat(paste("\t", paste(x$species[1:printlen], collapse = ", "),", ...\n", sep = ""))
  cat("Site names:\n")
  cat(paste("\t", paste(x$coords$sites[1:printlen], collapse = ", "),", ...\n\n", sep = ""))
  cat(paste("Species richness:  min", min(x$richness$richness), "max", max(x$richness$richness), "mean",mean(x$richness$richness),"\n"))
  cat(paste("Species occupancy:  min", min(x$occupancy, "max", max(x$occupancy), "mean",mean(x$occupancy,"\n"))))
}

summary.nodiv_data <- function(object, ...)
{
  ret <- summary.distrib_data(object, ...)
  ret$nodes <- nodenumbers(object$phylo)
  ret$node.label = object$phylo.node.label
  class(ret) <- "summary_nodiv_data"
  ret
}

print.summary_nodiv_data <- function(x, printlen = 4, ...)
{
  cat(paste("Data object with", x$type,"distributions and phylogenetic relationships of", length(x$species) ,"species in", length(x$sites),"sites\n\n"))
  cat("Species names:\n")
  cat(paste("\t", paste(x$species[1:printlen], collapse = ", "),", ...\n", sep = ""))
  cat("Site names:\n")
  cat(paste("\t", paste(x$coords$sites[1:printlen], collapse = ", "),", ...\n\n", sep = ""))
  cat(paste("Species richness:   min:", min(x$richness$richness), "\tmax:", max(x$richness$richness), "\tmean:", round(mean(x$richness$richness),2),"\n"))
  cat(paste("Species occupancy:  min:", min(x$occupancy), "\tmax:", max(x$occupancy), "\tmean:",round(mean(x$occupancy),2),"\n\n"))
  cat(paste("The phylogeny has", length(x$nodes), "internal nodes"))
  if (!is.null(x$node.label)) 
  {
    cat("Node labels:\n")
    if (x$nodes > printlen) 
    {
      cat(paste("\t", paste(x$node.label[1:printlen], collapse = ", "),", ...\n", sep = ""))
    } else print(x$node.label)
  }
}

add_shape <- function(distrib_data, shape)
{
  if(!inherits(distrib_data, "distrib_data"))
    stop("argument must be an object of types distrib_data, nodiv_data or nodiv_results")
  distrib_data$shape <- shape
  distrib_data
}

plot.distrib_data <- function(x, ...)
{
  if(is.null(x$shape)) shape <- NULL else shape <- x$shape
  if(x$type == "grid")
    plot_grid(richness(x), x$coords, shape = shape, ...) else
    plot_points(richness(x), x$coords, shape = shape, ...)
}  

plot_richness <- function(distrib_data, ...)
{
  if(!inherits(distrib_data, "distrib_data"))
    stop("argument must be an object of type distrib_data, nodiv_data or nodiv_results")
  plot.distrib_data(distrib_data, ...)
}

plot_node <- function(nodiv_data, node = basal_node(nodiv_data), ...)
{
  if(!inherits(nodiv_data, node))
    stop("argument must be an object of type nodiv_data or nodiv_result")
  node <- identify_node(node, nodiv_data)
  plot_richness(subsample(nodiv_data, node = node), ...)
}

plot.nodiv_data <- function(x,  col = rev(heat.colors(64)), ...)
{
  par(mfrow = c(1,2))
  plot.distrib_data(x, col = col, ...)
  plot(x$phylo, show.tip.label = isTRUE(Nspecies(x) < 40), cex = 0.7) #need to specify explicitly which
}

subsample<- function(x, ...) UseMethod("subsample")

subsample.distrib_data <- function(x, sites = NULL, species = NULL, ...)
{
  if(inherits(sites, "SpatialPoints")) sites <- as.character(sites@data)
  keep_sites <- F
  if(is.character(sites)) 
  {
    if(length(sites) == 1)
    {
      if(sites == "all") keep_sites <- T 
    } else 
      sites <- match(sites, x$coords$sites)
  }  
  if(is.null(sites) | keep_sites) sites <- 1:Nsites(x)
  
  keep_species <- F
  if(is.character(species))
  {
    if(length(species) == 1)
    {
      if(species == "all") keep_species <- T 
    } else
    species <- match(species, x$species)
  }
  
  if(is.null(species) | keep_species) species <- 1:Nspecies(x)
  
  ret <- x
  ret$comm <- ret$comm[sites, species]
  
  if(keep_sites) sites_keep <- rep_len(TRUE, nrow(ret$comm)) else sites_keep <- (rowSums(ret$comm, na.rm = T) > 0)
  
  if(keep_species) sites_keep <- rep_len(TRUE, nrow(ret$comm)) else species_keep <- (colSums(ret$comm, na.rm = T) > 0)

  ret$comm <- ret$comm[sites_keep, species_keep]
  
  ret$species <- colnames(ret$comm)
  ret$coords <- ret$coords[ret$coords$sites %in% rownames(ret$comm),]
  
  if(!is.null(x$shape)) ret$shape <- x$shape
  
  return(ret)
}

subsample.nodiv_data <- function(x, sites = NULL, species = NULL, node = NULL, ...)
{
#   if(sum(!is.null(species), !is.null(node), !is.null(sites)) > 1) stop("you can only specify one of sites, species or node")
#   if(sum(!is.null(species), !is.null(node), !is.null(sites)) == 0) stop("you must specify one of sites, species or node")
  
  ret_phylo <- x$phylo
  ret_phylo$node.label <- nodenumbers(x)  #this line
  
  if(!is.null(node))
  {
    node <- node[1]
    if(is.character(node))
    {
      if(node %in% x$phylo$node.label)
        node <- which(x$phylo$node.label == node) else stop("no nodelabels in phylo corresponds to that node")
    } 

    ret_phylo <- extract.clade(ret_phylo, node)
    species <- ret_phylo$tip.label
  } 
  
  ret <- subsample.distrib_data(x, sites, species)
  dat <- match.phylo.comm(ret_phylo, ret$comm)
  new_phylo <- dat$phy
  ret$phylo <- drop.tip(x$phylo, which(! x$species %in% new_phylo$tip.label))  #this line
  ret$species <- ret$phylo$tip.label
  
  ret$hcom <- subset(x$hcom, x$hcom$plot %in% ret$coords$sites & x$hcom$id %in% ret$species)
  ret$node_species <- x$node_species[, colnames(x$node_species) %in% ret$species]
  ret$node_species <- ret$node_species[rowSums(ret$node_species) > 0,]
  ret$node_species <- ret$node_species[as.numeric(rownames(ret$node_species)) %in% as.numeric(new_phylo$node.label),] # and this line are an ugly hack to make sure the node_species matrix does not get perverted
  rownames(ret$node_species) <- nodenumbers(ret$phylo)
  return(ret)
}

richness <- function(x)
{  
  if(!inherits(x, "distrib_data"))
    stop("x must be an object of type distrib_data")
  return(summary(x)$richness$richness)
}



