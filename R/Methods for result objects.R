
# 
# 
# 
# # TODO rename and make an S3 method of function
# plotnode_maps_new_parent <- function(node_number, par_rep_matrix = parent_rep_matrix, tree = htree, comm = hcom, coords = dat.LL, new.window = FALSE, plottype = c("map", "points"))
#   # An exploratory plotting function
#   # Given a node number, this draws four maps in the same window: 
#   #  - a map of the location of overrepresented sites for that node
#   #	- a map of the location of underrepresented sites for that node
#   #	- a map of the species richness of species descending from that node
#   #	- a map of the overall species richness
# {
#   # node_number   :  the internal ape node number (the number showed when calling nodelabels() on a plot of the phylogenetic tree)
#   # dispersion	: the representation matrix (resulting from calling measure_dispersion())
#   
#   require(ape)
#   
#   if(new.window) quartz()
#   par(mfrow = c(2,2))
#   
#   plotoverrep(node_number, par_rep_matrix, coords, tree, new.window = FALSE, plottype = plottype)
#   
#   plotnode(node_number, comm, tree, coords, new.window = FALSE, plottype = plottype, main = paste("richness of node", node_number))
#   
#   plotnode(Descendants(node_number, tree)[1], comm, tree, coords, new.window = FALSE, plottype = plottype, main = paste("richness of descendant 1", Descendants(node_number)[1]))
#   
#   plotnode(Descendants(node_number, tree)[2], comm, tree, coords, new.window = FALSE, plottype = plottype, main = paste("richness of descendant 2", Descendants(node_number)[2]))
#   
# }

print.nodiv_result <- function(x, printlen = 4, ...)
{
  cat(paste("Result of nodiv analysis on", x$type,"data\n"))
  cat(paste("Repeats",x$repeats,"\n"))
  cat(paste("Null model", x$method,"\n\n"))
  cat(paste("Species names: (n = ", length(x$species), "):\n", sep = ""))
cat(paste("\t", paste(x$species[1:printlen], collapse = ", "),", ...\n\n", sep = ""))
cat(paste("Site names (n = ", nrow(x$coords),"):\n", sep = ""))
cat(paste("\t", paste(x$coords$sites[1:printlen], collapse = ", "),", ...\n\n", sep = ""))
  cat(paste("\t", paste(x$coords$sites[1:printlen], collapse = ", "),", ...\n\n", sep = ""))
  cat(paste("GND and SOS calculated for", sum(!is.na(x$GND)), "nodes\n"))
}


summary.nodiv_result <- function(object, ...)
{
  ret <- summary.nodiv_data(object, ...)
  ret$GND <- object$GND
  ret$method <- object$method
  ret$repeats <- object$repeats
  if(sum(ret$GND > 0.6, na.rm = T) > 0) ret$sign <- which(ret$GND > 0.6) else ret$sign <- numeric()
  class(ret) <- "summary_nodiv_result"  
  ret
}

print.summary_nodiv_result <- function(x, printlen = 4, ...)
{
  cat(paste("Result of nodiv analysis on", x$type,"data\n"))
  cat(paste("Repeats:",x$repeats,"\n"))
  cat(paste("Null model:", x$method,"\n\n"))

  cat(paste("Species names: (n = ", length(x$species), "):\n", sep = ""))
  cat(paste("\t", paste(x$species[1:printlen], collapse = ", "),", ...\n\n", sep = ""))
  cat(paste("Site names (n = ", nrow(x$coords),"):\n", sep = ""))
  cat(paste("\t", paste(x$coords$sites[1:printlen], collapse = ", "),", ...\n\n", sep = ""))
  cat(paste("GND and SOS calculated for", sum(!is.na(x$GND)), "nodes\n"))
  cat(paste("GND values:  min", round(min(x$GND, na.rm = T),2), "max", round(max(x$GND, na.rm = T),2), "mean", round(mean(x$GND, na.rm = T),2),"\n\n"))
  cat(paste(length(x$sign),"nodes of", sum(!is.na(x$GND)), "were above 0.6"))  

  if (length(x$sign)>0) 
  {
    cat(":\n")
    for(i in seq_along(x$sign))
    {
      nodelab <- ifelse(!is.null(x$phylo$node.label), x$phylo$node.label[x$sign[i]], "")
      cat(paste("\tNode number ", x$nodes[x$sign[i]],": ", round(x$GND[x$sign[i]],2),"\t", nodelab,"\n", sep = ""))
    }
  } else(cat("\n"))
}

plot.nodiv_result <- function(x, label = nodenumbers(x), main = "", zlims = 0:1, ...)
{
  plot_nodes_phylo(x$GND, label = label, tree = x$phylo, main = main, zlims = zlims, show.legend = TRUE,...)
}

plotSOS <- function(nodiv_result, node, col = brewer.pal(11, "RdYlBu"), zlim, ...)
{

  if(!inherits(nodiv_result, "nodiv_result"))
    stop("nodiv_result must be the result object from running Node_analysis")
  
  node <- node[1]
  if(is.character(node))
  {
    if(is.null(nodiv_result$phylo$node.label))
      stop("node could not be matched, as the phylogeny does not have node labels")
    node <- match(node, nodiv_result$phylo$node.label)
    if(is.na(node))
      stop("the node could not be matched to the node labels")
  }
  
  if(node > Nspecies(nodiv_result))
    node <- node - Nspecies(nodiv_result)
  
  SOS <- nodiv_result$SOS[,node]
  
  if(missing(zlim)) 
  {
    maxabs <- max(abs(SOS))
    zlim <- c(-maxabs, maxabs)
  }
  
  if(is.null(nodiv_result$shape)) shape <- NULL else shape <- nodiv_result$shape
  if(nodiv_result$type == "grid")
    map_var(SOS, nodiv_result$coords, col = brewer.pal(11, "RdYlBu"), zlim = zlim, shape = shape, ...) else
    plot_points(SOS, nodiv_result$coords, col = brewer.pal(11, "RdYlBu"), zlim = zlim, shape = shape, ...)
}
