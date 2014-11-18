

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

plot.nodiv_result <- function(x, label = nodenumbers(x), zlim = 0:1, ...)
{
  plot_nodes_phylo(x$GND, tree = x$phylo, label = label, main = "",  zlim = zlim, show.legend = TRUE,...)
}

plotSOS <- function(nodiv_result, node, col = cm.colors(64), zlim, ...)
{
  sos <- SOS(nodiv_result, node)
  if(missing(zlim)) 
  {
    maxabs <- max(abs(sos), na.rm = T)
    zlim <- c(-maxabs, maxabs)
  }
  plot_sitestat(nodiv_result, sos, col = col, zlim = zlim, ...)
}

SOS <- function(nodiv_result, node)
{ 
  if(!inherits(nodiv_result, "nodiv_result"))
    stop("nodiv_result must be the result object from running Node_analysis")
  
  node <- identify_node(node, nodiv_result)
  
  if(node > Nspecies(nodiv_result))
    node <- node - Nspecies(nodiv_result)
  
  SOS <- nodiv_result$SOS[,node]
  SOS
}

