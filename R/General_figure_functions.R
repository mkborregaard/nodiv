



##########################################################
# Here comes a list of functions to use for the analysis. These definitions should all be loaded into R

create.cols <- function(vec, col = rev(rainbow(256,start=0.02,end=0.59)), zlims = c(min(vec,na.rm = T),max(vec,na.rm = T)))
{
  vec = vec - zlims[1]
  vec = floor(vec * (length(col)-1)/(zlims[2]- zlims[1]))+1
  return(col[vec])
}


#TODO should this be redone as a raster function?
# spatial grids data frame.

map.var = function(variable, Long, Lat, new.window = FALSE, xrange = range(Long), yrange = range(Lat), bordercolor = "darkgrey", cut.to.coast = FALSE, ...)
# a function that will create a pretty color-coded map of any variable with Long/Lat information
{
	# variable : the variable controlling the color
	# Long : longitude of each point
	# Lat : latitude of each point
	# new.window : should the plot be created in a new window? Set to FALSE, e.g. for saving to a postscript file
	# ... : extra arguments to be passed to the "image.plot" function

	require(fields)
	require(maptools)
	variable[is.na(variable)] <- -9999  #I am doing this trick to ensure that the whole range of Lat/Long values are printed, not just the ones with valid values of 'variable'
	cellsize_x <- min(diff(sort(unique(Long))),na.rm = T)
	cellsize_y <- min(diff(sort(unique(Lat))),na.rm = T)
	map <- as.image(variable, x = cbind(Long, Lat), grid = list(x = seq(xrange[1],xrange[2], by = cellsize_x), y = seq(yrange[1],yrange[2],by = cellsize_y)))
	map$z[map$z == -9999] <- NA

	if(new.window) quartz()
	image.plot(map, ...)
  sa <- readShapePoly('~/shapecoast/World_coast.shp')
	plot(sa,add=T, border = bordercolor)
	if(cut.to.coast)
	{
		sc <- readShapePoly('~/shapecoast/World_coast.shp')
		plot(sc, add = T, border = NULL, col = "white", usePolypath = TRUE)
		box()
	}
}



# TODO edit to give the output format used in the figure

plot_tree_nodes <- function(variable, label = variable, tree = htree, main = deparse(substitute(variable)),  new.window = FALSE, col = rev(rainbow(256,start=0.02,end=0.59)), show.legend = TRUE, sig.cutoff, nodes, roundoff= TRUE, genera.lines = FALSE, ...)
# plots a tree, where the colors of the nodes reflects the values of variable
{
	# variable  : the variable that controls the colors of nodes - given for ALL nodes, even though some nodes are not plotted!
	# tree 		: the phylogenetic tree
	# main 		: the title to give the plot
	# new.window: should the plot be made in a new window?
	# col		: the color palette to use for making the color plot

	require(ape)
	
	if(!length(variable) == Nnode(tree))
		stop("The length of the variable vector must be the same length as the number of nodes on the tree")
	
	plotvar <- variable
	if(roundoff & is.numeric(label)) label = round(label,2)
	
	if(missing(nodes)) node_index = rep(TRUE, Nnode(tree)) else {
		node_index <- rep(FALSE, Nnode(tree))
		node_index[nodes] <- TRUE
	}
	
	node_index[is.na(variable)] <- FALSE
	
	if(!missing(sig.cutoff))
		node_index[variable <= sig.cutoff] <- FALSE
	
	nodes <- nodenumbers(tree)[node_index]
	plotvar <- plotvar[node_index]
	label <- label[node_index]
	
	if(new.window) quartz(width = 12, height = 10)
	if(show.legend)
	{
		oldpar <- par()
		par(mar = c(5,4,4,6) + 0.1)
	}
	
	plot(tree, show.tip = F, ...)
	title(main)
	nodelabels(text = as.character(label), node = nodes, cex = 0.7, bg = create.cols(plotvar, col)) 
	if(genera.lines) Put_genera_lines()
	if(show.legend)
	{
		library(fields)
		par <- oldpar
		image.plot(zlim = range(variable, na.rm = T), col = col, legend.only = T)
	}
}


# TODO can probably just be deleted
# plot_richness_nodesize <- function(dats = dispersion, cutoff = 0.95)  #refactor
# # a diagnostic plot showing the species richness of sites against the mean species richness of overrepresented nodes
# {
# 	# dats   : the representation matrix (resulting from measure_dispersion()) used to calculate the measures
# 	# cutoff : the cutoff value to use for significant over/under representation	 
# 	
# 	mean_richness_of_overrepresented_nodes <- numeric()
# 	for (i in 1:length(siteresults$cell))
# 	{
# 		site <- siteresults$cell[i]
# 		represented_nodes <- which(dats[rownames(dats) == as.character(site),] > cutoff)
# 		mean_richness_of_overrepresented_nodes[i] <- mean(noderesults$nodesize[represented_nodes])
# 	}
# 	plot(siteresults$richness, mean_richness_of_overrepresented_nodes, pch = 16)
# 	abline(lm(mean_richness_of_overrepresented_nodes ~ siteresults$richness))
# 	abline(h = mean(noderesults$nodesize), lty = 2)
# 	print(summary(lm(mean_richness_of_overrepresented_nodes ~ siteresults$richness)))
# }


#TODO refactor
plotsite = function(site, dats = dispersion_corr, tree = htree, comm = hcom, genera_lines = TRUE) #update to be reasonable
# an exploratory plotting function 
# given a site number, this shows the location of the site, and a tree illustrating which nodes are overrepresented at the site
{
	# site			:  the number of the site, given as an index between 1 and N (the total number of sites). This number corresponds to the row number of the site in the community matrix
	# dats   : the representation matrix (resulting from measure_dispersion())
	
	require(ape)
	quartz(width=20,height=12)
	layout(matrix(c(1,2),ncol = 2),width = c(4,4))

	spec_index = which(tree$tip.label %in% Site.species(site, comm))
	colortip <- rep("transparent",Ntip(tree)) 
	colortip[spec_index] <- "black"

	plot(tree,show.tip = F)	
	tiplabels(tip = spec_index, frame="none", pch=16, col="black",bg=colortip,adj=1,cex=1)
	overreps = nodenumbers(tree)[dats[site,] > 1.96]
	underreps = nodenumbers(tree)[dats[site,] < -1.96] 
	
	if(length(overreps) > 0) nodelabels(node = overreps, cex=0.9, text = " + " ,font=2,bg="black",frame="rect",adj=c(1.1,0.0),col="white")
	if(length(underreps) > 0) nodelabels(node = underreps, cex=0.9, text = " - " ,font=2,bg="white",frame="rect",adj=c(-0.2,0.0),col="black")
	
	if(genera_lines) Put_genera_lines()
	
	plot(siteresults$Long, siteresults$Lat, type = "n")
	background()
	points(siteresults$Long[site], siteresults$Lat[site],col = "red", pch = 16, cex = 1.5)
}






plotnode <- function(node_number, comm, tree = htree, coords = dat.LL, plottype = c("map", "points"), ...)
{
	plottype = match.arg(plottype)
	if(plottype == "map")
	{
		richs <- Node.richness(node_number, comm, tree, coords)
		map.var(richs$richness, richs$Long, richs$Lat, ...)
	} 	else  {
		richs <- Node.richness(node_number, comm, tree, coords)
		oldpar <- par()
		par(mar = c(5,4,4,4) + 0.1)
    plot(richs$Long, richs$Lat, col = create.cols(richs$richness), pch = 16, cex = 1.5, ...)#sqrt(richs$richness))
		plot(sa, add = T)
		library(fields)
		par <- oldpar
		image.plot(zlim = range(richs$richness, na.rm = T), col = rev(rainbow(256,start=0.02,end=0.59)), legend.only = T)
	}
}

plotoverrep <- function(node_number, dispersion_eff = rep_matrix, coords = dat.LL, tree = htree, plottype = c("map", "points"), new.window = FALSE,...)
{
	plottype = match.arg(plottype)
	if(plottype == "map")
	{
	  overrep <- dispersion_eff[,which(colnames(dispersion_eff) == as.character(node_number))]
    max_z = max(abs(overrep[abs(overrep) != Inf]), na.rm = T)
    zlims = c(-max_z,max_z) #should maybe a global range across nodes, but keep in mind to control for extreme values, maybe by cutting all values down to mean +- 3 sd .
		map.var(overrep, coords$Long, coords$Lat, new.window, zlim = zlims, ...)
	} 	else  {
		overrep <- dispersion_eff[,which(nodenumbers(tree) == node_number)]
		oldpar <- par()
		par(mar = c(5,4,4,4) + 0.1)
		max_z = max(abs(overrep), na.rm = T)
		zlims = c(-max_z,max_z)
		cols = create.cols(overrep, zlims = zlims)
    plot(coords$Long, coords$Lat, col = cols, pch = 16, cex = 1.5, ...)#sqrt(abs(overrep)))#, ...)
		plot(sa, add = T)
		library(fields)
		par <- oldpar
		image.plot(zlim = zlims, col = rev(rainbow(256,start=0.02,end=0.59)), legend.only = T)
	}
	
	title(paste("representation eff size of node",node_number))
}


# TODO rename and make an S3 method of function
plotnode_maps_new_parent <- function(node_number, par_rep_matrix = parent_rep_matrix, tree = htree, comm = hcom, coords = dat.LL, new.window = FALSE, plottype = c("map", "points"))
  # An exploratory plotting function
  # Given a node number, this draws four maps in the same window: 
  #  - a map of the location of overrepresented sites for that node
  #	- a map of the location of underrepresented sites for that node
  #	- a map of the species richness of species descending from that node
  #	- a map of the overall species richness
{
  # node_number   :  the internal ape node number (the number showed when calling nodelabels() on a plot of the phylogenetic tree)
  # dispersion	: the representation matrix (resulting from calling measure_dispersion())
  
  require(ape)
  
  if(new.window) quartz()
  par(mfrow = c(2,2))
  
  plotoverrep(node_number, par_rep_matrix, coords, tree, new.window = FALSE, plottype = plottype)
  
  plotnode(node_number, comm, tree, coords, new.window = FALSE, plottype = plottype, main = paste("richness of node", node_number))
  
  plotnode(Descendants(node_number)[1], comm, tree, coords, new.window = FALSE, plottype = plottype, main = paste("richness of descendant 1", Descendants(node_number)[1]))
  
  plotnode(Descendants(node_number)[2], comm, tree, coords, new.window = FALSE, plottype = plottype, main = paste("richness of descendant 2", Descendants(node_number)[2]))
  
}

