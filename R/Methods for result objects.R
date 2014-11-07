







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
  
  spec_index = which(tree$tip.label %in% Site_species(site, comm))
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
  
  plotnode(Descendants(node_number, tree)[1], comm, tree, coords, new.window = FALSE, plottype = plottype, main = paste("richness of descendant 1", Descendants(node_number)[1]))
  
  plotnode(Descendants(node_number, tree)[2], comm, tree, coords, new.window = FALSE, plottype = plottype, main = paste("richness of descendant 2", Descendants(node_number)[2]))
  
}


