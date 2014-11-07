



##########################################################
# Here comes a list of functions to use for the analysis. These definitions should all be loaded into R

create.cols <- function(vec, col = rev(rainbow(256,start=0.02,end=0.59)), zlims = c(min(vec,na.rm = T),max(vec,na.rm = T)))
{
  vec = vec - zlims[1]
  vec = floor(vec * (length(col)-1)/(zlims[2]- zlims[1]))+1
  return(col[vec])
}


map_var <- function(x, coords, polygon = NULL, ...)
{
  if(inherits(x, "SpatialPixelsDataFrame"))
    rast <- raster(x) else
    if(missing(coords)) stop("coords must be defined if x is not a SpatialPixelsDataFrame") else
      if(inherits(coords, "SpatialPoints"))
      {
        if(!isGrid(coords)) stop("this function is only suitable for gridded data")
        coords <- SpatialPoints(coords)
        rast <- raster(SpatialPixelsDataFrame( coords, as.data.frame(x)))
      } else {
        if(is.matrix(coords)) 
          coords <- as.data.frame(coords)
        if(is.data.frame(coords))
        {
          if(!is.vector(x)) stop("x must be a vector or a SpatialPixelsDataFrame")
          if(!nrow(coords) == length(x)) stop("x and coords must have the same number of elements")
          if(!ncol(coords) == 2) stop("coords must have exactly two columns, for the x and y coordinates of data")
          if(!isGrid(coords)) stop("this function is only suitable for gridded data")
          coords <- SpatialPoints(coords)
          rast <- raster(SpatialPixelDataFrame(coords, as.data.frame(x)))
        } else stop("Undefined arguments")
      }
  
  plot(rast, ...)
  if(!is.null(polygon)) plot(polygon, add = T, col = "darkgrey")
  invisible(rast)
}

plot_nodes_phylo <- function(variable, label = variable, tree, main = deparse(substitute(variable)), zlims, col = rev(rainbow(256,start=0.02,end=0.59)), show.legend = TRUE, sig.cutoff, nodes, roundoff= TRUE, genera.lines = FALSE, ...)
  # plots a tree, where the colors of the nodes reflects the values of variable
{
  # variable  : the variable that controls the colors of nodes - given for ALL nodes, even though some nodes are not plotted!
  # tree   	: the phylogenetic tree
  # main 		: the title to give the plot
  # new.window: should the plot be made in a new window?
  # col		: the color palette to use for making the color plot
  
  if(!length(variable) == Nnode(tree))
    stop("The length of the variable vector must be the same length as the number of nodes on the tree")
  
  plotvar <- variable
  if(roundoff & is.numeric(label)) label = round(label,2)
  
  if(missing(zlims)) zlims <- c(min(plotvar, na.rm = T), max(plotvar, na.rm = T))
  
  sizes <- par("cex") * 3 * sqrt((plotvar - zlims[1])/zlims[2])
  
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
  sizes <- sizes[node_index]
  
  
  if(show.legend)
  {
    oldpar <- par()
    par(mar = c(5,4,4,6) + 0.1)
  }
  
  plot(tree, show.tip = F, ...)
  title(main)
  nodelabels(pch = 16, node = nodes, col = create.cols(plotvar, col, zlims), cex = sizes)
  nodelabels(text = as.character(label), node = nodes, cex = 0.6, frame = "none") 
  
  if(show.legend)
  {
    #     #use this instead of using 'fields' or use 'fields' for all the plotting
    #     for(i in seq_along(valcols)) segments(0, 0 + i/6, 0.03 * root_tip_length, 0 + i/6, col = valcols[i])
    #     text(0.035 * root_tip_length, 1, min(vals), adj = c(0,0.5))
    #     text(0.035 * root_tip_length, length(valcols)/6, round(max(vals),1), adj = c(0,0.5))
    
    library(fields)
    par <- oldpar
    image.plot(zlim = range(variable, na.rm = T), col = col, legend.only = T, zlim = zlims)
  }
}


plot_points <- function(x, coords, pch = 16, cols = rev(rainbow(256,start=0.02,end=0.59)),  ...)
{  
  if(inherits(x, "SpatialPointsDataFrame"))
  {
    x <- x@data
    coords <- SpatialPoints(x)
  } else {
      if(missing(coords)) stop("coords must be defined if x is not a SpatialPixelsDataFrame") else {
        if(!inherits(coords, "SpatialPoints"))
        {
          if(is.matrix(coords)) coords <- as.data.frame(coords) 
          if(is.data.frame(coords))
          {
            if(!is.vector(x)) stop("x must be a vector or a SpatialPointsDataFrame")
            if(!nrow(coords) == length(x)) stop("x and coords must have the same number of elements")
            if(!ncol(coords) == 2) stop("coords must have exactly two columns, for the x and y coordinates of data")
            
          } else stop("Undefined arguments")
        }
    }
  }
  
  coords <- SpatialPoints(coords)

  oldpar <- par()
  par(mar = c(5,4,4,4) + 0.1)
  plot(coords, col = create.cols(x, cols), pch = pch, ...)
  
  library(fields) #TODO replace imports
  par <- oldpar
  image.plot(zlim = range(x, na.rm = T), col = cols, legend.only = T)
}

# map.var = function(variable, Long, Lat, new.window = FALSE, xrange = range(Long), yrange = range(Lat), bordercolor = "darkgrey", cut.to.coast = FALSE, ...)
# # a function that will create a pretty color-coded map of any variable with Long/Lat information
# {
# 	# variable : the variable controlling the color
# 	# Long : longitude of each point
# 	# Lat : latitude of each point
# 	# new.window : should the plot be created in a new window? Set to FALSE, e.g. for saving to a postscript file
# 	# ... : extra arguments to be passed to the "image.plot" function

# 	require(fields)
# 	require(maptools)
# 	variable[is.na(variable)] <- -9999  #I am doing this trick to ensure that the whole range of Lat/Long values are printed, not just the ones with valid values of 'variable'
# 	cellsize_x <- min(diff(sort(unique(Long))),na.rm = T)
# 	cellsize_y <- min(diff(sort(unique(Lat))),na.rm = T)
# 	map <- as.image(variable, x = cbind(Long, Lat), grid = list(x = seq(xrange[1],xrange[2], by = cellsize_x), y = seq(yrange[1],yrange[2],by = cellsize_y)))
# 	map$z[map$z == -9999] <- NA

# 	if(new.window) quartz()
# 	image.plot(map, ...)
#   sa <- readShapePoly('~/shapecoast/World_coast.shp')
# 	plot(sa,add=T, border = bordercolor)
# 	if(cut.to.coast)
# 	{
# 		sc <- readShapePoly('~/shapecoast/World_coast.shp')
# 		plot(sc, add = T, border = NULL, col = "white", usePolypath = TRUE)
# 		box()
# 	}
# }



