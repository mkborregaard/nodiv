



##########################################################
# Here comes a list of functions to use for the analysis. These definitions should all be loaded into R

create.cols <- function(vec, col , zlim)
{
  if(missing(col)) 
    if(min(col) < 0) col <- cm.colors(64) else col <- rev(heat.colors(64))
  if(min(vec, na.rm = TRUE) == 0 & identical(vec, floor(vec))) col <- c("grey", col)
  if(missing(zlim)) zlim <- c(min(vec,na.rm = T),max(vec,na.rm = T))
  vec = vec - zlim[1]
  vec = floor(vec * (length(col)-1)/(zlim[2]- zlim[1]))+1
  return(col[vec])
}


plot_grid <- function(x, coords, col = rev(terrain.colors(64)), shape = NULL, shapefill = "grey", zlim = NULL, zoom_to_points = FALSE, ...)
{
  if(inherits(x, "SpatialPixelsDataFrame"))
    rast <- raster(x) else
    if(missing(coords)) stop("coords must be defined if x is not a SpatialPixelsDataFrame") else
      if(inherits(coords, "SpatialPoints"))
      {
        if(!isGrid(coords)) stop("this function is only suitable for gridded data")
        coords <- SpatialPoints(coords)
        suppressWarnings(rast <- raster(SpatialPixelsDataFrame( coords, as.data.frame(x)))) #TODO NB
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
          suppressWarnings(rast <- raster(SpatialPixelsDataFrame(coords, as.data.frame(x)))) #TODO NB
        } else stop("Undefined arguments")
      }
  
  
  if(is.null(zlim)) zlim <- c(min(getValues(rast),na.rm = T),max(getValues(rast),na.rm = T))
  if(min(zlim) == 0 & identical(getValues(rast), floor(getValues(rast)))) col <- c("grey", col) #TODO an experimental hack
  
  if(is.null(shape)) plot(rast, zlim = zlim, col = col, ...) else
  {
    if(zoom_to_points)
      plot(shape, col = shapefill,  xlim = bbox(coords)[1,], ylim = bbox(coords)[2,], ...) else plot(shape, col = shapefill, ...)
    plot(rast, add = T, zlim = zlim, col = col)
  }
  invisible(rast)
}


plot_points <- function(x, coords, col = rev(terrain.colors(64)), shape = NULL, shapefill = "grey", zlim= NULL,  zoom_to_points = FALSE, pch = 16, bg = par("bg"), ...)
{  
  if(inherits(x, "SpatialPointsDataFrame"))
  {
    coords <- SpatialPoints(x)
    x <- as.numeric(x@data[,1])
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
            coords <- SpatialPoints(coords)
            
          } else stop("Undefined arguments")
        }
    }
  }

  
  if(is.null(zlim)) zlim <- c(min(x,na.rm = T),max(x,na.rm = T))
  coords <- SpatialPoints(coords)

  #split.screen( rbind(c(0, .8,0,1), c(.8,1,0,1)))
  #screen(1)
  oldpar <- par()
  par(mar = c(5,4,4,6) + 0.1)
  
  plotcol <- create.cols(x, col, zlim = zlim)
  if(pch %in% 21:25)
  {
    bg <- plotcol
    plotcol <- 1
  }
    
  
  if(is.null(shape)) plot(coords, col = plotcol, pch = pch, bg = bg, ...) else
  {
    if(inherits(shape, "Raster")) legend <- FALSE else legend <- NULL
    if(inherits(shape, "SpatialPolygonsDataFrame"))
      border <- shapefill else border <- NULL
    if(zoom_to_points)
      plot(shape, col = shapefill,  xlim = bbox(coords)[1,], ylim = bbox(coords)[2,], legend = legend, border = border, ...) else  plot(shape, col = shapefill, legend = legend, border = border, ...)
    plot(coords, col = plotcol, pch = pch, bg = bg, add = T)
  } 

  #screen(2)
  par <- oldpar

  add_legend(col = col, zlim = zlim)
  #image.plot( zlim = zlim,legend.only=TRUE, smallplot=c(.85,.87, .38,.65), col=col)
  ret <- data.frame(x)
  names(ret) <- deparse(substitute(x))
  invisible(SpatialPointsDataFrame(coords, ret))
}


plot_nodes_phylo <- function(variable, tree, label = variable, main = deparse(substitute(variable)), zlim, col = rev(heat.colors(64)), show.legend = TRUE, sig.cutoff, nodes, roundoff= TRUE, show.tip.label = NULL, cex = NULL, ...)
{
    if(!length(variable) == Nnode(tree))
    stop("The length of the variable vector must be the same length as the number of nodes on the tree")
  
  if(is.null(show.tip.label)) show.tip.label <- isTRUE(Ntip(tree) < 50)
  
  plotvar <- variable
  if(roundoff & is.numeric(label)) label = round(label,2)
  
  if(missing(zlim)) zlim <- c(min(plotvar, na.rm = T), max(plotvar, na.rm = T))
  
  if(is.null(cex)) cex <- par("cex")
  sizes <- cex * 4 * sqrt((plotvar - zlim[1])/zlim[2])
  
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
  #  split.screen( rbind(c(0, .9,0,1), c(.9,1,0,1)))
  #  screen(1)
    oldpar <- par()
    par(mar = c(5,4,4,6) + 0.1)
  }
  
  plot(tree, show.tip.label = show.tip.label, ...)
  title(main)
  nodelabels(pch = 16, node = nodes, col = create.cols(plotvar, col, zlim), cex = sizes)
  nodelabels(text = as.character(label), node = nodes, cex = 0.6, frame = "none") 
  
  if(show.legend)
  {
    par <- oldpar #screen(2)
    add_legend(col = col, zlim = zlim)
  }
}

add_legend <- function (zlim, smallplot=c(.85,.866, .38,.65), col)
{
  old.par <- par()
  if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
    stop("plot region too small to add legend\n")
  }
  
  ix <- 1
  minz <- zlim[1]
  maxz <- zlim[2]
  nlevel <- length(col)
  binwidth <- (maxz - minz)/nlevel
  midpoints <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
  iy <- midpoints
  iz <- matrix(iy, nrow = 1, ncol = length(iy))
  
  par(new = TRUE, pty = "m", plt = smallplot, err = -1)
  
  axis.args <- list(side =  4, mgp = c(3, 1, 0), las = 2)
  image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = col)
  
  do.call("axis", axis.args)
  box()
  par(new = FALSE, pty = old.par$pty, plt = old.par$plt, err = old.par$err)
  invisible()
}

