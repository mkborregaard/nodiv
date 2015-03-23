



##########################################################
# Here comes a list of functions to use for the analysis. These definitions should all be loaded into R

parula <- function(){
  vals <- data.frame(r <- c(0.2081, 0.2116, 0.2123, 0.2081, 0.1959, 0.1707, 0.1253, 0.0591, 0.0117, 0.006, 0.0165, 0.0329, 0.0498, 0.0629, 0.0723, 0.0779, 0.0793, 0.0749, 0.0641, 0.0488, 0.0343, 0.0265, 0.0239, 0.0231, 0.0228, 0.0267, 0.0384, 0.059, 0.0843, 0.1133, 0.1453, 0.1801, 0.2178, 0.2586, 0.3022, 0.3482, 0.3953, 0.442, 0.4871, 0.53, 0.5709, 0.6099, 0.6473, 0.6834, 0.7184, 0.7525, 0.7858, 0.8185, 0.8507, 0.8824, 0.9139, 0.945, 0.9739, 0.9938, 0.999, 0.9955, 0.988, 0.9789, 0.9697, 0.9626, 0.9589, 0.9598, 0.9661, 0.9763), 
                     g <- c(0.1663, 0.1898, 0.2138, 0.2386, 0.2645, 0.2919, 0.3242, 0.3598, 0.3875, 0.4086, 0.4266, 0.443, 0.4586, 0.4737, 0.4887, 0.504, 0.52, 0.5375, 0.557, 0.5772, 0.5966, 0.6137, 0.6287, 0.6418, 0.6535, 0.6642, 0.6743, 0.6838, 0.6928, 0.7015, 0.7098, 0.7177, 0.725, 0.7317, 0.7376, 0.7424, 0.7459, 0.7481, 0.7491, 0.7491, 0.7485, 0.7473, 0.7456, 0.7435, 0.7411, 0.7384, 0.7356, 0.7327, 0.7299, 0.7274, 0.7258, 0.7261, 0.7314, 0.7455, 0.7653, 0.7861, 0.8066, 0.8271, 0.8481, 0.8705, 0.8949, 0.9218, 0.9514, 0.9831), 
                     b <- c(0.5292, 0.5777, 0.627, 0.6771, 0.7279, 0.7792, 0.8303, 0.8683, 0.882, 0.8828, 0.8786, 0.872, 0.8641, 0.8554, 0.8467, 0.8384, 0.8312, 0.8263, 0.824, 0.8228, 0.8199, 0.8135, 0.8038, 0.7913, 0.7768, 0.7607, 0.7436, 0.7254, 0.7062, 0.6859, 0.6646, 0.6424, 0.6193, 0.5954, 0.5712, 0.5473, 0.5244, 0.5033, 0.484, 0.4661, 0.4494, 0.4337, 0.4188, 0.4044, 0.3905, 0.3768, 0.3633, 0.3498, 0.336, 0.3217, 0.3063, 0.2886, 0.2666, 0.2403, 0.2164, 0.1967, 0.1794, 0.1633, 0.1475, 0.1309, 0.1132, 0.0948, 0.0755, 0.0538))
  
  ret <- rgb(vals)
  return(ret)
}

choose.colors <- function(vec, zlim = NULL, type = c("auto", "ramp", "monochrome", "divergent", "individual"))
{

  type = match.arg(type)

  if(is.null(zlim)){
    zlim <- range(vec, na.rm = TRUE)
    
    if(zlim[1] * zlim[2] < 0)
      zlim <- c(-max(abs(vec), na.rm = TRUE), max(abs(vec), na.rm = TRUE))
  }
  if(is.character(vec))
    vec <- as.factor(vec)
  if(type == "auto"){
    if(is.factor(vec))
      type <- "individual" else {
   
      if(zlim[1] * zlim[2] < 0){
        if(sum(zlim) == 0){
         type <- "divergent" 
        } else {
          warning("ramp color scheme chosen - if you want divergent colors use zlims symmetric around 0")
          type <- "ramp"
        }
         
      } else type <- "ramp" 
    }
  }
  
  if(type == "individual"){
    n <- length(unique(vec))
    if(requireNamespace("RColorBrewer")){
      if(n <= 9)
        ret <- RColorBrewer::brewer.pal(n, "Set1") else
          if(n <= 12)
            ret <- RColorBrewer::brewer.pal(n, "Set3") else
              if(n <= 21)
                ret <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(12, "Set3"))[1:n] else 
                  ret <- rainbow(n)
    } else {
      if(n <= 8)
        ret <- palette()[1:n] else
        {
          if(requireNamespace("colorspace"))
            ret <- colorspace::diverge_hcl(64) else 
              ret <- rainbow(n)
        }
          
    }
  }
  
  if(type == "ramp"){
    ret <- parula()
    #if(requireNamespace("fields"))
    #  ret <- fields::tim.colors(64) else
    #    ret <- rev(terrain.colors(64))
  }
  
  if(type == "divergent"){
    if(requireNamespace("RColorBrewer"))
      ret <- RColorBrewer::brewer.pal(11, "RdBu") else
        if(requireNamespace("colorspace"))
          ret <- colorspace::diverge_hcl(64) else
            ret <- cm.colors(64)
  }
  
  if(type == "monochrome"){
    if(requireNamespace("RColorBrewer"))
      ret <- RColorBrewer::brewer.pal(9, "YlOrRd") else
        if(requireNamespace("colorspace"))
          ret <- colorspace::sequential_hcl(64) else
            ret <- rev(heat.colors(64))
  }
  
  if(length(unique(na.omit(vec))) == 2)
    ret <- c("darkgrey", "red")
  
  attr(ret, "zlim") <- zlim
  
  return(ret)
}

create.cols <- function(vec, col, zlim, type = c("auto", "ramp", "monochrome", "divergent", "individual"))
{
  type = match.arg(type)
  
  if(missing(zlim)){
    zlim <- range(vec, na.rm = TRUE)
    if(zlim[1] * zlim[2] < 0)
      zlim <- c(-max(abs(vec), na.rm = TRUE), max(abs(vec), na.rm = TRUE))
  }
  
  if(missing(col)) 
    col <- choose.colors(vec, zlim, type)
  
  if(length(col) == length(unique(vec)))
    return(col[match(vec, sort(unique(vec)))])
    
  if(min(vec, na.rm = TRUE) == 0 & identical(vec, floor(vec)) & length(col) > 3) col <- c("grey", col)

  vec = vec - zlim[1]
  vec = floor(vec * (length(col)-1)/(zlim[2]- zlim[1]))+1
  return(col[vec])
}


plot_grid <- function(x, coords, col, shape = NULL, shapefill = "grey", zlim = NULL, zoom_to_points = FALSE, ...)
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
  
  if(missing(col)){
    if(is.character(x)) x <- as.factor(x)
    if(is.factor(x))
      type <- "individual" else type <- "auto"      
    col <- choose.colors(getValues(rast), zlim, type = type)
    zlim <- attr(col, "zlim")
  }
  
  if(is.null(zlim)) zlim <- c(min(getValues(rast),na.rm = T),max(getValues(rast),na.rm = T))
  if(min(zlim) == 0 & identical(getValues(rast), floor(getValues(rast)))) col <- c("grey", col) #TODO an experimental hack
  

  
  if(is.null(shape)) plot(rast, zlim = zlim, col = col, ...) else
  {
    if(inherits(shape, "SpatialPolygonsDataFrame"))
      border <- shapefill else border <- NULL
    if(zoom_to_points)
      plot(shape, col = shapefill, border = border, xlim = bbox(coords)[1,], ylim = bbox(coords)[2,], ...) else plot(shape, col = shapefill, border = border, ...)
    plot(rast, add = T, zlim = zlim, col = col)
  }
  invisible(rast)
}


plot_points <- function(x, coords, col , shape = NULL, shapefill = "grey", zlim= NULL,  zoom_to_points = FALSE, pch = 16, bg = par("bg"), ...)
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

  if(missing(col)){
    col <- choose.colors(x, zlim)
    zlim <- attr(col, "zlim")
  }
  
  if(is.null(zlim)) zlim <- c(min(x,na.rm = T),max(x,na.rm = T))
  
  
  coords <- SpatialPoints(coords)

  #split.screen( rbind(c(0, .8,0,1), c(.8,1,0,1)))
  #screen(1)
  oldpar <- par()
  #par(plt = c(0,0.8,0,1))
  par(mar = c(5,4,4,6) + 0.1)
  
  plotcol <- create.cols(x, col, zlim = zlim)
  if(pch %in% 21:25)
  {
    bg <- plotcol
    plotcol <- 1
  }
  coords <- coords[order(x),]
  plotcol <- plotcol[order(x)]  
  
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


plot_nodes_phylo <- function(variable, tree, label = variable, main = deparse(substitute(variable)), zlim = NULL, col , show.legend = TRUE, sig.cutoff, nodes, roundoff= TRUE, show.tip.label = NULL, cex = NULL, ...)
{
    if(!length(variable) == Nnode(tree))
    stop("The length of the variable vector must be the same length as the number of nodes on the tree")
  
  if(is.null(show.tip.label)) show.tip.label <- isTRUE(Ntip(tree) < 50)
  
  plotvar <- variable
  if(roundoff & is.numeric(label)) label = round(label,2)
  

  if(missing(col)){
    col <- choose.colors(plotvar, zlim)
    zlim <- attr(col, "zlim")
  }
  
  if(is.null(zlim)) zlim <- c(min(plotvar, na.rm = T), max(plotvar, na.rm = T))
  
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

