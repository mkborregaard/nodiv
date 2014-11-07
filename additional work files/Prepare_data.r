#TODO replace with import
library(picante)

# Definition of nodiv_data-class
# 
# spatialDataPoints coords : the points/grid
# phylo htree             : the tree
# commatrix commat        : the community matrix
# phylocom  hcom          : the phylocom dataset
# attribute(type = c("grid", "points")
#           #baseret på korteste afstand mellem unikke værdier
# 
# Definition of nodiv-class
# #samme som ovenfor (har det samme i sig) +
# 
# matrix parent_representation_matrix
# spatialGrid/PointsDataframe sitestatistics
# data.frame nodestatistics
# vector GND




nodiv_data <- function(phylo, commatrix, coords, proj4string = CRS(as.character(NA)), type = c("auto", "grid", "points"))
{
  type = match.arg(type)
  
  if(!(class(phylo) == "phylo")) stop ("phylo must be a phylogeny in the ape format") 
  
  if(!class(commatrix) == "distrib_data")
  {
    if(missing(coords)) stop("if commatrix is not an object of type distrib_data, coords must be specified")
    ret <- distrib_data(commatrix, coords, proj4string , type)
  } else ret <- commatrix
  
  cat("Comparing taxon names in phylogeny and communities with picante\n")
  dat <- match.phylo.comm(phylo, ret$commatrix)
  phylo <- dat$phy
  commatrix <- dat$comm
  if(!is.matrix(commatrix)) stop("The names do not match")
  
  ret <- match_commat_coords(commatrix, ret$coords)
  
  ret$phylo <- phylo
  ret$hcom <- matrix2sample(commat)
  class(ret) <- c("nodiv_data")
  return(ret)
}


distrib_data <- function(commatrix, coords, proj4string = CRS(as.character(NA)), type = c("auto", "grid", "points"))
{
  if(class(coords) == "SpatialPointsDataFrame" | class(coords) == "SpatialPixelsDataFrame")
    if(!proj4string == proj4string(coords))
    { 
      proj4string <- proj4string(coords)
      warning("specified proj4string overridden by the coords data")
    } 

  if(is.data.frame(commatrix) & ncol(commatrix) == 3 & !is.numeric(commatrix[,3])) commatrix <- sample2matrix(commatrix) #i.e. is the commatrix in phylocom format?
  
  if(is.data.frame(commatrix)) commatrix <- as.matrix(commatrix)
  if(!is.matrix(commatrix)) stop("commatrix must be a matrix of 0's and 1's, indicating presence or absence")
  if(!all.equal(sort(unique(as.numeric(commatrix))), 0:1)) stop("commatrix must be a matrix of 0's and 1's, indicating presence or absence")
  
  if(is.matrix(coords)) coords <- as.data.frame(coords)
  if(is.data.frame(coords)) coords <- toSpatialPoints(coords, commatrix, type, proj4string)
  if(!(class(coords) == "SpatialPointsDataFrame" | class(coords) == "SpatialPixelsDataFrame")) stop("coords must be a data.frame of coordinates or a SpatialPoints object")
  
  ret <- match_commat_coords(commatrix, coords)
  class(ret) <- "distrib_data"
  return(ret)

}


#internal functions

#TODO #much of this testing can be done with try-catch phrases
#use the testthat library to test everything


#TODO implement plot function for the result object



match_commat_coords <- function(commatrix, spatialPointData)
{
  cat("Comparing sites in community data and spatial points\n")
  
  if(!nrow(commatrix) == nrow(spatialPointData)) 
  {  
    if(is.null(rownames(commatrix)))
      stop("The number of sites in coords and the data matrix do not fit") else rownames(commatrix) <- spatialPointData$cell
    if(sum(spatialPointData$cell %in% rownames(commatrix)) < 2)
      stop("the coordinate names and the rownammes of the community matrix do not match")
    
    spatialPointData <- spatialPointData[spatialPointData$cell %in% rownames(commatrix),]
  }
    
  commatrix <- commatrix[match(spatialPointData$cell, rownames(commatrix)),]
  return(list(commatrix = commatrix, coords = spatialPointData)) 
}




toSpatialPoints <- function(coords, commatrix, type, proj4string)
{
    xcol <- 0
    ycol <- 0
    
    ret <- coords
    
    colnames(ret) <- tolower(colnames(ret))
    if('x' %in% colnames(ret) & 'y' %in% colnames(ret))
    {
      xcol <- which(colnames(ret) == 'x')
      ycol <- which(colnames(ret) == 'y')
      ret <- ret[,c(xcol, ycol)]
      
    } else if('lon' %in% substr(colnames(ret),1,3) & 'lat' %in% substr(colnames(ret), 1, 3)) {
      
      colnames(ret) <- substr(colnames(ret), 1, 3)
      xcol <- which(colnames(ret) == 'lon')[1]
      ycol <- which(colnames(ret) == 'lat')[1]
      ret <- ret[,c( xcol, ycol)]
    }
    
    if(!ncol(ret) == 2) stop("ret should be a data.frame or spatial data.frame with 2 columns, giving the x/longitude, and y/latitude of all sites")
    
    names(ret) <- c("Long", "Lat")
    
    sitenames <- character()
    if (ncol(coords)==3 & !(xcol + ycol == 0)) sitenames <- coords[,-c(xcol, ycol)] else
      if(nrow(coords) == nrow(commatrix) & !is.null(rownames(commatrix))) sitenames <- rownames(commatrix) else
        if(!is.null(rownames(coords))) sitenames <- rownames(coords) else 
          stop("There must be valid site names in the rownames of commatrix or in the coords data")
  
  
  type_auto <- ifelse(isGrid(ret), "grid", "points")
  
  if(type == "auto") type <- type_auto else 
    if(!type == type_auto)
      warning("The specified type of data (points or grid) seems to conflict with the automatic setting. This may cause problems")
  
  ret <- SpatialPoints(ret, proj4string)
  if(type == "grid") ret <- SpatialPixelsDataFrame(ret, data.frame(cell = sitenames)) else
    ret <- SpatialPointsDataFrame(ret, data.frame(cell = sitenames))
  
  return(ret)  
}

isGrid <- function(coords)
  return(isGridVar(coordinates(coords)[,1]) & isGridVar(coordinates(coords)[,2]))

isGridVar <- function(gridVar)
{
  dists <- diff(sort(unique(gridVar)))
  distab <- table(dists)
  smallest <- as.numeric(names(distab[1]))
  most_common <- as.numeric(names(distab))[distab == max(distab)]
  return(all.equal(dists/smallest, floor(dists/smallest)) & smallest %in% most_common)
  #if all differences are a multiplum of the smallest, and the smallest distance is the most common, it is probably a grid
}



# 
# 
# commat <- read.delim("../OCCUR_ALL.xls")
# hcom <- matrix2sample(commat)
# 
# htree <- read.tree("../phylogeny.txt")
# 
# commat = commat[, match(htree$tip.label, colnames(commat))]
# 
# source("Node_based_analysis.R")
# 
# node_species = Create_node_by_species_matrix(htree)
# 
# dat.LL <- NA
# 
# save(hcom, dat.LL, htree, commat, node_species, file = "GlobMammals_data.RData")
# 
