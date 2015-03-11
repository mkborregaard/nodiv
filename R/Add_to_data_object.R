add_sitestat <- function(distrib_data, sitestat, site = NULL){
  if(!inherits(distrib_data, "distrib_data"))
    stop("distrib_data must be an object of type distrib_data, nodiv_data or nodiv_result")

  if(is.matrix(sitestat))
    sitestat <- as.data.frame(sitestat, stringsAsFactors = FALSE)
  
  if(is.vector(sitestat)){
    nam <- deparse(substitute(sitestat))
    sitestat <- as.data.frame(sitestat, stringsAsFactors = FALSE)
    names(sitestat) <- nam
  } 

  
  if(is.null(site))
  {
    temp <- infer_sites(distrib_data, sitestat)
    site <- temp$site
    sitestat <- temp$sitestat
  }
  
  site <- identify_sites(site, distrib_data)
  
  mergeframe <- as.data.frame(sapply(sitestat, function(column) {
    ret <- rep(NA, Nsites(distrib_data))
    ret[site] <- column
    ret
  }), stringsAsFactors = FALSE)
  
  if(sum(names(mergeframe) %in% names(distrib_data$coords@data)) > 0){
    matches <- which(names(distrib_data$coords@data) %in% names(mergeframe))
    deleted <- names(distrib_data$coords@data)[matches]
    distrib_data$coords@data[,matches] <- NULL
    warning(paste("Some data in the original distrib_data overwritten:\n"), paste(deleted, sep = "\t"))
  }
  
  
  distrib_data$coords@data <- cbind(distrib_data$coords@data, mergeframe) 
  distrib_data
}
  
infer_sites <- function(distrib_data, sitestat) # a non-exported convenience function
{

  suppressWarnings(numsites <- as.numeric(sites(distrib_data)))
  if(sum(is.na(numsites)) < 0.2*Nsites(distrib_data)){
    if(all.equal(numsites, floor(numsites))) {      # if site names are just integers, matching is not attempted
      if(nrow(sitestat) == Nsites(distrib_data))
        return(list(site = sites(distrib_data), sitestat = sitestat)) else 
          warning("Site matching was done based on name matching, which is tricky when site names are integer values")
    }
  }
  
  potentials <- list()
  if(!is.null(sitestat$sites))
    potentials$sites <- sitestat$sites
  
  if(!is.null(sitestat$site))
    potentials$site <- sitestat$site
  
  if(!is.null(sitestat$ID))
    potentials$ID <- sitestat$ID
  
  if(!is.null(sitestat$Sites))
    potentials$sites <- sitestat$sites
  
  if(!is.null(sitestat$Site))
    potentials$site <- sitestat$site
  
  if(!is.null(sitestat$id))
    potentials$ID <- sitestat$ID

  if(!is.null(rownames(sitestat)))
    potentials$rownames <- rownames(sitestat)

  temp <- sapply(1:length(potentials), function(index){
    matches <- sum(potentials[[index]] %in% sites(distrib_data))
    return(matches/nrow(sitestat))
  })
  
  res <- which(temp == max(temp))[1]
  name <- names(potentials)[res]
  site <- potentials[[res]]
  
  if(temp[res] < 0.8){
    temp <- sapply(1:length(sitestat), function(index){
      matches <- sum(sitestat[[index]] %in% sites(distrib_data))
      return(matches/nrow(sitestat))
    })
    res <- which(temp == max(temp))[1]
    name <- names(sitestat)[res]
    site <- sitestat[[res]]
  }
    
  if(temp[res] < 0.8)
    stop("Sites could not be matched automatically, please supply the site argument explicitly")
  
  if(!name == "rownames")
    sitestat[[name]] <- NULL
  
  cat(paste("Matching sites by", name, "\n"))
    
  return(list(site = site, sitestat = sitestat))
}


