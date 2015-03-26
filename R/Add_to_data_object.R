# I made some major changes here and to infer_sites, that also need to be done for species!
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
  
  if(length(site) < Nsites(distrib_data))
    cat(paste(Nsites(distrib_data)- length(site), "sites were not found in", deparse(substitute(distrib_data)), "\n"))
  
  mergeframe <- as.data.frame(sapply(sitestat, function(column) {
    ret <- rep(NA, Nsites(distrib_data))
    ret[site] <- column
    ret
  }), stringsAsFactors = FALSE)
  
  if(sum(names(mergeframe) %in% names(distrib_data$coords@data)) > 0){
    matches <- which(names(distrib_data$coords@data) %in% names(mergeframe))
    deleted <- names(distrib_data$coords@data)[matches]
    distrib_data$coords <- distrib_data$coords[, -matches]
    warning(paste("Some data in the original distrib_data overwritten:\n"), paste(deleted, collapse = "\t"))
  }
  
  
  distrib_data$coords@data <- cbind(distrib_data$coords@data, mergeframe) 
  distrib_data
}

add_species_stat <- function(distrib_data, species_stat, specs = NULL){
  #restorepoint::restore.point("add_species_stat")
  
  
  if(!inherits(distrib_data, "distrib_data"))
    stop("distrib_data must be an object of type distrib_data, nodiv_data or nodiv_result")
  if(is.null(distrib_data$species_stats))
    stop("The distrib_data object is from an earlier version of nodiv. Please run update_object on the object before proceeding")
  
  if(is.matrix(species_stat))
    species_stat <- as.data.frame(species_stat, stringsAsFactors = FALSE)
  
  if(is.vector(species_stat)){
    nam <- deparse(substitute(species_stat))
    species_stat <- as.data.frame(species_stat, stringsAsFactors = FALSE)
    names(species_stat) <- nam
  } 
  
  if(is.factor(species_stat)){
    nam <- deparse(substitute(species_stat))
    species_stat <- as.data.frame(species_stat)
    names(species_stat) <- nam
  } 
  
  num <- nrow(species_stat)
  
  if(is.null(specs))
  {
    temp <- infer_species(distrib_data, species_stat)
    specs <- temp$species
    species_stat <- temp$species_stat
  }
    
  specs <- identify_species(specs, distrib_data)
  
  if(length(specs) < Nspecies(distrib_data))
    cat(paste(num - length(specs), "species were not found in", deparse(substitute(distrib_data)), "\n"))
  
  mergeframe <- as.data.frame(sapply(species_stat, function(column) {
    ret <- rep(NA, Nspecies(distrib_data))
    if(is.factor(column)) ret <- factor(ret, levels = levels(column))
    ret[specs] <- column
    ret
  }), stringsAsFactors = FALSE)
  
  if(sum(names(mergeframe) %in% names(distrib_data$species_stats)) > 0){
    matches <- which(names(distrib_data$species_stats) %in% names(mergeframe))
    deleted <- names(distrib_data$species_stats)[matches]
    distrib_data$species_stats[,matches] <- NULL
    warning(paste("Some data in the original distrib_data overwritten:\n"), paste(deleted, sep = "\t"))
  }
  
  
  distrib_data$species_stats <- cbind(distrib_data$species_stats, mergeframe) 
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
  
  ##### We need a matching function here to do the actual matching!
  
  sitestat <- subrow_data.frame(sitestat, which(!is.na(site)))
  site <- site[!is.na(site)]
  
  suppressWarnings(sitenames <- identify_sites(site, distrib_data, as.name = TRUE))
  matchsite <- match(site, sitenames)
  
  sitestat <- subrow_data.frame(sitestat, which(!is.na(matchsite)))
  site <- site[!is.na(matchsite)]

  cat(paste("Matching sites by", name, "\n"))
    
  return(list(site = site, sitestat = sitestat))
}

match_speciesnames <- function(reference_name, new_name){
  chars <- c(" ", "[.]", "_")
  ref_ll <- sapply(chars, function(x) length(grep(x, reference_name)))
  ref_char <- chars[which(ref_ll == max(ref_ll))[1]]
  new_ll <- sapply(chars, function(x) length(grep(x, new_name)))
  new_char <- chars[which(new_ll == max(new_ll))[1]]  
  
  ret <- match(gsub(new_char, ref_char, new_name), reference_name)
  ret
}

infer_species <- function(distrib_data, species_stat) # a non-exported convenience function
{
    species_stat$rownames <- rownames(species_stat)
    temp <- sapply(1:length(species_stat), function(index){
      matches <- match_speciesnames(species(distrib_data), species_stat[[index]])
      return(sum(!is.na(matches))/nrow(species_stat))
    })
    res <- which(temp == max(temp))[1]
    name <- names(species_stat)[res]
    spec <- species(distrib_data)[match_speciesnames(species(distrib_data), species_stat[[res]])]

  if(temp[res] < 0.5 & sum(spec %in% species(distrib_data)) < 0.5 * Nspecies(distrib_data))
    stop("Species could not be matched automatically, please supply the species argument explicitly")
  
  if(!name == "rownames")
    species_stat[[name]] <- NULL
  
  species_stat$rownames <- NULL
  

  species_stat <- subrow_data.frame(species_stat, which(!is.na(spec)))
  spec <- spec[!is.na(spec)]
  
  
  suppressWarnings(specnames <- identify_species(spec, distrib_data, as.name = TRUE))
  matchspec <- match(spec, specnames)
  
  species_stat <- subrow_data.frame(species_stat, which(!is.na(matchspec)))
  spec <- spec[!is.na(matchspec)]
  
  cat(paste("Matching species by", name, "\n"))
  
  return(list(species = spec, species_stat = species_stat))
}

