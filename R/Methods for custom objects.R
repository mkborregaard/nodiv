

Nspecies <- function(x)
{
  if (!inherits(x, "distrib_data")) 
    stop("object is not of class \"distrib_data\"")
  length(x$species)
}

Nsites<- function(x)
{
  if (!inherits(x, "distrib_data")) 
    stop("object is not of class \"distrib_data\"")
  length(x$sites)
}

print.distrib_data <- function(x, printlen = 5, ...)
{
  cat(paste("Data object with distributions of", Nspecies(x),"species in", Nsites(x),"sites\n\n"))
  cat("Species names:\n")
  cat(paste("\t", paste(x$species[1:printlen], collapse = ", "),"\n\n"))
  cat("Site names:\n")
  cat(paste("\t", paste(x$sites[1:printlen], collapse = ", "),"\n\n"))
}

