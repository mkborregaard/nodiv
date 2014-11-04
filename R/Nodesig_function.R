
## These functions are all internal

pval_sig <- function(pval) 
  2*(0.5-abs(pval-0.5))

logit <- function(p)
{
  if(!(sum(p >= 1, na.rm = T) + sum(p <= 0, na.rm = T) == 0))
    warning("logit is only defined between the interval 0 and 1")
  p[p == 0] <- NA
  p[p == 1] <- NA
  return(log(p/(1-p)))
}

inv_logit <- function(a)
{
  return(exp(a)/(1 + exp(a)))
}


quasiswap_nodesig <- function(simcom, node.sp, repeats, show = F)
{
  ll <- 0
  if (show) pb <- txtProgressBar(min = 0, max = repeats, style = 3)
  replicate(repeats-1,
  {
    simcom <- commsimulator(simcom, method = "quasiswap")
    rowSums(simcom[, node.sp])
    if (show) setTxtProgressBar(pb, ll <<- ll + 1)
    })
}

rdtable_nodesig <- function(simcom, node.sp, repeats, show = F)
{
  tempcom <- cbind(rowSums(simcom[,node.sp]), rowSums(simcom[,!node.sp]))
  r <- rowSums(tempcom)
  c <- colSums(tempcom)
  ll <- 0
  if (show) pb <- txtProgressBar(min = 0, max = repeats, style = 3)
  
  # Use the quasiswap algorithm to created random matrices, basing each new matrix on the previous
  nodereps <- replicate(repeats-1,
  {
    
    # Illustrate the progress in the chart window
    # perform the randomization
    simcom <- r2dtable(1, r, c)[[1]]
    # return the simulated species richness of sites
    if (show) setTxtProgressBar(pb, ll <<- ll + 1)
    simcom[, 1]
    
  })
  return(nodereps)
}


## TODO write a wrapper for this function, that will allow you to compare any two clades
Nodesig <- function(commat, node.sp, repeats = 100, method = c("quasiswap","rdtable"), show = F)
{
  if(sum(node.sp)== 1 | sum(!node.sp) == 1) return(NA) #if one of the descendant clades is a single species
  method = match.arg(method)
  # A global variable to count the number of repeats
  require(vegan)
  require(fields)
  simcom <- commat
  
  nodereps <- switch(method,
         quasiswap = quasiswap_nodesig(simcom, node.sp, repeats, show),
         rdtable = rdtable_nodesig(simcom, node.sp, repeats, show)
  )
  
  
  nodeemp <- rowSums(commat[, node.sp])
  nodereps <- cbind(nodeemp, nodereps)
  ord <- apply(nodereps, 1,  rank)[1,]
  ord[ord == repeats] <- repeats - 1
  ord[rowSums(commat) == 0] <- NA
  
  #calculates the effect size (take care of sd = 0)!
  nodemeans <- rowMeans(nodereps)
  nodesds <- apply(nodereps, 1, sd)
  
  # all sites occupied by the parent node get an SES value - the rest retain NA
  ses <- (nodeemp - nodemeans )/ nodesds
  pval <-  ord/repeats
  
  return(list(ses = ses,  pval = pval, nodeemp = nodeemp, nodemeans = nodemeans, nodesds = nodesds))
}
