

# #INTERNAL FUNCTION  
# nodiv_anal_old <- function(node, nodiv_data, repeats, method)
# {
#   # create the return vector of NA's. This means that all sites occupied by the parent node will retain the value NA in the final matrix
#   ret <- rep(NA, Nsites(nodiv_data))
#   ret_ses <- ret
#   ret_pval <- ret
#   # return a vector of NAs if the node is the basal node (as that has no parent node)
#   if (node == basal.node(nodiv_data$phylo))  
#     # forloop res_list <- append(res_list, NA) else {
#     return(NA)
  
#   parentmat <- nodiv_data$commatrix[,match(Node_species2(Parent(node, nodiv_data$phylo)), colnames(nodiv_data$commatrix))] #TODO or just nodiv$species
  
#   #TODO do I need this?
#   # remove empty columns (may be created if id is a factor in hcom)
#   ########################parentmat = parentmat[,colSums(parentmat)>0]
  
#   # a boolean vector indicating which of all sites are considered for this node
#   parNode_sites <- (sort(unique(nodiv_data$hcom$plot)) %in% Node_sites(Parent(node, nodiv_data$phylo),nodiv_data$hcom, nodiv_data$phylo)) #TODO make into an S3 method for nodiv_data
  
#   # a boolean vector indicating which of species descending from the parent node that descend from the focal node
#   Node_sp <- (colnames(parentmat) %in% Node_species2(node))  #TODO also just species?
  
#   # A global variable to count the number of repeats
#   res_object <- Nodesig(parentmat, Node_sp, repeats, method)
#   if(length(res_object) > 1) #i.e. if it is not an NA
#       res_object$sites <- parNode_sites
#   return(res_object)
# }
 
nodiv_anal <- function(node, nodiv_data, repeats, method)
{
  # return a vector of NAs if the node is the basal node (as that has no parent node)
  if (node == basal_node(nodiv_data))  
    # forloop res_list <- append(res_list, NA) else {
    return(data.frame(ses = rep(NA, Nsites(nodiv_data)),  pval = rep(NA, Nsites(nodiv_data)), nodeemp = rep(NA, Nsites(nodiv_data)), nodemeans = rep(NA, Nsites(nodiv_data)), nodesds = rep(NA, Nsites(nodiv_data))))
  
  parent_data <- subsample(nodiv_data, node = Parent(node, nodiv_data))
  
  parNode_sites <- match(parent_data$coords$sites, nodiv_data$coords$sites)
  # a vector indicating which of all sites are considered for this node
  
  # a boolean vector indicating which of species descending from the parent node that descend from the focal node
  Node_sp <- (parent_data$species %in% Node_species(node, nodiv_data))  #TODO also just species?
  
  res_object <- Nodesig(parent_data, Node_sp, repeats, method)

  res_object <- lapply(res_object, function(x) #ensure that dimensions are the same for merging into matrices
  {
    ret <- rep(NA, Nsites(nodiv_data))
    ret[parNode_sites] <- x
    ret
  })
  
  return(as.data.frame(res_object))
}

##TODO make nodenumbers a class function  
# Exported interface

Node_analysis <- function(nodiv_data, repeats = 100, method = c("rdtable", "quasiswap"), cores = 1)
{
  if(!inherits(nodiv_data, "nodiv_data")) stop("This function only works on objects of class nodiv_data")
  method <- match.arg(method)
 
  if(cores > 1)
  {
    stop("The parallel interface is currently unimplemented") 
#   require(parallel)
  # dedicate the desired number of cores to parallel processing
#   cl = makeCluster(cores)  
#   pb <- txtProgressBar(min = Nspecies(nodiv_data) + 1, max = Ntip(nodiv_data$phylo) + Nnode(nodiv_data$phylo), style = 3)
  
#   results <- parLapply(nodenumbers(nodiv_data), function(node)
#   {
#     setTxtProgressBar(pb, node)
#     ret <- nodiv_anal(node, nodiv_data, repeats, method)
#     ret[,1:2]
#   })

#   stopCluster(cl)
    
  } else {
    
    pb <- txtProgressBar(min = Nspecies(nodiv_data) + 1, max = Ntip(nodiv_data$phylo) + Nnode(nodiv_data$phylo), style = 3)
    results <- lapply(nodenumbers(nodiv_data), function(node)
    {
      setTxtProgressBar(pb, node)
      ret <- nodiv_anal(node, nodiv_data, repeats, method)
      ret[,1:2] #to keep object size down, I am currently only using these two variables for the big analysis
    })
    
  }
  
  nodiv_result <- nodiv_res(results, nodiv_data, repeats, method)
  
  return(nodiv_result)
}
  


