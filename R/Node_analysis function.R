

#INTERNAL FUNCTION  
nodiv_anal <- function(node, nodiv_data, repeats, method)
{
  # create the return vector of NA's. This means that all sites occupied by the parent node will retain the value NA in the final matrix
  ret = rep(NA, length(unique(nodiv_data$hcom$plot)))
  ret_ses <- ret
  ret_pval <- ret
  # return a vector of NAs if the node is the basal node (as that has no parent node)
  if (node == basal.node(nodiv_data$phylo))  
    # forloop res_list <- append(res_list, NA) else {
    return(NA)
  
  parentmat <- nodiv_data$commatrix[,match(Node_species2(Parent(node, nodiv_data$phylo)), colnames(nodiv_data$commatrix))] #TODO or just nodiv$species
  
  #TODO do I need this?
  # remove empty columns (may be created if id is a factor in hcom)
  ########################parentmat = parentmat[,colSums(parentmat)>0]
  
  # a boolean vector indicating which of all sites are considered for this node
  parNode_sites = (sort(unique(nodiv_data$hcom$plot)) %in% Node_sites(Parent(node, nodiv_data$phylo),nodiv_data$hcom, nodiv_data$phylo)) #TODO make into an S3 method for nodiv_data
  
  # a boolean vector indicating which of species descending from the parent node that descend from the focal node
  Node_sp = (colnames(parentmat) %in% Node_species2(node))  #TODO also just species?
  
  # A global variable to count the number of repeats
  res_object <- Nodesig(parentmat, Node_sp, repeats, method)
  if(length(res_object) > 1) #i.e. if it is not an NA
      res_object$sites <- parNode_sites
  return(res_object)
}
  
##TODO make nodenumbers a class function  
# Exported interface

Node_analysis <- function(nodiv_data, repeats = 100, cores = 1, method = c("quasiswap", "rdtable"))
{
  require(foreach)
  if(!"nodiv_data" %in% class(nodiv_data)) stop("This function only works on objects of class nodiv_data")
  
  if(cores > 1)
  {

    require(doSNOW)
  
  # dedicate the desired number of cores to parallel processing
    cl = makeCluster(cores)
    registerDoSNOW(cl)

  }
    pb <- txtProgressBar(min = Ntip(nodiv_data$phylo) + 1, max = Ntip(nodiv_data$phylo) + Nnode(nodiv_data$phylo), style = 3)
    
    results <- foreach (node = nodenumbers(nodiv_data$phylo), .packages = c("picante","vegan","ape"))  %dopar% 
    {
      setTxtProgressBar(pb, node)
      nodiv_anal(node, nodiv_data, repeats, method)
    } 
  
  return(results)
}
  






