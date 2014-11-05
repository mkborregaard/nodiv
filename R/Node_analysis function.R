

#INTERNAL FUNCTION  
nodiv_anal <- function(node, nodiv_data, repeats, method = c("quasiswap", "rdtable"))
{
  # create the return vector of NA's. This means that all sites occupied by the parent node will retain the value NA in the final matrix
  ret = rep(NA, length(unique(hcom$plot)))
  ret_ses <- ret
  ret_pval <- ret
  # return a vector of NAs if the node is the basal node (as that has no parent node)
  if (node == basal.node(htree))  
    # forloop res_list <- append(res_list, NA) else {
    return(NA)
  
  #TODO delete this
  # take the subset of species that descend from the parent node
  ########################parentcom = Node_comm(Parent(node, htree), hcom, htree)
  # create a species x site community matrix for those species
  ########################parentmat = sample2matrix(parentcom)
  # make sure the community matrix species are sorted in the same way as the tiplabels
  ########################parentmat = parentmat[,match(htree$tip.label, colnames(parentmat))] 
  # remove empty columns (may be created if Species is a factor in hcom)
  ########################parentmat = parentmat[,colSums(parentmat)>0]
  
  parentmat <- commat[,match(Node_species2(Parent(node, htree)), colnames(commat))]
  
  
  # a boolean vector indicating which of all sites are considered for this node
  parNode_sites = (sort(unique(hcom$plot)) %in% Node_sites(Parent(node, htree),hcom, htree))
  # a boolean vector indicating which of species descending from the parent node that descend from the focal node
  Node_sp = (colnames(parentmat) %in% Node_species2(node))  
  
  # A global variable to count the number of repeats
  res_object <- Nodesig(parentmat, Node_sp, repeats, method)
  if(!is.na(res_object)) res_object$sites <- parNode_sites
  return(res_object)
  
  
  
# Exported interface

Node_analysis <- function(nodiv_data, repeats = 100, cores = 1)
{
  if("nodiv_data" %in% class(nodiv_data) stop("This function only works on objects of class nodiv_data"))
  
  if(cores > 1)
  {
    require(foreach)
    require(doSNOW)
  
  # dedicate the desired number of cores to parallel processing
    cl = makeCluster(cores)
    registerDoSNOW(cl)
  }
  
  results <- foreach (node = nodenumbers(htree), .combine = cbind, .packages = c("picante","vegan","ape"))  %dopar% 
  {
    nodiv_anal(node, nodiv_data)
  }
  
  return(results)
}
  