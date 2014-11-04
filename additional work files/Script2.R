#FUNCTIONS


#TODO replace with imports
source("Node_based_analysis.R")
source("Nodesig_function.R")

  
# Set up parallel processing (recommended, as the analysis is computationally heavy for large matrices)
#TODO this is an alternative method.
require(foreach)
require(doSNOW)

# dedicate the desired number of cores to parallel processing
number_of_cores <- 10
cl = makeCluster(number_of_cores)
registerDoSNOW(cl)
repeats <- 100
  

#TODO put all this in a function

show <- T #TODO remove
representation_matrix <- foreach (node = nodenumbers(htree), .combine = cbind, .packages = c("picante","vegan","ape"))  %dopar%
{
  # create the return vector of NA's. This means that all sites occupied by the parent node will retain the value NA in the final matrix
  ret = rep(NA, length(unique(hcom$Plot)))
  ret_ses <- ret
  ret_pval <- ret
  # return a vector of NAs if the node is the basal node (as that has no parent node)
  if (node == basal.node(htree))  
    # forloop res_list <- append(res_list, NA) else {
    return(NA)
    
  # take the subset of species that descend from the parent node
  ########################parentcom = Node.comm(Parent(node), hcom, htree)
  # create a species x site community matrix for those species
  ########################parentmat = sample2matrix(parentcom)
  # make sure the community matrix species are sorted in the same way as the tiplabels
  ########################parentmat = parentmat[,match(htree$tip.label, colnames(parentmat))] 
  # remove empty columns (may be created if Species is a factor in hcom)
  ########################parentmat = parentmat[,colSums(parentmat)>0]
  
  parentmat <- commat[,match(Node.species2(Parent(node)), colnames(commat))]

  
  # a boolean vector indicating which of all sites are considered for this node
  parnode.sites = (sort(unique(hcom$Plot)) %in% Node.sites(Parent(node),hcom, htree))
  # a boolean vector indicating which of species descending from the parent node that descend from the focal node
  node.sp = (colnames(parentmat) %in% Node.species2(node))	
  
  # A global variable to count the number of repeats
  res_object <- Nodesig(parentmat, node.sp, repeats, "quasiswap")
  if(!is.na(res_object)) res_object$sites <- parnode.sites

  save(res_object, file = paste(node))
  return(res_object)
  #forloopres_list <- append(res_list, res_object) 
  #forloop print(node)} 
}






### After, we need to generate GNO values, the SOS dispersion matrix and possibly other things (e.g. the noderesults with GNO should also have ancestors and descendants)

pp <- pval_sig(result$pval)
pp[pp > 1-2/repeats] <- 1-2/repeats # because otherwise may give values == 1 where logit is undefined.
GNO <- inv_logit(mean(logit(pp), na.rm = T))



# Assign names to the matrix margins  
rownames(dispersion_eff_) <- sort(unique(hcom$Plot))
colnames(dispersion_eff_pval) <- nodenumbers(htree)
  
# save the resulting matrix
save(dispersion_eff, file = "dispersion_eff_newtree.RData")

# create a summary list of site statistics (richness, node numbers etc.)
sitestatistics <- Sitestats(hcom, htree)
sitestatistics <- merge(sitestatistics, dat.LL, by = "cell")

# summarize the dispersion_eff matrix to create SES values of each node etc.
datalist <- create.datalist(dispersion_eff, sitestatistics)
save(datalist, hcom, htree, dat.LL, file = "Results.RData")


# 'datalist' is the final result. It is a list with four elements: 
#   siteresults, which summarizes the analysis over different sites
#   noderesults, which summarizes the analysis over individual nodes
#   rep_matrix, which is really just the dispersion_eff matrix
#   parent_rep_matrix, which summarizes the results over the parent node, to create the results shown in the paper

# The figures in the article are based on mapping individual columns of datalist$parent_rep_matrix, and by
# coloring node.labels on phylogenetic plots by datalist$noderesults$parent_overdisp. The rest of the information
# in datalist is included for the purposes of exploratory data analysis


