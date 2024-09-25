#################################################################
##                  Loading packages and data                  ##
#################################################################
set.seed(4242)

# Loading required packages
library("parSim")
library("bootnet")
library("dplyr")
library("tidyr")
library("ggplot2")
library("psychonetrics")

# Loading data
wmat <- as.matrix(read.csv("./data/weightsmatrix.csv")) #loads weights matrix from TPESS study
rownames(wmat) <- colnames(wmat)

# Loading custom reorder function
reorder <- function(mat, p_reorder = 0){
  n_node <- ncol(mat)
  nodes <- seq_len(n_node)
  nodes_to_reorder <- which(sample(c(TRUE,FALSE),n_node,TRUE,prob=c(p_reorder,1-p_reorder)))
  if (length(nodes_to_reorder) > 0){
    for (i in nodes_to_reorder){
      change_to <- sample(1:n_node,1)
      n1 <- nodes[i]
      n2 <- nodes[change_to]
      nodes[i] <- n2
      nodes[change_to] <- n1
    }
    mat <- mat[nodes,nodes]
  }
  return(mat)
}

##################################################################
##                          Simulation                          ##
##################################################################

simres <- parSim(
  
  # Enter conditions to vary over here (a vector means you vary over those conditions):
  n_original = 907, # Sample size of original dataset
  n_replication = c(100, 200, 300, 400, 500), # Sample size of replication dataset
  p_reorder = c(0, 0.1, 0.2, 0.3), # Node reorder probability of generating network for replication dataset. 0 should mean replication should be successful
  
  # Other arguments:
  reps = 100,  # Number of simulation conditions
  nCores = 16,  # Number of computer cores to use (set to 1 to debut easier)
  progress = FALSE,
  export = c("reorder","wmat"), # Objects to export to the simulation code. Here I will use the estimated weights matrix as assumed true network structure
  
  # Enter the simulation code here (the condition arguments can be used as object names):
  expression = {
    # Needed packages in simulation:
    library("bootnet")
    library("psychonetrics")
    library("dplyr")
    library("qgraph")
    
    # Generate an original dataset:  
    data_original <- ggmGenerator()(n_original, wmat)
    
    # Little code to enforce proper network structure:
    n_test <- 0
    repeat{
      wmat2 <- reorder(wmat, p_reorder)
      ev <- eigen(diag(ncol(wmat2))-wmat2)$values
      if (all(ev > 0)) break
      if (n_test > 10) stop("No positive definite network to be made")
      n_test <- n_test + 1
    }
    
    # Generate a replication dataset:
    data_replication <- ggmGenerator()(n_replication, wmat2)
    
    # Fit an EBICglasso network to original data:
    net_original <- estimateNetwork(data_original, 
                                    default = "EBICglasso", 
                                    corMethod = "spearman",
                                    missing = "pairwise",
                                    sampleSize = "pairwise_average",
                                    tuning = 0.5)$graph
    
    # Obtain adjacency matrix from the "original" network:
    adj_original <- 1*(net_original!=0)
    
    # Fit confirmatory model:
    fit_replication <- ggm(data_replication, omega = adj_original) %>% runmodel
    
    # Get weights matrix with significant edges only:
    net_replication <- getmatrix(fit_replication, "omega", threshold = TRUE, alpha = 0.05)
    
    # Add also network comparisons:
    source("./scripts/compareNetworks.R")
    
    # Replication network for semi-confirmatory (TPESS only):
    fit_replication@fitmeasures <- c(fit_replication@fitmeasures,CompareNetworks(wmat[41,-41], net_replication[41,-41]))
    
    # Label the results:
    fit_replication@fitmeasures$type <- "semi"
    
    # Return replication fit measures:
      as.data.frame( fit_replication@fitmeasures)
  })

# Save results:
write.csv(simres,"./data/simres.csv")