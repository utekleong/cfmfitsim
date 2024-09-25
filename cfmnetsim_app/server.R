##################################################################
##                          Setting up                          ##
##################################################################
# Loading required packages

library(shiny)
library(parSim)
library(bootnet)
library(tidyverse)
library(psychonetrics)

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

#################################################################
##                         Server code                         ##
#################################################################
function(input, output, session) {
  
  #weights matrix upload
  wmat <- reactive({
    req(input$upload)
    
    ext <- tools::file_ext(input$upload$name)
    switch(ext,
           csv = vroom::vroom(input$upload$datapath, delim = ","),
           validate("Invalid file; Please upload a .csv file"))
  })
  
  #notification to allow users to check dimensions of matrix
  observeEvent(input$upload, {
    showNotification(paste0("Success! You have uploaded a ", nrow(wmat()), "x", ncol(wmat()), " weights matrix"), type = "warning")
  })
  
  #simulation triggered by button press
  
  simres <- reactiveValues()
  
  observeEvent(input$beginsim, {
    
    waiter <- waiter::Waiter$new()
    waiter$show()
    on.exit(waiter$hide())
    
    simres$x <- parSim(
      
      # Enter conditions to vary over here (a vector means you vary over those conditions):
      n_original = as.numeric(input$n_original), # Sample size of original dataset
      n_replication = eval(parse(text = input$n_replication)), # Sample size of replication dataset
      p_reorder = eval(parse(text = input$p_reorder)), # Node reorder probability of generating network for replication dataset. 0 should mean replication should be successful
      
      # Other arguments:
      reps = input$reps,  # Number of simulation conditions
      nCores = 1,  # Number of computer cores to use (set to 1 to debut easier)
      progress = FALSE,
      
      # Enter the simulation code here (the condition arguments can be used as object names):
      expression = {
        # Needed packages in simulation:
        library("bootnet")
        library("psychonetrics")
        library("dplyr")
        library("qgraph")
        
        # Generate an original dataset:  
        data_original <- ggmGenerator()(n_original, wmat())
        
        # Little code to enforce proper network structure:
        n_test <- 0
        repeat{
          wmat2 <- reorder(wmat(), p_reorder)
          ev <- eigen(diag(ncol(wmat2))-wmat2)$values
          if (all(ev > 0)) break
          if (n_test > 10) stop("No positive definite network to be made")
          n_test <- n_test + 1
        }
        
        # Generate a replication dataset:
        data_replication <- ggmGenerator()(n_replication, wmat2)
        
        # Fit an EBICglasso network to original data:
        net_original_raw <- eval(parse(input$networkcode))
        net_original <- net_original$graph
        
        # Obtain adjacency matrix from the "original" network:
        adj_original <- 1*(net_original!=0)
        
        # Fit confirmatory model:
        fit_replication <- ggm(data_replication, omega = adj_original) %>% runmodel
        
        # Get weights matrix with significant edges only:
        net_replication <- getmatrix(fit_replication, "omega", threshold = TRUE, alpha = 0.05)
        
        # Add also network comparisons:
        source("./scripts/compareNetworks.R")
        
        # Replication network for semi-confirmatory:
        fit_replication@fitmeasures <- c(fit_replication@fitmeasures,CompareNetworks(wmat()[41,-41], net_replication[41,-41]))
        
        # Label the results:
        fit_replication@fitmeasures$type <- "semi"
        
        # Return replication fit measures:
        as.data.frame( fit_replication@fitmeasures)
      })
  })
  
  output$downloadsimres <- downloadHandler(
    filename = function() {
      ("simresdata.csv")
    },
    content = function(file) {
      write.csv(simres$x, file)
    })
  
}