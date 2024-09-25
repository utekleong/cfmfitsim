
library(shiny)

# defining UI for application
fluidPage(
  
  waiter::use_waiter(),
  
  # Application title
  titlePanel("CFM-NET Fit Evaluation"),
  
  navlistPanel(
    id = "navlist",
    
    "Simulation",
    tabPanel("Import data",
             fileInput("upload", "Upload the .csv file of your weights matrix", multiple = FALSE, accept = ".csv")),
    tabPanel("Enter simulation parameters",
             numericInput("n_original", "Sample size of original dataset", value = 0, min = 0),
             textInput("n_replication", "Enter a vector of replication sample sizes to vary over", placeholder = "c(100, 200, 300)"),
             textInput("p_reorder", "Enter a vector of node reorder probabilities to vary over", placeholder = "c(0, 0.1, 0.2, 0.3)"),
             numericInput("reps", "Enter the number of simulation replications to perform", value = 10, min = 10, max = 1000)),
    tabPanel("Enter network code",
             textAreaInput("networkcode", "Enter the code to fit your ORIGINAL network here, replacing the name of the data variable with data_original", 
                       placeholder = "estimateNetwork(data_original, 
                                      default = 'EBICglasso', 
                                      corMethod = 'spearman',
                                      missing = 'pairwise',
                                      sampleSize = 'pairwise_average',
                                      tuning = 0.5)",
                       width = 600, height = 300),
             actionButton("beginsim", "Begin Simulation", class = "btn-block btn-success")
             ),
    
    "Results",
    tabPanel("Download simulation results",
             downloadButton("downloadsimres", "Download", class = "btn-lrg btn-success"))
    
  )
)