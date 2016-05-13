#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram ###
ui <- shinyUI(fluidPage(
   ## Application title ######
   titlePanel("Cluster Ensemble Consensus Solution"),
   
   ## Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         radioButtons("nClusts", "Number of Clusters:",
             c("2" = 2,
               "3" = 3,
               "4" = 4,
               "5" = 5)),
         h3('Documentation'),
         h4('1. Select the number of clusters above.'),
         h4('2. The application will create an ensemble of paritions.'),
         h4('3. From the ensemble, the clue package will construct a consensus partition.'),
         h4('4. The consensus partion will be profiled by:'),
         h5('  a. a simple barchart of cluster frequencies and'),
         h5('  b. a silhoutte chart based on original data distances')), 
         
   # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot"),
         plotOutput("silPlot") 
       )
      )
))

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
  library(clue); library(cluster)
  # Load data, keep only complete case obs, standardize attributes   
  data(rock)
  myRock0 <- na.omit(rock) # listwise deletion of missing
  myRock  <- scale(myRock0) # standardize variables (area, peri, shape,  perm)
  setCols=c("magenta","blue","red","green","orange")
  
  # K means K's and methods, all combos
  kmeans_k <- c(2, 4, 5, 6, 7, 8)
  kmeans_methods <- c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen")
  kmeans_results=list() 
  for (k in kmeans_k){
    kmeans_results <- c(kmeans_results,lapply(kmeans_methods, function(m) kmeans(myRock, k, algorithm = m, iter.max=100)))
  }
  km_names <- character() 
  for (i in 1:length(kmeans_k))
    km_names <- c(km_names,paste(kmeans_methods,kmeans_k[i]))
  names(kmeans_results) <- km_names 
  myDist <- daisy(rock)
  
  # Hierachical clustering, K's and methods, all combos
  hclust_k       <- c(2, 4, 5, 6, 7, 8)
  hclust_methods <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")			
  hclust_results=list() 
  for (i in 1:length(hclust_methods))
  {
    hc <- hclust(dist(myRock)^2, hclust_methods[i])
    hclust_results <- c(hclust_results,lapply(hclust_k, function(m) as.cl_partition(cutree(hc, k = m))))
  }
  hc_names <- character() 
  for (i in 1:length(hclust_methods))
    hc_names <- c(hc_names,paste(hclust_methods[i],hclust_k))
  names(hclust_results) <- hc_names 
  
  # Combine Results into ensemble
  All <- cl_ensemble(list = c(kmeans_results, hclust_results))
  
   output$distPlot <- renderPlot({
     ce_solution <- as.cl_hard_partition(cl_consensus(All, control=list(k=as.integer(input$nClusts))))[[1]]  # 2 Cluster Consensus Solution
     barplot(table(ce_solution), col=setCols)
   })
   
   output$silPlot <- renderPlot({
     ce_solution <- as.cl_hard_partition(cl_consensus(All, control=list(k=as.integer(input$nClusts))))[[1]]  # 2 Cluster Consensus Solution
     colFreqs <- table(unlist(ce_solution))
     myCols <- integer()
     for (i in 1:length(colFreqs))
       myCols <- c(myCols,rep(setCols[i],colFreqs[i]))
     plot(silhouette(ce_solution, myDist),main="Actual Clusters",col=myCols,border=NA)
   })
})

# Run the application 
shinyApp(ui = ui, server = server)

