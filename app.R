#####################################################################
####################### Popgraph in Shiny ###############################
#####################################################################
## Matt DeSaix
# July 19, 2018

library(leaflet)
library(shiny)
library(gstudio)
library(popgraph)
library(igraph)
library(ggmap)
library(ggplot2)
library(ggrepel)

data <- readRDS("df_AB_gstudio_2alleles.RDS")
fst <- readRDS("loc_fst_ranked_2alleles.RDS")
sites <- read.csv("Sampling_Sites.csv")
sites.breeding <- sites[1:11,]

ui <- fluidPage(
  titlePanel("Population Graphs"),
  sidebarLayout(
    sidebarPanel(
      # Use the Slider to adjust number of loci
        sliderInput("n","Number of Loci", value = 100, min = 0, max = 4000, step = 20),
        # 'Plot' button for initiating popgraph
        actionButton("go", label = "Plot"),
        # a couple of side tables
        tableOutput("fst"),
        tableOutput("pca_df")
    ),
    mainPanel(
      # Break up app into three different 
      tabsetPanel(
        tabPanel("PCA and Popgraph",
                 fluidRow(
                   column( width = 12, plotOutput("pca")),
                   column( width = 12, plotOutput("popgraph")))
        ),
        tabPanel("Popgraph on a map",
                 leafletOutput("leaf", height = 600) 
        ),
        tabPanel("Sites",
                 tableOutput("samples"))
      )
    )
  )
)

server <- function(input,output) {
  output$samples <- renderTable(
    sites.breeding
  )
  # Fst table of info, updates all the time
  output$fst <- renderTable({
    fst.sub <- fst[1:input$n,]
    fst.df <- data.frame((matrix(nrow = 4, ncol = 2)))
    fst.df[,1] <- c("n Loci", "Min Fst", "Mean Fst", "Max Fst")
    fst.df[,2] <- c(input$n, min(fst.sub[,2]), mean(fst.sub[,2]), max(fst.sub[,2]))
    colnames(fst.df) <- c(" ", " ")
    fst.df
  })
  ## Create table of just summary of PC1 and PC2
  output$pca_df <- renderTable({
    fit.pca <- prcomp(x())
    pca.df <- data.frame( matrix( nrow = 3, ncol = 3) )
    colnames(pca.df) <- c(" ", "PC1", "PC2")
    pca.df[,1] <- c("SD", "Proportion Var", "Cumulative Var")
    pca.df[1,2:3] <- fit.pca$sdev[1:2]
    cumul <- cumsum(fit.pca$sdev^2 / sum(fit.pca$sdev^2))
    pca.df[3,2:3] <- cumul[1:2]
    pca.df[2,2:3] <- c(cumul[1], (cumul[2] - cumul[1]))
    pca.df
    
  })
  ## The reactive input that makes the plots only update when the action button is clicked
  # saves the output from gstudio's 'to_mv()' function as x
  # NOTE: when this x is used in the future, it is represented by 'x()', because it is reactive
  x <- eventReactive(input$go, {
    fst.sub <- fst[1:input$n,]
    data.sub <- cbind( data[,1:3], data[,colnames(data) %in% fst.sub$Locus])
    to_mv( data.sub, drop.allele = FALSE)
  })
  
  ## Leaflet map
  output$leaf <- renderLeaflet({
    
    ## Making all the data happen
    graph <- popgraph( x = x(), groups = data$Pop)
    
    coords <- strata_coordinates(sites[1:11,], stratum = "Population", 
                                 longitude = "lon", 
                                 latitude = "lat")
    map_graph <- decorate_graph(graph, coords, stratum = "Stratum")
    
    ## the object 'map_graph' is a popgraph/igraph, not helpful for plotting in leaflet
    ## this for loop pulls out the names of all the connections listed in the popgraph
    # the output of this loop is a vector of vertices names, and the connections are listed next to each other
    # i.e. full.verts[1] and full.verts[2] are connected, full.verts[3] and full.verts[4] are connected, etc.
    full.verts <- c()
    for(i in 1:length(map_graph)){
      tmp.verts <- names(unlist(map_graph[[i]])) %>%
        strsplit(split = "[.]") %>%
        unlist()
      full.verts <- c(full.verts, tmp.verts)
    }
    ## convert full.verts to a df and add lat/lon
    df.verts <- as.data.frame( matrix( nrow = length(full.verts), ncol = 3))
    colnames(df.verts) <- c("verts", "lon", "lat")
    df.verts$verts <- full.verts
    df.verts$lat <- sites[match(df.verts$verts, sites$Population), "lat"]
    df.verts$lon <- sites[match(df.verts$verts, sites$Population), "lon"]
    
    ## get the nodes and sizes
    node.df <- data.frame("node" = V(map_graph)$name, 
                          "size" = V(map_graph)$size)
    node.df$node <- as.character(node.df$node)
    sites.breeding$Population <- as.character(sites.breeding$Population)
    sites.breeding$size <- node.df[match(sites.breeding$Population,node.df$node), "size"]
    
    
    ## Plot the leafmap
    leafmap <- leaflet(data = sites.breeding) %>%
                addProviderTiles(providers$Stamen.TerrainBackground)
    ## This loop is how to add the lines, because addPolylines will connected all lat/long listed in the df
    # Therefore, this loops through df.verts[1]&[2], [3] & [4], etc. and draws the connections one at a time
    for(i in 1:(nrow(df.verts)/2)){
      l1 <- 2*i-1
      leafmap <- addPolylines(leafmap,
                              lng = df.verts$lon[l1:(l1+1)],
                              lat = df.verts$lat[l1:(l1+1)],
                              color = "white",
                              opacity = 0.75)
    } 
    # formatting the circle markers
    leafmap <- addCircleMarkers(leafmap,
                                data = sites.breeding[c(1:6,9,11),],
                                lng = ~lon, 
                                lat = ~lat,
                                radius = ~size,
                                opacity = 0.6,
                                label = ~Population,
                                labelOptions = labelOptions(noHide = T, 
                                                            textOnly = T,
                                                            textsize = "14px",
                                                            direction = "right",
                                                            offset = c(15,0))) %>%
      addCircleMarkers(data = sites.breeding[c(7,8,10),],
                       lng = ~lon, 
                       lat = ~lat,
                       radius = ~size,
                       opacity = 0.6,
                       label = ~Population,
                       labelOptions = labelOptions(noHide = T, 
                                                   textOnly = T,
                                                   textsize = "14px",
                                                   direction = "left",
                                                   offset = c(-15,0)))
    leafmap
  })

  ## Popgraph
  output$popgraph <- renderPlot({
    ## Popgraph
    graph <- popgraph( x = x(), groups = data$Pop)
    coords <- strata_coordinates(sites.breeding, stratum = "Population", 
                                 longitude = "lon", 
                                 latitude = "lat")
    map_graph <- decorate_graph(graph, coords, stratum = "Stratum")
    ## use the fruchterman.reingold method, could set up interactive toggles here...
    c <- layout.fruchterman.reingold(map_graph)
    V(map_graph)$x <- c[,1]
    V(map_graph)$y <- c[,2]
    V(map_graph)$region <- c("Mississippi", "Mississippi", "Mississippi", 
                             "Mississippi", "Atlantic", "Mississippi", "Atlantic", 
                             "Atlantic", "Atlantic", "Atlantic", "Mississippi")
    
    frre <- data.frame("x" = c[,1], "y" = c[,2], "population" = V(map_graph)$name, "region" = V(map_graph)$region)
    
    
    q2 <- ggplot() + 
      geom_edgeset( aes(x, y), color = "darkgrey", map_graph) + 
      geom_nodeset( aes(x, y, size = size, color = region), alpha = 0.5, map_graph) +
      theme_empty() + guides(color = F, size = F) +
      geom_text_repel(aes(x, y, label = population), frre)
    q2
    })
  ## PCA
  # it'd be nice to color by pop or region
  output$pca <- renderPlot({
    fit.pca <- prcomp(x())
    # summary(fit.pca)
    pred <- predict(fit.pca)
    df <- data.frame(PC1 = pred[,1], PC2 = pred[,2], Pop = data$Pop)
    ggplot(df) + 
      geom_point(aes(x = PC1, y = PC2, color = Pop), size = 3, alpha = 0.75) +
      theme_bw()
    
  })
  
  
}

shinyApp(ui, server)