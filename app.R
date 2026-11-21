library(shiny)
library(ggplot2)
library(Biobase)
library(NMF)

# Read the data
tb_data <- readRDS("TBnanostring.rds")
tb_status <- factor(tb_data$TB_Status)

genes <- as.matrix(
  tb_data[, names(tb_data) != "TB_Status", drop = FALSE]
)

# Calculate the three reductions once when the app starts
pca_fit <- prcomp(genes, center = TRUE, scale. = FALSE)
pca_coordinates <- pca_fit$x[, 1:2]

set.seed(123)
umap_coordinates <- uwot::umap(
  genes,
  n_neighbors = 15,
  n_components = 2
)

shifted_genes <- sweep(
  genes, 2, apply(genes, 2, min), FUN = "-"
)

nmf_fit <- nmf(
  shifted_genes,
  rank = 2,
  nrun = 10,
  seed = 123
)

nmf_coordinates <- basis(nmf_fit)

# Store the coordinates under the method names
reductions <- list(
  PCA = pca_coordinates,
  UMAP = umap_coordinates,
  NMF = nmf_coordinates
)

# Build the interface
ui <- fluidPage(
  titlePanel("TB NanoString clustering"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "method",
        "Dimension reduction method",
        choices = c("PCA", "UMAP", "NMF"),
        selected = "PCA"
      ),
      
      numericInput(
        "clusters",
        "Number of clusters",
        value = 2,
        min = 1,
        max = 9,
        step = 1
      ),
      
      helpText("Colors show clusters. Shapes show TB status."),
      helpText("Black X marks show cluster centers.")
    ),
    
    mainPanel(
      plotOutput("cluster_plot", height = "420px"),
      h4("Clusters versus TB status"),
      tableOutput("cluster_table"),
      textOutput("ari")
    )
  )
)

# Define how the app responds to inputs
server <- function(input, output, session) {
  
  selected_data <- reactive({
    coordinates <- reductions[[input$method]]
    
    data.frame(
      Dimension1 = coordinates[, 1],
      Dimension2 = coordinates[, 2]
    )
  })
  
  clusters <- reactive({
    req(input$clusters)
    
    validate(
      need(
        input$clusters >= 1 &&
          input$clusters <= 9 &&
          input$clusters == floor(input$clusters),
        "Choose a whole number from 1 to 9."
      )
    )
    
    set.seed(123)
    
    kmeans(
      selected_data(),
      centers = input$clusters,
      nstart = 25
    )
  })
  
  output$cluster_plot <- renderPlot({
    plot_data <- selected_data()
    plot_data$Cluster <- factor(clusters()$cluster)
    plot_data$TB_Status <- tb_status
    
    centers <- as.data.frame(clusters()$centers)
    
    axis_labels <- switch(
      input$method,
      PCA = c("PC1", "PC2"),
      UMAP = c("UMAP1", "UMAP2"),
      NMF = c("Component 1 weight", "Component 2 weight")
    )
    
    ggplot(
      plot_data,
      aes(
        x = Dimension1,
        y = Dimension2,
        color = Cluster,
        shape = TB_Status
      )
    ) +
      geom_point(size = 2.5, alpha = 0.8) +
      geom_point(
        data = centers,
        aes(x = Dimension1, y = Dimension2),
        inherit.aes = FALSE,
        shape = 4,
        size = 5,
        stroke = 2,
        color = "black"
      ) +
      labs(
        title = paste(input$method, "with", input$clusters, "clusters"),
        x = axis_labels[1],
        y = axis_labels[2],
        color = "Cluster",
        shape = "TB status"
      ) +
      theme_minimal(base_size = 13)
  })
  
  output$cluster_table <- renderTable({
    counts <- table(
      Cluster = clusters()$cluster,
      TB_Status = tb_status
    )
    
    data.frame(
      Cluster = rownames(counts),
      as.data.frame.matrix(counts),
      row.names = NULL,
      check.names = FALSE
    )
  }, striped = TRUE, bordered = TRUE)
  
  output$ari <- renderText({
    score <- mclust::adjustedRandIndex(
      clusters()$cluster,
      tb_status
    )
    
    paste("Adjusted Rand index (ARI):", round(score, 3))
  })
}

shinyApp(ui = ui, server = server)