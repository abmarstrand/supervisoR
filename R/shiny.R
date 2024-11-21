#' Launch Pathway Visualization Shiny Application
#'
#' This function launches a Shiny application for visualizing pathway enrichment scores
#' and their relationships within a pathway graph. Users can interactively explore
#' top-level pathways and their sub-pathways, view glyphs representing enrichment scores
#' across multiple conditions, and reset the view to the initial state.
#'
#' @importFrom shiny shinyApp fluidPage sidebarLayout sidebarPanel mainPanel selectInput textInput actionButton plotOutput titlePanel h4 p br
#' @importFrom visNetwork visIgraphLayout visNodes visEdges visOptions visEvents visNetworkOutput renderVisNetwork visInteraction visNetwork
#' @importFrom ggplot2 ggplotGrob
#' @importFrom magrittr %>%
#' @importFrom igraph degree V add_vertices induced_subgraph `V<-` add_edges drl_defaults
#' @param enrichment_scores A data.frame with pathways as rows and conditions as columns,
#'   containing enrichment scores. Pathways should correspond to the labels in the pathway graph.
#' @param conditions A character vector specifying the conditions to visualize. These should match
#'   the column names in enrichment_scores.
#' @param mapping A data.frame containing mapping information between pathway names and their
#'   IDs. This is typically obtained from the load_default_data or
#'   load_and_preprocess_gene_sets functions.
#' @param g An igraph object representing the pathway graph. Nodes should have a
#'   label attribute corresponding to pathway names.
#' @param gene_sets A named list where each element corresponds to a pathway and contains
#'   the associated genes. Names should match the pathway labels in the graph.
#' @param enrichment_limits Optional numeric vector of length two specifying the minimum and
#'   maximum enrichment scores for visualization. If NULL, limits are computed from
#'   enrichment_scores.
#' @param glyph_size Optional numeric vector specifying the width and height of glyphs in
#'   pixels. Default is c(80, 60).
#' @param res Optional numeric value specifying the resolution of glyph images in DPI.
#'   Default is 96.
#' @param cache_dir Optional character string specifying the directory to store cached glyphs.
#'   Default is "glyph_cache".
#' @param layout_options_plot Optional list specifying layout options for ggraph.
#'   Default is list(layout = "fr", type = "full", physics = FALSE).
#'
#' @return Launches a Shiny application for interactive pathway visualization.
#' @export
run_pathway_shiny_app <- function(
    enrichment_scores,         # Data frame: pathways x conditions
    conditions,                # Character vector: conditions to visualize
    mapping,                   # Data frame: pathway to exactSource mapping
    g,                         # igraph object: pathway graph
    gene_sets,                 # Named list: pathway -> genes
    enrichment_limits = NULL,  # Numeric vector: c(min, max) for enrichment
    glyph_size = c(80, 60),    # Numeric vector: width and height of glyphs in pixels
    res = 96,                  # Numeric: resolution of glyph images in DPI
    cache_dir = file.path(tempdir(),"glyph_cache"), # Character: directory to store cached glyphs
    layout_options_plot = list(
      layout = "layout_with_dh",
      physics = FALSE,
      smooth = T,
      type = "full",
      randomSeed = 42,
      options = NULL
    ) # List: Options for plot layout
) {

  # -----------------------------
  # 1. Define Shiny UI
  # -----------------------------

  ui <- fluidPage(
    titlePanel("Pathway Visualization Shiny App"),

    sidebarLayout(
      sidebarPanel(
        h4("Instructions"),
        p("Use the search box below to find specific pathways. Click on a top-level pathway glyph to view its sub-pathways. Use the 'Reset View' button to return to the top-level pathways."),

        # Dropdown for selecting comparisons
        selectInput(
          inputId = "selected_comparisons",
          label = "Select Conditions:",
          choices = conditions,
          selected = conditions,  # Default is all comparisons
          multiple = TRUE
        ),

        # Dropdown for selecting reference comparison
        selectInput(
          inputId = "reference_comparison",
          label = "Select Reference Condition:",
          choices = c("None", conditions),
          selected = "None",  # Default to "None"
          multiple = FALSE
        ),

        textInput("search_pathway", "Search Pathway:", value = ""),
        actionButton("reset", "Reset View"),
        br(),
        br(),

        h4("Legend"),
        p("Glyphs represent enrichment scores across conditions."),

        # Add the dynamic ggplot legend below the legend text
        plotOutput("legend_plot", height = "250px"),  # Adjust height as needed

        br(),
        br(),

        # Help button
        actionButton("help", "Help")
      ),

      mainPanel(
        visNetworkOutput("pathway_network", height = "90vh")
      )
    )
  )

  # -----------------------------
  # 2. Define Shiny Server
  # -----------------------------

  server <- function(input, output, session) {
    # Identify top-level pathways: pathways with no incoming edges
    top_level_pathways <- V(g)[degree(g, mode = "in") == 0]$name

    # Add a root node if not already present
    if (!"root" %in% V(g)$name) {
      g <<- add_vertices(g, 1, name = "root", label = "Root", color = "red")
    }

    # Connect the root node to all top-level pathways
    existing_edges <- igraph::as_data_frame(g, what = "edges")
    # Prepare edges from root to top-level pathways
    new_edges <- c()
    for (pathway in top_level_pathways) {
      # Check if the edge already exists to prevent duplication
      if (!any(existing_edges$from == "root" & existing_edges$to == pathway)) {
        new_edges <- c(new_edges, "root", pathway)
      }
    }
    if (length(new_edges) > 0) {
      # Add edges using 'name' attributes
      g <<- add_edges(g, new_edges)
    }

    # Refresh top-level pathways to include root
    top_level_pathways <- c("root", top_level_pathways)

    # Create a subgraph with only the root and top-level pathways
    initial_graph <- induced_subgraph(g, vids = top_level_pathways)

    # Reactive value to store the current graph state
    current_graph <- reactiveVal(initial_graph)

    # Reactive value to store glyph images
    glyph_images <- reactiveVal(list())

    # Reactive expressions for dropdown selections
    selected_comparisons <- reactive({
      req(input$selected_comparisons)
      # Exclude the reference comparison if it's not "None"
      if (input$reference_comparison != "None") {
        setdiff(input$selected_comparisons, input$reference_comparison)
      } else {
        input$selected_comparisons
      }
    })

    reference_comparison <- reactive({
      req(input$reference_comparison)
      input$reference_comparison
    })

    # Reactive expression to filter enrichment scores based on selected comparisons
    filtered_enrichment_scores <- reactive({
      comparisons <- selected_comparisons()

      # If no comparisons are selected, return NULL
      if (length(comparisons) == 0) {
        return(NULL)
      }

      # Subset the enrichment_scores based on selected comparisons
      enriched_subset <- enrichment_scores[, comparisons, drop = FALSE]

      # If a reference is selected and it's not "None", adjust the enrichment scores
      ref <- reference_comparison()
      if (!is.null(ref) && ref != "None" && ref %in% conditions) {
        # Access the reference scores from the original enrichment_scores
        ref_scores <- enrichment_scores[, ref, drop = FALSE]
        # Subtract the reference comparison scores from all selected comparisons
        enriched_subset <- sweep(enriched_subset, 1, ref_scores[, 1], FUN = "-")
      }

      return(enriched_subset)
    })

    # Reactive expression to determine new enrichment limits based on filtered scores
    new_enrichment_limits <- reactive({
      scores <- filtered_enrichment_scores()
      if (is.null(scores)) {
        return(NULL)
      }
      c(min(scores, na.rm = TRUE), max(scores, na.rm = TRUE))
    })

    # Pre-generate glyph images with caching and progress indication
    observe({
      withProgress(message = 'Generating Glyph Images...', value = 0, {
        g_current <- current_graph()
        data_hash <- digest::digest(list(enrichment_scores = filtered_enrichment_scores(),
                                         conditions = selected_comparisons()), algo = "md5")
        
        # Generate glyph images with caching based on filtered enrichment scores
        glyphs <- generate_glyph_images_cached(
          g = g_current,
          gene_sets = gene_sets,
          enrichment_scores = filtered_enrichment_scores(),
          conditions = selected_comparisons(),
          mapping = mapping,
          enrichment_limits = new_enrichment_limits(),
          glyph_size = glyph_size,
          res = res,
          cache_dir = file.path(cache_dir,data_hash),
          progress = Progress$new(),
          force_regenerate = FALSE  # Set to TRUE to force cache invalidation
        )
        glyph_images(glyphs)
      })
    })

    # Sanity Check: Ensure all pathway nodes have mappings
    observe({
      graph_names <- V(current_graph())$name
      mapping_names <- mapping$exactSource
      missing_mappings <- setdiff(graph_names, mapping_names)
      if (length(missing_mappings) > 0) {
        warning(paste("Missing mappings for:", paste(missing_mappings, collapse = ", ")))
      }
    })

    # Reactive expression to filter nodes based on search input
    filtered_nodes <- reactive({
      search_term <- tolower(trimws(input$search_pathway))
      if (search_term == "") {
        # No filtering; show all nodes
        V(current_graph())$name
      } else {
        # Filter nodes that contain the search term in their labels
        matched_nodes <- V(current_graph())[grepl(search_term, tolower(V(current_graph())$label))]
        matched_nodes$name
      }
    })

    # Render the ggplot legend
    output$legend_plot <- renderPlot({
      req(selected_comparisons(), new_enrichment_limits())

      # Generate the legend plot
      legend_plot <- generate_legend_plot(
        conditions = selected_comparisons(),
        enrichment_limits = new_enrichment_limits(),
        reference = reference_comparison()
      )

      # Check if the plot was generated successfully
      if (is.null(legend_plot)) {
        # Create an empty plot with a message
        ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "Legend could not be generated.", size = 5) +
          theme_void()
      } else {
        legend_plot
      }
    }, height = 250)  # Adjust height as needed

    # Render the visNetwork graph
    output$pathway_network <- renderVisNetwork({
      g_current <- current_graph()

      # Retrieve pre-generated glyph images
      glyphs <- glyph_images()

      # Prepare nodes data frame for visNetwork
      nodes <- data.frame(
        id = V(g_current)$name,
        label = V(g_current)$label,  # Pathway names as labels
        title = ifelse(
          V(g_current)$name == "root",
          "The Root node serves as a common parent for all top-level pathways.",
          paste(
            "Pathway:", V(g_current)$label, "<br>",
            "Genes:", vapply(V(g_current)$label, function(x) {
              if (!is.null(gene_sets[[x]])) {
                paste(as.character(gene_sets[[x]]), collapse = ", ")
              } else {
                "No genes available"
              }
            }, FUN.VALUE = character(1)),
            sep = ""
          )
        ),
        stringsAsFactors = FALSE
      )

      # Assign images to nodes if available
      nodes$image <- sapply(V(g_current)$label, function(x) {
        if (x == "Root") {
          ""  # Root node doesn't have an image
        } else {
          img <- glyphs[[x]]
          if (!is.null(img) && !is.na(img) && img != "") {
            img
          } else {
            ""  # Use default image or keep empty
          }
        }
      })

      # Define node shapes and sizes
      nodes$shape <- ifelse(nodes$image != "", "image",
                            ifelse(nodes$id == "root", "dot", "dot"))
      nodes$size <- ifelse(nodes$id == "root", 30,
                           ifelse(nodes$shape == "dot", 10, 15))

      # Define label properties for better visibility
      nodes$font.size <- ifelse(nodes$id == "root", 20, 14)  # Larger font for root
      nodes$font.face <- "Tahoma"  # Consistent font face

      # Highlight nodes based on search
      if (input$search_pathway != "") {
        nodes$color.background <- ifelse(nodes$id %in% filtered_nodes(), "orange", "lightblue")
      } else {
        nodes$color.background <- "lightblue"
      }

      # Prepare edges data frame for visNetwork
      edges <- igraph::as_data_frame(g_current, what = "edges")
      colnames(edges) <- c("from", "to")

      # Create visNetwork object with layout
      visNetwork(nodes, edges, height = "800px", width = "100%") %>%
        visNodes(
          shapeProperties = list(useImageSize = TRUE),
          shadow = FALSE,
          font = list(color = "black", size = nodes$font.size, face = nodes$font.face, bold = TRUE)
        ) %>%
        visEdges(arrows = "to") %>%  # Add arrows to edges for directionality
        visIgraphLayout(
          layout = layout_options_plot$layout,
          type = layout_options_plot$type,
          physics = layout_options_plot$physics
        ) %>%
        visOptions(
          highlightNearest = list(enabled = TRUE, degree = 1, hover = FALSE),
          nodesIdSelection = TRUE
        ) %>%
        visInteraction(
          navigationButtons = TRUE,
          zoomView = TRUE,
          dragView = TRUE,
          multiselect = FALSE
        ) %>%
        visEvents(
          click = "function(nodes) {
                    Shiny.setInputValue('pathway_network_selected', nodes.nodes, {priority: 'event'});
                  }"
        )
    })

    # Observe clicks on the network
    observeEvent(input$pathway_network_selected, {
      selected_node <- input$pathway_network_selected

      # Prevent triggering when no node is selected
      if (is.null(selected_node) || length(selected_node) == 0 || selected_node == "") return()

      selected_node <- selected_node[1]  # Take the first selected node

      # Prevent actions on the "root" node
      if (selected_node == "root") {
        showModal(modalDialog(
          title = "Root Node",
          "The Root node serves as a common parent for all top-level pathways.",
          easyClose = TRUE,
          footer = NULL
        ))
        return()
      }

      # Get the pathway name from the selected node using existing mapping
      selected_pathway <- mapping[which(mapping$exactSource == selected_node),"processed_name"]

      # Check if the selected node exists in the graph
      if (!(selected_node %in% V(g)$name)) {
        showModal(modalDialog(
          title = "Pathway Not Found",
          paste("The selected pathway", selected_pathway, "is not found in the graph."),
          easyClose = TRUE,
          footer = NULL
        ))
        return()
      }

      # Use the existing get_child_pathways function to find sub-pathways
      child_pathways <- get_child_pathways(
        parent_geneset_name = selected_pathway,
        mapping = mapping,
        g = g
      )

      if (length(child_pathways) == 0) {
        showModal(modalDialog(
          title = "No Sub-Pathways",
          paste("The selected pathway", selected_pathway, "has no sub-pathways."),
          easyClose = TRUE,
          footer = NULL
        ))
        return()
      }

      # Map child pathway names to exactSource IDs
      child_exactSource <- mapping$exactSource[match(child_pathways, mapping$processed_name)]

      # Remove any NAs resulting from unmatched pathways
      child_exactSource <- child_exactSource[!is.na(child_exactSource)]

      if (length(child_exactSource) == 0) {
        showModal(modalDialog(
          title = "No Valid Sub-Pathways",
          paste("The selected pathway", selected_pathway, "has no valid sub-pathways with matching exactSource IDs."),
          easyClose = TRUE,
          footer = NULL
        ))
        return()
      }

      # Induce subgraph with selected pathway and its children
      sub_g <- induced_subgraph(g, vids = c(selected_node, child_exactSource))

      # Update the current graph
      current_graph(sub_g)
    })

    # Handle the reset button to show the initial top-level pathways
    observeEvent(input$reset, {
      current_graph(initial_graph)
    })

    # Help Modal
    observeEvent(input$help, {
      showModal(modalDialog(
        title = "Help",
        "Use the search box to find specific pathways. Click on a top-level pathway glyph to view its sub-pathways. Use the 'Reset View' button to return to the top-level pathways.",
        easyClose = TRUE,
        footer = NULL
      ))
    })
  }

  # -----------------------------
  # 3. Launch the Shiny App
  # -----------------------------

  shinyApp(ui = ui, server = server)
}
