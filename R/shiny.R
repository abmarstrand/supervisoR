#' Launch Pathway Visualization Shiny Application
#'
#' This function launches a Shiny application for visualizing pathway enrichment scores
#' and their relationships within a pathway graph. Users can interactively explore
#' top-level pathways and their sub-pathways, view glyphs representing enrichment scores
#' across multiple conditions, and reset the view to the initial state.
#'
#' @importFrom shiny shinyApp renderUI fluidPage reactiveValues modalButton uiOutput tags HTML sidebarLayout sidebarPanel mainPanel selectInput selectizeInput updateSelectizeInput actionButton plotOutput titlePanel h4 p br reactiveVal reactive icon req observe withProgress Progress renderPlot showModal modalDialog observeEvent
#' @importFrom shinyWidgets searchInput updateSearchInput 
#' @importFrom visNetwork visIgraphLayout visNodes visEdges visOptions visEvents visNetworkOutput renderVisNetwork visInteraction visNetwork
#' @importFrom ggplot2 ggplotGrob annotate
#' @importFrom stringr str_to_upper str_trim
#' @importFrom dplyr arrange
#' @importFrom magrittr %>%
#' @importFrom cowplot plot_grid
#' @importFrom igraph degree V add_vertices induced_subgraph `V<-` add_edges drl_defaults as_data_frame
#' @importFrom utils head tail 
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
  mapping,                   # Data frame: pathway to exact_source mapping
  g,                         # igraph object: pathway graph
  gene_sets,                 # Named list: pathway -> genes
  enrichment_limits = NULL,  # Numeric vector: c(min, max) for enrichment
  glyph_size = c(80, 60),    # Numeric vector: width and height of glyphs in pixels
  res = 96,                  # Numeric: resolution of glyph images in DPI
  cache_dir = file.path(tempdir(), "glyph_cache"), # Character: directory to store cached glyphs
  layout_options_plot = list(
    layout = "layout_with_dh",
    physics = FALSE,
    smooth = TRUE,
    type = "full",
    randomSeed = 42,
    options = NULL
  ) # List: Options for plot layout
) {
  # Prepare list of all pathways for search
  all_pathways <- unique(mapping$processed_name)
  names(all_pathways) <- all_pathways
  gene_list <- unlist(gene_sets)
  pathway_list <- rep(names(gene_sets), times = lengths(gene_sets))
  genes_upper <- toupper(gene_list)  # Convert genes to uppercase for consistency
  # Now create a mapping: genes -> pathways
  gene_to_pathways <- split(pathway_list, genes_upper)
  all_genes <- sort(unique(names(gene_to_pathways)))  # all uppercase gene names
  enriched_pathways <- all_pathways[all_pathways %in% rownames(enrichment_scores)]

  # 1. Define Shiny UI
  # -----------------------------

  ui <- fluidPage(
    titlePanel("supervisoR Shiny App"),

    sidebarLayout(
      sidebarPanel(
        h4("Instructions"),
        p("Use the search box below to find specific pathways. Click on a top-level pathway glyph to view its sub-pathways. Use the 'Reset View' button to return to the top-level pathways. Use the 'Back' button to return to the previous view."),
        # Dropdown for selecting comparisons
        selectInput(
          inputId = "selected_comparisons",
          label = "Select Data:",
          choices = conditions,
          selected = conditions,  # Default is all comparisons
          multiple = TRUE
        ),
        # Dropdown for selecting reference comparison
        selectInput(
          inputId = "reference_comparison",
          label = "Select Data as Reference:",
          choices = c("None", conditions),
          selected = "None",  # Default to "None"
          multiple = FALSE
        ),
        # Dropdown for finding pathways
        selectizeInput(
          inputId = "search_pathway",
          label = "Search Pathway:",
          choices = NULL,   # Choices will be populated in the server
          selected = "",  # Ensure no default selection
          multiple = FALSE,
          options = list(
            placeholder = 'Type to search pathways',
            maxOptions = 5000  # Adjust as needed
          )
        ),
        # Dropdown for finding pathways which contain a given gene
        selectizeInput(
          inputId = "gene_input",
          label = "Search Gene:",
          choices = NULL,             # Provide all genes to select from
          selected = "",
          multiple = FALSE,
          options = list(
            placeholder = 'Type to search for a gene',
            maxOptions = 5000,              # Adjust as needed
            create = FALSE,                 # Disallow creating new entries
            highlight = TRUE
          )
        ),
        actionButton("gene_input_search", "Search Gene"),
        uiOutput("pathway_select_ui"), 
        actionButton("reset", "Reset View"),
        actionButton("back", "Back"),
        br(), br(),

        h4("Legend"),
        p("Glyphs represent enrichment scores across conditions."),

        # Add the dynamic ggplot legend below the legend text
        plotOutput("legend_plot", height = "250px"),  # Adjust height as needed

        br(), br(),
        
        # Help button
        actionButton("help", "Help")
      ),

      mainPanel(
        visNetworkOutput("pathway_network", height = "90vh")
      )
    ),
    tags$script(HTML("
      $(window).on('resize', function() {
        var width = $(window).width();
        var cutoff = Math.floor(width / 300); // Adjust based on desired width per gene
        Shiny.setInputValue('cutoff_update', cutoff, {priority: 'event'});
      });
      $(document).ready(function() {
        $(window).trigger('resize'); // Trigger on load to set initial value
      });
    ")),
    uiOutput("gene_modal")
  )

  # -----------------------------
  # 2. Define Shiny Server
  # -----------------------------

  server <- function(input, output, session) {
    # Add a reactive cutoff value based on window size
    cutoff <- reactiveVal(4) # Default to 4
    last_selected_node <- reactiveVal(NULL)
    
    
    # Identify top-level pathways: pathways with no incoming edges
    top_level_pathways <- V(g)[degree(g, mode = "in") == 0]$name
    
    # Add a root node if not already present
    if (!"root" %in% V(g)$name) {
      g <<- add_vertices(g, 1, name = "root", label = "Root", color = "red")
    }

    # Connect the root node to all top-level pathways
    existing_edges <- as_data_frame(g, what = "edges")
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
      g <<- add_edges(g, new_edges, relation = NA_character_)
    }

    # Refresh top-level pathways to include root
    top_level_pathways <- c("root", top_level_pathways)
    
    # Create a subgraph with only the root and top-level pathways
    initial_graph <- induced_subgraph(g, vids = top_level_pathways)
    
    # Initialize history stack
    history <- reactiveValues(stack = list(), current = initial_graph)

    # Reactive value to store the current graph state
    current_graph <- reactiveVal(initial_graph)

    # Reactive value to store glyph images
    glyph_images <- reactiveVal(list())

    # Reactive expressions for dropdown selections
    selected_comparisons <- reactive({
      req(input$selected_comparisons)
      comps <- input$selected_comparisons
      if (input$reference_comparison != "None") {
        comps <- setdiff(comps, input$reference_comparison)
      }
      comps
    })

    reference_comparison <- reactive({
      req(input$reference_comparison)
      input$reference_comparison
    })

    # Reactive expression to filter enrichment scores based on selected comparisons
    filtered_enrichment_scores <- reactive({
      comps <- selected_comparisons()

      # If no comparisons are selected, return NULL
      if (length(comps) == 0) return(NULL)

      # Subset the enrichment_scores based on selected comparisons
      enriched_subset <- enrichment_scores[, comps, drop = FALSE]

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
      if (is.null(scores)) return(NULL)
      c(min(scores, na.rm = TRUE), max(scores, na.rm = TRUE))
    })

    # Pre-generate glyph images with caching and progress indication
    observe({
      req(filtered_enrichment_scores())
      withProgress(message = "Generating Glyph Images...", value = 0, {
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
          cache_dir = file.path(cache_dir, data_hash),
          progress = Progress$new(),
          force_regenerate = FALSE  # Set to TRUE to force cache invalidation
        )
        glyph_images(glyphs)
      })
    })

    # Sanity Check: Ensure all pathway nodes have mappings
    observe({
      graph_names <- V(current_graph())$name
      mapping_names <- mapping$exact_source
      missing_mappings <- setdiff(graph_names, mapping_names)
      if (length(missing_mappings) > 0) {
        warning(paste("Missing mappings for:", paste(missing_mappings, collapse = ", ")))
      }
    })

    # Reactive expression to filter nodes based on search input
    filtered_nodes <- reactive({
      search_term_input <- input$search_pathway
      
      # If the input is NULL or has length 0, set search_term to an empty string
      if (is.null(search_term_input) || length(search_term_input) == 0) {
        search_term <- ""
      } else {
        # Trim whitespace
        search_term <- trimws(search_term_input)
        # If trimming leads to zero-length, default to empty string
        if (length(search_term) == 0) search_term <- ""
      }
      
      # Convert to lowercase
      search_term <- tolower(search_term)
      
      labels <- V(current_graph())$label
      labels[is.na(labels)] <- ""
      
      # Now it's safe to compare search_term to ""
      if (search_term == "") {
        V(current_graph())$name
      } else {
        matched_nodes <- V(current_graph())[grepl(search_term, tolower(labels))]
        matched_nodes$name
      }
    })
    
    # Update the selectizeInput choices with server-side processing
    observe({
      updateSelectizeInput(session, 
                           choices = all_genes, 
                           "gene_input",
                           selected = "",
                           server = TRUE)
      updateSelectizeInput(
        session,
        "search_pathway",
        choices = enriched_pathways,
        selected = "",  # Ensure no selection is made
        server = TRUE             # Enable server-side selectize
      )
    })
    
    # Handle search input selection
    observeEvent(input$search_pathway, {
      selected_pathway <- input$search_pathway
      if (is.null(selected_pathway) || selected_pathway == "") return()
      
      # Get the exact_source ID of the selected pathway
      selected_node <- mapping$exact_source[match(selected_pathway, mapping$processed_name)]
      if (is.na(selected_node)) {
        showModal(modalDialog(
          title = "Pathway Not Found",
          paste("The selected pathway", selected_pathway, "is not found in the mapping."),
          easyClose = TRUE
        ))
        return()
      }
      
      # Use the existing get_child_pathways function to find sub-pathways
      child_pathways <- get_child_pathways(
        parent_geneset_name = selected_pathway,
        mapping = mapping,
        g = g
      )
      
      # Map child pathway names to exact_source IDs
      child_exact_source <- mapping$exact_source[match(child_pathways, mapping$processed_name)]
      child_exact_source <- child_exact_source[!is.na(child_exact_source)]
      
      # Include the selected node itself
      vids <- c(selected_node, child_exact_source)
      
      if (length(vids) == 0) {
        showModal(modalDialog(
          title = "No Sub-Pathways",
          paste("The selected pathway", selected_pathway, "has no sub-pathways."),
          easyClose = TRUE
        ))
        return()
      }
      
      # Induce subgraph with selected pathway and its children
      sub_g <- induced_subgraph(g, vids = vids)
      
      # Add previous view to history stack
      history$stack <- c(history$stack, list(current_graph()))
      
      # Update the current graph
      current_graph(sub_g)
      last_selected_node(NULL)
      
      # Clear the search input
      updateSelectizeInput(session, "search_pathway", selected = character(0))
    })
    
    # Reactive value to store pathways for the entered gene
    pathways_for_gene <- reactiveVal(NULL)
    
    # Handle the gene search button click
    # Handle the gene search button click
    observeEvent(input$gene_input_search, {
      ### CHANGES START
      # Convert gene to uppercase and trim whitespace
      gene_input_value <- input$gene_input
      if (is.null(gene_input_value) || length(gene_input_value) == 0) {
        gene <- ""
      } else {
        gene <- toupper(trimws(gene_input_value))
        if (length(gene) == 0) gene <- ""
      }
      
      if (gene == "") return()
      
      # Look up pathways containing the gene
      pathways <- gene_to_pathways[[gene]]
      
      if (is.null(pathways) || length(pathways) == 0) {
        showModal(modalDialog(
          title = "Gene Not Found",
          paste("No pathways containing the gene", gene, "were found."),
          easyClose = TRUE
        ))
        return()
      }
      
      # Filter pathways so that only those present in the original graph are shown
      # Map processed_name to exact_source, then check if exact_source is in V(g)$name
      valid_nodes <- V(g)$name
      exact_sources <- mapping$exact_source[match(pathways, mapping$processed_name)]
      valid_pathways <- pathways[exact_sources %in% valid_nodes]
      
      if (length(valid_pathways) == 0) {
        showModal(modalDialog(
          title = "No Valid Pathways",
          paste("No pathways containing the gene", gene, "exist in the original graph."),
          easyClose = TRUE
        ))
        return()
      }
      
      # Store only valid pathways
      pathways_for_gene(valid_pathways)
      
      # Update the UI to show the pathway selection input
      output$pathway_select_ui <- renderUI({
        selectInput(
          inputId = "selected_pathway_from_gene",
          label = "Pathways containing selected gene:",
          choices = c(" " = "", valid_pathways),
          selected = " "
        )
      })
      ### CHANGES END
    })
    
    # Handle the pathway selection from gene search
    observeEvent(input$selected_pathway_from_gene, {
      selected_pathway <- input$selected_pathway_from_gene
      if (is.null(selected_pathway) || selected_pathway == "") return()
      
      # Get the exact_source ID of the selected pathway
      selected_node <- mapping$exact_source[match(selected_pathway, mapping$processed_name)]
      if (is.na(selected_node)) {
        showModal(modalDialog(
          title = "Pathway Not Found",
          paste("The selected pathway", selected_pathway, "is not found in the mapping."),
          easyClose = TRUE
        ))
        return()
      }
      
      # Use get_child_pathways to find sub-pathways
      child_pathways <- get_child_pathways(
        parent_geneset_name = selected_pathway,
        mapping = mapping,
        g = g
      )
      
      child_exact_source <- mapping$exact_source[match(child_pathways, mapping$processed_name)]
      child_exact_source <- child_exact_source[!is.na(child_exact_source)]
      
      vids <- c(selected_node, child_exact_source)
      
      if (length(vids) == 0) {
        showModal(modalDialog(
          title = "No Sub-Pathways",
          paste("The selected pathway", selected_pathway, "has no sub-pathways."),
          easyClose = TRUE,
          footer = NULL
        ))
        return()
      }
      
      sub_g <- induced_subgraph(g, vids = vids)
      
      history$stack <- c(history$stack, list(current_graph()))
      current_graph(sub_g)
      last_selected_node(NULL)
      
      # Clear the gene input and pathway selection UI
      updateSelectizeInput(session, "gene_input", selected = "", server = TRUE)
      pathways_for_gene(NULL)
      output$pathway_select_ui <- renderUI(NULL)
    })
  
    # Render the ggplot legend
    output$legend_plot <- renderPlot({
      req(selected_comparisons(), new_enrichment_limits())
      
      g_current <- current_graph()
      edges_df <- as_data_frame(g_current, what = "edges")
      # Identify which relations are present (excluding NA)
      
      relation_colors <- c(
        "is_a" = "black",
        "part_of" = "darkcyan",
        "regulates" = "orange",
        "positively_regulates" = "steelblue",
        "negatively_regulates" = "salmon"
      )
      
      relation_linetypes <- c(
        "is_a" = "solid",
        "part_of" = "dashed",
        "regulates" = "dotted",
        "positively_regulates" = "dotdash",
        "negatively_regulates" = "longdash"
      )
      
      if ("relation" %in% colnames(edges_df)) {
        relations <- unique(edges_df$relation)
      } else {
        relations <- character(0)
      }
      
      relations_present <- relations[!is.na(relations)]
      relations_present <- relations_present[relations_present %in% names(relation_colors)]
      
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
        if (length(relations_present) > 0) {
          relation_legend <- generate_relation_legend(
            relations_present = relations_present,
            relation_colors = relation_colors,
            relation_linetypes = relation_linetypes
          )
          # Combine the two legends vertically
          combined_plot <- plot_grid(legend_plot, relation_legend, ncol = 1, rel_heights = c(3,1))
          combined_plot
        } else {
          legend_plot
          
        }
      }
    }, height = 250)  # Adjust height as needed

    # Render the visNetwork graph
    output$pathway_network <- renderVisNetwork({
      g_current <- current_graph()

      # Retrieve pre-generated glyph images
      glyphs <- glyph_images()
      
      # Safely handle search_value
      search_value <- input$search_pathway
      if (is.null(search_value) || length(search_value) == 0) {
        search_value <- ""
      } else {
        search_value <- trimws(search_value)
        if (length(search_value) == 0) search_value <- ""
      }
      

      # Prepare nodes data frame for visNetwork
      nodes <- data.frame(
        id = V(g_current)$name,
        label = V(g_current)$label,  # Pathway names as labels
        title = ifelse(
          V(g_current)$label == "root",
          "The Root node serves as a common parent for all top-level pathways.",
          paste(
            "Pathway:", V(g_current)$label, "<br>",
            "Genes:", sapply(V(g_current)$label, function(x) {
              if (!is.null(gene_sets[[x]])) {
                gene_list <- gene_sets[[x]]
                if (length(gene_list) > 10) {
                  # Truncate the gene list and add a clickable link
                  paste(
                    paste(gene_list[1:cutoff()], collapse = ", "),
                    "... <a href='#' class='show-full-genes' data-node='", x, "'>Show more</a>"
                  )
                } else {
                  paste(gene_list, collapse = ", ")
                }
              } else {
                "No genes available"
              }
            }),
            sep = ""
          )
        ),
        stringsAsFactors = FALSE
      )

      # Assign images to nodes if available
      nodes$image <- sapply(V(g_current)$label, function(x) {
        img <- glyphs[[x]]
        if (!is.null(img) && !is.na(img) && img != "") {
          img
        } else {
          ""
        }
      })

      # Define node shapes and sizes
      nodes$shape <- ifelse(nodes$image != "", "image",
                            ifelse(nodes$id == "root", "dot", "dot"))
      nodes$size <- ifelse(nodes$id == "root", 30,
                           ifelse(nodes$shape == "dot", 10, 15))

      # Define label properties for better visibility
      nodes$font.size <- ifelse(nodes$id == "root", 20, 14)  # Larger font for root
      nodes$font.face <- "Arial"  # Consistent font face

      # Highlight nodes based on search
      if (is.null(input$search_pathway) || input$search_pathway == "") {
        # No search term
        nodes$color.background <- "lightblue"
      } else {
        # Highlight filtered nodes
        nodes$color.background <- ifelse(nodes$id %in% filtered_nodes(), "orange", "lightblue")
      }
      nodes <- nodes %>% arrange(label)

      # Prepare edges data frame for visNetwork
      # Prepare edges data frame for visNetwork
      edges <- as_data_frame(g_current, what = "edges")
      
      # Define line types mapping from relations
      relation_values <- c(
        "is_a" = "solid",
        "part_of" = "dashed",
        "regulates" = "dotted",
        "positively_regulates" = "dotdash",
        "negatively_regulates" = "longdash"
      )
      
      # Define dash patterns to simulate these line types.
      # The array values represent lengths of dashes and gaps.
      line_patterns <- list(
        "solid" = FALSE,        # No dashes, a solid line
        "dashed" = c(5, 5),     # Even dash-gap pattern
        "dotted" = c(1, 5),     # Short dash, longer gap
        "dotdash" = c(1, 5, 5, 5), # Dot then dash pattern
        "longdash" = c(10, 5)   # Longer dash than 'dashed'
      )
      
      # Define colors for each relation
      relation_colors <- c(
        "is_a" = "black",
        "part_of" = "darkcyan",
        "regulates" = "orange",
        "positively_regulates" = "steelblue",
        "negatively_regulates" = "salmon"
      )
      
      # In your renderVisNetwork block:
      edges <- as_data_frame(g_current, what = "edges")
      
      # Determine the line type based on relation:
      # If NA or "is_a", then "solid"
      edges$linetype <- ifelse(is.na(edges$relation) | edges$relation == "is_a",
                               "solid",
                               relation_values[edges$relation])
      #Make dotted lines a bit thicker
      edges$width <- ifelse(is.na(edges$relation) | edges$relation == "is_a",
                               1,
                               2)
      
      # Assign dash patterns based on the determined linetype
      edges$dashes <- sapply(edges$linetype, function(lt) line_patterns[[lt]])
      
      # Assign colors:
      # For NA relation, treat it like "is_a" (solid)
      edges$color <- ifelse(is.na(edges$relation),
                            relation_colors["is_a"],
                            relation_colors[edges$relation])
      
      # Now keep only needed columns for visNetwork
      edges <- edges[, c("from", "to", "color", "dashes", "width"), drop = FALSE]
      
      arrow_dir <- ifelse(grepl("GO",V(g)[1]$name), "from", "to")

      # Create visNetwork object with layout
      visNetwork(nodes, edges, height = "800px", width = "100%") %>%
        visEdges(arrows = arrow_dir,
                 smooth = c(enabled = T,
                            type = "curvedCCW",
                            roundness = 0.8)) %>%  # Add arrows to edges for directionality
        visNodes(
          shapeProperties = list(useImageSize = TRUE),
          shadow = FALSE,
          font = list(color = "black", size = nodes$font.size, face = nodes$font.face, bold = TRUE)
        ) %>%
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
          multiselect = FALSE,
          hoverConnectedEdges = TRUE,
          tooltipStay = 500
        ) %>%
        visEvents(
          click = "function(nodes) {
                    Shiny.setInputValue('pathway_network_selected', nodes.nodes, {priority: 'event'});
                  }",
          afterDrawing = "function() {
                    // Delegate event handling to the document
                    document.addEventListener('click', function(e) {
                      if (e.target && e.target.classList.contains('show-full-genes')) {
                        e.preventDefault(); // Prevent default anchor behavior
                        var nodeName = e.target.getAttribute('data-node').trim(); // Get the data-node attribute
                        Shiny.setInputValue('tooltip_click', nodeName, {priority: 'event'}); // Send to Shiny
                      }
                    });
                  }"
        )
    })
    
    
    # Observe clicks on the network
    observeEvent(input$pathway_network_selected, {
      selected_node <- input$pathway_network_selected

      # Prevent triggering when no node is selected
      if (is.null(selected_node) || length(selected_node) == 0 || selected_node == "") return()

      selected_node <- selected_node[1]  # Take the first selected node
      
      if (!is.null(last_selected_node()) && selected_node == last_selected_node()) return()
      
      last_selected_node(selected_node)

      # Prevent actions on the "root" node
      if (selected_node == "root") {
        showModal(modalDialog(
          title = "Root Node",
          "The Root node serves as a common parent for all top-level pathways.",
          easyClose = TRUE
        ))
        return()
      }

      # Get the pathway name from the selected node using existing mapping
      selected_pathway <- mapping[which(mapping$exact_source == selected_node), "processed_name"]

      # Check if the selected node exists in the graph
      if (!(selected_node %in% V(g)$name)) {
        showModal(modalDialog(
          title = "Pathway Not Found",
          paste("The selected pathway", selected_pathway, "is not found in the graph."),
          easyClose = TRUE
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

      # Map child pathway names to exact_source IDs
      child_exact_source <- mapping$exact_source[match(child_pathways, mapping$processed_name)]

      # Remove any NAs resulting from unmatched pathways
      child_exact_source <- child_exact_source[!is.na(child_exact_source)]

      if (length(child_exact_source) == 0) {
        showModal(modalDialog(
          title = "No Valid Sub-Pathways",
          paste("The selected pathway", selected_pathway, "has no valid sub-pathways with matching exact_source IDs."),
          easyClose = TRUE
        ))
        return()
      }

      # Induce subgraph with selected pathway and its children
      sub_g <- induced_subgraph(g, vids = c(selected_node, child_exact_source))
      
      # Add previous view to history stack
      history$stack <- c(history$stack, list(current_graph()))
      
      # Update the current graph
      current_graph(sub_g)
      last_selected_node(NULL)
    }, ignoreInit = TRUE)
    
    # Update cutoff based on window size
    session$onFlushed(function() {
      session$sendCustomMessage(type = "updateCutoff", message = NULL)
    })
    
    observeEvent(input$back, {
      if (length(history$stack) > 0) {
        # Pop the last graph from the stack
        last_graph <- tail(history$stack, 2)[[1]]
        history$stack <- head(history$stack, -2)
        current_graph(last_graph)
        # Reset last_selected_node
        last_selected_node(NULL)
      } else {
        # Stack is empty, can't go back further
        showModal(modalDialog(
          title = "No Previous View",
          "There is no previous view to go back to.",
          easyClose = TRUE
        ))
      }
    })
    
    observeEvent(input$cutoff_update, {
      cutoff(input$cutoff_update)
    })
    
    observeEvent(input$tooltip_click, {
      node_label <- input$tooltip_click
      if (!is.null(node_label)) {
        genes <- gene_sets[[node_label]]
        showModal(
          modalDialog(
            title = paste("All genes in", node_label),
            if (!is.null(genes)) {
              paste(genes, collapse = ", ")
            } else {
              "No genes available for this pathway."
            },
            easyClose = TRUE,
            footer = modalButton("Close")
          )
        )
      }
    })

    # Handle the reset button to show the initial top-level pathways
    observeEvent(input$reset, {
      history$stack <- list()  # Clear the history stack
      current_graph(initial_graph)
      last_selected_node(NULL) # Reset last_selected_node
    })

    # Help Modal
    observeEvent(input$help, {
      showModal(modalDialog(
        title = "Help",
        "Use the search box to find specific pathways. Click on a top-level pathway glyph to view its sub-pathways.
        Use the 'Reset View' button to return to the top-level pathways.
        Use the 'Back' button to returnt to the previous view.",
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
