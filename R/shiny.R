#' Launch Pathway Visualization Shiny Application
#'
#' @importFrom shiny shinyApp renderUI fluidPage reactiveValues modalButton uiOutput tags HTML sidebarLayout sidebarPanel mainPanel selectInput selectizeInput updateSelectizeInput actionButton plotOutput titlePanel h4 p br reactiveVal reactive icon req observe withProgress Progress renderPlot showModal modalDialog observeEvent
#' @importFrom shinyWidgets searchInput updateSearchInput
#' @importFrom visNetwork visIgraphLayout visNodes visEdges visOptions visEvents visNetworkOutput renderVisNetwork visInteraction visNetwork
#' @importFrom ggplot2 ggplotGrob annotate
#' @importFrom dplyr arrange
#' @importFrom magrittr %>%
#' @importFrom cowplot plot_grid
#' @importFrom igraph degree V add_vertices induced_subgraph add_edges as_data_frame
#' @importFrom utils head tail
#' @param enrichment_scores A data.frame with pathways as rows and conditions as columns.
#' @param conditions A character vector specifying the conditions to visualize.
#' @param pathway_column Character, integer or NULL. Defaults to 1L.
#' @param species Species string. Default "Homo sapiens".
#' @param database Database string. Default "reactome".
#' @param mapping Optional mapping data frame.
#' @param g Optional igraph object.
#' @param gene_sets Optional named list.
#' @param enrichment_limits Optional numeric vector of length two.
#' @param glyph_size Numeric vector. Default c(80, 60).
#' @param res Numeric resolution. Default 96.
#' @param cache_dir Character cache directory path.
#' @param layout_options_plot List of layout options.
#' @param hide_nodes_without_enrichment Logical. Default TRUE.
#' @param lollipop_plot Logical. Default FALSE.
#' @param lollipop_colors Optional named list.
#' @return Launches a Shiny application.
#' @export
run_pathway_shiny_app <- function(
    enrichment_scores, conditions,
    pathway_column = 1L, species = "Homo sapiens", database = "reactome",
    mapping = NULL, g = NULL, gene_sets = NULL,
    enrichment_limits = NULL, glyph_size = c(80, 60), res = 96,
    lollipop_plot = FALSE, lollipop_colors = NULL,
    cache_dir = file.path(tempdir(), "glyph_cache"),
    hide_nodes_without_enrichment = TRUE,
    layout_options_plot = list(layout = "layout_with_dh", physics = FALSE, smooth = TRUE,
                               type = "full", randomSeed = 42, options = NULL)
) {
  
  if (is.null(mapping) || is.null(g) || is.null(gene_sets)) {
    data <- load_default_data(database = database, species = species)
    gs <- data$gene_sets
    rel <- data$relations
    trans <- data.frame(gene_sets_name = data$mapping$processed_name,
                        relation_name = data$mapping$exact_source, stringsAsFactors = FALSE)
    data <- load_and_preprocess_gene_sets(gene_sets = gs, pathways_relation = rel, translation_layer = trans)
    if (is.null(mapping)) mapping <- data$mapping
    if (is.null(g)) g <- data$pathway_graph
    if (is.null(gene_sets)) gene_sets <- data$gene_sets
  }
  
  if (!is.null(pathway_column)) {
    if (is.character(enrichment_scores[, pathway_column])) {
      rownames(enrichment_scores) <- enrichment_scores[, pathway_column]
      enrichment_scores <- enrichment_scores[, -pathway_column, drop = FALSE]
    } else {
      stop("Pathway column does not contain characters.")
    }
  }
  
  if (hide_nodes_without_enrichment) {
    node_labels <- V(g)$label
    node_labels[is.na(node_labels)] <- ""
    valid_labels <- node_labels[node_labels %in% rownames(enrichment_scores)]
    if (length(valid_labels) > 0) {
      scores <- enrichment_scores[valid_labels, conditions, drop = FALSE]
      missing_scores <- rowSums(is.na(scores)) == length(conditions)
      missing_scores[is.na(missing_scores)] <- TRUE
      keep_labels <- names(missing_scores)[!missing_scores]
      keep_nodes <- which(
        (node_labels %in% keep_labels) |
          (V(g)$name == "root") |
          grepl("^ellipsis_for_", V(g)$name)
      )
      g <- induced_subgraph(g, vids = keep_nodes)
    } else {
      keep_nodes <- which(V(g)$name == "root" | grepl("^ellipsis_for_", V(g)$name))
      g <- induced_subgraph(g, vids = keep_nodes)
    }
  }
  
  all_pathways <- unique(mapping$processed_name)
  names(all_pathways) <- all_pathways
  gene_list <- unlist(gene_sets)
  pathway_list <- rep(names(gene_sets), times = lengths(gene_sets))
  genes_upper <- toupper(gene_list)
  gene_to_pathways <- split(pathway_list, genes_upper)
  all_genes <- sort(unique(names(gene_to_pathways)))
  enriched_pathways <- all_pathways[all_pathways %in% rownames(enrichment_scores)]
  
  ui <- fluidPage(
    titlePanel("supervisoR Shiny App"),
    sidebarLayout(
      sidebarPanel(
        h4("Instructions"),
        p("Use the search box below to find specific pathways. Click on a top-level pathway glyph to view its sub-pathways. Use the 'Reset View' button to return to the top-level pathways. Use the 'Back' button to return to the previous view."),
        selectInput("selected_comparisons", "Select Data:", choices = conditions, selected = conditions, multiple = TRUE),
        selectInput("reference_comparison", "Select Data as Reference:", choices = c("None", conditions), selected = "None", multiple = FALSE),
        selectizeInput("search_pathway", "Search Pathway:", choices = NULL, selected = "", multiple = FALSE,
                              options = list(placeholder = "Type to search pathways", maxOptions = 5000)),
        selectizeInput("gene_input", "Search Gene:", choices = NULL, selected = "", multiple = FALSE,
                              options = list(placeholder = "Type to search for a gene", maxOptions = 5000, create = FALSE, highlight = TRUE)),
        actionButton("gene_input_search", "Search Gene"),
        uiOutput("pathway_select_ui"),
        actionButton("reset", "Reset View"),
        actionButton("back", "Back"),
        br(), br(),
        h4("Legend"),
        p("Glyphs represent enrichment scores across conditions."),
        plotOutput("legend_plot", height = "250px"),
        br(), br(),
        actionButton("help", "Help")
      ),
      mainPanel(
        visNetworkOutput("pathway_network", height = "90vh")
      )
    ),
    tags$script(HTML("
      $(window).on('resize', function() {
        var width = $(window).width();
        var cutoff = Math.floor(width / 300);
        Shiny.setInputValue('cutoff_update', cutoff, {priority: 'event'});
      });
      $(document).ready(function() { $(window).trigger('resize'); });
    ")),
    uiOutput("gene_modal")
  )
  
  server <- function(input, output, session) {
    cutoff <- reactiveVal(4)
    last_selected_node <- reactiveVal(NULL)
    
    g_local <- g
    top_level_pathways <- V(g_local)[degree(g_local, mode = "in") == 0]$name
    
    if (!"root" %in% V(g_local)$name) {
      g_local <- add_vertices(g_local, 1, name = "root", label = "Root", color = "red")
    }
    existing_edges <- as_data_frame(g_local, what = "edges")
    new_edges <- c()
    for (pathway in top_level_pathways) {
      if (!any(existing_edges$from == "root" & existing_edges$to == pathway)) {
        new_edges <- c(new_edges, "root", pathway)
      }
    }
    if (length(new_edges) > 0) {
      g_local <- add_edges(g_local, new_edges, relation = NA_character_)
    }
    top_level_pathways <- c("root", top_level_pathways)
    initial_graph <- induced_subgraph(g_local, vids = top_level_pathways)
    
    history <- reactiveValues(stack = list(), current = initial_graph)
    current_graph <- reactiveVal(initial_graph)
    glyph_images <- reactiveVal(list())
    
    selected_comparisons <- reactive({
      req(input$selected_comparisons)
      comps <- input$selected_comparisons
      if (input$reference_comparison != "None") comps <- setdiff(comps, input$reference_comparison)
      comps
    })
    
    reference_comparison <- reactive({ req(input$reference_comparison); input$reference_comparison })
    
    filtered_enrichment_scores <- reactive({
      comps <- selected_comparisons()
      if (length(comps) == 0) return(NULL)
      enriched_subset <- enrichment_scores[, comps, drop = FALSE]
      ref <- reference_comparison()
      if (!is.null(ref) && ref != "None" && ref %in% conditions) {
        ref_scores <- enrichment_scores[, ref, drop = FALSE]
        enriched_subset <- sweep(enriched_subset, 1, ref_scores[, 1], FUN = "-")
      }
      return(enriched_subset)
    })
    
    new_enrichment_limits <- reactive({
      scores <- filtered_enrichment_scores()
      if (is.null(scores)) return(NULL)
      c(min(scores, na.rm = TRUE), max(scores, na.rm = TRUE))
    })
    
    observe({
      req(filtered_enrichment_scores())
      withProgress(message = "Generating Glyph Images...", value = 0, {
        g_current <- current_graph()
        data_hash <- digest(list(enrichment_scores = filtered_enrichment_scores(),
                                         conditions = selected_comparisons()), algo = "spookyhash")
        glyphs <- generate_glyph_images_cached(
          g = g_current, gene_sets = gene_sets,
          enrichment_scores = filtered_enrichment_scores(), conditions = selected_comparisons(),
          mapping = mapping, enrichment_limits = new_enrichment_limits(),
          glyph_size = glyph_size, res = res,
          cache_dir = file.path(cache_dir, data_hash),
          progress = Progress$new(), force_regenerate = FALSE,
          lollipop_plot = lollipop_plot, lollipop_colors = lollipop_colors
        )
        glyph_images(glyphs)
      })
    })
    
    observe({
      graph_names <- V(current_graph())$name
      missing_mappings <- setdiff(graph_names, mapping$exact_source)
      if (length(missing_mappings) > 0) warning(paste("Missing mappings for:", paste(missing_mappings, collapse = ", ")))
    })
    
    filtered_nodes <- reactive({
      search_term_input <- input$search_pathway
      search_term <- if (is.null(search_term_input) || length(search_term_input) == 0) "" else trimws(search_term_input)
      if (length(search_term) == 0) search_term <- ""
      search_term <- tolower(search_term)
      labels <- V(current_graph())$label
      labels[is.na(labels)] <- ""
      if (search_term == "") V(current_graph())$name
      else V(current_graph())[grepl(search_term, tolower(labels))]$name
    })
    
    observe({
      updateSelectizeInput(session, choices = all_genes, "gene_input", selected = "", server = TRUE)
      updateSelectizeInput(session, "search_pathway", choices = enriched_pathways, selected = "", server = TRUE)
    })
    
    observeEvent(input$search_pathway, {
      selected_pathway <- input$search_pathway
      if (is.null(selected_pathway) || selected_pathway == "") return()
      selected_node <- mapping$exact_source[match(selected_pathway, mapping$processed_name)]
      if (is.na(selected_node)) { showModal(modalDialog(title = "Pathway Not Found", paste("Not found:", selected_pathway), easyClose = TRUE)); return() }
      child_pathways <- get_child_pathways(parent_geneset_name = selected_pathway, mapping = mapping, g = g_local)
      child_exact_source <- mapping$exact_source[match(child_pathways, mapping$processed_name)]
      child_exact_source <- child_exact_source[!is.na(child_exact_source)]
      vids <- c(selected_node, child_exact_source)
      if (length(vids) == 0) { showModal(modalDialog(title = "No Sub-Pathways", paste(selected_pathway, "has no sub-pathways."), easyClose = TRUE)); return() }
      sub_g <- induced_subgraph(g_local, vids = vids)
      history$stack <- c(history$stack, list(current_graph()))
      current_graph(sub_g)
      last_selected_node(NULL)
      updateSelectizeInput(session, "search_pathway", selected = character(0))
    })
    
    pathways_for_gene <- reactiveVal(NULL)
    
    observeEvent(input$gene_input_search, {
      gene_input_value <- input$gene_input
      gene <- if (is.null(gene_input_value) || length(gene_input_value) == 0) "" else toupper(trimws(gene_input_value))
      if (length(gene) == 0 || gene == "") return()
      pathways <- gene_to_pathways[[gene]]
      if (is.null(pathways) || length(pathways) == 0) { showModal(modalDialog(title = "Gene Not Found", paste("No pathways for gene", gene), easyClose = TRUE)); return() }
      valid_nodes <- V(g_local)$name
      exact_sources <- mapping$exact_source[match(pathways, mapping$processed_name)]
      valid_pathways <- pathways[exact_sources %in% valid_nodes]
      if (length(valid_pathways) == 0) { showModal(modalDialog(title = "No Valid Pathways", paste("No valid pathways for gene", gene), easyClose = TRUE)); return() }
      pathways_for_gene(valid_pathways)
      output$pathway_select_ui <- renderUI({
        selectInput("selected_pathway_from_gene", "Pathways containing selected gene:", choices = c(" " = "", valid_pathways), selected = " ")
      })
    })
    
    observeEvent(input$selected_pathway_from_gene, {
      selected_pathway <- input$selected_pathway_from_gene
      if (is.null(selected_pathway) || selected_pathway == "") return()
      selected_node <- mapping$exact_source[match(selected_pathway, mapping$processed_name)]
      if (is.na(selected_node)) { showModal(modalDialog(title = "Pathway Not Found", paste("Not found:", selected_pathway), easyClose = TRUE)); return() }
      child_pathways <- get_child_pathways(parent_geneset_name = selected_pathway, mapping = mapping, g = g_local)
      child_exact_source <- mapping$exact_source[match(child_pathways, mapping$processed_name)]
      child_exact_source <- child_exact_source[!is.na(child_exact_source)]
      vids <- c(selected_node, child_exact_source)
      if (length(vids) == 0) { showModal(modalDialog(title = "No Sub-Pathways", paste(selected_pathway, "has no sub-pathways."), easyClose = TRUE, footer = NULL)); return() }
      sub_g <- induced_subgraph(g_local, vids = vids)
      history$stack <- c(history$stack, list(current_graph()))
      current_graph(sub_g)
      last_selected_node(NULL)
      updateSelectizeInput(session, "gene_input", selected = "", server = TRUE)
      pathways_for_gene(NULL)
      output$pathway_select_ui <- renderUI(NULL)
    })
    
    output$legend_plot <- renderPlot({
      req(selected_comparisons(), new_enrichment_limits())
      g_current <- current_graph()
      edges_df <- as_data_frame(g_current, what = "edges")
      relation_colors <- c("is_a" = "black", "part_of" = "darkcyan", "regulates" = "orange",
                           "positively_regulates" = "steelblue", "negatively_regulates" = "salmon")
      relation_linetypes <- c("is_a" = "solid", "part_of" = "dashed", "regulates" = "dotted",
                              "positively_regulates" = "dotdash", "negatively_regulates" = "longdash")
      relations <- if ("relation" %in% colnames(edges_df)) unique(edges_df$relation) else character(0)
      relations_present <- relations[!is.na(relations)]
      relations_present <- relations_present[relations_present %in% names(relation_colors)]
      
      legend_plot <- generate_legend_plot(conditions = selected_comparisons(), enrichment_limits = new_enrichment_limits(),
                                          reference = reference_comparison(), lollipop_plot = lollipop_plot, lollipop_colors = lollipop_colors)
      if (is.null(legend_plot)) {
        ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Legend could not be generated.", size = 5) + theme_void()
      } else if (length(relations_present) > 0) {
        relation_legend <- generate_relation_legend(relations_present = relations_present, relation_colors = relation_colors, relation_linetypes = relation_linetypes)
        plot_grid(legend_plot, relation_legend, ncol = 1, rel_heights = c(3, 1))
      } else {
        legend_plot
      }
    }, height = 250)
    
    output$pathway_network <- renderVisNetwork({
      g_current <- current_graph()
      glyphs <- glyph_images()
      search_value <- input$search_pathway
      if (is.null(search_value) || length(search_value) == 0) search_value <- ""
      else { search_value <- trimws(search_value); if (length(search_value) == 0) search_value <- "" }
      
      node_labels <- V(g_current)$label
      node_names <- V(g_current)$name
      cutoff_val <- cutoff()
      
      gene_summaries <- vapply(node_labels, function(x) {
        if (!is.null(gene_sets[[x]])) {
          gl <- gene_sets[[x]]
          if (length(gl) > 10) paste(paste(gl[seq_len(cutoff_val)], collapse = ", "),
                                     "... <a href='#' class='show-full-genes' data-node='", x, "'>Show more</a>")
          else paste(gl, collapse = ", ")
        } else "No genes available"
      }, FUN.VALUE = character(1))
      
      nodes <- data.frame(
        id = node_names, label = node_labels,
        title = ifelse(node_labels == "root", "The Root node serves as a common parent for all top-level pathways.",
                       paste("Pathway:", node_labels, "<br>", "Genes:", gene_summaries, sep = "")),
        stringsAsFactors = FALSE
      )
      nodes$image <- vapply(node_labels, function(x) {
        img <- glyphs[[x]]; if (!is.null(img) && !is.na(img) && img != "") img else ""
      }, FUN.VALUE = character(1))
      nodes$shape <- ifelse(nodes$image != "", "image", "dot")
      nodes$size <- ifelse(nodes$id == "root", 30, ifelse(nodes$shape == "dot", 10, 15))
      nodes$font.size <- ifelse(nodes$id == "root", 20, 14)
      nodes$font.face <- "Arial"
      nodes$color.background <- if (is.null(input$search_pathway) || input$search_pathway == "") "lightblue"
      else ifelse(nodes$id %in% filtered_nodes(), "orange", "lightblue")
      nodes <- nodes[order(nodes$label), ]
      
      edges <- as_data_frame(g_current, what = "edges")
      relation_values <- c("is_a" = "solid", "part_of" = "dashed", "regulates" = "dotted",
                           "positively_regulates" = "dotdash", "negatively_regulates" = "longdash")
      line_patterns <- list("solid" = FALSE, "dashed" = c(5, 5), "dotted" = c(1, 5),
                            "dotdash" = c(1, 5, 5, 5), "longdash" = c(10, 5))
      rel_colors <- c("is_a" = "black", "part_of" = "darkcyan", "regulates" = "orange",
                      "positively_regulates" = "steelblue", "negatively_regulates" = "salmon")
      
      edges$linetype <- ifelse(is.na(edges$relation) | edges$relation == "is_a", "solid", relation_values[edges$relation])
      edges$width <- ifelse(is.na(edges$relation) | edges$relation == "is_a", 1, 2)
      edges$dashes <- lapply(edges$linetype, function(lt) line_patterns[[lt]])
      edges$color <- ifelse(is.na(edges$relation), rel_colors["is_a"], rel_colors[edges$relation])
      edges <- edges[, c("from", "to", "color", "dashes", "width"), drop = FALSE]
      arrow_dir <- ifelse(grepl("GO", V(g_local)[1]$name), "from", "to")
      
      visNetwork(nodes, edges, height = "800px", width = "100%") %>%
        visEdges(arrows = arrow_dir, smooth = c(enabled = TRUE, type = "curvedCCW", roundness = 0.8)) %>%
        visNodes(shapeProperties = list(useImageSize = TRUE), shadow = FALSE,
                             font = list(color = "black", size = nodes$font.size, face = nodes$font.face, bold = TRUE)) %>%
        visIgraphLayout(layout = layout_options_plot$layout, type = layout_options_plot$type, physics = layout_options_plot$physics) %>%
        visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = FALSE), nodesIdSelection = TRUE) %>%
        visInteraction(navigationButtons = TRUE, zoomView = TRUE, dragView = TRUE, multiselect = FALSE,
                                   hoverConnectedEdges = TRUE, tooltipStay = 500) %>%
        visEvents(
          type = "on",
          click = "function(nodes) { Shiny.setInputValue('pathway_network_selected', nodes.nodes, {priority: 'event'}); }",
          afterDrawing = "function() {
            document.addEventListener('click', function(e) {
              if (e.target && e.target.classList.contains('show-full-genes')) {
                e.preventDefault();
                var nodeName = e.target.getAttribute('data-node').trim();
                Shiny.setInputValue('tooltip_click', nodeName, {priority: 'event'});
              }
            });
          }"
        )
    })
    
    observeEvent(input$pathway_network_selected, {
      selected_node <- input$pathway_network_selected
      if (is.null(selected_node) || length(selected_node) == 0 || selected_node == "") return()
      selected_node <- selected_node[1]
      if (!is.null(last_selected_node()) && selected_node == last_selected_node()) return()
      last_selected_node(selected_node)
      if (selected_node == "root") { showModal(modalDialog(title = "Root Node", "Common parent for all top-level pathways.", easyClose = TRUE)); return() }
      selected_pathway <- mapping[which(mapping$exact_source == selected_node), "processed_name"]
      if (!(selected_node %in% V(g_local)$name)) { showModal(modalDialog(title = "Not Found", paste(selected_pathway, "not found."), easyClose = TRUE)); return() }
      child_pathways <- get_child_pathways(parent_geneset_name = selected_pathway, mapping = mapping, g = g_local)
      if (length(child_pathways) == 0) { showModal(modalDialog(title = "No Sub-Pathways", paste(selected_pathway, "has no sub-pathways."), easyClose = TRUE, footer = NULL)); return() }
      child_exact_source <- mapping$exact_source[match(child_pathways, mapping$processed_name)]
      child_exact_source <- child_exact_source[!is.na(child_exact_source)]
      if (length(child_exact_source) == 0) { showModal(modalDialog(title = "No Valid Sub-Pathways", paste(selected_pathway, "has no valid sub-pathways."), easyClose = TRUE)); return() }
      sub_g <- induced_subgraph(g_local, vids = c(selected_node, child_exact_source))
      history$stack <- c(history$stack, list(current_graph()))
      current_graph(sub_g)
      last_selected_node(NULL)
    }, ignoreInit = TRUE)
    
    session$onFlushed(function() { session$sendCustomMessage(type = "updateCutoff", message = NULL) })
    
    observeEvent(input$back, {
      if (length(history$stack) > 0) {
        last_graph <- tail(history$stack, 1)[[1]]
        history$stack <- head(history$stack, -1)
        current_graph(last_graph)
        last_selected_node(NULL)
      } else {
        showModal(modalDialog(title = "No Previous View", "No previous view to go back to.", easyClose = TRUE))
      }
    })
    
    observeEvent(input$cutoff_update, { cutoff(input$cutoff_update) })
    
    observeEvent(input$tooltip_click, {
      node_label <- input$tooltip_click
      if (!is.null(node_label)) {
        genes <- gene_sets[[node_label]]
        showModal(modalDialog(
          title = paste("All genes in", node_label),
          if (!is.null(genes)) paste(genes, collapse = ", ") else "No genes available for this pathway.",
          easyClose = TRUE, footer = modalButton("Close")
        ))
      }
    })
    
    observeEvent(input$reset, { history$stack <- list(); current_graph(initial_graph); last_selected_node(NULL) })
    
    observeEvent(input$help, {
      showModal(modalDialog(
        title = "Help",
        "Use the search box to find specific pathways. Click on a top-level pathway glyph to view its sub-pathways. Use 'Reset View' to return to the top-level pathways. Use 'Back' to return to the previous view.",
        easyClose = TRUE, footer = NULL
      ))
    })
  }
  shinyApp(ui = ui, server = server)
}