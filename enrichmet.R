enrichmet <- function(inputMetabolites, PathwayVsMetabolites, example_data, top_n = 100, p_value_cutoff = 1) {
  
  # ----------------- Convert Adjacency Matrix to List -------------------------
  matrix_to_list <- function(pws) {
    pws.l <- list()
    for (pw in colnames(pws)) {
      pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
    }
    return(pws.l)
  }
  
  # ----------------- Prepare GMT File from PathwayVsMetabolites -----------------
  prepare_gmt <- function(gmt_file, metabolites_in_data, savefile = FALSE) {
    gmt <- gmtPathways(gmt_file)
    hidden <- unique(unlist(gmt))
    
    mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                  nrow = length(hidden), ncol = length(gmt))
    for (i in 1:dim(mat)[2]){
      mat[,i] <- as.numeric(hidden %in% gmt[[i]])
    }
    
    hidden1 <- intersect(metabolites_in_data, hidden)
    mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,]) > 5)]] 
    
    final_list <- matrix_to_list(mat)
    
    if(savefile){
      saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
    }
    
    return(final_list)
  }
  
  # ----------------- Convert PathwayVsMetabolites to GMT Format -----------------
  PathwayVsMetabolites$description <- "https://www.genome.jp/kegg/pathway.html#metabolism"
  
  convert_to_gmt <- function(pathway, description, metabolites) {
    pathway_underscore <- gsub(" ", "_", pathway)
    metabolites_vector <- unlist(strsplit(metabolites, ","))
    gmt_line <- paste(pathway_underscore, description, paste(metabolites_vector, collapse = "\t"), sep = "\t")
    return(gmt_line)
  }
  
  gmt_data <- mapply(convert_to_gmt, PathwayVsMetabolites$Pathway, PathwayVsMetabolites$description, PathwayVsMetabolites$Metabolites, SIMPLIFY = TRUE)
  
  gmt_file <- "output.gmt"
  writeLines(gmt_data, gmt_file)
  
  # ----------------- Prepare Data -----------------
  data <- PathwayVsMetabolites %>%
    mutate(Metabolites = strsplit(Metabolites, ",")) %>%
    unnest(Metabolites)
  
  allMetabolitesSet <- unique(data$Metabolites)
  
  # ----------------- Compute Betweenness Centrality -----------------
  edge_list_pathways <- data.frame(from = rep(data$Pathway, lengths(data$Metabolites)), to = unlist(data$Metabolites))
  g_pathways <- graph_from_data_frame(d = edge_list_pathways, directed = FALSE)
  
  edge_list_metabolites <- data.frame(from = unlist(data$Metabolites), to = rep(data$Pathway, lengths(data$Metabolites)))
  g_metabolites <- graph_from_data_frame(d = edge_list_metabolites, directed = FALSE)
  
  betweenness_metabolites <- betweenness(g_metabolites, directed = FALSE, normalized = TRUE)
  metabolite_centrality <- data.frame(Metabolite = names(betweenness_metabolites), RBC_Metabolite = betweenness_metabolites)
  
  input_metabolite_centrality <- metabolite_centrality %>%
    filter(Metabolite %in% inputMetabolites) %>%
    arrange(desc(RBC_Metabolite)) %>%
    filter(RBC_Metabolite > 0)
  
  # ----------------- Perform Fisher’s Exact Test -----------------
  results <- list()
  
  for (i in 1:nrow(PathwayVsMetabolites)) {
    row <- PathwayVsMetabolites[i, ]
    pathway <- row$Pathway
    pathwayMetabolites <- unlist(strsplit(row$Metabolites, ","))
    
    a <- length(intersect(pathwayMetabolites, inputMetabolites))
    b <- length(setdiff(inputMetabolites, pathwayMetabolites))
    c <- length(setdiff(pathwayMetabolites, inputMetabolites))
    d <- length(setdiff(allMetabolitesSet, union(inputMetabolites, pathwayMetabolites)))
    
    contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    
    fisher_test_result <- fisher.test(contingency_table, alternative = "two.sided")
    
    results[[i]] <- list(Pathway = pathway, P_value = fisher_test_result$p.value, Log_P_value = -log10(fisher_test_result$p.value))
  }
  
  results_df <- do.call(rbind, lapply(results, as.data.frame))
  results_df$Adjusted_P_value <- p.adjust(results_df$P_value, method = "BH")
  
  significant_results_df <- results_df %>%
    filter(Adjusted_P_value < p_value_cutoff) %>%
    arrange(desc(Log_P_value))
  
  if (!is.null(top_n)) {
    significant_results_df <- head(significant_results_df, top_n)
  }
  
  # ----------------- Pathway Enrichment Plot -----------------
  pathway_plot <- ggplot(significant_results_df, aes(x = reorder(Pathway, Log_P_value), y = Log_P_value, fill = Adjusted_P_value)) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient(low = "red", high = "blue") +
    labs(x = "Pathway", y = "-log10(P-value)", fill = "Adj P-value") +
    theme_minimal()
  
  # ----------------- GSEA Preparation -----------------
  # Prepare example data for GSEA
  example_data = read_xlsx("Documents/Neha/example_data.xlsx")
  example_filtered_data <- example_data %>%
    arrange(pval) %>%           # Sort by p-value
    distinct(met_id, .keep_all = TRUE) 
  
  example_cleaned <- example_filtered_data %>%
    filter(met_id != "No Metabolites found")
  
  meta = example_cleaned$met_id
  bg_metabolites = prepare_gmt(gmt_file, meta, savefile = FALSE)
  
  # Prepare ranked list of metabolites for GSEA
  example_cleaned$pval = as.numeric(example_cleaned$pval)
  rankings = sign(example_cleaned$log2fc) * (-log10(example_cleaned$pval))
  names(rankings) <- example_cleaned$met_id  # metabolites as names
  rankings <- sort(rankings, decreasing = TRUE)
  
  # Run GSEA
  MSEAres <- fgsea(pathways = bg_metabolites, 
                   stats = rankings,
                   scoreType = 'std',
                   minSize = 10,
                   maxSize = 500,
                   nproc = 1)
  
  # ----------------- GSEA Plot -----------------
  # Ensure leadingEdge is in character format and then split by commas
  MSEAres$input_count <- sapply(strsplit(as.character(MSEAres$leadingEdge), ","), length)
  
  # Sorting based on p-value
  MSEAres <- MSEAres %>%
    arrange(pval) 
  
  # Convert 'pathway' to a factor with levels in the order they appear in the file
  MSEAres$pathway <- factor(MSEAres$pathway, levels = MSEAres$pathway)
  
  # GSEA dot plot
  gsea_plot <- ggplot(MSEAres, aes(x = -log10(pval), y = pathway, size = input_count, color = NES)) +
    geom_point() +
    labs(x = "-log10(p-value)", y = "Pathway", size = "Metabolite count", color = "NES") +
    scale_size_continuous(labels = scales::number_format(accuracy = 1)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
          axis.ticks.y = element_line(color = "black"))
  
  # ----------------- Network Plot -----------------
  set.seed(42)  # For reproducibility
  edges <- data.frame(
    from = sample(input_metabolite_centrality$Metabolite, 40, replace = TRUE),
    to = sample(input_metabolite_centrality$Metabolite, 40, replace = TRUE)
  )
  
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  V(g)$RBC <- input_metabolite_centrality$RBC_Metabolite[match(V(g)$name, input_metabolite_centrality$Metabolite)]
  
  network_plot <- ggraph(g, layout = "fr") +
    geom_edge_link(alpha = 0.3, color = "black") +
    geom_node_point(aes(size = RBC, color = RBC, fill = RBC), shape = 21, stroke = 0.5) +
    geom_node_text(aes(label = name), repel = TRUE, size = 4) +
    scale_color_gradient(low = "blue", high = "red") +
    scale_fill_gradient(low = "blue", high = "red") +
    labs(title = "Metabolite Interaction Network") +
    theme_void() +
    theme(legend.position = "none")
  
  return(list(pathway_plot = pathway_plot, gsea_plot = gsea_plot, network_plot = network_plot))
}

results <- enrichmet(inputMetabolites, PathwayVsMetabolites,example_data, top_n = 10, p_value_cutoff = 1)

# Display plots
results$pathway_plot
results$network_plot
results$gsea_plot

