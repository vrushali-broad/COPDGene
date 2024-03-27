# This script provides a comprehensive workflow for conducting enrichment analysis and visualization
# on differential expression data, specifically tailored for exploring the genetic underpinnings of
# diseases like COPD (Chronic Obstructive Pulmonary Disease) and ACO (Asthma-COPD Overlap). Utilizing
# a suite of Bioconductor packages, the script performs KEGG pathway and GO term enrichment analysis,
# generates chord and tree diagrams for visual representation, and logs actions with timestamps for
# traceability.

# Key features include:
# - Data loading and preprocessing, including mapping gene IDs and filtering based on expression levels.
# - Enrichment analysis using `clusterProfiler` for KEGG pathways and GO terms.
# - Visualization of enrichment results with `GOplot`, `enrichplot`, and custom plotting functions.
# - Logging of actions and package versions to ensure reproducibility.

# Prerequisites:
# The script requires R with the following packages installed: `limma`, `edgeR`, `clusterProfiler`,
# `org.Hs.eg.db`, `enrichplot`, `GOplot`, `data.table`, and `dplyr`. Additional dependencies may include
# `AnnotationDbi`, `biomaRt`, `RColorBrewer`, `EnhancedVolcano`, `DOSE`, `tibble`, `GOSemSim`, `ReactomePA`.

# Usage:
# Modify the `conditions` vector to include the conditions of interest (e.g., COPD, ACO) and set the
# appropriate file paths for input data. Ensure differential expression data is in the expected format
# and location. Execute the script in an R environment, and results will be generated and saved in
# the specified directories.

# Example:
# To perform analysis for COPD and ACO, ensure data files are prepared, and run:
# conditions <- c('COPD', 'ACO')
# main(conditions)

# Author: Vrushali D. Fangal, 
#         Channing Division of Network Medicine, 
#         Harvard Medical School, 
#         Boston, MA.
# Date: March 22, 2024
# License: MIT license


# Loading required libraries
library('limma')
library("edgeR")
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("GOplot")
library(data.table)
library(dplyr)


# library('Biobase')
# library('limma')
# library("biomaRt")
# library("edgeR")
# library('RColorBrewer')
# library('EnhancedVolcano')
# library("clusterProfiler")
# library("DOSE")
# library("tibble")
# library("AnnotationDbi")
# library("org.Hs.eg.db")
# library("enrichplot")
# library("GOSemSim")
# library("GOplot")
# library("biomaRt")
# library('ReactomePA')
# library(data.table)

# Logging function
log_with_timestamp <- function(message) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(paste("[", timestamp, "] ", message, "\n", sep=""))
}

# Log package versions
log_with_timestamp(paste("limma version:", packageVersion("limma")))
log_with_timestamp(paste("edgeR version:", packageVersion("edgeR")))
log_with_timestamp(paste("clusterProfiler version:", packageVersion("clusterProfiler")))
log_with_timestamp(paste("org.Hs.eg.db version:", packageVersion("org.Hs.eg.db")))
log_with_timestamp(paste("GOplot version:", packageVersion("GOplot")))
log_with_timestamp(paste("enrichplot version:", packageVersion("enrichplot")))
log_with_timestamp(paste("data.table version:", packageVersion("data.table")))
log_with_timestamp(paste("dplyr version:", packageVersion("dplyr")))

load_data <- function(PATH, condition = "COPD", lfc_threshold = 0) {
    log_with_timestamp("Loading annotations...")
    annot_file <- file.path(PATH, 'Processed_Data', 'Annotations.tsv')
    gene_names_file <- file.path(PATH, 'Processed_Data', "gene_names.txt")
    de_genes_file <- file.path(PATH, 'Processed_Data', paste('DE_', condition, '.csv', sep = ''))

    # Check if necessary files exist
    if (!file.exists(annot_file) || !file.exists(gene_names_file) || !file.exists(de_genes_file)) {
        log_with_timestamp(paste("One or more required files are missing. Skipping data loading for condition:", condition))
        return(NULL)
    }

    annot <- read.csv(annot_file, sep = '\t', header = TRUE)
    annot <- annot[!duplicated(annot$ensembl_gene_id), ]
    rownames(annot) <- annot$ensembl_gene_id
    log_with_timestamp(paste("Dimensions of annotations:", paste(dim(annot), collapse = " x ")))

    log_with_timestamp("Loading gene names...")
    lcpm <- readLines(gene_names_file)
    universe <- annot[lcpm, 'hgnc_symbol']
    log_with_timestamp(paste("Number of genes in universe:", length(universe)))
    
    log_with_timestamp(paste("Loading DE analysis results for", condition, "..."))
    de_genes <- read.csv(de_genes_file, header = TRUE, row.names = 1)
    log_with_timestamp(paste("Dimensions of DE analysis results:", paste(dim(de_genes), collapse = " x ")))

    # Filter for upregulated and downregulated genes
    fcvals_up <- de_genes[de_genes$logFC > lfc_threshold, ]
    fcvals_down <- de_genes[de_genes$logFC < -lfc_threshold, ]

    return(list(annot = annot, 
                universe = universe, 
                de_genes = de_genes,
                fcvals_up = fcvals_up,
                fcvals_down = fcvals_down))
}


# KEGG Enrichment Analysis Function
perform_kegg_analysis <- function(de_genes, annot, condition) {
    
    log_with_timestamp("Performing KEGG Enrichment Analysis...")
    de_genes_list <- de_genes
    # de_genes_list <- de_genes[de_genes$logFC > 0, ]
    gene_list <- de_genes_list$logFC
    names(gene_list) <- de_genes_list$entrezgene_id
    gene_list <- na.omit(gene_list)
    gene_list <- sort(gene_list, decreasing = TRUE)

    kegg <- enrichKEGG(gene = names(gene_list),
                       organism = 'hsa',
                       pvalueCutoff = 0.05)
    kegg <- as.data.frame(setReadable(kegg, 'org.Hs.eg.db', keyType="ENTREZID"))
    
    
    # Calculate enrichment scores
    enrichment_scores <- kegg %>%
      mutate(
        BgRatioTotal = as.numeric(gsub("/.*", "", BgRatio)), # Extract total gene count from BgRatio
        ExpectedCount = (as.numeric(gsub(".*/", "", GeneRatio)) / BgRatioTotal) * Count, # Calculate expected count
        EnrichmentScore = Count / ExpectedCount # Calculate enrichment score
      ) %>%
      select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count, EnrichmentScore)

    # View the results
    # print(enrichment_scores)
    
    # # Sort the data based on EnrichmentScore
    # sorted_kegg <- enrichment_scores %>%
    #   arrange(desc(EnrichmentScore))
    
    # # Sort the data based on Count, pvalue, and EnrichmentScore
    # sorted_kegg <- enrichment_scores %>%
    #   arrange(desc(Count), p.adjust, desc(EnrichmentScore))
    
    # Sort the data based on Count, pvalue, and EnrichmentScore
    sorted_kegg <- enrichment_scores %>%
      arrange(p.adjust, desc(Count), desc(EnrichmentScore))
    

    # Directory for saving results
    results_dir <- file.path("Figures", "Enrichment_analysis")
    if (!dir.exists(results_dir)) {
        dir.create(results_dir, recursive = TRUE)
    }

    # Create a file name using the condition
    file_name <- file.path(results_dir, paste0(condition, "_KEGG_Enrichment_Analysis_Results.csv"))
    # # Save the results to a CSV file
    # write.csv(kegg, file_name, row.names = FALSE)
    # Save the sorted results to a CSV file
    write.csv(sorted_kegg, file_name, row.names = FALSE)
    log_with_timestamp(paste("KEGG Enrichment Analysis results saved in", file_name))
    
    # return(kegg)
}

# GO Enrichment Analysis Function
perform_go_analysis <- function(fcvals_up, annot, universe, lfc_threshold = 0, condition) {
    log_with_timestamp("Performing GO Enrichment Analysis...")

    # Filter for upregulated genes
    # fcvals_up <- de_genes[de_genes$logFC > lfc_threshold, ]

    # Process upregulated genes (you can replicate this for downregulated genes if needed)
    degenes <- subset(fcvals_up, select=c('hgnc_symbol', 'logFC'))
    degenes$logFC <- degenes$logFC 
    degenes <- degenes[!duplicated(degenes$hgnc_symbol), ]
    colnames(degenes) <- c('ID', 'logFC')

    genelist <- fcvals_up[,c('hgnc_symbol', 'logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')]
    names(genelist)[1] <- "ID"

    gene_list <- genelist$logFC
    names(gene_list) <- annot[rownames(genelist),]$hgnc_symbol
    gene_list <- na.omit(gene_list)
    gene_list <- sort(gene_list, decreasing = TRUE)

    go <- enrichGO(genelist$ID, 
                   OrgDb = "org.Hs.eg.db", 
                   keyType = 'SYMBOL',
                   ont = "all",
                   minGSSize = 10, 
                   maxGSSize = 100, 
                   universe = universe, 
                   readable = TRUE)
    
    # Directory for saving results
    results_dir <- file.path("Figures", "Enrichment_analysis")
    if (!dir.exists(results_dir)) {
        dir.create(results_dir, recursive = TRUE)
    }

    # Save the GO results to a CSV file
    go_file_name <- file.path(results_dir, paste0("GO_Enrichment_", condition, ".csv"))
    write.csv(as.data.frame(go), go_file_name, row.names = FALSE)
    log_with_timestamp(paste("GO Enrichment Analysis results saved in", go_file_name))

    return(list(go = go, genelist = genelist, degenes = degenes))
}

generate_chord_diagram <- function(go_results, genelist, degenes, condition) {
    # Filter GO results for Biological Processes
    go <- filter(go_results, ONTOLOGY == 'BP')
    print(dim(go))

    # Create a data frame for the enrichment results
    EC <- data.frame(Category = go$ONTOLOGY,
                     ID = go$ID, 
                     Term = go$Description,
                     Genes = gsub("/", ", ", go$geneID),
                     adj_pval = go$p.adjust,
                     count = go$Count)

    # Process the data for chord diagram
    circ <- circle_dat(EC, genelist)
    reduced_circ <- reduce_overlap(circ, overlap = 0.9)
    
    # Create chord data for the reduced overlapping terms
    chord <- chord_dat(data = circ, genes = degenes, process = reduced_circ[1:5, 'term'])

    # Plotting the chord diagram with specified parameters
    p <- GOChord(chord, 
                 space = 0.02, 
                 gene.order = 'logFC', 
                 lfc.col = c('#001ecc', '#ffffff', '#43529E'),
                 gene.space = 0.25,
                 ribbon.col = c("#ffe599", "#ea9990", "#d9d2e9", "#d9ead3", "#f9cb9c"),
                 border.size = 0,
                 limit = c(0, 0),
                 gene.size = 2.7,
                 process.label = 5.5)
    
    # Directory for saving results
    results_dir <- file.path("Figures", "Enrichment_analysis")
    if (!dir.exists(results_dir)) {
        dir.create(results_dir, recursive = TRUE)
    }
    
        
    # Save the data used to generate chord plot
    file_name <- file.path(results_dir, paste0('chord_', condition, '.csv'))
    write.csv(EC[EC$Term %in% reduced_circ[1:5, 'term'],], file_name, row.names = FALSE)
    
    
    # Construct the file name and save the plot
    file_name <- file.path(results_dir, paste0('chord_', condition, '.pdf'))
    # Save the plot
    ggsave(p, file = file_name, width = 5.7, height = 7.2)
    # Log the location of the saved file
    log_with_timestamp(paste("Chord diagram saved in", file_name))

}

plot_tree_diagram <- function(de_genes, annot, fcvals_up, fcvals_down, midpoint, condition) {
    log_with_timestamp("Preparing data for Tree Diagram...")

    # Combine up and downregulated genes
    de_genes_list <- rbind(fcvals_down, fcvals_up)
    de_genes_list['hgnc_symbol'] <- annot[rownames(de_genes_list), 'hgnc_symbol']
    
    # Filter out rows with missing values
    # de_genes_list <- na.omit(de_genes_list)
    
    # Create a gene list with necessary columns
    genelist <- de_genes_list[, c('hgnc_symbol', 'logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')]
    names(genelist)[1] <- "ID"
    
    # Ensure logFC is a continuous variable for proper plotting
    genelist$logFC <- as.numeric(as.character(genelist$logFC))
    
    # Create a sorted gene list for enrichment
    gene_list <- genelist$logFC
    names(gene_list) <- annot[rownames(genelist), ]$hgnc_symbol
    gene_list <- na.omit(gene_list)
    gene_list <- sort(gene_list, decreasing = TRUE)

    # Perform GO enrichment analysis
    log_with_timestamp("Performing GO Enrichment Analysis for Tree Diagram...")
    go_all <- enrichGO(genelist$ID, 
                       OrgDb = "org.Hs.eg.db", 
                       keyType = 'SYMBOL',
                       ont = "all",
                       readable = TRUE)
    
    # Filter for Biological Process and compute pairwise term similarity
    ego <- go_all #filter(go_all, ONTOLOGY == 'BP') 
    
    # Step 1: Extract the data frame from the go_all object
    ego_data <- ego@result

    # Step 2: Filter for 'BP' ontology and then keep unique geneIDs
    ego_filtered_data <- ego_data %>%
      filter(ONTOLOGY == 'BP') %>%
      { .[!duplicated(.$geneID), ] }
    ego@result <- ego_filtered_data
    
    # ego <- filter(go_all, ONTOLOGY == 'BP') 
    edox <- pairwise_termsim(ego)
    edox <- mutate(edox, log10.p.adjust = -log10(as.numeric(p.adjust)))
    edox2 <- pairwise_termsim(edox)
    

    # Plot the tree diagram
    log_with_timestamp("Plotting Tree Diagram...")
    p <- treeplot(edox2,
                   cluster.params = list(label_words_n = 5,
                                         method = 'ward.D',
                                         n = 3),
                   color = 'log10.p.adjust',
                   showCategory = 10, fontsize = 4)
    
    # Comment out the conflicting color scale if treeplot already sets one
    p <- p + scale_color_gradient2(name = expression("-Log"[10]*"FDR"),
                                   low = '#ffe2d2',
                                   mid = "#ffcc8c", 
                                   high = "#e5661c", 
                                   midpoint = midpoint,
                                   guide = "colourbar")
    
    # Save the plot
    # Directory for saving results
    results_dir <- file.path("Figures", "Enrichment_analysis")
    if (!dir.exists(results_dir)) {
        dir.create(results_dir, recursive = TRUE)
    }
    
    
    # Save the data used to generate tree plot
    file_name <- file.path(results_dir, paste0('tree_', condition, '.csv'))
    write.csv(head(edox2, 10), file_name, row.names = FALSE)

    
    # Construct the file name based on the condition
    file_name <- file.path(results_dir, paste0('tree_', condition, '.pdf'))
    # file_name <- paste0("Figures/tree_", condition, ".pdf")

    # Save the plot with the dynamically generated file name
    ggsave(p, file = file_name, width = 8, height = 8)

    log_with_timestamp(paste("Tree diagram saved in", file_name))

}


# Main Execution Function
main <- function(conditions) {
    log_with_timestamp("Starting the enrichment analysis pipeline...")

    for (condition in conditions) {
        # Load data
        log_message <- paste("\n##################################################",
                             "\n        Loading data for condition:", condition, 
                             "       \n##################################################\n")
        log_with_timestamp(log_message)
        
        # Setting up paths
        PATH <- getwd()
        PATH <- paste(PATH, 'Data', sep = .Platform$file.sep)
        data <- load_data(PATH, condition = condition)
        log_with_timestamp("Data loading complete.")
        
        # Check if data is NULL, and stop execution if it is
        if (is.null(data)) {
            log_with_timestamp(paste("Data loading failed for condition:", condition, ". Skipping to next condition."))
            next  # Skip to the next iteration of the loop
        }
        
        # Perform KEGG Enrichment Analysis
        log_with_timestamp("Starting KEGG Enrichment Analysis...")
        kegg_results <- perform_kegg_analysis(data$de_genes, data$annot, condition)
        log_with_timestamp("KEGG Enrichment Analysis complete.")

        # Perform GO Enrichment Analysis
        log_with_timestamp("Starting GO Enrichment Analysis...")
        results <- perform_go_analysis(data$fcvals_up, data$annot, data$universe,lfc_threshold = 0, condition)
        go_results <- results$go
        genelist = results$genelist
        degenes <- results$degenes
        log_with_timestamp("GO Enrichment Analysis complete.")

        # Plot Chord Diagram
        # log_with_timestamp("Plotting Chord Diagram...")
        generate_chord_diagram(go_results, genelist, degenes, condition)
        log_with_timestamp("Chord Diagram plotting complete.")

        # Plot Tree Diagram
        log_with_timestamp("Plotting Tree Diagram...")
        # Set the midpoint based on condition
        if (condition == "COPD") {
            midpoint_value <- 6.25
        } else if (condition == "ACO") {
            midpoint_value <- 2.83
        } else {
            midpoint_value <- 7.5  # Default value or handle as needed
        }
        plot_tree_diagram(data$de_genes, data$annot, data$fcvals_up, data$fcvals_down, midpoint_value, condition)
        log_with_timestamp("Tree Diagram plotting complete.")

        log_with_timestamp(paste("Enrichment analysis completed for condition:", condition))
    }
}


# Conditions to run the analysis for
conditions <- c('COPD', 'ACO')  #'Asthma',
# conditiona <- c('ACO')

# Run the main function for the specified conditions
main(conditions)

# sessionInfo()


















# # Load Data Function with Condition Parameter
# load_data <- function(PATH, condition = "COPD", lfc_threshold = 0) {
#     log_with_timestamp("Loading annotations...")
#     annot <- read.csv(file.path(PATH, 'Processed_Data', 'Annotations.tsv'), sep = '\t', header = TRUE)
#     annot <- annot[!duplicated(annot$ensembl_gene_id), ]
#     rownames(annot) <- annot$ensembl_gene_id
#     # log_with_timestamp(paste("Dimensions of annotations:", dim(annot)))
#     log_with_timestamp(paste("Dimensions of annotations:", paste(dim(annot), collapse = " x ")))

#     log_with_timestamp("Loading gene names...")
#     lcpm <- readLines(file.path(PATH, 'Processed_Data', "gene_names.txt"))
#     universe <- annot[lcpm, 'hgnc_symbol']
#     log_with_timestamp(paste("Number of genes in universe:", length(universe)))
    
#     log_with_timestamp(paste("Loading DE analysis results for", condition, "..."))
#     de_genes_file <- paste('DE_', condition, '.csv', sep = '')
#     de_genes <- read.csv(file.path(PATH, 'Processed_Data', de_genes_file), header = TRUE, row.names = 1)
#     log_with_timestamp(paste("Dimensions of DE analysis results:", paste(dim(de_genes), collapse = " x ")))

#     # Filter for upregulated genes
#     fcvals_up <- de_genes[de_genes$logFC > lfc_threshold, ]
#     # Filter for downregulated genes
#     fcvals_down <- de_genes[de_genes$logFC < -lfc_threshold, ]

#     return(list(annot = annot, 
#                 universe = universe, 
#                 de_genes = de_genes,
#                 fcvals_up = fcvals_up,
#                 fcvals_down = fcvals_down))
# }

# # Chord Diagram Plotting Function
# plot_chord_diagram <- function(go, genelist, degenes) {
#     log_with_timestamp("Plotting Chord Diagram...")
    
#     ## Extract biological processes ontology
#     go <- filter(go, ONTOLOGY == 'BP')
    
#     # Set the GO results to a readable format
#     # go <- setReadable(go, 'org.Hs.eg.db', keyType = "ENSEMBL")
    
#     # Assuming 'circle_dat' and 'chord_dat' are custom functions or part of a specific package
#     EC <- data.frame(Category = go$ONTOLOGY,
#                      ID = go$ID, 
#                      Term = go$Description,
#                      Genes = gsub("/", ", ", go$geneID),
#                      adj_pval = go$p.adjust,
#                      count = go$Count)
    
#     # Use circle_dat to process the data
#     circ <- circle_dat(EC, genelist)
    
#     # Apply reduce_overlap to minimize process overlap
#     reduced_circ <- reduce_overlap(circ, overlap = 0.9)
    
#     log_with_timestamp(reduced_circ)
    
#     chord <- chord_dat(data = circ, 
#                        genes = degenes,
#                        process = reduced_circ[1:5, 'term'])

#     # Save the plot to a PDF
#     # ch_name <- paste0('Figures/chord_', dis, '.pdf')
#     # pdf(file = ch_name, bg = "transparent", width = 5.7, height = 7.2)
#     pdf(file = "Figures/chord_diagram.pdf", bg = "transparent", width = 8, height = 6)
    
#     # Plotting the chord diagram with specified parameters
#     GOChord(chord, 
#                  space = 0.02, 
#                  gene.order = 'logFC',
#                  lfc.col = c('#001ecc', '#ffffff', '#43529E'),
#                  gene.space = 0.25,
#                  ribbon.col = c("#ffe599", "#ea9990", "#d9d2e9", "#d9ead3", "#f9cb9c"),
#                  border.size = 0,
#                  limit = c(0, 0),
#                  gene.size = 2.7,
#                  process.label = 5.5)

#     dev.off()
# }

# # Tree Diagram Plotting Function
# plot_tree_diagram <- function(go_all) {
#     log_with_timestamp("Plotting Tree Diagram...")
#     ego <- filter(go_all, ONTOLOGY == 'BP')
#     edox <- pairwise_termsim(ego)
#     edox <- mutate(edox, log10.p.adjust = -log10(as.numeric(p.adjust)))
#     edox2 <- pairwise_termsim(edox)

#     p2 <- treeplot(edox2, color = 'log10.p.adjust')
#     ggsave(p2, file = "tree_diagram.pdf", width = 8, height = 8)
# }


# # Load annotations
# log_with_timestamp("Loading annotations")
# annot <- read.csv(file.path(PATH_DATA, 'Annotations.tsv'), sep = '\t', header = TRUE)
# annot <- annot[!duplicated(annot$ensembl_gene_id), ]
# rownames(annot) <- annot$ensembl_gene_id

# # Load gene names
# log_with_timestamp("Loading gene names")
# lcpm <- readLines("row_names.txt")
# universe = annot[lcpm, 'hgnc_symbol']

# # Load DE analysis results
# log_with_timestamp("Loading DE analysis results")
# de_genes <- read.csv(file.path(PATH_DATA, 'Processed_Data', 'DE_COPD.csv'), header = TRUE, row.names = 1)

# # Set log fold change threshold and filter genes
# lfc_threshold = 0
# fcvals_up <- de_genes[de_genes$logFC > lfc_threshold, ]
# fcvals_down <- de_genes[de_genes$logFC < -lfc_threshold, ]

# # Load Data Function
# load_data <- function() {
#     log_with_timestamp("Loading annotations...")
#     annot <- read.csv(file.path(PATH_DATA, 'Annotations.tsv'), sep = '\t', header = TRUE)
#     annot <- annot[!duplicated(annot$ensembl_gene_id), ]
#     rownames(annot) <- annot$ensembl_gene_id

#     log_with_timestamp("Loading gene names...")
#     lcpm <- readLines(file.path(PATH_DATA, "row_names.txt"))
#     universe <- annot[lcpm, 'hgnc_symbol']

#     log_with_timestamp("Loading DE analysis results...")
#     de_genes <- read.csv(file.path(PATH_DATA, 'Processed_Data', 'DE_COPD.csv'), header = TRUE, row.names = 1)
    
#     return(list(annot = annot, universe = universe, de_genes = de_genes))
# }

# # Function for KEGG Enrichment Analysis
# perform_kegg_analysis <- function(gene_list, organism = 'hsa') {
#     log_with_timestamp("Performing KEGG enrichment analysis")
#     kegg <- enrichKEGG(gene = names(gene_list),
#                        organism = organism,
#                        keyType = "ENTREZID")
#     kegg <- as.data.frame(setReadable(kegg, 'org.Hs.eg.db', keyType="ENTREZID"))
#     return(kegg)
# }
# # Function for GO Enrichment Analysis
# perform_go_analysis <- function(gene_list, universe, ontology = "all") {
#     log_with_timestamp("Performing GO enrichment analysis")
#     go_all <- enrichGO(gene = names(gene_list),
#                        OrgDb = "org.Hs.eg.db",
#                        keyType = 'SYMBOL',
#                        ont = ontology,
#                        universe = universe,
#                        readable = TRUE)
#     return(go_all)
# }

# # Main Execution
# main <- function() {
#     data <- load_data()
#     kegg <- perform_kegg_analysis(data$de_genes, data$annot)
#     go <- perform_go_analysis(data$de_genes, data$annot)

#     plot_chord_diagram(go, data$de_genes)
#     plot_tree_diagram(go)
# }

# # Run the main function
# main()

# # GO Enrichment Analysis Function
# perform_go_analysis <- function(de_genes, annot, universe, lfc_threshold = 0) {
#     log_with_timestamp("Performing GO Enrichment Analysis...")
    
#     # Filter for upregulated genes
#     de_genes <- de_genes[de_genes$logFC > lfc_threshold, ]
#     degenes <- subset(de_genes, select=c('hgnc_symbol','logFC'))
#     degenes$logFC <- degenes$logFC 
#     degenes <- degenes[!duplicated(degenes$hgnc_symbol), ]
#     colnames(degenes) <- c('ID', 'logFC')

#     genelist <- de_genes[,c('hgnc_symbol','logFC', 'AveExpr', 't','P.Value', 'adj.P.Val', 'B')]
#     names(genelist)[1] <- "ID"

#     gene_list <- genelist$logFC
#     names(gene_list) <- annot[rownames(genelist),]$hgnc_symbol
#     gene_list <- na.omit(gene_list)
#     gene_list <- sort(gene_list, decreasing = TRUE)
    
#     go <- enrichGO(genelist$ID, 
#                     OrgDb = "org.Hs.eg.db", 
#                     keyType = 'SYMBOL',
#                     ont = "all",
#                     minGSSize = 10, 
#                     maxGSSize = 100, 
#                     universe = universe, 
#                     readable = TRUE)
    
#     # go <- enrichGO(genelist$ID, 
#     #                OrgDb = "org.Hs.eg.db", 
#     #                keyType = 'SYMBOL',
#     #                ont = "all",
#     #                pvalueCutoff = 0.05,
#     #                minGSSize = 10, 
#     #                maxGSSize = 500,
#     #                universe = universe,
#     #                readable = TRUE)
    
#     # return(go)
#     return(list(go = go, degenes = degenes))
# }

# process_go_results <- function(go, genelist) {
#     # Set the GO results to a readable format
#     go <- setReadable(go, 'org.Hs.eg.db', keyType = "ENSEMBL")
    
#     # Create a data frame for the enrichment results
#     EC <- data.frame(Category = go$ONTOLOGY,
#                      ID = go$ID, 
#                      Term = go$Description,
#                      Genes = gsub("/", ", ", go$geneID),
#                      adj_pval = go$p.adjust,
#                      count = go$Count)

#     # Use circle_dat to process the data
#     circ <- circle_dat(EC, genelist)

#     # Apply reduce_overlap to minimize process overlap
#     reduced_circ <- reduce_overlap(circ, overlap = 0.9)

#     return(circ, reduced_circ)
# }
