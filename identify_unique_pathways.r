# Load necessary libraries
library(VennDiagram)
library(grid)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

# Logging function
log_with_timestamp <- function(message) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(paste("[", timestamp, "] ", message, "\n", sep=""))
}

# Initialize lists to store DE genes
genes_Asthma <- c()
genes_COPD <- c()
genes_ACO <- c()
genes_Asthma_up <- c()
genes_COPD_up <- c()
genes_ACO_up <- c()
genes_Asthma_down <- c()
genes_COPD_down <- c()
genes_ACO_down <- c()

# Helper function to read and filter DE genes
read_and_filter_genes <- function(filepath, lfc_threshold) {
  if (file.exists(filepath)) {
    DE_genes <- read.csv(filepath)
    DE_genes_up <- na.omit(DE_genes$hgnc_symbol[DE_genes$logFC > lfc_threshold])  # Upregulated
    DE_genes_down <- na.omit(DE_genes$hgnc_symbol[DE_genes$logFC < -lfc_threshold])  # Downregulated
    DE_genes_all <- na.omit(DE_genes$hgnc_symbol)  # All DE genes
    DE_gene_list <- read.csv(filepath, header = TRUE, row.names = 1)
    return(list(all = DE_genes_all, up = DE_genes_up, down = DE_genes_down, de_genes = DE_gene_list))
  } else {
    message(paste("DE genes file for", filepath, "does not exist."))
    return(list(all = c(), up = c(), down = c()))
  }
}

# Set log fold change threshold
lfc_threshold <- 0

# Read and filter DE genes for each condition
genes_Asthma <- read_and_filter_genes("Data/Processed_Data/DE_Asthma.csv", lfc_threshold)
genes_Asthma_all <- genes_Asthma$all
genes_Asthma_up <- genes_Asthma$up
genes_Asthma_down <- genes_Asthma$down
genes_Asthma_full <- genes_Asthma$de_genes

genes_COPD <- read_and_filter_genes("Data/Processed_Data/DE_COPD.csv", lfc_threshold)
genes_COPD_all <- genes_COPD$all
genes_COPD_up <- genes_COPD$up
genes_COPD_down <- genes_COPD$down
genes_COPD_full <- genes_COPD$de_genes

genes_ACO <- read_and_filter_genes("Data/Processed_Data/DE_ACO.csv", lfc_threshold)
genes_ACO_all <- genes_ACO$all
genes_ACO_up <- genes_ACO$up
genes_ACO_down <- genes_ACO$down
genes_ACO_full <- genes_ACO$de_genes

# Function to create and save Venn diagram
create_venn_diagram <- function(gene_sets, filename) {
  if (length(gene_sets) >= 2) {
    num_sets <- length(gene_sets)
    
    venn.plot <- venn.diagram(
      x = gene_sets,
      category.names = names(gene_sets),
      output = TRUE,
      filename = NULL,  # To prevent direct file output
      imagetype = "pdf",
      height = 2,  # Adjusted height
      width = 2,   # Adjusted width
      resolution = 300,  # Adjusted resolution
      col = "transparent",
      fill = c('#AFE4DE','#FAA07A','#FCFBCB')[1:num_sets], 
      alpha = 0.4,
      cex = 0.8,
      fontface = "bold",
      fontfamily = "sans",
      cat.cex = 1, # title - ACO and COPD
      cat.fontface = "bold",
      cat.fontfamily = "sans",
      cat.pos = c(-22, 25),  # Adjusted category position
      cat.dist = c(0.06, 0.06),  # Adjusted category distance
      margin = 0.1,  # Adjusted margin
      log = NULL  
    )
    
    # Save the Venn diagram as a PDF
    pdf(filename, width = 2.5, height = 2.5)
    
    grid.draw(venn.plot)
    dev.off()
  } else {
    message("Not enough non-empty gene sets to generate a Venn diagram.")
  }
}

# Create Venn diagrams for all DE genes
gene_sets_all <- list()
if (length(genes_Asthma_all) > 0) {
  gene_sets_all$Asthma <- genes_Asthma_all
}
if (length(genes_COPD_all) > 0) {
  gene_sets_all$COPD <- genes_COPD_all
}
if (length(genes_ACO_all) > 0) {
  gene_sets_all$ACO <- genes_ACO_all
}
create_venn_diagram(gene_sets_all, "Figures/Enrichment_analysis/Venn_All_DE_Genes.pdf")

# Create Venn diagrams for upregulated genes
gene_sets_up <- list()
if (length(genes_Asthma_up) > 0) {
  gene_sets_up$Asthma <- genes_Asthma_up
}
if (length(genes_COPD_up) > 0) {
  gene_sets_up$COPD <- genes_COPD_up
}
if (length(genes_ACO_up) > 0) {
  gene_sets_up$ACO <- genes_ACO_up
}
create_venn_diagram(gene_sets_up, "Figures/Enrichment_analysis/Venn_Upregulated_Genes.pdf")

# Create Venn diagrams for downregulated genes
gene_sets_down <- list()
if (length(genes_Asthma_down) > 0) {
  gene_sets_down$Asthma <- genes_Asthma_down
}
if (length(genes_COPD_down) > 0) {
  gene_sets_down$COPD <- genes_COPD_down
}
if (length(genes_ACO_down) > 0) {
  gene_sets_down$ACO <- genes_ACO_down
}
create_venn_diagram(gene_sets_down, "Figures/Enrichment_analysis/Venn_Downregulated_Genes.pdf")

######## Find Unique KEGG Pathways in ACO and COPD ##########
# Function to perform KEGG pathway enrichment analysis
perform_kegg_analysis <- function(genes_COPD, genes_ACO, species = "hsa") {
  # Perform KEGG pathway enrichment analysis
  enrich_kegg <- function(de_genes_list) {
    gene_list <- de_genes_list$logFC
    names(gene_list) <- de_genes_list$entrezgene_id
    gene_list <- na.omit(gene_list)
    gene_list <- sort(gene_list, decreasing = TRUE)
    kegg <- enrichKEGG(gene = names(gene_list), 
                       organism = species, 
                       keyType = "kegg", 
                       pvalueCutoff = 0.05)
    kegg <- as.data.frame(setReadable(kegg, 'org.Hs.eg.db', keyType = "ENTREZID"))
    return(kegg)
  }

  # Enrichment for COPD DE genes
  kegg_COPD <- enrich_kegg(genes_COPD_full)
  
  # Enrichment for ACO DE genes
  kegg_ACO <- enrich_kegg(genes_ACO_full)
  
  # Function to annotate pathways with upregulated and downregulated genes
  annotate_pathways <- function(kegg_df, de_genes_list) {
    kegg_df$Upregulated_Genes <- sapply(kegg_df$geneID, function(genes) {
      gene_ids <- unlist(strsplit(genes, "/"))
      up_genes <- de_genes_list$hgnc_symbol[de_genes_list$hgnc_symbol %in% gene_ids & de_genes_list$logFC > 0]
      paste(up_genes, collapse = ", ")
    })
    
    kegg_df$Downregulated_Genes <- sapply(kegg_df$geneID, function(genes) {
      gene_ids <- unlist(strsplit(genes, "/"))
      down_genes <- de_genes_list$hgnc_symbol[de_genes_list$hgnc_symbol %in% gene_ids & de_genes_list$logFC < 0]
      paste(down_genes, collapse = ", ")
    })
    
    return(kegg_df)
  }
  
  # Annotate pathways with upregulated and downregulated genes
  kegg_COPD_annotated <- annotate_pathways(kegg_COPD, genes_COPD_full)
  kegg_ACO_annotated <- annotate_pathways(kegg_ACO, genes_ACO_full)
  
  # Identify unique KEGG pathways
  unique_pathways_ACO <- setdiff(kegg_ACO_annotated$Description, kegg_COPD_annotated$Description)
  unique_pathways_COPD <- setdiff(kegg_COPD_annotated$Description, kegg_ACO_annotated$Description)
    
  # # Identify unique KEGG pathways
  # unique_pathways_ACO <- kegg_ACO_annotated[!kegg_ACO_annotated$Description %in% kegg_COPD_annotated$Description, ]
  # unique_pathways_COPD <- kegg_COPD_annotated[!kegg_COPD_annotated$Description %in% kegg_ACO_annotated$Description, ]
    
  # Identify unique KEGG pathways
  unique_pathways_ACO <- kegg_ACO_annotated[!kegg_ACO_annotated$Description %in% kegg_COPD_annotated$Description, 
                                            c("ID", "Description", "geneID", "Upregulated_Genes", "Downregulated_Genes")]
  unique_pathways_COPD <- kegg_COPD_annotated[!kegg_COPD_annotated$Description %in% kegg_ACO_annotated$Description, 
                                              c("ID", "Description", "geneID", "Upregulated_Genes", "Downregulated_Genes")]

  # Return the results as a list
  return(list(
    kegg_COPD = kegg_COPD_annotated,
    kegg_ACO = kegg_ACO_annotated,
    unique_ACO = unique_pathways_ACO,
    unique_COPD = unique_pathways_COPD
  ))
}

# Perform KEGG pathway analysis
kegg_results <- perform_kegg_analysis(genes_COPD_all, genes_ACO_all)

# Save the KEGG pathway results
write.csv(kegg_results$kegg_ACO, "Figures/Enrichment_analysis/KEGG_Pathways_alldegs_ACO.csv", row.names = FALSE)
write.csv(kegg_results$kegg_COPD, "Figures/Enrichment_analysis/KEGG_Pathways_alldegs_COPD.csv", row.names = FALSE)

# Save the unique pathways
write.csv(kegg_results$unique_ACO, "Figures/Enrichment_analysis/Unique_KEGG_Pathways_ACO.csv", row.names = FALSE)
write.csv(kegg_results$unique_COPD, "Figures/Enrichment_analysis/Unique_KEGG_Pathways_COPD.csv", row.names = FALSE)

# Display unique pathways
print("==============Unique KEGG Pathways in ACO:==============")
print(kegg_results$unique_ACO)

print("==============Unique KEGG Pathways in COPD:==============")
print(kegg_results$unique_COPD)

################ Unique GO pathways ########################
# Function to perform GO enrichment analysis
perform_go_analysis <- function(de_genes_list, species = "hsa") {
    
  # fcvals_up <- de_genes_list[de_genes_list$logFC > 0, ]
  fcvals_up <- de_genes_list[de_genes_list$logFC > 0, ]

  # Process upregulated genes (you can replicate this for downregulated genes if needed)
  degenes <- subset(fcvals_up, select=c('hgnc_symbol', 'logFC'))
  degenes$logFC <- degenes$logFC 
  degenes <- degenes[!duplicated(degenes$hgnc_symbol), ]
  colnames(degenes) <- c('ID', 'logFC')

  genelist <- fcvals_up[,c('hgnc_symbol', 'logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')]
  names(genelist)[1] <- "ID"
    
  annot_file <- file.path('Processed_Data', 'Annotations.tsv')
  annot <- read.csv(annot_file, sep = '\t', header = TRUE)
  annot <- annot[!duplicated(annot$ensembl_gene_id), ]
  rownames(annot) <- annot$ensembl_gene_id
  log_with_timestamp(paste("Dimensions of annotations:", paste(dim(annot), collapse = " x ")))
    
  lcpm <- readLines(file.path('Data/Processed_Data', "gene_names.txt"))
  universe <- annot[lcpm, 'hgnc_symbol']
    
  go <- enrichGO(genelist$ID,  
                 OrgDb = "org.Hs.eg.db",
                 keyType = 'SYMBOL',
                 ont = "all",
                 minGSSize = 10, 
                 maxGSSize = 100, 
                 universe = universe, 
                 pvalueCutoff = 0.05,
                readable = TRUE)

  return(go)
}

# Perform GO enrichment analysis for Biological Process (BP)
go_COPD <- perform_go_analysis(genes_COPD_full)
go_ACO <- perform_go_analysis(genes_ACO_full)

# Convert enrichment results to data frames for easier handling
go_COPD_df <- as.data.frame(go_COPD)
go_ACO_df <- as.data.frame(go_ACO)

# # Identify unique GO pathways
# unique_go_ACO <- setdiff(go_ACO_df$Description, go_COPD_df$Description)
# unique_go_COPD <- setdiff(go_COPD_df$Description, go_ACO_df$Description)

# Identify unique GO pathways
unique_go_ACO <- go_ACO_df[!go_ACO_df$Description %in% go_COPD_df$Description, c( "Description", "ONTOLOGY")]
unique_go_COPD <- go_COPD_df[!go_COPD_df$Description %in% go_ACO_df$Description, c( "Description", "ONTOLOGY")]

# Return the results as a list
go_results <- list(
  go_COPD = go_COPD_df,
  go_ACO = go_ACO_df,
  unique_ACO = unique_go_ACO,
  unique_COPD = unique_go_COPD
)

# # Save the GO pathway results
# write.csv(go_results$go_ACO, "Processed_Data/GO_Pathways_ACO.csv", row.names = FALSE)
# write.csv(go_results$go_COPD, "Processed_Data/GO_Pathways_COPD.csv", row.names = FALSE)

# Save the unique pathways
write.csv(go_results$unique_ACO, "Figures/Enrichment_analysis/Unique_GO_Pathways_ACO.csv", row.names = FALSE)
write.csv(go_results$unique_COPD, "Figures/Enrichment_analysis/Unique_GO_Pathways_COPD.csv", row.names = FALSE)

# Display unique pathways
print("==============Unique GO Pathways in ACO:==============")
print(go_results$unique_ACO)

print("==============Unique GO Pathways in COPD:==============")
print(go_results$unique_COPD)

















######## Find Unique KEGG Pathways in ACO and COPD ##########
# # Function to perform KEGG pathway enrichment analysis
# perform_kegg_analysis <- function(genes_COPD, genes_ACO, species = "hsa") {
#   # Perform KEGG pathway enrichment analysis
#   enrich_kegg <- function(de_genes_list) {
#         gene_list <- de_genes_list$logFC
#         names(gene_list) <- de_genes_list$entrezgene_id
#         gene_list <- na.omit(gene_list)
#         gene_list <- sort(gene_list, decreasing = TRUE)
#         kegg <- enrichKEGG(gene = names(gene_list), #gene_list, 
#                                     organism = species, 
#                                     keyType = "kegg", 
#                                     pvalueCutoff = 0.05)
#         kegg <- as.data.frame(setReadable(kegg, 'org.Hs.eg.db', keyType="ENTREZID"))
#         return(kegg)
#   }

#   # Enrichment for COPD DE genes
#   kegg_COPD <- enrich_kegg(genes_COPD_full)
  
#   # Enrichment for ACO DE genes
#   kegg_ACO <- enrich_kegg(genes_ACO_full)
  
#   # Convert enrichment results to data frames for easier handling
#   kegg_COPD_df <- as.data.frame(kegg_COPD)
#   kegg_ACO_df <- as.data.frame(kegg_ACO)
  
#   # Identify unique KEGG pathways
#   unique_pathways_ACO <- setdiff(kegg_ACO_df$Description, kegg_COPD_df$Description)
#   unique_pathways_COPD <- setdiff(kegg_COPD_df$Description, kegg_ACO_df$Description)
    
#   # Identify unique KEGG pathways
#   unique_pathways_ACO <- kegg_ACO_df[!kegg_ACO_df$Description %in% kegg_COPD_df$Description, ]
#   unique_pathways_COPD <- kegg_COPD_df[!kegg_COPD_df$Description %in% kegg_ACO_df$Description, ]

#   # Return the results as a list
#   return(list(
#     kegg_COPD = kegg_COPD_df,
#     kegg_ACO = kegg_ACO_df,
#     unique_ACO = unique_pathways_ACO,
#     unique_COPD = unique_pathways_COPD
#   ))
# }

# # Perform KEGG pathway analysis
# kegg_results <- perform_kegg_analysis(genes_COPD_all, genes_ACO_all)


# # Save the KEGG pathway results
# write.csv(kegg_results$kegg_ACO, "Processed_Data/KEGG_Pathways_ACO.csv", row.names = FALSE)
# write.csv(kegg_results$kegg_COPD, "Processed_Data/KEGG_Pathways_COPD.csv", row.names = FALSE)

# # Save the unique pathways
# write.csv(kegg_results$unique_ACO, "Processed_Data/Unique_KEGG_Pathways_ACO.csv", row.names = FALSE)
# write.csv(kegg_results$unique_COPD, "Processed_Data/Unique_KEGG_Pathways_COPD.csv", row.names = FALSE)

# # Display unique pathways
# print("==============Unique KEGG Pathways in ACO:==============")
# print(kegg_results$unique_ACO)

# print("==============Unique KEGG Pathways in COPD:==============")
# print(kegg_results$unique_COPD)
