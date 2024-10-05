# This script performs KEGG pathway analysis and visualization for differential expression data
# in the context of COPD and ACO using the gage and pathview packages. It includes data preparation,
# mapping of ENSEMBL gene IDs to Entrez IDs, execution of GAGE pathway analysis, and generation of
# KEGG pathway plots with customized color settings for visual differentiation of gene expression.
#
# Prerequisites:
# - AnnotationDbi: Interface to Bioconductor annotation databases and resources.
# - org.Hs.eg.db: Genome wide annotation for Human, primarily based on mapping using Entrez Gene IDs.
# - pathview: A tool for pathway-based data integration and visualization.
# - gage: Generally Applicable Gene-set Enrichment for pathway analysis.
# - gageData: Data package for gage, includes sample data for demonstration.
#
# The script first loads necessary libraries and datasets for KEGG pathway sets and significant metabolic
# indices. It defines a logging function to timestamp script actions and logs the versions of the used
# packages. The main functionality is encapsulated in the `plot_kegg_pathway` function, which takes a data file
# containing differential expression results, a KEGG pathway ID, and color settings for the plot. It processes
# the data, maps gene IDs, conducts GAGE analysis, and plots the KEGG pathway highlighting differential
# expression in the context of COPD and ACO. Finally, it executes the pathway plotting for specified
# pathways using predefined color schemes for each condition.
#
# Usage:
# - Define color settings for specific conditions (COPD, ACO).
# - Call `plot_kegg_pathway` with the differential expression data file, KEGG pathway ID, and color settings.
#
# Note: Ensure the differential expression data file is correctly formatted and accessible at the specified path.
#
# Author: Vrushali D. Fangal, 
#         Channing Division of Network Medicine, 
#         Harvard Medical School, 
#         Boston, MA.
# Date: March 22, 2024
# License: MIT license

# Load necessary libraries
library("AnnotationDbi")
library("org.Hs.eg.db")
library("pathview")
library("gage")
library(gageData)

data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]

# Define the logging function
log_with_timestamp <- function(message) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(paste("[", timestamp, "] ", message, "\n", sep=""))
}

# Logging versions of the packages used
log_with_timestamp(paste("AnnotationDbi version:", packageVersion("AnnotationDbi")))
log_with_timestamp(paste("org.Hs.eg.db version:", packageVersion("org.Hs.eg.db")))
log_with_timestamp(paste("pathview version:", packageVersion("pathview")))
log_with_timestamp(paste("gage version:", packageVersion("gage")))
log_with_timestamp(paste("gage data version:", packageVersion("gageData")))

# Function to perform KEGG pathway analysis and plotting
plot_kegg_pathway <- function(data_file, kegg_id, colors) {
    
    # Logging data loading
    log_with_timestamp(paste("Loading differential expression results from", data_file))
    
    
    # Load differential expression results
    res <- read.csv(data_file, row.names = 1)
    
    # Map ENSEMBL gene IDs to Entrez IDs
    log_with_timestamp("Mapping ENSEMBL IDs to Entrez IDs")
    res$entrez <- mapIds(org.Hs.eg.db,
                         keys = row.names(res), 
                         column = "ENTREZID",
                         keytype = "ENSEMBL",
                         multiVals = "first")

    # Prepare fold changes for GAGE analysis
    foldchanges <- res$logFC
    names(foldchanges) <- res$entrez
    
    # GAGE analysis
    log_with_timestamp("Performing GAGE analysis")
    keggres <- gage(foldchanges, gsets = kegg.sets.hs, same.dir = TRUE)
    
    # Pathview analysis
    log_with_timestamp(paste("Plotting KEGG pathway:", kegg_id))
    pathview(gene.data = foldchanges*1000, pathway.id = kegg_id, species = "hsa",
             high = list(gene = colors$high, cpd = colors$cpd_high),
             low = list(gene = colors$low, cpd = colors$cpd_low),
             mid = list(gene = colors$mid, cpd = colors$cpd_mid),
             kegg.native = T)
    
    log_with_timestamp("KEGG pathway plot completed")
}

# Define color settings for COPD
colors_copd <- list(high = '#A6CDC7', cpd_high = 'blue',
                    low = '#e89fb4', cpd_low = 'yellow',
                    mid = 'gray88', cpd_mid = 'yellow')

# Define color settings for ACO
colors_aco <- colors_copd 
            # list(high = '#a5d732', cpd_high = 'green',
            #        low = '#e49999', cpd_low = 'pink',
            #        mid = 'gray88', cpd_mid = 'yellow')

# Plot KEGG pathway for COPD
# Necroptosis --> hsa04217
# Neutrophil extracellular trap - hsa04613
# plot_kegg_pathway(data_file = 'Data/Processed_Data/DE_COPD.csv', kegg_id = 'hsa04064', colors = colors_copd)
plot_kegg_pathway(data_file = 'Data/Processed_Data/DE_COPD.csv', kegg_id = 'hsa04217', colors = colors_copd)
plot_kegg_pathway(data_file = 'Data/Processed_Data/DE_COPD.csv', kegg_id = 'hsa04613', colors = colors_copd)

## Plot ACO
plot_kegg_pathway(data_file = 'Data/Processed_Data/DE_ACO.csv', kegg_id = 'hsa04066', colors = colors_copd) # HIF-1 signaling pathway
# plot_kegg_pathway(data_file = 'Data/Processed_Data/DE_ACO.csv', kegg_id = 'hsa00051', colors = colors_copd) # Fructose and mannose metabolism
plot_kegg_pathway(data_file = 'Data/Processed_Data/DE_ACO.csv', kegg_id = 'hsa00520', colors = colors_copd) # Amino sugar and nucleotide sugar metabolism
plot_kegg_pathway(data_file = 'Data/Processed_Data/DE_ACO.csv', kegg_id = 'hsa01250', colors = colors_copd) # Biosynthesis of nucleotide sugars






