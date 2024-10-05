# This script is designed for differential expression analysis and visualization, specifically targeting 
# diseases like Chronic Obstructive Pulmonary Disease (COPD) and Asthma-COPD Overlap (ACO). It integrates 
# several R packages including 'limma', 'edgeR', 'org.Hs.eg.db', 'data.table', and 'biomaRt' to facilitate 
# comprehensive gene expression analysis, from data preparation and normalization to differential expression 
# testing, and further enrichment analysis. The script supports gene annotation, logging actions with timestamps 
# for reproducibility, and offers functionality for both batch effect correction and visualization of results 
# through various plotting methods.

# Features:
# - Data loading and preprocessing, with support for gene annotation using 'org.Hs.eg.db'.
# - Differential expression analysis utilizing 'limma' and 'edgeR'.
# - Batch effect correction to adjust for non-biological variations among samples.
# - Enrichment analysis to identify significantly impacted biological pathways and processes.
# - Custom logging function for tracking script execution progress and package versions.
# - Visualization tools including 'biomaRt' and 'GOplot' for enhanced interpretation of results.

# Usage Notes:
# - Ensure all required libraries are installed and loaded before running the script.
# - Set the `PATH` variable to the directory containing your data files.
# - Customize the `log_with_timestamp` function calls to suit your logging preferences.
# - Adjust the `load_data`, `perform_kegg_analysis`, and `perform_go_analysis` functions according to your 
#   data format and analysis needs.
# - Consider the parameters within each function, such as `lfc_threshold` and `condition`, for tailored analysis.
# - Execute the script within an R environment capable of handling the computational demands of genomic data analysis.

# Dependencies: 
# R version 3.6.0 or later, and the Bioconductor packages: 'limma', 'edgeR', 'org.Hs.eg.db', 'data.table', 'biomaRt', 
# and potentially others based on the specific needs of your analysis.

# Author: Vrushali D. Fangal, 
#         Channing Division of Network Medicine, 
#         Harvard Medical School, 
#         Boston, MA.
# Date: March 22, 2024
# License: MIT license


# Load necessary libraries
library('limma')         # For linear models and statistical tests
library("edgeR")         # For differential expression analysis
library("org.Hs.eg.db")  # For gene annotation
library(data.table)      # For data manipulation
library("biomaRt")

# Define the logging function
log_with_timestamp <- function(message) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(paste("[", timestamp, "] ", message, "\n", sep=""))
}

# Log package versions
log_with_timestamp(paste("limma version:", packageVersion("limma")))
log_with_timestamp(paste("edgeR version:", packageVersion("edgeR")))
log_with_timestamp(paste("org.Hs.eg.db version:", packageVersion("org.Hs.eg.db")))
log_with_timestamp(paste("data.table version:", packageVersion("data.table")))
log_with_timestamp(paste("biomaRt version:", packageVersion("biomaRt")))

# Set paths
# message("Setting up paths...")
log_with_timestamp("Setting up paths...")
PATH <- getwd()

# # Save human annotations 
# log_with_timestamp("Retrieving and saving human annotations...")
# ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# annot <- getBM(attributes =  c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "chromosome_name", "strand", "start_position", "end_position","gene_biotype"),
#               mart=ensembl)
# write.table(annot,"Processed_Data/Annotations.tsv",sep="\t",row.names=FALSE)
# log_with_timestamp(paste("Annotations saved. Number of genes:", dim(annot)[1]))

# Load annotations
log_with_timestamp("Loading annotations...")
annot <- read.csv(paste('Processed_Data', 'Annotations.tsv', sep = .Platform$file.sep), sep = '\t', header = TRUE)
annot <- annot[!duplicated(annot$ensembl_gene_id), ]
rownames(annot) <- annot$ensembl_gene_id
# Concatenate the message and dimensions into a single string
message <- paste("Annotations loaded. Number of genes after removing duplicates:", dim(annot)[1], "\nDimensions:", dim(annot)[1], "*", dim(annot)[2])
log_with_timestamp(message)

# Load expression data
log_with_timestamp("Loading expression data...")
counts <- as.matrix(fread(paste(PATH, "Processed_Data/counts.csv", sep = .Platform$file.sep)), rownames = 1)
message <- paste("Expression data loaded. Dimensions:", dim(counts)[1], "*", dim(counts)[2])
log_with_timestamp(message)

# Load metadata
log_with_timestamp("Loading metadata...")
meta <- read.delim(paste(PATH, "Processed_Data/pheno.tsv", sep = .Platform$file.sep), row.names = 1)
metadata_message <- paste("Metadata loaded. Dimensions:", dim(meta)[1], "*", dim(meta)[2])
log_with_timestamp(metadata_message)

# Subset counts to match metadata
log_with_timestamp("Subsetting counts to match metadata...")
# meta <- meta[meta$smoking_status_P2 == 2, ] ## Subset to Current smokers
counts <- counts[, rownames(meta)]
annot <- annot[rownames(counts),]

# Find the different diseases in this subset
unique_diseases <- unique(meta$Disease)
log_with_timestamp(paste("Unique diseases in former smokers:", paste(unique_diseases, collapse = ", ")))


# Prepare DGEList object
log_with_timestamp("Preparing DGEList object...")
x <- DGEList(counts)
samplenames <- rownames(meta)
group <- factor(meta$Disease)
Batch <- factor(meta$Lc.Batch)
Batch <- gsub("-", "_", Batch)
x$genes <- annot
x$samples$Batch <- Batch

## Use control as reference
group <- relevel(group, ref = "Control")

# Filter lowly expressed genes
log_with_timestamp("Filtering lowly expressed genes...")
# keep.exprs <- filterByExpr(x, min.count = 20, group=group)
keep.exprs <- filterByExpr(x, group=group, min.count = 20, min.total.count = length(samplenames)%/%2 , min.prop = 0.5)
x <- x[keep.exprs,, keep.lib.sizes = FALSE]
message <- paste("Number of genes after filtering:", dim(x)[1])
log_with_timestamp(message)

# Normalize and batch correct data
log_with_timestamp("Normalizing data...") # and batch correcting data...")
x <- calcNormFactors(x, method = "TMM")
# lcpm_bc <- removeBatchEffect(cpm(x, log = TRUE), batch = factor(meta$Lc.Batch))

# Define design matrix for DE analysis
log_with_timestamp("Defining design matrix for DE analysis...")
Gender <- factor(meta$gender)
Currentsmoking <- factor(meta$smoking_status_P2)
ATS_PackYears <- meta$ATS_PackYears_P2
Age_Enroll <- meta$Age_P2
design <- model.matrix(~0 + group + Batch + Age_Enroll + Gender + Currentsmoking + ATS_PackYears)
# design <- model.matrix(~0 + group + Batch + Age_Enroll + Gender + ATS_PackYears)

# Define the design matrix for batch correction
# Log-transform the counts
log_with_timestamp("Log-transforming the counts...")
logCPM <- cpm(x, log = TRUE, prior.count=5) # 'prior.count' adds a small number to avoid log(0)
log_with_timestamp("Log-transformation completed.")

# Batch correction using removeBatchEffect
log_with_timestamp("Performing batch correction with removeBatchEffect...")
bc_design <- model.matrix(~0 + group + Age_Enroll + Gender + Currentsmoking + ATS_PackYears)  # Exclude Batch if it's part of the design matrix for DE analysis
# bc_design <- model.matrix(~0 + group + Age_Enroll + Gender +  ATS_PackYears)  # Exclude Batch if it's part of the design matrix for DE analysis
lcpm_bc <- removeBatchEffect(logCPM, batch = x$samples$Batch, design = bc_design)
write.csv(lcpm_bc, paste(PATH, 'Data/Processed_Data/normalized_batch_corrected_gene_counts.csv', sep = .Platform$file.sep))
log_with_timestamp("Batch correction process completed.")

# ## Batch correction using ComBat
# library(sva)
# # Performing batch correction with ComBat
# log_with_timestamp("Performing batch correction with ComBat...")
# bc_design <- model.matrix(~0 + group + Age_Enroll + Gender + Currentsmoking + ATS_PackYears)  # Exclude Batch if it's part of the design matrix for DE analysis
# # lcpm_bc <- ComBat(dat = logCPM, batch = x$samples$Batch, mod = NULL)
# lcpm_bc <- ComBat(dat = logCPM, batch = x$samples$Batch, mod = bc_design)
# # Saving the batch-corrected data
# log_with_timestamp("Saving batch-corrected data...")
# write.csv(lcpm_bc, paste(PATH, 'Data/Processed_Data/normalized_batch_corrected_gene_counts.csv', sep = .Platform$file.sep))
# log_with_timestamp("Batch correction process completed.")

# Define contrasts
log_with_timestamp("Defining contrasts for DE analysis...")
cm <- makeContrasts(
    COPDvsControl = groupCOPD - groupControl,
    ACOvsControl = groupACO - groupControl,
    AsthmavsControl = groupAsthma - groupControl,
    ACOvsCOPD = groupACO - groupCOPD,
    ACOvsAsthma = groupACO - groupAsthma,
    AsthmavsCOPD = groupAsthma - groupCOPD,
    # CombinedEffect = (groupAsthma - groupControl) + (groupACO - groupControl) - (groupCOPD - groupControl),
    # Effect = groupAsthma + groupACO - groupCOPD - groupControl,
    
#     interaction_term = groupACO - groupAsthma:groupCOPD,
#     ACO.Asthma = (groupACO:groupAsthma),
#     ACO.COPD = (groupACO:groupCOPD),
#     ACOinAsthmaContext = (groupACO + groupACO:groupAsthma) - groupAsthma,  #unique effect of the ACO condition in the presence of asthma
#     ACOinCOPDContext = (groupACO + groupACO:groupCOPD) - groupCOPD,
#     AsthmaWithoutACO = groupAsthma - (groupACO + groupACO:groupAsthma),
#     COPDWithoutACO = groupCOPD - (groupACO + groupACO:groupAsthma),
#     check = groupACO - groupACO:groupCOPD,
    
    # ACOvsAsthmaCOPDInteraction = AsthmaCOPDInteractionACO - AsthmaCOPDInteractionOther,
    # AsthmaCOPDInteraction = groupACO - AsthmaCOPDInteraction,
    # ACOvsAsthmaCOPDInteraction = groupACO - (groupAsthma + groupCOPD)/2,
    # ACOvsAsthmaCOPDInteraction = interaction_term,
    # InteractionEffect = (groupACO - groupAsthma) - (groupCOPD - groupControl),
    # CombinedEffect = (groupCOPD - groupControl) + (groupACO - groupControl) - (groupAsthma - groupControl),
    # ACOvsAsthmaCOPDInteraction = groupACO - (groupAsthma + groupCOPD - groupAsthma:groupCOPD),
    # Alleviation = groupACO - (groupCOPD + groupAsthma - groupControl),
    # interaction_term = groupACO - groupAsthma:groupCOPD,
    # InteractionEffect = groupACO - (groupAsthma + groupCOPD)/2, 
    # AsthmaCOPDInteraction = (groupAsthma + groupCOPD + groupAsthma:groupCOPD) - groupControl,
    # ACO_vs_Asthma_COPD_Interaction = (groupACO - (groupAsthma + groupCOPD) + groupControl),

    levels = design
)

# Perform voom transformation and linear modeling
log_with_timestamp("Performing voom transformation and linear modeling...")
v <- voom(x, design, plot = FALSE)
vfit <- lmFit(v, design)
# lmFitResults <- vfit
vfit <- contrasts.fit(vfit, contrasts = cm)
efit <- eBayes(vfit)

# # Identify batch coefficients (assuming they are the first coefficients)
# batch_coef_indices <- grep("Batch", colnames(design))
# lmFitResults$coefficients[, batch_coef_indices] <- 0
# # Recompute the fitted values without the batch effect
# batch_corrected <- fitted(lmFitResults)
# colnames(batch_corrected) <- colnames(x) 
# # head(batch_corrected)
# write.csv(batch_corrected, paste(PATH, 'Data/Processed_Data/normalized_batch_corrected_gene_counts.csv', sep = .Platform$file.sep))

# Summarize results
log_with_timestamp("Summarizing DE analysis results...")
pvalue = 0.01
logFC_cutoff = 0
# Log fold change cutoff and p-value used
log_with_timestamp(paste("Fold change cutoff used:", logFC_cutoff))  # Replace 0 with your actual fold change cutoff if different
log_with_timestamp(paste("P-value cutoff used:", pvalue))  # Adjust if you use a different p-value cutoff

results <- decideTests(efit, p.value=pvalue, lfc=logFC_cutoff)
summary(results)

# Process each disease compared to control
diseases <- c('Asthma', 'COPD', 'ACO')

for (dis in diseases) {
    # Concatenate the message with the disease name
    message <- paste("Processing", dis, "vs Control")
    log_with_timestamp(message)
    
    cond <- paste(dis, 'vsControl', sep ='')

    down_genes <- summary(results)['Down',cond]
    up_genes <- summary(results)['Up',cond]
    
    res <- topTable(efit, coef=cond, number=Inf, sort.by='p')
    de_genes <- head(res, up_genes + down_genes)
    
    ## Save DE analysis results
    if (nrow(de_genes) > 0) {
        log_with_timestamp(paste("Preparing to save DE genes for", dis))
        name = paste('Data/Processed_Data/DE_', dis, sep ='')
        file_name = paste(name, '.csv', sep ='')
        write.csv(de_genes[,c('entrezgene_id','hgnc_symbol','gene_biotype','logFC','AveExpr','t','P.Value','adj.P.Val','B')], file_name)
        log_with_timestamp(paste("DE genes saved for", dis, "to", file_name))
    } else {
        log_with_timestamp(paste("No DE genes found for", dis, ". File not saved."))
    }
    
    lfc <- 0
    fcvals_up <- de_genes[de_genes$logFC > lfc, ]
    fcvals_down <- de_genes[de_genes$logFC < -lfc, ]
    
    # Log dimensions of upregulated and downregulated genes
    up_gene_count <- paste("Number of upregulated genes for", dis, ":", dim(fcvals_up)[1])
    down_gene_count <- paste("Number of downregulated genes for", dis, ":", dim(fcvals_down)[1])
    log_with_timestamp(up_gene_count)
    log_with_timestamp(down_gene_count)

    de_plot_genes <- lcpm_bc[c(rownames(fcvals_up), rownames(fcvals_down)),]
    # de_plot_genes <- batch_corrected[c(rownames(fcvals_up), rownames(fcvals_down)),]
    # Check if there are rows in de_plot_genes
    if (nrow(de_plot_genes) > 0) {
        log_with_timestamp(paste("Preparing to save expression of DE genes for", dis))
        file_name <- paste('Data/Processed_Data/DE_expression_', dis, '_Control_p.01_fc0.csv', sep = '')
        write.csv(de_plot_genes[, rownames(meta[meta$Disease %in% c('Control', dis), ])], file_name)
        log_with_timestamp(paste("Expression data of DE genes saved for", dis, "to", file_name))
    } else {
        log_with_timestamp(paste("No expression data found for DE genes for", dis, ". File not saved."))
    }
}

# Save gene names from lcpm_bc matrix
gene_names <- rownames(lcpm_bc)
# gene_names <- rownames(batch_corrected)
write.table(gene_names, file = paste(PATH, "Data/Processed_Data/gene_names.txt", sep = .Platform$file.sep), row.names = FALSE, col.names = FALSE, quote = FALSE)

log_with_timestamp("DE analysis completed.")


