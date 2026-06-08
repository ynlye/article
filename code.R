# ============================================================================
# Integrated Analysis Pipeline: Differential Expression, WGCNA, Enrichment,
# Mendelian Randomization, Machine Learning, and Immune Correlation
# ============================================================================

# Load required libraries
library(DESeq2)
library(WGCNA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(TwoSampleMR)
library(glmnet)
library(e1071)
library(caret)
library(GSVA)
library(GSEABase)
library(limma)
library(dplyr)

# ============================================================================
# 1. Differential Expression Analysis (DESeq2)
# ============================================================================
run_deseq2 <- function(count_matrix, col_data, group_col, ref_level, treat_level) {
  # count_matrix: gene x sample count matrix
  # col_data: data frame with sample info, rownames = sample names
  # group_col: column name in col_data containing group variable
  # ref_level: control group name
  # treat_level: treatment group name
  
  col_data$group <- factor(col_data[[group_col]], levels = c(ref_level, treat_level))
  dds <- DESeqDataSetFromMatrix(countData = round(count_matrix), colData = col_data, design = ~ group)
  dds <- dds[rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("group", treat_level, ref_level))
  res <- res[order(res$pvalue), ]
  
  # Add change column
  logFC_cutoff <- 0.5
  res_df <- as.data.frame(res)
  res_df$change <- as.factor(ifelse(res_df$pvalue < 0.05 & abs(res_df$log2FoldChange) > logFC_cutoff,
                                    ifelse(res_df$log2FoldChange > logFC_cutoff, 'UP', 'DOWN'), 'NOT'))
  sig_diff <- subset(res_df, pvalue < 0.05 & abs(log2FoldChange) >= logFC_cutoff)
  
  return(list(all_results = res_df, sig_diff = sig_diff))
}

# ============================================================================
# 2. WGCNA: Network Construction and Module Detection
# ============================================================================
run_wgcna <- function(expr_matrix, clinical_data, power_range = 1:20, min_module_size = 50, mediss_thres = 0.3) {
  # expr_matrix: samples x genes matrix (transposed from usual)
  # clinical_data: data frame with sample traits, rownames matching expr_matrix
  
  # Choose soft threshold
  sft <- pickSoftThreshold(expr_matrix, powerVector = power_range, verbose = 0)
  softPower <- sft$powerEstimate
  if (is.na(softPower)) softPower <- 6  # default if no estimate
  
  # Adjacency and TOM
  adj <- adjacency(expr_matrix, power = softPower)
  TOM <- TOMsimilarity(adj)
  dissTOM <- 1 - TOM
  geneTree <- hclust(as.dist(dissTOM), method = "average")
  
  # Dynamic module cutting
  dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2,
                               pamRespectsDendro = FALSE, minClusterSize = min_module_size)
  dynamicColors <- labels2colors(dynamicMods)
  
  # Merge similar modules
  MEList <- moduleEigengenes(expr_matrix, colors = dynamicColors)
  MEs <- MEList$eigengenes
  MEDiss <- 1 - cor(MEs)
  METree <- hclust(as.dist(MEDiss), method = "average")
  merge <- mergeCloseModules(expr_matrix, dynamicColors, cutHeight = mediss_thres, verbose = 0)
  mergedColors <- merge$colors
  mergedMEs <- merge$newMEs
  
  # Module-trait correlation
  moduleTraitCor <- cor(mergedMEs, clinical_data, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(expr_matrix))
  
  # Export module genes
  module_genes <- list()
  unique_modules <- unique(mergedColors)
  for (mod in unique_modules) {
    genes <- colnames(expr_matrix)[mergedColors == mod]
    module_genes[[mod]] <- genes
  }
  
  return(list(module_colors = mergedColors, module_eigengenes = mergedMEs,
              module_trait_cor = moduleTraitCor, module_trait_pval = moduleTraitPvalue,
              module_genes = module_genes))
}

# ============================================================================
# 3. Gene Set Enrichment Analysis (GO and KEGG)
# ============================================================================
run_enrichment <- function(gene_list, org_db = org.Hs.eg.db, ont = "all", pvalue_cutoff = 0.05) {
  # gene_list: vector of gene symbols
  # Convert symbols to Entrez IDs
  entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db)
  genes_entrez <- entrez$ENTREZID
  
  # GO enrichment
  go <- enrichGO(gene = genes_entrez, OrgDb = org_db, ont = ont,
                 pvalueCutoff = pvalue_cutoff, qvalueCutoff = 1, readable = TRUE)
  go_df <- as.data.frame(go)
  
  # KEGG enrichment
  kegg <- enrichKEGG(gene = genes_entrez, organism = "hsa",
                     pvalueCutoff = pvalue_cutoff, qvalueCutoff = 1)
  kegg_df <- as.data.frame(kegg)
  # Map back gene symbols
  if (nrow(kegg_df) > 0) {
    kegg_df$geneID_symbol <- sapply(strsplit(kegg_df$geneID, "/"), function(x) {
      paste(entrez$SYMBOL[match(x, entrez$ENTREZID)], collapse = "/")
    })
  }
  
  return(list(go = go_df, kegg = kegg_df))
}

# ============================================================================
# 4. Mendelian Randomization (TwoSampleMR)
# ============================================================================
run_mr_analysis <- function(exposure_ids, outcome_id, pval_thresh = 5e-6, r2_thresh = 0.001, kb = 10000) {
  # exposure_ids: vector of exposure IDs (e.g., eQTL IDs)
  # outcome_id: outcome GWAS ID
  # Extract instruments for each exposure
  exposure_dat_list <- list()
  for (id in exposure_ids) {
    dat <- tryCatch({
      extract_instruments(outcomes = id, clump = TRUE, p1 = pval_thresh, r2 = r2_thresh, kb = kb)
    }, error = function(e) NULL)
    if (!is.null(dat)) exposure_dat_list[[id]] <- dat
  }
  exposure_dat <- bind_rows(exposure_dat_list)
  
  # Extract outcome data
  outcome_dat <- extract_outcome_data(snps = exposure_dat$SNP, outcomes = outcome_id)
  
  # Harmonise
  harmonised <- harmonise_data(exposure_dat, outcome_dat, action = 2)
  
  # MR analysis
  mr_results <- mr(harmonised)
  
  # Heterogeneity
  heterogeneity <- mr_heterogeneity(harmonised)
  
  # Pleiotropy (MR-Egger intercept)
  pleiotropy <- mr_pleiotropy_test(harmonised)
  
  # Steiger test for directionality
  steiger <- directionality_test(harmonised)
  
  # Leave-one-out sensitivity (optional, not plotting)
  leaveoneout <- mr_leaveoneout(harmonised)
  
  return(list(mr_results = mr_results, heterogeneity = heterogeneity,
              pleiotropy = pleiotropy, steiger = steiger, leaveoneout = leaveoneout,
              harmonised_data = harmonised))
}

# ============================================================================
# 5. Machine Learning: LASSO and SVM-RFE for Feature Selection
# ============================================================================
run_feature_selection <- function(x_matrix, y_factor, nfolds = 10, svm_sizes = NULL) {
  # x_matrix: feature matrix (samples x features)
  # y_factor: binary factor outcome
  
  # LASSO with cross-validation
  cv_lasso <- cv.glmnet(x_matrix, y_factor, family = "binomial", alpha = 1, nfolds = nfolds)
  lasso_coef <- coef(cv_lasso, s = "lambda.min")
  lasso_genes <- rownames(lasso_coef)[which(lasso_coef != 0)][-1]  # remove intercept
  
  # SVM-RFE
  if (is.null(svm_sizes)) svm_sizes <- c(1:ncol(x_matrix))
  ctrl <- rfeControl(functions = caretFuncs, method = "cv", number = nfolds)
  svm_rfe <- rfe(x_matrix, y_factor, sizes = svm_sizes, rfeControl = ctrl)
  svm_genes <- svm_rfe$optVariables
  
  # Intersection
  common_genes <- intersect(lasso_genes, svm_genes)
  
  return(list(lasso_genes = lasso_genes, svm_genes = svm_genes, common_genes = common_genes,
              lasso_cv = cv_lasso, svm_rfe = svm_rfe))
}

# ============================================================================
# 6. ssGSEA for Immune Infiltration Scoring
# ============================================================================
run_ssgsea <- function(expr_matrix, gene_set_file) {
  # expr_matrix: genes x samples expression matrix (log-normalized recommended)
  # gene_set_file: GMT file containing gene sets (e.g., immune signatures)
  
  geneSet <- getGmt(gene_set_file, geneIdType = SymbolIdentifier())
  ssgseaPar <- ssgseaParam(expr_matrix, geneSet)
  ssgseaScores <- gsva(ssgseaPar)
  # Normalize to [0,1]
  normalize <- function(x) (x - min(x)) / (max(x) - min(x))
  scores_norm <- apply(ssgseaScores, 1, normalize)
  return(t(scores_norm))
}

# ============================================================================
# 7. Correlation between Immune Scores and Target Genes
# ============================================================================
compute_immune_correlations <- function(gene_expression_matrix, immune_scores_matrix, method = "spearman") {
  # gene_expression_matrix: samples x genes
  # immune_scores_matrix: samples x immune cell types
  common_samples <- intersect(rownames(gene_expression_matrix), rownames(immune_scores_matrix))
  genes <- gene_expression_matrix[common_samples, , drop = FALSE]
  immune <- immune_scores_matrix[common_samples, , drop = FALSE]
  
  results <- list()
  for (gene in colnames(genes)) {
    for (cell in colnames(immune)) {
      ct <- cor.test(genes[, gene], immune[, cell], method = method)
      results <- c(results, list(data.frame(
        Gene = gene, Immune_Cell = cell,
        Correlation = ct$estimate, P_value = ct$p.value
      )))
    }
  }
  cor_df <- do.call(rbind, results)
  return(cor_df)
}

# ============================================================================
# 8. GSEA on Ranked Gene List (Correlation-based)
# ============================================================================
run_gsea <- function(ranked_gene_list, gmt_file, pvalue_cutoff = 0.05) {
  # ranked_gene_list: named numeric vector (gene symbol -> score)
  # gmt_file: GMT file for pathways
  gmt <- read.gmt(gmt_file)
  set.seed(1)
  gsea_res <- GSEA(ranked_gene_list, TERM2GENE = gmt, pvalueCutoff = pvalue_cutoff)
  return(as.data.frame(gsea_res))
}

# ============================================================================
# Example usage (commented out)
# ============================================================================
# # Differential expression
# deg_results <- run_deseq2(count_matrix = my_counts, col_data = sample_info,
#                           group_col = "condition", ref_level = "Control", treat_level = "Treat")
# sig_genes <- rownames(deg_results$sig_diff)
# 
# # WGCNA
# wgcna_results <- run_wgcna(expr_matrix = t(normalized_expr), clinical_data = traits)
# 
# # Enrichment
# enrich_results <- run_enrichment(gene_list = sig_genes)
# 
# # MR
# mr_results <- run_mr_analysis(exposure_ids = eqtl_ids, outcome_id = "ebi-a-GCST90029022")
# 
# # Feature selection
# fs_results <- run_feature_selection(x_matrix = train_x, y_factor = train_y)
# 
# # ssGSEA
# immune_scores <- run_ssgsea(expr_matrix = log2expr, gene_set_file = "immune.gmt")
# 
# # Correlation
# immune_cor <- compute_immune_correlations(gene_expression_matrix = t(gene_expr),
#                                           immune_scores_matrix = immune_scores)