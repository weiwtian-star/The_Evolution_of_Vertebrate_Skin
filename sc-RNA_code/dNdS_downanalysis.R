library(Seurat)
library(dplyr)
library(tibble)
library(phytools)
library(ape)
library(ggplot2)

###########################################################
#################expression-weighted mean dN/dS######################
sce<-readRDS('skin.rds')
dNdS_table<-read.table('result_df.txt')

##### function to calculate weighted average per cell #####
#Create function to calculate for each cell the weighted average across all genes expressed
function_average_within_cell <- function(Seurat_obj, ident_of_Seurat_obj, table, value_column_of_table, geneID_column_of_table){
  
  #Create matrix of data and table of metadata for Seurat object
  data <- GetAssayData(Seurat_obj,layer = 'data')
  meta <- Seurat_obj@meta.data
  
  #Remove any potential NA from the table
  table <- table[!is.na(table[, value_column_of_table]), c(geneID_column_of_table, value_column_of_table)]
  rownames(table) <- table[ ,geneID_column_of_table]
  
  #select genes expressed in data object
  data_genes_in_common <- data[intersect(rownames(data), table[ ,geneID_column_of_table]),]
  
  #Get unique values for ident
  meta[,ident_of_Seurat_obj] <- factor(meta[,ident_of_Seurat_obj])
  uniq_ident_of_Seurat_obj <- unique(meta[,ident_of_Seurat_obj])
  
  
  value_by_ident <- lapply(uniq_ident_of_Seurat_obj, function(x){
    y <- rownames(meta[meta[,ident_of_Seurat_obj] == x,]) #get cell IDs that belong to the ident
    matrix <- data_genes_in_common[,y] #subset the data object to these cells
    
    #Calculate the weighted average value per cell
    weighted_average_per_cell <- lapply(colnames(matrix), function(z){
      gene_exp <- matrix[,z]
      gene_exp <- gene_exp[gene_exp !=0] #only select genes expressed in the cell
      gene_exp <- lapply(gene_exp, function(a){return(expm1(a))}) 
      gene_exp_value <- lapply(names(gene_exp), function(n) {return(table[n, 'dnds']*gene_exp[[n]])})
      gene_exp_value <- unlist(gene_exp_value)
      gene_exp_value_ave <- sum(gene_exp_value)/sum(unlist(gene_exp)) 
      return(gene_exp_value_ave)
    })
    
    values_per_ident <- data.frame("weightAve"=unlist(weighted_average_per_cell),"ident"=x)
    return(values_per_ident)
  })
  
  #Create a dataframe containing all values across all idents
  all_values <- Reduce(rbind, value_by_ident) 
  return(all_values)
}

colnames(mam_dnds_omega )<-c('gene','dnds')

##### calculate values per cell using the function #####
dNdS <- function_average_within_cell(sce, "celltype", dNdS_table, "dnds", "gene")




required_packages <- c("SeuratObject", "Matrix", "data.table", "splines")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Please install the following packages: ",
    paste(missing_packages, collapse = ", ")
  )
}


# ---- Internal helper functions ------------------------------------------------

# Read a dN/dS table from a data.frame, matrix, named vector, RDS file, or text
# file. The final internal table always contains the columns `gene` and `dnds`.
.read_dnds_simple <- function(
    dnds,
    gene_col = NULL,
    dnds_col = NULL,
    text_header = TRUE,
    text_sep = ""
) {
  if (is.character(dnds) && length(dnds) == 1L && file.exists(dnds)) {
    if (grepl("\\.rds$", dnds, ignore.case = TRUE)) {
      dnds <- readRDS(dnds)
    } else {
      dnds <- utils::read.table(
        dnds,
        header = text_header,
        sep = text_sep,
        check.names = FALSE,
        stringsAsFactors = FALSE,
        comment.char = "",
        quote = "\""
      )
    }
  }
  
  if (is.numeric(dnds) && !is.null(names(dnds))) {
    dnds <- data.frame(
      gene = names(dnds),
      dnds = as.numeric(dnds),
      stringsAsFactors = FALSE
    )
  }
  
  if (is.matrix(dnds)) {
    dnds <- as.data.frame(dnds, stringsAsFactors = FALSE)
  }
  
  if (!is.data.frame(dnds)) {
    stop("`dnds` must be a data.frame, matrix, named numeric vector, or file path.")
  }
  
  gene_candidates <- c("gene", "Gene", "gene_id", "geneID", "GeneID")
  value_candidates <- c("dnds", "dN_dS", "dN.dS", "omega", "Omega")
  
  if (is.null(gene_col)) {
    gene_col <- intersect(gene_candidates, colnames(dnds))[1L]
  }
  if (is.null(dnds_col)) {
    dnds_col <- intersect(value_candidates, colnames(dnds))[1L]
  }
  
  if (is.null(gene_col) || is.na(gene_col) || !gene_col %in% colnames(dnds)) {
    stop("Gene column not found. Please specify `gene_col`.")
  }
  if (is.null(dnds_col) || is.na(dnds_col) || !dnds_col %in% colnames(dnds)) {
    stop("dN/dS column not found. Please specify `dnds_col`.")
  }
  
  table <- data.table::data.table(
    gene = trimws(as.character(dnds[[gene_col]])),
    dnds = suppressWarnings(as.numeric(dnds[[dnds_col]]))
  )
  
  table <- table[
    !is.na(gene) & gene != "" & is.finite(dnds)
  ]
  
  if (nrow(table) == 0L) {
    stop("No valid dN/dS values remain after filtering.")
  }
  
  # If a gene occurs more than once, use its median dN/dS.
  table <- table[, .(dnds = stats::median(dnds)), by = gene]
  table[]
}


# Match the expression matrix and dN/dS table and return the common data needed
# by both public functions.
.prepare_dnds_simple <- function(
    seurat_object,
    dnds,
    celltype_col,
    assay,
    layer,
    undo_log1p,
    gene_col,
    dnds_col,
    minimum_genes
) {
  if (!inherits(seurat_object, "Seurat")) {
    stop("`seurat_object` must be a Seurat object.")
  }
  
  metadata <- seurat_object[[]]
  if (!celltype_col %in% colnames(metadata)) {
    stop("Metadata column `", celltype_col, "` was not found.")
  }
  
  expression_matrix <- SeuratObject::GetAssayData(
    object = seurat_object,
    assay = assay,
    layer = layer
  )
  expression_matrix <- methods::as(expression_matrix, "dgCMatrix")
  
  # Seurat's normalized `data` layer is normally log1p transformed.
  # Applying expm1 returns it to the linear scale used as the weight.
  if (undo_log1p) {
    expression_matrix@x <- expm1(expression_matrix@x)
  }
  
  if (any(!is.finite(expression_matrix@x)) || any(expression_matrix@x < 0)) {
    stop("The selected expression layer contains invalid or negative values.")
  }
  
  metadata <- metadata[colnames(expression_matrix), , drop = FALSE]
  valid_cells <- !is.na(metadata[[celltype_col]]) & metadata[[celltype_col]] != ""
  
  expression_matrix <- expression_matrix[, valid_cells, drop = FALSE]
  ident <- as.character(metadata[[celltype_col]][valid_cells])
  
  dnds_table <- .read_dnds_simple(
    dnds = dnds,
    gene_col = gene_col,
    dnds_col = dnds_col
  )
  
  common_genes <- intersect(rownames(expression_matrix), dnds_table$gene)
  if (length(common_genes) < minimum_genes) {
    stop(
      "Only ", length(common_genes),
      " genes were shared between the expression matrix and dN/dS table."
    )
  }
  
  expression_matrix <- expression_matrix[common_genes, , drop = FALSE]
  gene_dnds <- dnds_table$dnds[match(common_genes, dnds_table$gene)]
  mean_expression <- Matrix::rowMeans(expression_matrix)
  
  keep_gene <- is.finite(gene_dnds) & is.finite(mean_expression) & mean_expression > 0
  expression_matrix <- expression_matrix[keep_gene, , drop = FALSE]
  gene_dnds <- gene_dnds[keep_gene]
  mean_expression <- mean_expression[keep_gene]
  
  if (nrow(expression_matrix) < minimum_genes) {
    stop("Too few matched and expressed genes remain after filtering.")
  }
  
  list(
    expression = expression_matrix,
    gene_dnds = gene_dnds,
    mean_expression = mean_expression,
    ident = ident
  )
}


# Calculate an expression-weighted mean gene value for every cell.
.weighted_score_simple <- function(expression_matrix, gene_value) {
  denominator <- Matrix::colSums(expression_matrix)
  numerator <- as.numeric(Matrix::crossprod(gene_value, expression_matrix))
  
  score <- numerator / denominator
  score[denominator <= 0] <- NA_real_
  score
}


# ---- Public function 1: expression-bin correction ----------------------------

calculate_dnds_expression_bins <- function(
    seurat_object,
    dnds,
    celltype_col = "celltype",
    assay = "RNA",
    layer = "data",
    undo_log1p = TRUE,
    gene_col = NULL,
    dnds_col = NULL,
    n_bins = 5L,
    require_all_bins = TRUE,
    minimum_genes = 20L
) {
  prepared <- .prepare_dnds_simple(
    seurat_object = seurat_object,
    dnds = dnds,
    celltype_col = celltype_col,
    assay = assay,
    layer = layer,
    undo_log1p = undo_log1p,
    gene_col = gene_col,
    dnds_col = dnds_col,
    minimum_genes = minimum_genes
  )
  
  n_bins <- as.integer(n_bins)
  n_genes <- nrow(prepared$expression)
  if (n_bins < 2L || n_genes < n_bins) {
    stop("`n_bins` must be at least 2 and no larger than the number of genes.")
  }
  
  # Rank genes by their mean expression. Gene ID is used to break ties so that
  # the bin assignment is reproducible.
  gene_order <- order(prepared$mean_expression, rownames(prepared$expression))
  bin_number <- integer(n_genes)
  
  # This exact rank formula is used by the original analysis. It assigns ranks
  # 1,...,n to n approximately equal-sized bins while keeping the first-ranked
  # gene in Q1 and the last-ranked gene in Qn.
  bin_number[gene_order] <-
    floor((seq_along(gene_order) - 1L) * n_bins / n_genes) + 1L
  
  # Calculate one expression-weighted dN/dS score per cell within each bin.
  score_by_bin <- vapply(
    seq_len(n_bins),
    FUN.VALUE = numeric(ncol(prepared$expression)),
    FUN = function(current_bin) {
      keep <- bin_number == current_bin
      .weighted_score_simple(
        expression_matrix = prepared$expression[keep, , drop = FALSE],
        gene_value = prepared$gene_dnds[keep]
      )
    }
  )
  
  # Give every expression bin equal weight in the final score.
  if (require_all_bins) {
    final_score <- rowMeans(score_by_bin)
  } else {
    final_score <- rowMeans(score_by_bin, na.rm = TRUE)
    final_score[!is.finite(final_score)] <- NA_real_
  }
  
  data.frame(
    weightAve = as.numeric(final_score),
    ident = prepared$ident,
    stringsAsFactors = FALSE
  )
}


# ---- Public function 2: expression-covariate regression ----------------------

calculate_dnds_expression_covariate <- function(
    seurat_object,
    dnds,
    celltype_col = "celltype",
    assay = "RNA",
    layer = "data",
    undo_log1p = TRUE,
    gene_col = NULL,
    dnds_col = NULL,
    spline_df = 3L,
    minimum_genes = 20L
) {
  prepared <- .prepare_dnds_simple(
    seurat_object = seurat_object,
    dnds = dnds,
    celltype_col = celltype_col,
    assay = assay,
    layer = layer,
    undo_log1p = undo_log1p,
    gene_col = gene_col,
    dnds_col = dnds_col,
    minimum_genes = minimum_genes
  )
  
  # Estimate the dN/dS normally expected from each gene's mean expression.
  expression_model <- stats::lm(
    prepared$gene_dnds ~ splines::ns(
      log1p(prepared$mean_expression),
      df = spline_df
    )
  )
  
  # Keep the part of dN/dS not predicted by expression. Adding the overall mean
  # places the adjusted values back on an interpretable dN/dS scale.
  adjusted_gene_dnds <- mean(prepared$gene_dnds) +
    stats::residuals(expression_model)
  
  final_score <- .weighted_score_simple(
    expression_matrix = prepared$expression,
    gene_value = adjusted_gene_dnds
  )
  
  data.frame(
    weightAve = as.numeric(final_score),
    ident = prepared$ident,
    stringsAsFactors = FALSE
  )
}


# ---- Example ------------------------------------------------------------------
#
# source("two_expression_corrected_dnds_functions.R")
# sce <- readRDS("mam.rds")
# dnds <- readRDS("mam_dnds_omega.rds")
#
# # Method 1: final dN/dS values after expression-bin correction
# dnds_bins <- calculate_dnds_expression_bins(
#   seurat_object = sce,
#   dnds = dnds,
#   celltype_col = "celltype"
# )
#
# # Method 2: final dN/dS values after expression-covariate correction
# dnds_covariate <- calculate_dnds_expression_covariate(
#   seurat_object = sce,
#   dnds = dnds,
#   celltype_col = "celltype"
# )
#
# utils::write.csv(dnds_bins, "dNdS_mammalia_bins.csv", row.names = FALSE)
# utils::write.csv(
#   dnds_covariate,
#   "dNdS_mammalia_covariate.csv",
#   row.names = FALSE
# )
#
# # Boxplots use ident on the x-axis and weightAve on the y-axis.
# boxplot(weightAve ~ ident, data = dnds_bins, las = 2)
# boxplot(weightAve ~ ident, data = dnds_covariate, las = 2)











