#=============================================================================
#
#    ExIR
#
#=============================================================================

#' Experimental data-based Integrated Ranking
#'
#' This function runs the Experimental data-based Integrated Ranking (ExIR)
#' model for the classification and ranking of top candidate features. The input
#' data could come from any type of experiment such as transcriptomics and proteomics.
#' A shiny app has also been developed for Running the ExIR model, visualization of its results as well as computational
#' simulation of knockout and/or up-regulation of its top candidate outputs, which is accessible using
#' the `influential::runShinyApp("ExIR")` command.
#' You can also access the shiny app online at https://influential.erc.monash.edu/.
#' @param Desired_list (Optional) A character vector of your desired features. This vector could be, for
#' instance, a list of features obtained from cluster analysis, time-course analysis,
#' or a list of dysregulated features with a specific sign.
#' @param Diff_data A dataframe of all significant differential/regression data and their
#' statistical significance values (p-value/adjusted p-value). Note that the differential data
#' should be in the log fold-change (log2FC) format.
#' You may have selected a proportion of the differential data as the significant ones according
#' to your desired thresholds. A function, named \code{\link[influential]{diff_data.assembly}}, has also been
#' provided for the convenient assembling of the Diff_data dataframe.
#' @param Diff_value An integer vector containing the column number(s) of the differential
#' data in the Diff_data dataframe. The differential data could result from any type of
#' differential data analysis. One example could be the fold changes (FCs) obtained from differential
#' expression analyses. The user may provide as many differential data as he/she wish.
#' @param Regr_value (Optional) An integer vector containing the column number(s) of the regression
#' data in the Diff_data dataframe. The regression data could result from any type of regression
#' data analysis or other analyses such as time-course data analyses that are based on regression models.
#' @param Sig_value An integer vector containing the column number(s) of the significance values (p-value/adjusted p-value) of
#' both differential and regression data (if provided). Providing significance values for the regression data is optional.
#' @param Exptl_data Experimental data used by the ExIR model. This can be a data frame, tibble,
#' matrix, sparse matrix such as a \code{dgCMatrix}, or a Seurat object. For non-Seurat inputs,
#' the expected orientation is controlled by \code{Exptl_data_orientation}. By default,
#' features/genes are expected to be in rows and samples/cells in columns, which is the usual
#' omics layout. Internally, ExIR converts the data to its required analysis format, with
#' samples/cells in rows and features in columns.
#' @param Exptl_data_type Character string specifying the experimental data type. One of
#' \code{"bulk"} or \code{"sc"}. This is used for data-type checks, optional normalization,
#' and pseudo-sampling. For \code{"bulk"}, the input expression data may be either already
#' normalized/log-transformed or raw count-like data when \code{normalize = TRUE}. For
#' \code{"sc"}, raw counts are recommended and raw counts are required when
#' \code{pseudo_sample = TRUE}.
#' @param condition A character string or character/factor vector specifying the sample/cell
#' conditions. If a single character string is supplied, it is interpreted as the name of the
#' condition column/row in \code{Exptl_data}, or as the name of a metadata column when
#' \code{Exptl_data} is a Seurat object. If a vector is supplied, it must have the same length
#' and order as the samples/cells in \code{Exptl_data}. Default is \code{"condition"}.
#' @param Exptl_data_orientation Character string specifying the orientation of non-Seurat
#' \code{Exptl_data}. One of \code{"features_rows"} or \code{"samples_rows"}. If
#' \code{"features_rows"}, features are rows and samples/cells are columns. If
#' \code{"samples_rows"}, samples/cells are rows and features are columns. Default is
#' \code{"features_rows"}.
#' @param assay Character string specifying the assay to use when \code{Exptl_data} is a Seurat
#' object. Default is \code{"RNA"}.
#' @param layer Character string specifying the assay layer to use when \code{Exptl_data} is a
#' Seurat object. For pseudo-sampling of single-cell data, this should usually be a raw-count
#' layer such as \code{"counts"}. Default is \code{"counts"}.
#' @param normalize Logical; whether to normalize count-like input data using TMM normalization
#' followed by logCPM transformation with \pkg{edgeR}. Default is \code{FALSE}. For
#' \code{Exptl_data_type = "bulk"}, this can be used when raw bulk RNA-seq count-like data are
#' supplied. This normalization strategy is appropriate for many bulk omics count datasets,
#' especially bulk RNA-seq, but users should confirm that TMM/logCPM normalization is suitable
#' for their specific data modality. If the data modality requires a different normalization
#' strategy, users should pre-normalize their data and set \code{normalize = FALSE}. For
#' \code{Exptl_data_type = "sc"}, normalization is automatically applied after pseudo-bulking
#' when \code{pseudo_sample = TRUE}. If \code{Exptl_data_type = "sc"} and
#' \code{pseudo_sample = FALSE}, users should provide pre-normalized single-cell data and keep
#' \code{normalize = FALSE}.
#' @param pseudo_sample Logical; whether to perform pseudo-sampling before running ExIR.
#' Pseudo-sampling is recommended when the number of cells/samples is greater than 500 or when
#' computational resources are limited. For bulk data, pseudo-sampling averages normalized
#' log-expression values within non-overlapping condition-specific groups. For single-cell data,
#' pseudo-sampling sums raw counts within non-overlapping condition-specific groups, followed by
#' TMM normalization and logCPM transformation using \pkg{edgeR}. Default is \code{FALSE}.
#' @param pseudo_samples_per_group Integer specifying the target number of pseudo-samples to
#' generate per condition group when \code{pseudo_sample = TRUE}. For example, if one condition
#' contains 500 cells/samples and \code{pseudo_samples_per_group = 100}, each pseudo-sample will
#' contain 5 cells/samples. If another condition contains 536 cells/samples, 64 pseudo-samples
#' will contain 5 cells/samples and 36 pseudo-samples will contain 6 cells/samples. Default is
#' \code{100}.
#' @param Exptl_data_size_check Logical; whether to check the number of input samples/cells and,
#' in interactive sessions, prompt the user to consider pseudo-sampling when more than 500
#' samples/cells are provided and \code{pseudo_sample = FALSE}. In non-interactive sessions,
#' a message is shown and the function continues. Default is \code{TRUE}.  
#' @param feature_filter Logical; whether to apply conservative feature filtering before
#' running RF, PCA and correlation analysis. This filter is not a highly variable gene
#' filter. It removes only features with insufficient expression/prevalence or essentially
#' zero variance, which are unlikely to produce reliable correlations. Default is \code{TRUE}.
#' @param min_feature_prevalence Integer or \code{NULL}; minimum number of samples/cells/
#' pseudo-samples in which a feature must be non-zero to be retained. If \code{NULL}, an
#' adaptive conservative threshold is used based on \code{Exptl_data_type}, pseudo-sampling,
#' and sample size.
#' @param min_feature_total Numeric or \code{NULL}; minimum total abundance/count/expression
#' support required for a feature to be retained. If \code{NULL}, this criterion is not used.
#' For raw count-like data, values such as 10 or 20 may be useful. For normalized/log-scale
#' data, users should usually leave this as \code{NULL}.
#' @param min_feature_variance Numeric; minimum variance required for a feature to be retained.
#' This is intended only to remove zero-variance or near-zero-variance features, not to perform
#' HVG selection. Default is \code{1e-12}.
#' @param always_keep_diff_features Logical; whether to always retain features present in
#' \code{Diff_data} and \code{Desired_list}, even if they fail the conservative expression/
#' prevalence filters. This helps preserve candidate differential features while still reducing
#' uninformative non-DE background features. Default is \code{TRUE}.
#' @param cor_thresh_method A character string indicating the method for filtering the correlation results, either
#' "mr" (default; Mutual Rank) or "cor.coefficient".
#' @param mr An integer determining the threshold of mutual rank for the selection of correlated features (default is 20). Note that
#' higher mr values considerably increase the computation time.
#' @param r The threshold of Spearman correlation coefficient for the selection of correlated features (default is 0.5).
#' @param max.connections The maximum number of connections to be included in the association network.
#' Higher max.connections might increase the computation time, cost, and accuracy of the results (default is 50,000).
#' @param alpha The threshold of the statistical significance (p-value) used throughout the entire model (default is 0.05)
#' @param num_trees Number of trees to be used for the random forests classification (supervised machine learning). Default is set to 500.
#' @param mtry Number of features to possibly split at in each node. Default is the (rounded down) square root of the
#' number of variables. Alternatively, a single argument function returning an integer, given the number of independent variables.
#' @param num_permutations Number of permutations to be used for computation of the statistical significance (p-values) of
#' the importance scores resulted from random forests classification (default is 50).
#' @param inf_const The constant value to be multiplied by the maximum absolute value of differential (logFC)
#' values for the substitution with infinite differential values. This results in noticeably high biomarker values for features
#' with infinite differential values compared with other features. Having said that, the user can still use the
#' biomarker rank to compare all of the features. This parameter is ignored if no infinite value
#' is present within Diff_data. However, this is used in the case of sc-seq experiments where some genes are uniquely
#' expressed in a specific cell-type and consequently get infinite differential values. Note that the sign of differential
#' value is preserved (default is 10^10).
#' @param ncores Integer; the number of cores to be used for parallel processing. If ncores == "default" (default), the number of 
#' cores to be used will be the max(number of available cores) - 1. We recommend leaving ncores argument as is (ncores = "default").
#' @param seed The seed to be used for all of the random processes throughout the model (default is 1234).
#' @param verbose Logical; whether to display formatted progress messages and a progress bar
#' using \pkg{cli}. If \code{TRUE}, ExIR reports the major analysis stages, selected
#' warnings, and a final output summary. Default is \code{TRUE}.
#' @return A list of one graph and one to four tables including:
#'
#' - Driver table: Top candidate drivers
#'
#' - DE-mediator table: Top candidate differentially expressed/abundant mediators
#'
#' - nonDE-mediator table: Top candidate non-differentially expressed/abundant mediators
#'
#' - Biomarker table: Top candidate biomarkers
#'
#' The number of returned tables depends on the input data and specified arguments.
#' @aliases ExIR
#' @keywords exir
#' @family integrative ranking functions
#' @seealso \code{\link[influential]{exir.vis}},
#' \code{\link[influential]{diff_data.assembly}},
#' \code{\link[stats]{prcomp}},
#' \code{\link[ranger]{ranger}},
#' \code{\link[ranger]{importance_pvalues}}
#' @export exir
#' @examples
#' \dontrun{
#' MyDesired_list <- Desiredlist
#' MyDiff_data <- Diffdata
#' Diff_value <- c(1,3,5)
#' Regr_value <- 7
#' Sig_value <- c(2,4,6,8)
#' MyExptl_data <- Exptldata
#' condition <- "condition"
#' My.exir <- exir(Desired_list = MyDesired_list,
#'                Diff_data = MyDiff_data, Diff_value = Diff_value,
#'                Regr_value = Regr_value, Sig_value = Sig_value,
#'                Exptl_data = MyExptl_data, condition = condition)
#' }
exir <- function(Desired_list = NULL,
                 Diff_data, Diff_value, Regr_value = NULL, Sig_value,
                 Exptl_data,
                 Exptl_data_type = c("bulk", "sc"),
                 condition = "condition",
                 Exptl_data_orientation = c("features_rows", "samples_rows"),
                 assay = "RNA",
                 layer = "counts",
                 normalize = FALSE,
                 pseudo_sample = FALSE,
                 pseudo_samples_per_group = 100,
                 Exptl_data_size_check = TRUE,
                 feature_filter = TRUE,
                 min_feature_prevalence = NULL,
                 min_feature_total = NULL,
                 min_feature_variance = 1e-12,
                 always_keep_diff_features = TRUE,
                 cor_thresh_method = "mr", r = 0.5, mr = 20,
                 max.connections = 50000, alpha = 0.05,
                 num_trees = 500, mtry = NULL, num_permutations = 50,
                 inf_const = 10^10, ncores = "default", seed = 1234, verbose = TRUE) {
  
  # -------------------------------------------------------------------------
  # cli-based verbose helpers
  # -------------------------------------------------------------------------
  
  .exir_stage_id <- 0L
  .exir_n_stages <- 15L
  .exir_current_stage <- NULL
  
  .exir_start <- function() {
    if(isTRUE(verbose)) {
      cli::cli_h1("Running ExIR")
    }
  }
  
  .exir_step <- function(stage, value = NULL) {
    if(isTRUE(verbose)) {
      
      .exir_stage_id <<- .exir_stage_id + 1L
      .exir_current_stage <<- stage
      
      cli::cli_h2(
        paste0(
          "[", .exir_stage_id, "/", .exir_n_stages, "] ",
          stage
        )
      )
    }
  }
  
  .exir_progress <- function(value = NULL, status = .exir_current_stage) {
    if(isTRUE(verbose) && !is.null(status)) {
      cli::cli_alert_success(status)
    }
  }
  
  .exir_info <- function(...) {
    if(isTRUE(verbose)) {
      cli::cli_alert_info(paste0(...))
    }
  }
  
  .exir_success <- function(...) {
    if(isTRUE(verbose)) {
      cli::cli_alert_success(paste0(...))
    }
  }
  
  .exir_warn <- function(...) {
    cli::cli_alert_warning(paste0(...))
  }
  
  .exir_done <- function() {
    if(isTRUE(verbose)) {
      cli::cli_h2("ExIR completed")
      cli::cli_alert_success("ExIR completed successfully.")
    }
  }
  
  .exir_start()
  .exir_step("Preparing the input data")
  
  # Checking NAs in Diff_data
  if(any(is.na(Diff_data))) {
    cli::cli_abort(
      c(
        "NA values found in {.arg Diff_data}.",
        "i" = "Please remove unnecessary columns and make sure there are no missing values in Diff_data.",
        "i" = "Also check {.arg Diff_value} and {.arg Sig_value} if you modified Diff_data columns."
      )
    )
  }
  
  # Match argument choices
  Exptl_data_type <- match.arg(Exptl_data_type)
  Exptl_data_orientation <- match.arg(Exptl_data_orientation)
  
  # make sure Diff_data is of data frame class
  Diff_data <- as.data.frame(Diff_data)
  
  # -------------------------------------------------------------------------
  # Helper functions for Exptl_data preparation
  # -------------------------------------------------------------------------
  
  .is_sparse_matrix <- function(x) {
    inherits(x, "sparseMatrix")
  }
  
  .replace_na_with_zero <- function(x) {
    if(.is_sparse_matrix(x)) {
      x@x[is.na(x@x)] <- 0
      return(x)
    }
    x[is.na(x)] <- 0
    x
  }
  
  .get_zero_fraction <- function(x) {
    if(.is_sparse_matrix(x)) {
      return(1 - (length(x@x) / (nrow(x) * ncol(x))))
    }
    
    n_total <- length(x)
    if(n_total == 0) return(NA_real_)
    
    # Avoid expensive full scans for very large dense matrices
    if(n_total > 1e6) {
      set.seed(seed)
      sampled <- sample(as.vector(x), size = min(1e6, n_total))
      return(mean(sampled == 0, na.rm = TRUE))
    }
    
    mean(x == 0, na.rm = TRUE)
  }
  
  .looks_integer_like <- function(x, tolerance = 1e-8) {
    if(.is_sparse_matrix(x)) {
      values <- x@x
    } else {
      values <- as.vector(x)
    }
    
    values <- values[is.finite(values)]
    if(length(values) == 0) return(FALSE)
    
    if(length(values) > 1e6) {
      set.seed(seed)
      values <- sample(values, size = 1e6)
    }
    
    all(values >= 0) && all(abs(values - round(values)) < tolerance)
  }
  
  .has_cell_barcode_like_colnames <- function(x) {
    cn <- colnames(x)
    if(is.null(cn) || length(cn) == 0) return(FALSE)
    
    # Common 10x-like patterns, e.g. AAACCTGAGAAACCAT-1
    barcode_hits <- grepl("^[ACGTN]{8,}[-_][0-9]+$", cn) |
      grepl("^[ACGTN]{12,}$", cn)
    
    mean(barcode_hits) > 0.3
  }
  
  .check_exptl_data_type_clues <- function(expr_features_by_samples,
                                           Exptl_data_type,
                                           from_seurat = FALSE,
                                           layer = NULL,
                                           verbose = TRUE) {
    
    n_features <- nrow(expr_features_by_samples)
    n_samples <- ncol(expr_features_by_samples)
    zero_fraction <- .get_zero_fraction(expr_features_by_samples)
    barcode_like <- .has_cell_barcode_like_colnames(expr_features_by_samples)
    integer_like <- .looks_integer_like(expr_features_by_samples)
    
    sc_score <- 0
    bulk_score <- 0
    
    if(n_samples >= 1000) sc_score <- sc_score + 1
    if(!is.na(zero_fraction) && zero_fraction > 0.5) sc_score <- sc_score + 1
    if(barcode_like) sc_score <- sc_score + 1
    if(from_seurat) sc_score <- sc_score + 1
    if(!is.null(layer) && layer %in% c("counts", "data", "scale.data")) sc_score <- sc_score + 1
    
    if(n_samples <= 500) bulk_score <- bulk_score + 1
    if(!is.na(zero_fraction) && zero_fraction < 0.3) bulk_score <- bulk_score + 1
    if(!barcode_like) bulk_score <- bulk_score + 1
    if(!integer_like) bulk_score <- bulk_score + 1
    
    if(Exptl_data_type == "bulk" && sc_score >= 2) {
      .exir_warn(
        "The input was declared as bulk, but some properties look single-cell-like ",
        "(samples/cells = ", n_samples,
        ", zero fraction = ", round(zero_fraction, 3),
        ", barcode-like sample names = ", barcode_like, "). ",
        "Please make sure Exptl_data_type = 'bulk' is correct."
      )
    }
    
    if(Exptl_data_type == "sc" && bulk_score >= 3 && !from_seurat) {
      .exir_warn(
        "The input was declared as single-cell, but some properties look bulk-like ",
        "(samples/cells = ", n_samples,
        ", zero fraction = ", round(zero_fraction, 3),
        ", barcode-like sample names = ", barcode_like, "). ",
        "Please make sure Exptl_data_type = 'sc' is correct."
      )
    }
    
    invisible(list(
      n_features = n_features,
      n_samples = n_samples,
      zero_fraction = zero_fraction,
      barcode_like = barcode_like,
      integer_like = integer_like,
      sc_score = sc_score,
      bulk_score = bulk_score
    ))
  }
  
  .extract_seurat_expression_and_condition <- function(Exptl_data,
                                                       condition,
                                                       assay,
                                                       layer) {
    
    if(!requireNamespace("SeuratObject", quietly = TRUE)) {
      cli::cli_abort(
        "A Seurat object was provided, but the {.pkg SeuratObject} package is not available."
      )
    }
    
    available_assays <- SeuratObject::Assays(Exptl_data)
    
    if(!(assay %in% available_assays)) {
      cli::cli_abort(
        c(
          "The specified {.arg assay} is not among the assays of the Seurat object.",
          "i" = "Available assays: {.val {available_assays}}.",
          "i" = "You can also check them using {.code SeuratObject::Assays(Exptl_data)}."
        )
      )
    }
    
    available_layers <- tryCatch(
      SeuratObject::Layers(Exptl_data[[assay]]),
      error = function(e) character(0)
    )
    
    if(length(available_layers) > 0 && !(layer %in% available_layers)) {
      cli::cli_abort(
        c(
          "The specified {.arg layer} is not among the layers of the selected assay.",
          "i" = "Available layers: {.val {available_layers}}."
        )
      )
    }
    
    expr <- tryCatch(
      SeuratObject::GetAssayData(
        object = Exptl_data,
        assay = assay,
        layer = layer
      ),
      error = function(e) {
        cli::cli_abort(
          c(
            "Could not extract expression data from the Seurat object.",
            "i" = "Please check {.arg assay} and {.arg layer}.",
            "x" = conditionMessage(e)
          )
        )
      }
    )
    
    if(nrow(expr) == 0 || ncol(expr) == 0) {
      cli::cli_abort("The selected Seurat assay/layer contains no data.")
    }
    
    meta <- Exptl_data[[]]
    
    if(!(length(condition) == 1 && condition %in% colnames(meta))) {
      cli::cli_abort(
        c(
          "For Seurat input, {.arg condition} must be the name of a metadata column.",
          "i" = "Available metadata columns include: {.val {utils::head(colnames(meta), 20)}}."
        )
      )
    }
    
    sample_condition <- as.character(meta[[condition]])
    names(sample_condition) <- rownames(meta)
    
    if(!all(colnames(expr) %in% names(sample_condition))) {
      cli::cli_abort("Some cells/samples in the selected assay/layer are missing from the Seurat metadata.")
    }
    
    sample_condition <- sample_condition[colnames(expr)]
    
    list(
      expr_features_by_samples = expr,
      condition = sample_condition,
      from_seurat = TRUE
    )
  }
  
  .extract_matrix_expression_and_condition <- function(Exptl_data,
                                                       condition,
                                                       Exptl_data_orientation) {
    
    x <- Exptl_data
    
    if(tibble::is_tibble(x)) {
      x <- as.data.frame(x, check.names = FALSE)
    }
    
    # If a feature-row data.frame has a first column containing feature names,
    # automatically use it as rownames when rownames are uninformative.
    if(is.data.frame(x) && Exptl_data_orientation == "features_rows") {
      rn_uninformative <- is.null(rownames(x)) ||
        anyDuplicated(rownames(x)) ||
        identical(rownames(x), as.character(seq_len(nrow(x))))
      
      first_col <- x[[1]]
      first_col_can_be_rownames <- is.character(first_col) &&
        !anyDuplicated(first_col) &&
        !(length(condition) == 1 && names(x)[1] == condition)
      
      if(rn_uninformative && first_col_can_be_rownames) {
        rownames(x) <- first_col
        x <- x[, -1, drop = FALSE]
      }
    }
    
    if(is.data.frame(x)) {
      
      if(Exptl_data_orientation == "samples_rows") {
        
        if(length(condition) == 1 && condition %in% colnames(x)) {
          sample_condition <- as.character(x[[condition]])
          expr <- x[, setdiff(colnames(x), condition), drop = FALSE]
        } else if(length(condition) == nrow(x)) {
          sample_condition <- as.character(condition)
          expr <- x
        } else {
          cli::cli_abort(
            "For {.code Exptl_data_orientation = 'samples_rows'}, {.arg condition} must be either a column name of {.arg Exptl_data} or a vector with length equal to the number of samples/cells."
          )
        }
        
        expr <- data.frame(lapply(expr, function(z) as.numeric(as.character(z))),
                           check.names = FALSE)
        expr <- as.matrix(expr)
        rownames(expr) <- rownames(x)
        
        # convert to features x samples
        expr <- t(expr)
        
      } else {
        
        if(length(condition) == 1 && condition %in% rownames(x)) {
          sample_condition <- as.character(unlist(x[condition, , drop = TRUE]))
          expr <- x[setdiff(rownames(x), condition), , drop = FALSE]
        } else if(length(condition) == ncol(x)) {
          sample_condition <- as.character(condition)
          expr <- x
        } else {
          cli::cli_abort(
            "For {.code Exptl_data_orientation = 'features_rows'}, {.arg condition} must be either a row name of {.arg Exptl_data} or a vector with length equal to the number of samples/cells."
          )
        }
        
        expr <- data.frame(lapply(expr, function(z) as.numeric(as.character(z))),
                           check.names = FALSE)
        expr <- as.matrix(expr)
      }
      
    } else if(is.matrix(x) || .is_sparse_matrix(x)) {
      
      if(Exptl_data_orientation == "samples_rows") {
        
        if(length(condition) != nrow(x)) {
          cli::cli_abort(
            "For matrix/sparse matrix input with {.code Exptl_data_orientation = 'samples_rows'}, {.arg condition} must be a vector with length equal to the number of rows/samples."
          )
        }
        
        sample_condition <- as.character(condition)
        expr <- Matrix::t(x)
        
      } else {
        
        if(length(condition) == 1 && !is.null(rownames(x)) && condition %in% rownames(x)) {
          sample_condition <- as.character(x[condition, ])
          expr <- x[setdiff(rownames(x), condition), , drop = FALSE]
        } else if(length(condition) == ncol(x)) {
          sample_condition <- as.character(condition)
          expr <- x
        } else {
          cli::cli_abort(
            "For matrix/sparse matrix input with {.code Exptl_data_orientation = 'features_rows'}, {.arg condition} must be either a row name of {.arg Exptl_data} or a vector with length equal to the number of columns/samples."
          )
        }
      }
      
      if(!.is_sparse_matrix(expr)) {
        storage.mode(expr) <- "double"
      }
      
    } else {
      cli::cli_abort(
        "{.arg Exptl_data} must be a data frame, tibble, matrix, sparse matrix, or Seurat object."
      )
    }
    
    if(is.null(rownames(expr))) {
      cli::cli_abort("Feature names are missing. Please provide feature names as rownames of {.arg Exptl_data}.")
    }
    
    if(is.null(colnames(expr))) {
      colnames(expr) <- paste0("Sample_", seq_len(ncol(expr)))
    }
    
    list(
      expr_features_by_samples = expr,
      condition = sample_condition,
      from_seurat = FALSE
    )
  }
  
  .normalize_with_edger <- function(counts_features_by_samples) {
    
    if(!requireNamespace("edgeR", quietly = TRUE)) {
      cli::cli_abort("The {.pkg edgeR} package is required for TMM/logCPM normalization.")
    }
    
    if(!.looks_integer_like(counts_features_by_samples)) {
      cli::cli_abort(
        paste0(
          "TMM/logCPM normalization requires raw count-like non-negative integer data. ",
          "If your data are already normalized or are not suitable for edgeR normalization, ",
          "set normalize = FALSE and provide pre-normalized data."
        )
      )
    }
    
    counts_features_by_samples <- as.matrix(counts_features_by_samples)
    
    dge <- edgeR::DGEList(counts = counts_features_by_samples)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    edgeR::cpm(dge, log = TRUE, prior.count = 1)
  }
  
  .make_pseudo_samples <- function(expr_features_by_samples,
                                   sample_condition,
                                   Exptl_data_type,
                                   pseudo_samples_per_group,
                                   seed) {
    
    if(!requireNamespace("Matrix", quietly = TRUE)) {
      cli::cli_abort("The {.pkg Matrix} package is required for pseudo-sampling.")
    }
    
    if(Exptl_data_type == "sc" && !requireNamespace("edgeR", quietly = TRUE)) {
      cli::cli_abort("The {.pkg edgeR} package is required for single-cell pseudo-bulk normalization.")
    }
    
    if(!is.numeric(pseudo_samples_per_group) ||
       length(pseudo_samples_per_group) != 1 ||
       is.na(pseudo_samples_per_group) ||
       pseudo_samples_per_group < 1) {
      cli::cli_abort("{.arg pseudo_samples_per_group} must be a positive integer.")
    }
    
    pseudo_samples_per_group <- as.integer(pseudo_samples_per_group)
    sample_condition <- factor(sample_condition)
    
    condition_sizes <- table(sample_condition)
    
    if(any(condition_sizes < pseudo_samples_per_group)) {
      cli::cli_abort(
        c(
          "{.arg pseudo_samples_per_group} is larger than the number of samples/cells in at least one condition.",
          "i" = "Condition sizes: {.val {paste(names(condition_sizes), condition_sizes, sep = '=', collapse = ', ')}}.",
          "i" = "Use a smaller {.arg pseudo_samples_per_group}, or disable pseudo-sampling."
        )
      )
    }
    
    if(Exptl_data_type == "sc" && !.looks_integer_like(expr_features_by_samples)) {
      cli::cli_abort(
        paste0(
          "For Exptl_data_type = 'sc' with pseudo_sample = TRUE, the input expression data must be ",
          "raw counts or count-like non-negative integers. If your single-cell data are already normalized, ",
          "set pseudo_sample = FALSE and normalize = FALSE."
        )
      )
    }
    
    set.seed(seed)
    
    pseudo_expr_list <- list()
    pseudo_condition <- character(0)
    pseudo_names <- character(0)
    
    for(cond in levels(sample_condition)) {
      
      sample_idx <- which(sample_condition == cond)
      n_cond <- length(sample_idx)
      
      # Randomize sample/cell order within condition
      sample_idx <- sample(sample_idx, n_cond, replace = FALSE)
      
      # We generate exactly pseudo_samples_per_group groups per condition.
      # Group sizes differ by at most 1.
      base_group_size <- n_cond %/% pseudo_samples_per_group
      remainder <- n_cond %% pseudo_samples_per_group
      
      group_sizes <- rep(base_group_size, pseudo_samples_per_group)
      
      if(remainder > 0) {
        group_sizes[seq_len(remainder)] <- group_sizes[seq_len(remainder)] + 1
      }
      
      group_end <- cumsum(group_sizes)
      group_start <- c(1, utils::head(group_end, -1) + 1)
      
      groups <- Map(function(start, end) {
        sample_idx[start:end]
      }, group_start, group_end)
      
      for(i in seq_along(groups)) {
        
        group_idx <- groups[[i]]
        pseudo_name <- paste0(cond, "_pseudo_", i)
        
        if(Exptl_data_type == "bulk") {
          
          # For bulk data, input should already be normalized/log-transformed
          # before reaching this function.
          pseudo_vec <- Matrix::rowMeans(expr_features_by_samples[, group_idx, drop = FALSE])
          
        } else {
          
          # For scRNA-seq, raw counts are summed first.
          pseudo_vec <- Matrix::rowSums(expr_features_by_samples[, group_idx, drop = FALSE])
        }
        
        pseudo_expr_list[[pseudo_name]] <- pseudo_vec
        pseudo_condition <- c(pseudo_condition, as.character(cond))
        pseudo_names <- c(pseudo_names, pseudo_name)
      }
    }
    
    pseudo_expr <- do.call(cbind, pseudo_expr_list)
    rownames(pseudo_expr) <- rownames(expr_features_by_samples)
    colnames(pseudo_expr) <- pseudo_names
    
    if(Exptl_data_type == "sc") {
      
      # scRNA-seq pseudo-bulked counts are normalized after aggregation.
      pseudo_expr <- .normalize_with_edger(pseudo_expr)
    }
    
    list(
      expr_features_by_samples = pseudo_expr,
      condition = pseudo_condition
    )
  }
  
  .make_internal_exir_dataframe <- function(expr_features_by_samples,
                                            sample_condition) {
    
    # Convert from features x samples to samples x features
    expr_samples_by_features <- Matrix::t(expr_features_by_samples)
    
    # The current ExIR downstream code, ranger, PCA, and fcor expect dense data.
    # Delaying this conversion until after optional pseudo-sampling minimizes memory use.
    if(.is_sparse_matrix(expr_samples_by_features)) {
      expr_samples_by_features <- as.matrix(expr_samples_by_features)
    }
    
    Exptl_data_internal <- as.data.frame(
      expr_samples_by_features,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    
    Exptl_data_internal <- .replace_na_with_zero(Exptl_data_internal)
    
    Exptl_data_internal$condition <- as.factor(sample_condition)
    
    Exptl_data_internal
  }
  
  .filter_exptl_features <- function(expr_features_by_samples,
                                     sample_condition,
                                     Diff_data,
                                     Desired_list,
                                     Exptl_data_type,
                                     pseudo_sample,
                                     feature_filter,
                                     min_feature_prevalence,
                                     min_feature_total,
                                     min_feature_variance,
                                     always_keep_diff_features,
                                     verbose = TRUE) {
    
    if(!feature_filter) {
      return(expr_features_by_samples)
    }
    
    n_features_before <- nrow(expr_features_by_samples)
    n_samples <- ncol(expr_features_by_samples)
    
    if(is.null(rownames(expr_features_by_samples))) {
      cli::cli_abort("Feature names are missing from Exptl_data.")
    }
    
    # -----------------------------------------------------------------------
    # Adaptive conservative prevalence threshold
    # -----------------------------------------------------------------------
    
    if(is.null(min_feature_prevalence)) {
      
      if(Exptl_data_type == "bulk") {
        min_feature_prevalence <- max(2L, ceiling(0.10 * n_samples))
      } else {
        if(pseudo_sample) {
          min_feature_prevalence <- max(2L, ceiling(0.05 * n_samples))
        } else {
          min_feature_prevalence <- max(10L, ceiling(0.001 * n_samples))
        }
      }
      
      min_feature_prevalence <- min(min_feature_prevalence, n_samples)
    }
    
    # -----------------------------------------------------------------------
    # Prevalence / total / variance filters
    # -----------------------------------------------------------------------
    
    if(inherits(expr_features_by_samples, "sparseMatrix")) {
      
      prevalence <- Matrix::rowSums(expr_features_by_samples != 0)
      feature_total <- Matrix::rowSums(expr_features_by_samples)
      
      # Variance for sparse matrices.
      feature_mean <- Matrix::rowMeans(expr_features_by_samples)
      feature_mean_sq <- Matrix::rowMeans(expr_features_by_samples ^ 2)
      feature_var <- feature_mean_sq - feature_mean ^ 2
      
    } else {
      
      prevalence <- rowSums(expr_features_by_samples != 0, na.rm = TRUE)
      feature_total <- rowSums(expr_features_by_samples, na.rm = TRUE)
      feature_var <- apply(expr_features_by_samples, 1, stats::var, na.rm = TRUE)
    }
    
    keep_prevalence <- prevalence >= min_feature_prevalence
    keep_variance <- is.finite(feature_var) & feature_var > min_feature_variance
    
    if(is.null(min_feature_total)) {
      keep_total <- rep(TRUE, n_features_before)
    } else {
      keep_total <- feature_total >= min_feature_total
    }
    
    keep <- keep_prevalence & keep_variance & keep_total
    
    # -----------------------------------------------------------------------
    # Always retain candidate/differential features if requested
    # -----------------------------------------------------------------------
    
    if(always_keep_diff_features) {
      
      force_keep <- unique(c(
        rownames(Diff_data),
        if(!is.null(Desired_list)) Desired_list else character(0)
      ))
      
      force_keep <- force_keep[force_keep %in% rownames(expr_features_by_samples)]
      
      keep[match(force_keep, rownames(expr_features_by_samples))] <- TRUE
    }
    
    expr_filtered <- expr_features_by_samples[keep, , drop = FALSE]
    
    n_features_after <- nrow(expr_filtered)
    
    .exir_info(
      "Feature filtering retained ",
      n_features_after,
      " of ",
      n_features_before,
      " features before ExIR correlation analysis."
    )
    
    .exir_info(
      "Feature filtering thresholds: min_feature_prevalence = ",
      min_feature_prevalence,
      ", min_feature_total = ",
      ifelse(is.null(min_feature_total), "NULL", min_feature_total),
      ", min_feature_variance = ",
      min_feature_variance,
      "."
    )
    
    if(n_features_after < 2) {
      cli::cli_abort(
        "Feature filtering retained fewer than two features. Please relax the filtering thresholds."
      )
    }
    
    if(!always_keep_diff_features) {
      
      removed_diff_features <- setdiff(
        intersect(rownames(Diff_data), rownames(expr_features_by_samples)),
        rownames(expr_filtered)
      )
      
      if(length(removed_diff_features) > 0) {
        .exir_warn(
          length(removed_diff_features),
          " Diff_data features were removed by feature filtering. ",
          "Set always_keep_diff_features = TRUE to force their retention."
        )
      }
    }
    
    expr_filtered
  }
  
  .prepare_exptl_data <- function(Exptl_data,
                                  Exptl_data_type,
                                  condition,
                                  Exptl_data_orientation,
                                  assay,
                                  layer,
                                  normalize,
                                  pseudo_sample,
                                  pseudo_samples_per_group,
                                  Exptl_data_size_check,
                                  feature_filter,
                                  min_feature_prevalence,
                                  min_feature_total,
                                  min_feature_variance,
                                  always_keep_diff_features,
                                  Diff_data,
                                  Desired_list,
                                  seed,
                                  verbose) {
    
    from_seurat <- inherits(Exptl_data, "Seurat")
    
    if(from_seurat) {
      extracted <- .extract_seurat_expression_and_condition(
        Exptl_data = Exptl_data,
        condition = condition,
        assay = assay,
        layer = layer
      )
    } else {
      extracted <- .extract_matrix_expression_and_condition(
        Exptl_data = Exptl_data,
        condition = condition,
        Exptl_data_orientation = Exptl_data_orientation
      )
    }
    
    expr_features_by_samples <- extracted$expr_features_by_samples
    sample_condition <- extracted$condition
    
    if(length(sample_condition) != ncol(expr_features_by_samples)) {
      cli::cli_abort(
        "The length of {.arg condition} does not match the number of samples/cells in {.arg Exptl_data}."
      )
    }
    
    expr_features_by_samples <- .replace_na_with_zero(expr_features_by_samples)
    
    data_clues <- .check_exptl_data_type_clues(
      expr_features_by_samples = expr_features_by_samples,
      Exptl_data_type = Exptl_data_type,
      from_seurat = from_seurat,
      layer = layer,
      verbose = verbose
    )
    
    n_samples <- ncol(expr_features_by_samples)
    
    # -------------------------------------------------------------------------
    # Normalization handling
    # -------------------------------------------------------------------------
    
    if(normalize) {
      
      if(Exptl_data_type == "bulk") {
        
        .exir_info("Normalizing bulk count-like data using edgeR TMM normalization followed by logCPM transformation.")
        
        expr_features_by_samples <- .normalize_with_edger(expr_features_by_samples)
        
      } else if(Exptl_data_type == "sc" && !pseudo_sample) {
        
        cli::cli_abort(
          paste0(
            "normalize = TRUE is not supported for single-cell input when pseudo_sample = FALSE. ",
            "For single-cell data without pseudo-sampling, please provide pre-normalized data and set normalize = FALSE. ",
            "Alternatively, provide raw counts and set pseudo_sample = TRUE so that ExIR can pseudo-bulk the data ",
            "and then apply TMM/logCPM normalization."
          )
        )
        
      } else if(Exptl_data_type == "sc" && pseudo_sample) {
        
        .exir_info("Single-cell normalization will be applied after pseudo-bulk aggregation using edgeR TMM/logCPM.")
      }
    }
    
    if(Exptl_data_size_check && !pseudo_sample && n_samples > 500) {
      
      msg <- paste0(
        "The input Exptl_data contains ", n_samples, " samples/cells. ",
        "This may substantially increase computation time and memory use. ",
        "Pseudo-sampling is recommended when the number of samples/cells is greater than 500 ",
        "or when computational resources are limited."
      )
      
      if(interactive()) {
        
        do_pseudo <- utils::askYesNo(
          paste0(msg, "\n\nWould you like to perform pseudo-sampling?")
        )
        
        if(isTRUE(do_pseudo)) {
          
          pseudo_sample <- TRUE
          
          pseudo_input <- readline(
            prompt = paste0(
              "Number of pseudo-samples per condition group [default: ",
              pseudo_samples_per_group,
              "]: "
            )
          )
          
          if(nzchar(pseudo_input)) {
            pseudo_samples_per_group <- as.integer(pseudo_input)
          }
          
        } else if(isFALSE(do_pseudo)) {
          .exir_info("Continuing without pseudo-sampling.")
        } else if(is.na(do_pseudo)) {
          cli::cli_abort("ExIR execution was cancelled by the user.")
        }
        
      } else {
        .exir_warn(
          msg,
          " To enable pseudo-sampling, set pseudo_sample = TRUE and ",
          "pseudo_samples_per_group = 100, or another suitable value."
        )
      }
    }
    
    if(pseudo_sample) {
      
      .exir_info(
        "Performing pseudo-sampling with ",
        pseudo_samples_per_group,
        " pseudo-samples per condition group."
      )
      
      pseudo_res <- .make_pseudo_samples(
        expr_features_by_samples = expr_features_by_samples,
        sample_condition = sample_condition,
        Exptl_data_type = Exptl_data_type,
        pseudo_samples_per_group = pseudo_samples_per_group,
        seed = seed
      )
      
      expr_features_by_samples <- pseudo_res$expr_features_by_samples
      sample_condition <- pseudo_res$condition
    }
    
    if(Exptl_data_type == "bulk" && !normalize && .looks_integer_like(expr_features_by_samples)) {
      .exir_warn(
        "The input was declared as bulk and normalize = FALSE, but the values look count-like. ",
        "Bulk ExIR input should usually be normalized and log-transformed before running the model. ",
        "If these are raw bulk RNA-seq counts, consider setting normalize = TRUE."
      )
    }
    
    if(Exptl_data_type == "sc" && !pseudo_sample) {
      .exir_warn(
        "Single-cell input is being used without pseudo-sampling. ",
        "The current ExIR downstream workflow will eventually require dense operations for RF, PCA, and fcor, ",
        "which may be memory-intensive for large datasets."
      )
    }
    
    # -------------------------------------------------------------------------
    # Conservative feature filtering before dense conversion and correlation
    # -------------------------------------------------------------------------
    
    expr_features_by_samples <- .filter_exptl_features(
      expr_features_by_samples = expr_features_by_samples,
      sample_condition = sample_condition,
      Diff_data = Diff_data,
      Desired_list = Desired_list,
      Exptl_data_type = Exptl_data_type,
      pseudo_sample = pseudo_sample,
      feature_filter = feature_filter,
      min_feature_prevalence = min_feature_prevalence,
      min_feature_total = min_feature_total,
      min_feature_variance = min_feature_variance,
      always_keep_diff_features = always_keep_diff_features,
      verbose = verbose
    )
    
    Exptl_data_internal <- .make_internal_exir_dataframe(
      expr_features_by_samples = expr_features_by_samples,
      sample_condition = sample_condition
    )
    
    list(
      Exptl_data = Exptl_data_internal,
      n_samples_original = n_samples,
      pseudo_sample = pseudo_sample,
      pseudo_samples_per_group = pseudo_samples_per_group,
      data_clues = data_clues
    )
  }
  
  # -------------------------------------------------------------------------
  # Prepare Exptl_data
  # -------------------------------------------------------------------------
  
  prepared_exptl <- .prepare_exptl_data(
    Exptl_data = Exptl_data,
    Exptl_data_type = Exptl_data_type,
    condition = condition,
    Exptl_data_orientation = Exptl_data_orientation,
    assay = assay,
    layer = layer,
    normalize = normalize,
    pseudo_sample = pseudo_sample,
    pseudo_samples_per_group = pseudo_samples_per_group,
    Exptl_data_size_check = Exptl_data_size_check,
    feature_filter = feature_filter,
    min_feature_prevalence = min_feature_prevalence,
    min_feature_total = min_feature_total,
    min_feature_variance = min_feature_variance,
    always_keep_diff_features = always_keep_diff_features,
    Diff_data = Diff_data,
    Desired_list = Desired_list,
    seed = seed,
    verbose = verbose
  )
  
  Exptl_data <- prepared_exptl$Exptl_data
  
  # Internal condition column used by ExIR
  condition.index <- match("condition", colnames(Exptl_data))
  
  if(is.na(condition.index)) {
    cli::cli_abort("Internal error: the prepared Exptl_data does not contain a {.code condition} column.")
  }
  
  #change the colnames of Diff_data
  base::colnames(Diff_data) <- base::paste("source",
                                           base::colnames(Diff_data),
                                           sep = ".")
  
  #change the Inf/-Inf diff values (applicable to sc-Data)
  for(i in 1:base::length(Diff_value)) {
    
    if(any(base::is.infinite(Diff_data[, Diff_value[i]]))) {
      
      temp.max.abs.diff.value <-
        base::max(base::abs(Diff_data[, Diff_value[i]][!base::is.infinite(Diff_data[, Diff_value[i]])]))
      
      temp.inf.index <- base::which(base::is.infinite(Diff_data[, Diff_value[i]]))
      
      Diff_data[temp.inf.index, Diff_value[i]] <-
        base::ifelse(base::unlist(Diff_data[temp.inf.index, Diff_value[i]]) > 0,
                     temp.max.abs.diff.value * inf_const,
                     -1 * temp.max.abs.diff.value * inf_const)
    }
  }
  
  .exir_progress(5, "Input data prepared")
  
  .exir_step("Calculating the differential score", 5)
  
  #1 Calculate differential score
  Diff_data$sum.Diff_value <- base::abs(base::apply(Diff_data[,Diff_value, drop = FALSE],1,sum))
  #range normalize the differential score
  Diff_data$sum.Diff_value <- 1+(((Diff_data$sum.Diff_value-min(Diff_data$sum.Diff_value))*(100-1))/
                                   (max(Diff_data$sum.Diff_value)-min(Diff_data$sum.Diff_value)))
  
  .exir_progress(10, "Differential score calculated")
  
  if(!is.null(Regr_value)) {
    .exir_step("Calculating the regression/time-course R-squared score", 10)
  }
  
  #2 Calculate regression/time-course R-squared score (if provided)
  if (!is.null(Regr_value)) {
    Diff_data$sum.Regr_value <- base::apply(Diff_data[,Regr_value, drop = FALSE],1,sum)
    #range normalize the R-squared score
    Diff_data$sum.Regr_value <- 1+(((Diff_data$sum.Regr_value-min(Diff_data$sum.Regr_value))*(100-1))/
                                     (max(Diff_data$sum.Regr_value)-min(Diff_data$sum.Regr_value)))
  }
  
  if(!is.null(Regr_value)) {
    .exir_progress(15, "Regression/time-course score calculated")
  }
  
  .exir_step("Calculating the collective statistical significance score", 15)
  
  #3 Calculate statistical significance of differential/regression factors
  if (max(Diff_data[,Sig_value]) > 1 | min(Diff_data[,Sig_value]) < 0) {
    cli::cli_abort("Input significance values must all be in the range 0 to 1.")
  }
  
  for(m in 1:length(Sig_value)) {
    
    if(min(Diff_data[,Sig_value[m]])==0) {
      
      #range normalize the primitive Sig_value
      temp.min_Sig_value <- base::sort(base::unique(Diff_data[,Sig_value[m]]))[2]
      
      Diff_data[,Sig_value[m]] <- temp.min_Sig_value+
        (((Diff_data[,Sig_value[m]]-min(Diff_data[,Sig_value[m]]))*(max(Diff_data[,Sig_value[m]])-temp.min_Sig_value))/
           (max(Diff_data[,Sig_value[m]])-min(Diff_data[,Sig_value[m]])))
    }
  }
  
  Diff_data$sum.Sig_value <- base::apply(-log10(Diff_data[,Sig_value, drop = FALSE]),1,sum)
  #range normalize the statistical significance
  Diff_data$sum.Sig_value <- 1+(((Diff_data$sum.Sig_value-min(Diff_data$sum.Sig_value))*(100-1))/
                                  (max(Diff_data$sum.Sig_value)-min(Diff_data$sum.Sig_value)))
  
  .exir_progress(20, "Statistical significance score calculated")
  
  .exir_step("Performing random forest classification", 20)
  
  #4 Calculation of the Integrated Value of Influence (IVI)
  
  #a Separate a transcriptomic profile of diff features
  if(!is.null(Desired_list)) {
    sig.diff.index <- stats::na.omit(base::unique(base::match(Desired_list,
                                                              colnames(Exptl_data))))
  } else {
    sig.diff.index <- stats::na.omit(base::unique(base::match(rownames(Diff_data),
                                                              colnames(Exptl_data))))
  }
  
  exptl.for.super.learn <- Exptl_data[,sig.diff.index]
  exptl.for.super.learn$condition <- Exptl_data[,condition.index]
  
  #correct the names of features
  #first preserve a copy of original names
  features.exptl.for.super.learn <- colnames(exptl.for.super.learn)[-ncol(exptl.for.super.learn)]
  
  colnames(exptl.for.super.learn) <- janitor::make_clean_names(colnames(exptl.for.super.learn))
  
  # Defining the num.threads
  if(ncores == "default") {
    num.threads <- parallel::detectCores() - 1
  } else {
    num.threads <- ncores
  }
  
  #b Perform random forests classification
  base::set.seed(seed = seed)
  rf.diff.exptl <- tryCatch(
    
    {
      ranger::ranger(
        formula = condition ~ .,
        data = exptl.for.super.learn,
        num.trees = num_trees,
        mtry = mtry,
        importance = "impurity_corrected",
        write.forest = FALSE,
        num.threads = num.threads,
        seed = seed
      )
    },
    
    error = function(e) {
      
      if (grepl("protection stack overflow", conditionMessage(e))) {
        
        n_features <- ncol(exptl.for.super.learn) - 1
        
        stop(
          paste0(
            "Supervised machine learning failed due to extremely high feature dimensionality.\n\n",
            "The input data contains ", n_features, " features (e.g. genes/proteins), which exceeds ",
            "what the formula interface in R can safely handle on this system.\n\n",
            "Please reduce the number of features prior to model fitting, for example by:\n",
            "  - Applying more stringent filtering (e.g. edgeR::filterByExpr)\n",
            "  - Selecting only highly variable features\n",
            "  - Removing low-expression or low-variance features"
          ),
          call. = FALSE
        )
        
      } else {
        stop(e)  # re-throw unrelated errors
      }
    }
  )
  
  .exir_info("Estimating random forest feature-importance p-values")
  
  base::set.seed(seed = seed)
  rf.diff.exptl.pvalue <- as.data.frame(ranger::importance_pvalues(x = rf.diff.exptl,
                                                                   formula = condition ~ .,
                                                                   num.permutations = num_permutations,
                                                                   data = exptl.for.super.learn,
                                                                   method = "altmann",
                                                                   num.threads = num.threads,
                                                                   seed = seed))
  
  #replace feature names (rownames) with their original names
  rownames(rf.diff.exptl.pvalue) <- features.exptl.for.super.learn
  
  if(any(is.na(rf.diff.exptl.pvalue[,"pvalue"])) |
     any(is.nan(rf.diff.exptl.pvalue[,"pvalue"]))) {
    rf.diff.exptl.pvalue[c(which(is.na(rf.diff.exptl.pvalue[,"pvalue"])),
                           which(is.nan(rf.diff.exptl.pvalue[,"pvalue"]))),
                         "pvalue"] <- 1
  }
  
  # filtering the RF output data
  select.number <- ifelse(!is.null(Desired_list),
                          round(length(Desired_list)/2),
                          100)
  
  if(length(which(rf.diff.exptl.pvalue[,"pvalue"] <alpha)) >= select.number) {
    rf.diff.exptl.pvalue <- base::subset(rf.diff.exptl.pvalue, rf.diff.exptl.pvalue$pvalue < alpha)
    
  } else {
    if(length(which(rf.diff.exptl.pvalue[, "pvalue"] < alpha)) >= 10) {
      rf.pval.select <- which(rf.diff.exptl.pvalue[, "pvalue"] < alpha)
    } else {
      rf.pval.select <- which(order(rf.diff.exptl.pvalue[, "pvalue"]) <= 10)
    }
    rf.nonSig <- seq(nrow(rf.diff.exptl.pvalue))[-rf.pval.select]
    required.pos.importance <- select.number - length(rf.pval.select)
    
    temp.rf.diff.exptl.pvalue <- rf.diff.exptl.pvalue[rf.nonSig,]
    
    rf.importance.select <- utils::tail(order(temp.rf.diff.exptl.pvalue[,"importance"]),
                                        n = required.pos.importance)
    
    temp.rf.diff.exptl.pvalue <- temp.rf.diff.exptl.pvalue[rf.importance.select,]
    
    #combine pvalue-based and importance-based tables
    rf.diff.exptl.pvalue <- rbind(rf.diff.exptl.pvalue[rf.pval.select,],
                                  temp.rf.diff.exptl.pvalue)
  }
  
  # negative importance values could be considered as 0
  if(any(rf.diff.exptl.pvalue[,"importance"] < 0)) {
    rf.diff.exptl.pvalue[which(rf.diff.exptl.pvalue[,"importance"] < 0),
                         "importance"] <- 0
  }
  
  # taking care of zero p-values
  if(min(rf.diff.exptl.pvalue[,"pvalue"])==0) {
    
    #range normalize the primitive pvalue
    temp.min_pvalue <- base::sort(base::unique(rf.diff.exptl.pvalue[,"pvalue"]))[2]
    
    rf.diff.exptl.pvalue[,"pvalue"] <- temp.min_pvalue+
      (((rf.diff.exptl.pvalue[,"pvalue"]-min(rf.diff.exptl.pvalue[,"pvalue"]))*(max(rf.diff.exptl.pvalue[,"pvalue"])-temp.min_pvalue))/
         (max(rf.diff.exptl.pvalue[,"pvalue"])-min(rf.diff.exptl.pvalue[,"pvalue"])))
  }
  
  .exir_progress(35, "Random forest classification completed")
  
  
  .exir_step("Performing PCA", 35)
  
  #5 Unsupervised machine learning (PCA)
  
  Exptl_data.for.PCA.index <- stats::na.omit(base::match(base::rownames(rf.diff.exptl.pvalue),
                                                         base::colnames(Exptl_data)))
  temp.Exptl_data.for.PCA <- Exptl_data[,Exptl_data.for.PCA.index]
  
  set.seed(seed)
  temp.PCA <- irlba::prcomp_irlba(
    as.matrix(temp.Exptl_data.for.PCA),
    n = 1,
    center = TRUE,
    scale. = FALSE
  )
  temp.PCA.r <- base::abs(temp.PCA$rotation[,1])
  names(temp.PCA.r) <- colnames(temp.Exptl_data.for.PCA)
  
  #range normalize the rotation values
  temp.PCA.r <- 1+(((temp.PCA.r-min(temp.PCA.r))*(100-1))/
                     (max(temp.PCA.r)-min(temp.PCA.r)))
  
  .exir_progress(40, "PCA completed")
  
  .exir_step("Performing first-round association analysis", 40)
  .exir_info("Calculating full feature-feature Spearman correlation network")
  
  #c Performing correlation analysis
  temp.corr <- fcor(data = Exptl_data[, -condition.index],
                    method = "spearman",
                    mutualRank = ifelse(cor_thresh_method == "mr", TRUE, FALSE))
  
  # Define selected / differential features once
  selected.features <- rownames(rf.diff.exptl.pvalue)
  
  # -------------------------------------------------------------------------
  # First round of association analysis
  # -------------------------------------------------------------------------
  
  # Filter corr data for only those correlations between selected features
  # and themselves/others.
  filter.corr.index <- stats::na.omit(base::unique(c(
    base::which(temp.corr$row %in% selected.features),
    base::which(temp.corr$column %in% selected.features)
  )))
  
  temp.corr.round1 <- temp.corr[filter.corr.index, , drop = FALSE]
  
  # filtering low level correlations
  cor.thresh <- r
  mr.thresh <- mr
  
  if(cor_thresh_method == "mr") {
    
    temp.corr.round1 <- base::subset(temp.corr.round1,
                                     temp.corr.round1[, 4] < mr.thresh)
    
    if(nrow(temp.corr.round1) > (max.connections * 0.95)) {
      
      temp.corr.select.index <- utils::head(order(temp.corr.round1$mr),
                                            n = round(max.connections * 0.95))
      
      temp.corr.round1 <- temp.corr.round1[temp.corr.select.index, , drop = FALSE]
    }
    
  } else if(cor_thresh_method == "cor.coefficient") {
    
    temp.corr.round1 <- base::subset(temp.corr.round1,
                                     base::abs(temp.corr.round1[, 3]) > cor.thresh)
    
    if(nrow(temp.corr.round1) > (max.connections * 0.95)) {
      
      temp.corr.select.index <- utils::tail(order(temp.corr.round1$cor),
                                            n = round(max.connections * 0.95))
      
      temp.corr.round1 <- temp.corr.round1[temp.corr.select.index, , drop = FALSE]
    }
  }
  
  diff.only.temp.corr <- temp.corr.round1
  rm(temp.corr.round1)
  
  # Getting the list of selected features and their correlated features
  diff.plus.corr.features <- base::unique(c(
    base::as.character(diff.only.temp.corr[, 1]),
    base::as.character(diff.only.temp.corr[, 2])
  ))
  
  # Find selected features among diff.plus.corr.features
  diff.only.features.index <- stats::na.omit(base::unique(base::match(
    selected.features,
    diff.plus.corr.features
  )))
  
  non.diff.only.features <- diff.plus.corr.features[-diff.only.features.index]
  
  .exir_progress(50, "First-round association analysis completed")
  
  .exir_step("Performing second-round association analysis", 50)
  
  # -------------------------------------------------------------------------
  # Second round of association analysis
  # -------------------------------------------------------------------------
  
  if(base::length(non.diff.only.features) > 0) {
    
    # Filter the original full temp.corr directly.
    filter.corr.index <- stats::na.omit(base::unique(c(
      base::which(temp.corr$row %in% non.diff.only.features),
      base::which(temp.corr$column %in% non.diff.only.features)
    )))
    
    temp.corr.round2 <- temp.corr[filter.corr.index, , drop = FALSE]
    
    # filtering low level correlations
    cor.thresh <- r
    mr.thresh <- mr
    
    if(cor_thresh_method == "mr") {
      
      temp.corr.round2 <- base::subset(temp.corr.round2,
                                       temp.corr.round2[, 4] < mr.thresh)
      
    } else if(cor_thresh_method == "cor.coefficient") {
      
      temp.corr.round2 <- base::subset(temp.corr.round2,
                                       base::abs(temp.corr.round2[, 3]) > cor.thresh)
    }
    
    # Remove correlations that involve selected/differential features.
    temp.corr.diff.only.index <- stats::na.omit(base::unique(c(
      base::which(temp.corr.round2$row %in% selected.features),
      base::which(temp.corr.round2$column %in% selected.features)
    )))
    
    if(base::length(temp.corr.diff.only.index) > 0) {
      temp.corr.round2 <- temp.corr.round2[-temp.corr.diff.only.index, , drop = FALSE]
    }
    
    if(nrow(temp.corr.round2) > (max.connections - nrow(diff.only.temp.corr))) {
      
      if(cor_thresh_method == "mr") {
        
        temp.corr.select.index <- utils::head(order(temp.corr.round2$mr),
                                              n = (max.connections - nrow(diff.only.temp.corr)))
        
      } else if(cor_thresh_method == "cor.coefficient") {
        
        temp.corr.select.index <- utils::tail(order(temp.corr.round2$cor),
                                              n = (max.connections - nrow(diff.only.temp.corr)))
      }
      
      temp.corr.round2 <- temp.corr.round2[temp.corr.select.index, , drop = FALSE]
    }
    
    # Now the full correlation table is no longer needed.
    # Overwrite temp.corr with the final edge table.
    temp.corr <- base::rbind(temp.corr.round2, diff.only.temp.corr)
    
    rm(temp.corr.round2, diff.only.temp.corr)
    
  } else {
    
    # No second-round features; final edge table is only first-round edges.
    temp.corr <- diff.only.temp.corr
    rm(diff.only.temp.corr)
  }
  
  # release memory after the full correlation table has been overwritten
  gc()
  
  .exir_progress(60, "Second-round association analysis completed")
  
  .exir_step("Reconstructing the association network", 60)
  
  #d Graph reconstruction
  temp.corr.graph <- igraph::graph_from_data_frame(temp.corr[, c(1:2)])
  
  # Remove self-loops, if any.
  # Multiple edges are retained to avoid unintentionally changing graph structure.
  temp.corr.graph <- igraph::simplify(
    temp.corr.graph,
    remove.multiple = FALSE,
    remove.loops = TRUE
  )
  
  .exir_progress(65, "Association network reconstructed")
  
  .exir_step("Calculating integrated value of influence (IVI)", 65)
  
  #e Calculation of IVI
  temp.corr.ivi <- ivi(temp.corr.graph, ncores = ncores)
  
  .exir_progress(70, "IVI calculated")
  
  .exir_step("Calculating the primitive driver score", 70)
  
  ## Driver score and ranking
  
  #a calculate first level driver score based on #3 and #4
  
  Diff_data.IVI.index <- stats::na.omit(match(names(temp.corr.ivi),
                                              rownames(Diff_data)))
  
  if(length(Diff_data.IVI.index) > 0) {
    
    Diff_data$IVI <- 0
    
    temp.corr.ivi.for.Diff_data.IVI.index <- stats::na.omit(match(rownames(Diff_data)[Diff_data.IVI.index],
                                                                  names(temp.corr.ivi)))
    
    Diff_data$IVI[Diff_data.IVI.index] <- temp.corr.ivi[temp.corr.ivi.for.Diff_data.IVI.index]
    
    #range normalize the IVI
    Diff_data$IVI <- 1+(((Diff_data$IVI-min(Diff_data$IVI))*(100-1))/
                          (max(Diff_data$IVI)-min(Diff_data$IVI)))
    
  } else {
    Diff_data$IVI <- 1
  }
  
  diff_vals <- as.matrix(Diff_data[, Diff_value, drop = FALSE])
  
  has_neg <- rowSums(diff_vals < 0) > 0
  has_pos <- rowSums(diff_vals > 0) > 0
  mixed_sign <- has_neg & has_pos
  
  Diff_data$first.Driver.Rank <- ifelse(
    mixed_sign,
    0,
    Diff_data$sum.Sig_value * Diff_data$IVI
  )
  
  #range normalize the first driver rank
  if(any(Diff_data$first.Driver.Rank == 0)) {
    Diff_data$first.Driver.Rank <- 0+(((Diff_data$first.Driver.Rank-min(Diff_data$first.Driver.Rank))*(100-0))/
                                        (max(Diff_data$first.Driver.Rank)-min(Diff_data$first.Driver.Rank)))
  } else {
    Diff_data$first.Driver.Rank <- 1+(((Diff_data$first.Driver.Rank-min(Diff_data$first.Driver.Rank))*(100-1))/
                                        (max(Diff_data$first.Driver.Rank)-min(Diff_data$first.Driver.Rank)))
  }
  
  .exir_progress(75, "Primitive driver score calculated")
  
  .exir_step("Calculating the neighbourhood driver score", 75)
  
  #b (#6) calculate neighborhood score
  
  #get the list of network nodes
  network.nodes <- base::as.character(igraph::as_ids(igraph::V(temp.corr.graph)))
  
  # Build a named vector of first.Driver.Rank for all graph nodes.
  # Nodes not present in Diff_data receive score 0, matching the original behaviour.
  rank_vec <- stats::setNames(
    rep(0, length(network.nodes)),
    network.nodes
  )
  
  rank_match <- match(rownames(Diff_data), network.nodes)
  valid_rank_match <- !is.na(rank_match)
  
  rank_vec[rank_match[valid_rank_match]] <- Diff_data$first.Driver.Rank[valid_rank_match]
  
  # Build sparse adjacency matrix.
  # Since self-loops were removed after graph construction, self-correlations do not
  # contribute to the neighbourhood score.
  adj_mat <- igraph::as_adjacency_matrix(
    temp.corr.graph,
    type = "both",
    sparse = TRUE
  )
  
  # Ensure rank_vec follows the same node order as the adjacency matrix.
  rank_vec <- rank_vec[rownames(adj_mat)]
  
  # Neighbourhood score = sum of first.Driver.Rank across direct neighbours.
  N_scores <- as.numeric(adj_mat %*% rank_vec)
  
  neighborehood.score.table <- data.frame(
    node = rownames(adj_mat),
    N.score = N_scores,
    stringsAsFactors = FALSE
  )
  
  .exir_progress(80, "Neighbourhood driver score calculated")
  
  .exir_step("Preparing the driver table", 80)
  
  # Map N.score back to Diff_data
  Diff_data$N.score <- 0
  
  nscore_match <- match(rownames(Diff_data), neighborehood.score.table$node)
  valid_nscore_match <- !is.na(nscore_match)
  
  Diff_data$N.score[valid_nscore_match] <- neighborehood.score.table$N.score[nscore_match[valid_nscore_match]]
  
  #range normalize (1,100) the neighborhood score
  Diff_data$N.score <- ifelse(
    sum(Diff_data$N.score) == 0,
    1,
    1 + (((Diff_data$N.score - min(Diff_data$N.score)) * (100 - 1)) /
           (max(Diff_data$N.score) - min(Diff_data$N.score)))
  )
  
  #c calculate the final driver score
  
  Diff_data$final.Driver.score <- (Diff_data$first.Driver.Rank)*(Diff_data$N.score)
  Diff_data$final.Driver.score[Diff_data$final.Driver.score==0] <- NA
  
  # Create the Drivers table
  
  Driver.table <- Diff_data
  
  #remove the rows/features with NA in the final driver score
  Driver.table <- Driver.table[stats::complete.cases(Driver.table),]
  
  #filter the driver table by the desired list (if provided)
  if(!is.null(Desired_list)) {
    Driver.table.row.index <- stats::na.omit(match(Desired_list,
                                                   rownames(Driver.table)))
    Driver.table <- Driver.table[Driver.table.row.index,]
  }
  
  if(nrow(as.data.frame(Driver.table))==0) {Driver.table <- NULL} else {
    
    #add Z.score
    Driver.table$Z.score <- base::scale(Driver.table$final.Driver.score)
    
    #range normalize final driver score
    ifelse(length(unique(Driver.table$final.Driver.score)) > 1,
           Driver.table$final.Driver.score <- 1+(((Driver.table$final.Driver.score-min(Driver.table$final.Driver.score))*(100-1))/
                                                   (max(Driver.table$final.Driver.score)-min(Driver.table$final.Driver.score))),
           Driver.table$final.Driver.score <- 1)
    
    #add driver rank
    Driver.table$rank <- rank(-Driver.table$final.Driver.score,
                              ties.method = "min")
    
    #add P-value
    Driver.table$p.value <- stats::pnorm(Driver.table$Z.score,
                                         lower.tail = FALSE)
    
    #add adjusted pvalue
    Driver.table$padj <- stats::p.adjust(p = Driver.table$p.value,
                                         method = "BH")
    
    #add driver type
    driver_diff_sum <- rowSums(as.matrix(Driver.table[, Diff_value, drop = FALSE]))
    
    Driver.table$driver.type <- ifelse(
      driver_diff_sum < 0,
      "Decelerator",
      ifelse(
        driver_diff_sum > 0,
        "Accelerator",
        NA_character_
      )
    )
    
    Driver.table <- Driver.table[stats::complete.cases(Driver.table),]
    
    #remove redundent columns
    Driver.table <- Driver.table[,c("final.Driver.score",
                                    "Z.score",
                                    "rank",
                                    "p.value",
                                    "padj",
                                    "driver.type")]
    
    #rename column names
    colnames(Driver.table) <- c("Score", "Z.score",
                                "Rank", "P.value",
                                "P.adj", "Type")
    
    #filtering redundant (NaN) results
    Driver.table <- Driver.table[stats::complete.cases(Driver.table),]
    
    Driver.table <- cbind("Driver" = rownames(Driver.table), Driver.table)
    Driver.table$Z.score <- as.numeric(Driver.table$Z.score)
    Driver.table$P.value <- as.numeric(Driver.table$P.value)
    
    if(nrow(as.data.frame(Driver.table))==0) {Driver.table <- NULL}
    
  }
  
  .exir_progress(85, "Driver table prepared")
  
  .exir_step("Preparing the biomarker table", 85)
  
  # Create the Biomarker table
  
  Biomarker.table <- Diff_data
  
  #remove the rows/features with NA in the final driver score
  Biomarker.table <- Biomarker.table[stats::complete.cases(Biomarker.table),]
  
  #filter the biomarker table by the desired list (if provided)
  if(!is.null(Desired_list)) {
    Biomarker.table.row.index <- stats::na.omit(match(Desired_list,
                                                      rownames(Biomarker.table)))
    Biomarker.table <- Biomarker.table[Biomarker.table.row.index,]
  }
  
  if(nrow(as.data.frame(Biomarker.table))==0) {Biomarker.table <- NULL} else {
    
    #add RF importance score and p-value
    Biomarker.table$rf.importance <- 0
    Biomarker.table$rf.pvalue <- 1
    
    Biomarker.table.rf.index <- stats::na.omit(match(rownames(rf.diff.exptl.pvalue),
                                                     rownames(Biomarker.table)))
    
    rf.for.Biomarker.table <- stats::na.omit(match(rownames(Biomarker.table)[Biomarker.table.rf.index],
                                                   rownames(rf.diff.exptl.pvalue)))
    
    Biomarker.table$rf.importance[Biomarker.table.rf.index] <- rf.diff.exptl.pvalue$importance[rf.for.Biomarker.table]
    Biomarker.table$rf.pvalue[Biomarker.table.rf.index] <- rf.diff.exptl.pvalue$pvalue[rf.for.Biomarker.table]
    
    #range normalize rf.importance and rf.pvalue
    Biomarker.table$rf.importance <- 1+(((Biomarker.table$rf.importance-min(Biomarker.table$rf.importance))*(100-1))/
                                          (max(Biomarker.table$rf.importance)-min(Biomarker.table$rf.importance)))
    
    Biomarker.table$rf.pvalue <- -log10(Biomarker.table$rf.pvalue)
    Biomarker.table$rf.pvalue <- 1+(((Biomarker.table$rf.pvalue-min(Biomarker.table$rf.pvalue))*(100-1))/
                                      (max(Biomarker.table$rf.pvalue)-min(Biomarker.table$rf.pvalue)))
    
    #add rotation values
    Biomarker.table$rotation <- 0
    Biomarker.table.for.rotation <- stats::na.omit(match(names(temp.PCA.r),
                                                         rownames(Biomarker.table)))
    
    Biomarker.table.rotation.index <- stats::na.omit(match(rownames(Biomarker.table)[Biomarker.table.for.rotation],
                                                           names(temp.PCA.r)))
    
    Biomarker.table$rotation[Biomarker.table.for.rotation] <- temp.PCA.r[Biomarker.table.rotation.index]
    
    #range normalize rotation values
    Biomarker.table$rotation <- 1+(((Biomarker.table$rotation-min(Biomarker.table$rotation))*(100-1))/
                                     (max(Biomarker.table$rotation)-min(Biomarker.table$rotation)))
    
    #calculate biomarker score
    if(!is.null(Regr_value)) {
      Biomarker.table$final.biomarker.score <- (Biomarker.table$sum.Diff_value)*
        (Biomarker.table$sum.Regr_value)*(Biomarker.table$sum.Sig_value)*
        (Biomarker.table$rf.pvalue)*(Biomarker.table$rf.importance)*
        (Biomarker.table$rotation)
    } else {
      Biomarker.table$final.biomarker.score <- (Biomarker.table$sum.Diff_value)*
        (Biomarker.table$sum.Sig_value)*
        (Biomarker.table$rf.pvalue)*(Biomarker.table$rf.importance)*
        (Biomarker.table$rotation)
    }
    
    #add biomarker Z.score
    Biomarker.table$Z.score <- base::scale(Biomarker.table$final.biomarker.score)
    
    #range normalize biomarker score
    ifelse(length(unique(Biomarker.table$final.biomarker.score)) > 1,
           Biomarker.table$final.biomarker.score <- 1+(((Biomarker.table$final.biomarker.score-min(Biomarker.table$final.biomarker.score))*(100-1))/
                                                         (max(Biomarker.table$final.biomarker.score)-min(Biomarker.table$final.biomarker.score))),
           Biomarker.table$final.biomarker.score <- 1)
    
    #add biomarker rank
    Biomarker.table$rank <- rank(-Biomarker.table$final.biomarker.score, ties.method = "min")
    
    #add biomarker P-value
    Biomarker.table$P.value <- stats::pnorm(Biomarker.table$Z.score,
                                            lower.tail = FALSE)
    
    #add biomarker adjusted p-value
    Biomarker.table$padj <- stats::p.adjust(p = Biomarker.table$P.value,
                                            method = "BH")
    
    #add biomarker type
    biomarker_diff_sum <- rowSums(as.matrix(Biomarker.table[, Diff_value, drop = FALSE]))
    
    Biomarker.table$type <- ifelse(
      biomarker_diff_sum < 0,
      "Down-regulated",
      ifelse(
        biomarker_diff_sum > 0,
        "Up-regulated",
        NA_character_
      )
    )
    
    #remove redundent columns
    Biomarker.table <- Biomarker.table[,c("final.biomarker.score",
                                          "Z.score", "rank",
                                          "P.value", "padj", "type")]
    
    #rename column names
    colnames(Biomarker.table) <- c("Score", "Z.score",
                                   "Rank", "P.value",
                                   "P.adj", "Type")
    
    #filtering redundant (NaN) results
    Biomarker.table <- Biomarker.table[stats::complete.cases(Biomarker.table),]
    
    Biomarker.table <- cbind("Biomarker" = rownames(Biomarker.table), Biomarker.table)
    Biomarker.table$Z.score <- as.numeric(Biomarker.table$Z.score)
    Biomarker.table$P.value <- as.numeric(Biomarker.table$P.value)
    
    if(nrow(as.data.frame(Biomarker.table))==0) {Biomarker.table <- NULL}
    
  }
  
  .exir_progress(90, "Biomarker table prepared")
  
  #_________________________________
  
  get_associated_drivers <- function(nodes, graph, Driver.table) {
    
    if(length(nodes) == 0 || is.null(Driver.table) || nrow(Driver.table) == 0) {
      return(list(
        first_order = character(0),
        second_order = character(0)
      ))
    }
    
    driver_names <- rownames(Driver.table)
    
    # Compute all first- and second-order neighbourhoods in batch.
    # This preserves the same igraph::neighborhood() semantics as the original code.
    all_nbrs_1 <- igraph::neighborhood(
      graph = graph,
      order = 1,
      nodes = nodes,
      mode = "all"
    )
    
    all_nbrs_2 <- igraph::neighborhood(
      graph = graph,
      order = 2,
      nodes = nodes,
      mode = "all"
    )
    
    # First-order associated drivers
    first_order_drivers <- lapply(all_nbrs_1, function(nbrs) {
      
      nbr_names <- igraph::as_ids(nbrs)
      
      assoc_drivers <- Driver.table[
        driver_names %in% nbr_names,
        ,
        drop = FALSE
      ]
      
      assoc_drivers <- rownames(assoc_drivers[order(assoc_drivers$Rank), , drop = FALSE])
      
      assoc_drivers
    })
    
    # Second-order associated drivers
    second_order_drivers <- lapply(seq_along(all_nbrs_2), function(i) {
      
      nbr_names <- igraph::as_ids(all_nbrs_2[[i]])
      
      assoc_drivers <- Driver.table[
        driver_names %in% nbr_names,
        ,
        drop = FALSE
      ]
      
      # Preserve original behaviour:
      # remove first-order drivers from the second-order driver list.
      assoc_drivers <- assoc_drivers[
        !(rownames(assoc_drivers) %in% first_order_drivers[[i]]),
        ,
        drop = FALSE
      ]
      
      assoc_drivers <- rownames(assoc_drivers[order(assoc_drivers$Rank), , drop = FALSE])
      
      assoc_drivers
    })
    
    first_order_drivers <- vapply(
      first_order_drivers,
      function(x) paste0(x, collapse = ", "),
      character(1)
    )
    
    second_order_drivers <- vapply(
      second_order_drivers,
      function(x) paste0(x, collapse = ", "),
      character(1)
    )
    
    list(
      first_order = first_order_drivers,
      second_order = second_order_drivers
    )
  }
  
  #_________________________________
  
  .exir_step("Preparing the DE-mediator table", 90)
  
  # Create the DE mediators table
  
  DE.mediator.table <- Diff_data
  
  #include only rows/features with NA in the final driver score (which are mediators)
  DE.mediator.index <- which(is.na(DE.mediator.table$final.Driver.score))
  
  DE.mediator.table <- DE.mediator.table[DE.mediator.index,]
  
  DE.mediator.table$DE.mediator.score <- DE.mediator.table$sum.Sig_value*
    DE.mediator.table$IVI*
    DE.mediator.table$N.score
  
  #filter the DE mediators table by either the desired list or the list of network nodes
  
  DE.mediator.row.index <- stats::na.omit(match(network.nodes,
                                                rownames(DE.mediator.table)))
  
  if(!is.null(Desired_list)) {
    desired.DE.mediator.row.index <- stats::na.omit(match(Desired_list,
                                                          rownames(DE.mediator.table)[DE.mediator.row.index]))
    DE.mediator.row.index <- DE.mediator.row.index[desired.DE.mediator.row.index]
  }
  
  DE.mediator.table <- DE.mediator.table[DE.mediator.row.index,]
  if(nrow(as.data.frame(DE.mediator.table))==0) {DE.mediator.table <- NULL} else {
    
    #add DE mediators Z score
    DE.mediator.table$Z.score <- base::scale(DE.mediator.table$DE.mediator.score)
    
    #range normalize DE mediators score
    ifelse(length(unique(DE.mediator.table$DE.mediator.score)) > 1,
           DE.mediator.table$DE.mediator.score <- 1+(((DE.mediator.table$DE.mediator.score-min(DE.mediator.table$DE.mediator.score))*(100-1))/
                                                       (max(DE.mediator.table$DE.mediator.score)-min(DE.mediator.table$DE.mediator.score))),
           DE.mediator.table$DE.mediator.score <- 1)
    
    #add DE mediators rank
    DE.mediator.table$rank <- rank(-DE.mediator.table$DE.mediator.score, ties.method = "min")
    
    #add DE mediators P-value
    DE.mediator.table$P.value <- stats::pnorm(DE.mediator.table$Z.score,
                                              lower.tail = FALSE)
    
    #add DE mediators adjusted P-value
    DE.mediator.table$padj <- stats::p.adjust(p = DE.mediator.table$P.value,
                                              method = "BH")
    
    #remove redundent columns
    DE.mediator.table <- DE.mediator.table[,c("DE.mediator.score", "Z.score",
                                              "rank", "P.value", "padj")]
    
    #rename column names
    colnames(DE.mediator.table) <- c("Score", "Z.score",
                                     "Rank", "P.value",
                                     "P.adj")
    
    #filtering redundant (NaN) results
    DE.mediator.table <- DE.mediator.table[stats::complete.cases(DE.mediator.table),]
    
    if(nrow(as.data.frame(DE.mediator.table)) == 0) {
      
      DE.mediator.table <- NULL
      
    } else {
      
      # Adding the associated drivers to the table
      associated_drivers <- get_associated_drivers(
        nodes = rownames(DE.mediator.table),
        graph = temp.corr.graph,
        Driver.table = Driver.table
      )
      
      DE.mediator.table$First.order.Drivers <- associated_drivers$first_order
      DE.mediator.table$Second.order.Drivers <- associated_drivers$second_order
      
      DE.mediator.table <- cbind("DE.mediator" = rownames(DE.mediator.table), DE.mediator.table)
      DE.mediator.table$Z.score <- as.numeric(DE.mediator.table$Z.score)
      DE.mediator.table$P.value <- as.numeric(DE.mediator.table$P.value)
    }
  }
  
  .exir_progress(95, "DE-mediator table prepared")
  
  .exir_step("Preparing the non-DE-mediator table", 95)
  
  # Create the non-DE mediators table
  non.DE.mediators.index <- stats::na.omit(unique(match(rownames(Diff_data),
                                                        neighborehood.score.table$node)))
  
  non.DE.mediators.table <- neighborehood.score.table[-c(non.DE.mediators.index),]
  if(nrow(as.data.frame(non.DE.mediators.table))==0) {non.DE.mediators.table <- NULL} else {
    
    #filter the non-DE mediators table by either the desired list
    if(!is.null(Desired_list)) {
      non.DE.mediator.row.index <- stats::na.omit(match(Desired_list,
                                                        non.DE.mediators.table$node))
      non.DE.mediators.table <- non.DE.mediators.table[non.DE.mediator.row.index,]
    }
  }
  
  if(nrow(as.data.frame(non.DE.mediators.table))==0) {non.DE.mediators.table <- NULL} else {
    
    rownames(non.DE.mediators.table) <- non.DE.mediators.table$node
    
    non.DE.mediators.ivi.index <- stats::na.omit(match(rownames(non.DE.mediators.table),
                                                       names(temp.corr.ivi)))
    
    non.DE.mediators.table$ivi <- temp.corr.ivi[non.DE.mediators.ivi.index]
    
    non.DE.mediators.table$non.DE.mediator.score <- non.DE.mediators.table$N.score*non.DE.mediators.table$ivi
    
    #add non-DE mediators Z.score
    non.DE.mediators.table$Z.score <- base::scale(non.DE.mediators.table$non.DE.mediator.score)
    
    #range normalize nonDE mediators score
    ifelse(length(unique(non.DE.mediators.table$non.DE.mediator.score)) > 1,
           non.DE.mediators.table$non.DE.mediator.score <- 1+(((non.DE.mediators.table$non.DE.mediator.score-min(non.DE.mediators.table$non.DE.mediator.score))*(100-1))/
                                                                (max(non.DE.mediators.table$non.DE.mediator.score)-min(non.DE.mediators.table$non.DE.mediator.score))),
           non.DE.mediators.table$non.DE.mediator.score <- 1)
    
    #add non-DE mediators P-value
    non.DE.mediators.table$P.value <- stats::pnorm(non.DE.mediators.table$Z.score,
                                                   lower.tail = FALSE)
    
    #add non-DE mediators adjusted p-value
    non.DE.mediators.table$padj <- stats::p.adjust(p = non.DE.mediators.table$P.value,
                                                   method = "BH")
    
    #add non-DE mediators rank
    non.DE.mediators.table$rank <- rank(-non.DE.mediators.table$non.DE.mediator.score, ties.method = "min")
    
    #remove redundent columns
    non.DE.mediators.table <- non.DE.mediators.table[,c("non.DE.mediator.score",
                                                        "Z.score", "rank",
                                                        "P.value", "padj")]
    
    #rename column names
    colnames(non.DE.mediators.table) <- c("Score", "Z.score",
                                          "Rank", "P.value",
                                          "P.adj")
    
    #filtering redundant (NaN) results
    non.DE.mediators.table <- non.DE.mediators.table[stats::complete.cases(non.DE.mediators.table),]
    
    if(nrow(as.data.frame(non.DE.mediators.table)) == 0) {
      
      non.DE.mediators.table <- NULL
      
    } else {
      
      # Adding the associated drivers to the table
      associated_drivers <- get_associated_drivers(
        nodes = rownames(non.DE.mediators.table),
        graph = temp.corr.graph,
        Driver.table = Driver.table
      )
      
      non.DE.mediators.table$First.order.Drivers <- associated_drivers$first_order
      non.DE.mediators.table$Second.order.Drivers <- associated_drivers$second_order
      
      non.DE.mediators.table <- cbind("non.DE.mediator" = rownames(non.DE.mediators.table), non.DE.mediators.table)
      non.DE.mediators.table$Z.score <- as.numeric(non.DE.mediators.table$Z.score)
      non.DE.mediators.table$P.value <- as.numeric(non.DE.mediators.table$P.value)
    }
  }
  
  .exir_progress(100, "non-DE-mediator table prepared")
  
  Results <- list("Driver table" = Driver.table,
                  "DE-mediator table" = DE.mediator.table,
                  "nonDE-mediator table" = non.DE.mediators.table,
                  "Biomarker table" = Biomarker.table,
                  "Graph" = temp.corr.graph)
  
  Results.Non.Null.vector <- vector()
  
  for (i in 1:length(Results)) {
    Results.Non.Null.vector[i] <- !is.null(Results[[i]])
  }
  
  Results <- Results[Results.Non.Null.vector]
  
  # set the class of Results
  base::class(Results) <- "ExIR_Result"
  
  if(isTRUE(verbose)) {
    cli::cli_h2("ExIR output summary")
    cli::cli_ul(c(
      paste0("Driver table: ", ifelse(is.null(Driver.table), "not returned", paste0(nrow(Driver.table), " rows"))),
      paste0("DE-mediator table: ", ifelse(is.null(DE.mediator.table), "not returned", paste0(nrow(DE.mediator.table), " rows"))),
      paste0("non-DE-mediator table: ", ifelse(is.null(non.DE.mediators.table), "not returned", paste0(nrow(non.DE.mediators.table), " rows"))),
      paste0("Biomarker table: ", ifelse(is.null(Biomarker.table), "not returned", paste0(nrow(Biomarker.table), " rows"))),
      paste0("Graph nodes: ", igraph::vcount(temp.corr.graph)),
      paste0("Graph edges: ", igraph::ecount(temp.corr.graph))
    ))
  }
  
  .exir_done()
  
  return(Results)
}
