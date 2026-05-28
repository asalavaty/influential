#=============================================================================
#
#    Fast correlation and mutual rank analysis
#
#=============================================================================

#' Fast correlation and mutual rank analysis
#'
#' This function calculates Pearson/Spearman correlations between all pairs of features in a matrix/dataframe much faster than the base R cor function.
#' It can also calculate correlations between all pairs of features from two input matrices/dataframes when `data2` is provided.
#' It is also possible to simultaneously calculate mutual rank (MR) of correlations as well as their p-values and adjusted p-values.
#' Additionally, this function can automatically combine and flatten the result matrices. Selecting correlated features using an MR-based threshold
#' rather than based on their correlation coefficients or an arbitrary p-value is more efficient and accurate in inferring
#' functional associations in systems, for example in gene regulatory networks.
#'
#' @param data a numeric dataframe/matrix with features on columns and samples/observations on rows. If `data2` is not provided,
#' correlations are calculated between all pairs of features in `data`.
#' @param data2 an optional numeric dataframe/matrix with features on columns and samples/observations on rows. If provided,
#' correlations are calculated between all features in `data` and all features in `data2`. `data` and `data2` must have the same
#' number of rows, and the rows must correspond to the same samples/observations in the same order. Default is `NULL`.
#' @param na_to_zero logical, whether to convert NAs to 0 in the output (default) or not.
#' @param method a character string indicating which correlation coefficient is to be computed. One of `"pearson"` or `"spearman"` (default).
#' @param mutualRank logical, whether to calculate mutual ranks of correlations or not.
#' @param mutualRank_mode a character string indicating whether to rank based on `"signed"` or `"unsigned"` (default) correlation values.
#' In the `"unsigned"` mode, only the level of a correlation value is important and not its sign; therefore, the function ranks the absolute
#' values of correlations. Options are `"unsigned"` and `"signed"`.
#' @param pvalue logical, whether to calculate p-values of correlations or not.
#' @param adjust p-value correction method when `pvalue = TRUE`, a character string including any of `"BH"` (default),
#' `"bonferroni"`, `"holm"`, `"hochberg"`, `"hommel"`, or `"none"`.
#' @param flat logical, whether to combine and flatten the result matrices or not.
#' @param remove_self logical, whether to remove self-correlations from the flattened output when `data2` is provided. This is useful when
#' `data2` contains some or all of the same features as `data`. Default is `TRUE`.
#' @param remove_duplicate_pairs logical, whether to remove duplicate undirected feature pairs from the flattened output when `data2` is provided.
#' This is useful when `data2` contains the same features as `data`, because pairs such as `geneA-geneB` and `geneB-geneA` may otherwise both
#' be returned. Default is `TRUE`.
#'
#' @return Depending on the input data and the value of `flat`, a dataframe or list including `cor` correlation coefficients,
#' `mr` mutual ranks of correlation coefficients, `p` p-values of correlation coefficients, and `p.adj` adjusted p-values.
#' If `data2` is not provided and `flat = TRUE`, the flattened output contains the upper triangle of the all-pairs correlation matrix.
#' If `data2` is provided and `flat = TRUE`, the flattened output contains feature pairs between `data` and `data2`.
#'
#' @details
#' When `data2 = NULL`, the function performs the standard all-pairs correlation analysis among the features of `data`.
#' When `data2` is provided, the function performs a rectangular correlation analysis between the features of `data` and the features of `data2`.
#'
#' For Spearman correlation with `data2`, the two input matrices are internally combined before rank transformation so that feature-wise ranks
#' are calculated consistently across the same samples/observations.
#'
#' When `mutualRank = TRUE` and `data2` is provided, the calculated MR values are based on the rectangular correlation space between `data`
#' and `data2`. Therefore, these MR values are not necessarily identical to MR values obtained from a full all-pairs correlation matrix followed
#' by post hoc filtering.
#'
#' @keywords fcor
#' @seealso \code{\link[stats]{p.adjust}} and \code{\link[influential]{graph_from_data_frame}}
#' @export fcor
#'
#' @examples
#' \dontrun{
#' set.seed(1234)
#'
#' # All-pairs correlation among features
#' data <- datasets::attitude
#' cor <- fcor(data = data)
#'
#' # Correlation between two sets of features
#' data1 <- mtcars[, 1:4]
#' data2 <- mtcars[, 5:11]
#' cor_rect <- fcor(data = data1, data2 = data2)
#'
#' # Correlation between selected features and all features
#' selected_data <- mtcars[, 1:4]
#' all_data <- mtcars
#' cor_selected_all <- fcor(data = selected_data, data2 = all_data)
#' }
fcor <- function(data,
                 data2 = NULL,
                 na_to_zero = TRUE,
                 method = "spearman",
                 mutualRank = TRUE,
                 mutualRank_mode = "unsigned",
                 pvalue = FALSE,
                 adjust = "BH",
                 flat = TRUE,
                 remove_self = TRUE,
                 remove_duplicate_pairs = TRUE) {
  
  #________________________________________
  # Dealing with warnings
  
  old_warn <- getOption("warn")
  options(warn = -1)
  on.exit(options(warn = old_warn), add = TRUE)
  
  #________________________________________
  # Input checks
  
  method <- match.arg(method, choices = c("pearson", "spearman"))
  mutualRank_mode <- match.arg(mutualRank_mode, choices = c("unsigned", "signed"))
  
  rectangular <- !is.null(data2)
  
  data <- as.matrix(data)
  
  if(!is.numeric(data)) {
    storage.mode(data) <- "double"
  }
  
  var_names1 <- colnames(data)
  obs_names1 <- rownames(data)
  
  if(is.null(var_names1)) {
    var_names1 <- paste0("V", seq_len(ncol(data)))
    colnames(data) <- var_names1
  }
  
  if(rectangular) {
    
    data2 <- as.matrix(data2)
    
    if(!is.numeric(data2)) {
      storage.mode(data2) <- "double"
    }
    
    var_names2 <- colnames(data2)
    obs_names2 <- rownames(data2)
    
    if(is.null(var_names2)) {
      var_names2 <- paste0("V", seq_len(ncol(data2)))
      colnames(data2) <- var_names2
    }
    
    if(nrow(data) != nrow(data2)) {
      cli::cli_abort("{.arg data} and {.arg data2} must have the same number of rows/samples.")
    }
    
    if(!is.null(obs_names1) && !is.null(obs_names2)) {
      if(!identical(obs_names1, obs_names2)) {
        cli::cli_abort("{.arg data} and {.arg data2} must have samples/cells in the same order.")
      }
    }
  }
  
  #________________________________________
  # Spearman ranking
  
  if(method == "spearman") {
    
    if(rectangular) {
      
      # Rank-transform the combined matrix so that tied ranks are handled
      # consistently across data and data2.
      combined_data <- cbind(data, data2)
      combined_names <- colnames(combined_data)
      combined_obs <- rownames(combined_data)
      
      ranked_t <- rank_matrix_cpp(
        t(combined_data),
        descending = FALSE,
        use_abs = FALSE
      )
      
      combined_ranked <- t(ranked_t)
      colnames(combined_ranked) <- combined_names
      rownames(combined_ranked) <- combined_obs
      
      n_data_cols <- ncol(data)
      
      data <- combined_ranked[, seq_len(n_data_cols), drop = FALSE]
      data2 <- combined_ranked[, (n_data_cols + 1):ncol(combined_ranked), drop = FALSE]
      
    } else {
      
      data_t <- t(data)
      
      ranked_t <- rank_matrix_cpp(
        data_t,
        descending = FALSE,
        use_abs = FALSE
      )
      
      data <- t(ranked_t)
      colnames(data) <- var_names1
      rownames(data) <- obs_names1
    }
  }
  
  #________________________________________
  # Helper: normalize matrix columns as feature rows
  
  normalize_for_cor <- function(x) {
    
    variableNames <- colnames(x)
    totalVariables <- ncol(x)
    
    validVariablesIdx <- which(apply(x, 2, stats::var) > 0)
    filteredData <- x[, validVariablesIdx, drop = FALSE]
    
    if(length(validVariablesIdx) == 0L) {
      return(
        list(
          normalized = NULL,
          valid_idx = validVariablesIdx,
          variable_names = variableNames,
          total_variables = totalVariables
        )
      )
    }
    
    centeredMatrix <- t(filteredData)
    centeredMatrix <- centeredMatrix - rowMeans(centeredMatrix)
    
    l2Norms <- sqrt(rowSums(centeredMatrix^2))
    l2Norms[l2Norms == 0] <- NA_real_
    
    normalizedMatrix <- centeredMatrix / l2Norms
    
    list(
      normalized = normalizedMatrix,
      valid_idx = validVariablesIdx,
      variable_names = variableNames,
      total_variables = totalVariables
    )
  }
  
  #________________________________________
  # Correlation calculation
  
  r <- NULL
  mutR <- NULL
  p <- NULL
  pa <- NULL
  
  if(!rectangular) {
    
    norm1 <- normalize_for_cor(data)
    
    if(length(norm1$valid_idx) == 0L) {
      
      r <- matrix(
        NA_real_,
        nrow = norm1$total_variables,
        ncol = norm1$total_variables,
        dimnames = list(norm1$variable_names, norm1$variable_names)
      )
      
      if(na_to_zero) {
        r[!is.finite(r)] <- 0
      }
      
    } else {
      
      correlationMatrix <- tcrossprod(norm1$normalized)
      
      r <- matrix(
        NA_real_,
        nrow = norm1$total_variables,
        ncol = norm1$total_variables
      )
      
      r[norm1$valid_idx, norm1$valid_idx] <- correlationMatrix
      rownames(r) <- norm1$variable_names
      colnames(r) <- norm1$variable_names
    }
    
  } else {
    
    norm1 <- normalize_for_cor(data)
    norm2 <- normalize_for_cor(data2)
    
    r <- matrix(
      NA_real_,
      nrow = norm1$total_variables,
      ncol = norm2$total_variables,
      dimnames = list(norm1$variable_names, norm2$variable_names)
    )
    
    if(length(norm1$valid_idx) > 0L && length(norm2$valid_idx) > 0L) {
      
      correlationMatrix <- tcrossprod(
        norm1$normalized,
        norm2$normalized
      )
      
      r[norm1$valid_idx, norm2$valid_idx] <- correlationMatrix
    }
  }
  
  if(na_to_zero) {
    r[!is.finite(r)] <- 0
  }
  
  #________________________________________
  # P-values
  
  if(pvalue) {
    
    n <- nrow(data)
    
    tstat <- (r * sqrt(n - 2)) / sqrt(1 - r^2)
    
    p <- -2 * expm1(
      stats::pt(
        abs(tstat),
        df = n - 2,
        log.p = TRUE
      )
    )
    
    p[p > 1 | is.nan(p)] <- 1
    
    if(adjust != "none") {
      pa <- matrix(
        stats::p.adjust(as.vector(p), method = adjust),
        nrow = nrow(p),
        ncol = ncol(p),
        dimnames = dimnames(p)
      )
    }
  }
  
  #________________________________________
  # Mutual rank
  
  if(mutualRank) {
    
    use_abs <- (mutualRank_mode == "unsigned")
    
    if(!rectangular) {
      
      r_rank <- rank_matrix_cpp(
        r,
        descending = TRUE,
        use_abs = use_abs
      )
      
      dimnames(r_rank) <- dimnames(r)
      mutR <- sqrt(r_rank * t(r_rank))
      
    } else {
      
      # Rectangular MR-like score.
      # Row rank: rank each data feature against all data2 features.
      # Column rank: rank each data2 feature against all data features.
      #
      # IMPORTANT:
      # This is NOT identical to full-matrix MR unless data and data2
      # contain the same full feature set.
      
      row_rank <- rank_matrix_cpp(
        r,
        descending = TRUE,
        use_abs = use_abs
      )
      
      col_rank <- t(
        rank_matrix_cpp(
          t(r),
          descending = TRUE,
          use_abs = use_abs
        )
      )
      
      dimnames(row_rank) <- dimnames(r)
      dimnames(col_rank) <- dimnames(r)
      
      mutR <- sqrt(row_rank * col_rank)
    }
  }
  
  #________________________________________
  # Flatten
  
  if(flat) {
    
    if(!rectangular) {
      
      flt_list <- flatten_cor_matrix_cpp(
        r,
        if(mutualRank) mutR else NULL,
        if(pvalue) p else NULL,
        if(pvalue && adjust != "none") pa else NULL,
        rownames(r)
      )
      
      result <- data.frame(
        flt_list,
        stringsAsFactors = FALSE,
        row.names = NULL
      )
      
    } else {
      
      flt_list <- flatten_rect_cor_matrix_cpp(
        r,
        if(mutualRank) mutR else NULL,
        if(pvalue) p else NULL,
        if(pvalue && adjust != "none") pa else NULL,
        rownames(r),
        colnames(r)
      )
      
      result <- data.frame(
        flt_list,
        stringsAsFactors = FALSE,
        row.names = NULL
      )
      
      if(remove_self) {
        result <- result[result$row != result$column, , drop = FALSE]
      }
      
      if(remove_duplicate_pairs) {
        
        pair_id <- ifelse(
          result$row < result$column,
          paste(result$row, result$column, sep = "___"),
          paste(result$column, result$row, sep = "___")
        )
        
        result <- result[!duplicated(pair_id), , drop = FALSE]
      }
    }
    
  } else {
    
    result <- list(
      r = r,
      mr = mutR,
      p = p,
      p.adj = pa
    )
  }
  
  class(result) <- c(class(result), "fcor", "influential")
  
  return(result)
}
