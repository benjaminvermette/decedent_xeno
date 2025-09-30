#' Format Adaptive/immuneACCESS Data for Immunarch
#'
#' Adapted from format_immunarch() (cdr3tools, Chris Parks) for Adpative V2.
#' Accepts a list of data frames (e.g., `data_list`) or a single data frame.
#'
#' Required/assumed columns in each df (vs. v1):
#'   nucleotide (-> CDR3.nt, Sequence)
#'   aminoAcid  (-> CDR3.aa)
#'   count (templates/reads) (-> Clones)
#'   vMaxResolved / dMaxResolved / jMaxResolved (-> V/D/J.name)
#'   vIndex, n1Index, dIndex, n2Index, jIndex (-> V.end, D.start, D.end, J.start)
#'   n1Insertion, n2Insertion (-> VD.ins, DJ.ins)
#'
#' @param data A list of data frames or a single data frame (Adaptive exports).
#' @param try_imgt Logical. If TRUE and cdr3tools is available, convert gene names to IMGT.
#' @returns A list or a single data frame in Immunarch format.
#' @export
format_immunarch_adaptive <- function(data, try_imgt = TRUE) {
  suppressPackageStartupMessages({
    requireNamespace("dplyr", quietly = TRUE)
    requireNamespace("purrr", quietly = TRUE)
    requireNamespace("tidyr", quietly = TRUE)
    requireNamespace("rlang", quietly = TRUE)
    requireNamespace("tidyselect", quietly = TRUE)
  })
  
  # If a single data.frame, wrap into a list
  if (!inherits(data, "list") && inherits(data, "data.frame")) {
    data <- list(data)
  }
  if (!inherits(data, "list")) {
    rlang::abort("`data` must be a list or a single data frame.")
  }
  
  # Helper: check cols exist
  .need <- function(df, cols) {
    miss <- setdiff(cols, names(df))
    if (length(miss)) {
      rlang::abort(paste0("Missing required columns: ", paste(miss, collapse = ", ")))
    }
  }
  
  # Optional IMGT conversion safely
  .maybe_imgt <- function(x) {
    if (try_imgt && requireNamespace("cdr3tools", quietly = TRUE)) {
      cdr3tools::imgt_format_gene_names(x)
    } else {
      x
    }
  }
  
  out <- purrr::map(data, function(df) {
    # --- Required columns in your schema ---
    required_cols <- c(
      "nucleotide", "aminoAcid", "count (templates/reads)",
      "vMaxResolved", "dMaxResolved", "jMaxResolved",
      "vIndex", "n1Index", "dIndex", "n2Index", "jIndex",
      "n1Insertion", "n2Insertion"
    )
    .need(df, required_cols)
    
    # Work in dplyr
    df_fmt <- dplyr::as_tibble(df)
    
    # Compute Immunarch boundary indices relative to V start
    df_fmt <- df_fmt %>%
      dplyr::mutate(
        V.end  = .data$n1Index - .data$vIndex,
        D.start= .data$dIndex  - .data$vIndex,
        D.end  = .data$n2Index - .data$vIndex,
        J.start= .data$jIndex  - .data$vIndex
      ) %>%
      dplyr::mutate(
        dplyr::across(
          .cols = c(.data$V.end, .data$D.start, .data$D.end, .data$J.start),
          .fns = ~ dplyr::if_else(.x < 0, as.numeric(NA), as.numeric(.x))
        )
      )
    
    # Select and rename to Immunarch format
    df_fmt <- df_fmt %>%
      dplyr::select(
        Clones      = dplyr::all_of("count (templates/reads)"),
        # Duplicates = tidyselect::any_of("duplicates"), # (not present, leave out)
        CDR3.nt     = .data$nucleotide,
        CDR3.aa     = .data$aminoAcid,
        V.name      = .data$vMaxResolved,
        D.name      = .data$dMaxResolved,
        J.name      = .data$jMaxResolved,
        .data$V.end,
        .data$D.start,
        .data$D.end,
        .data$J.start,
        VD.ins      = .data$n1Insertion,
        DJ.ins      = .data$n2Insertion,
        Sequence    = .data$nucleotide
      ) %>%
      dplyr::mutate(
        Proportion = .data$Clones / sum(.data$Clones, na.rm = TRUE),
        .after = .data$Clones
      ) %>%
      # TRB: VJ.ins is NA by definition
      dplyr::mutate(VJ.ins = NA_real_, .after = .data$J.start) %>%
      # Optional IMGT normalization for gene names
      dplyr::mutate(V.name = .maybe_imgt(.data$V.name),
                    D.name = .maybe_imgt(.data$D.name),
                    J.name = .maybe_imgt(.data$J.name)) %>%
      dplyr::arrange(dplyr::desc(.data$Proportion))
    
    attr(df_fmt, "repertoire_data_format") <- "adaptive_immunarch"
    df_fmt
  })
  
  if (length(out) == 1) out <- out[[1]]
  out
}
