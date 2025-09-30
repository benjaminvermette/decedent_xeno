library (readr)

#### Download TCRb data from Adaptive Biotech immuneACCESS Database in Combined Rearrangements format -----

CombinedRearrangements = read_tsv("CombinedRearrangements.tsv")
CombinedRearrangements = as.data.frame(CombinedRearrangements)

#### Filter out known lab contaminants ---------

contaminants.sequences <- c(
  "CASSAETQYF", #20D11
  "CASSIEGPTGELFF", #D222D2
  "CASSSFWGSDTGELFF", #Clone 5
  "CASSLAGGANSPLHF", # #BRI164 R164 
  "CASSLSFGTEAFF",#DMF5 MART1:26-35
  "CASSSTGLPYGYTF",#HAI.7
  "CASSTRLAGDEQFF", #1.C8
  "CASSLSASGGATDTQYF",#A3.10
  "CASWTVSYNEQFF", #1F3 
  "CASSLGPGQRETQYF", #A5.5
  "CASSLERDGYTF" , #A1.9
  "CASSSIKSGGSTDTQYF", #DH1.1
  "CASSPGSGLVETQYF" , #DH17.3
  "CARGLAGQETQYF" , #DR4
  "CASSPYSGDTYEQYF", #CTL2
  "CASNQGYTEAFF" , #IFH1
  "CASSFNAGGSDTQYF" #ID3
)

CombinedRearrangements = CombinedRearrangements[!(CombinedRearrangements$`Amino Acid` %in% contaminants.sequences),]

#### Resolve ambiguous clones, ie. clones present in both CD4 & CD8s -----
# to solve for sorting errors
# need to open 'resolveambiguous' R script from Github and save functions

# --- select columns by pattern ---
cols_cd4_all <- grep("^(Stim-|Unstim-).*-CD4$", colnames(CombinedRearrangements), value = TRUE)
cols_cd8_all <- grep("^(Stim-|Unstim-).*-CD8$", colnames(CombinedRearrangements), value = TRUE)
cols_cd4_all <- setdiff(cols_cd4_all, "CD4")
cols_cd8_all <- setdiff(cols_cd8_all, "CD8")

# --- 1) depth-normalize each column independently ---
CombinedRearrangements[, cols_cd4_all] <- normalize(as.matrix(CombinedRearrangements[, cols_cd4_all, drop = FALSE]))
CombinedRearrangements[, cols_cd8_all] <- normalize(as.matrix(CombinedRearrangements[, cols_cd8_all, drop = FALSE]))

# --- 2) pool via rowMeans ---
CombinedRearrangements$CD4 <- rowMeans(CombinedRearrangements[, cols_cd4_all, drop = FALSE], na.rm = TRUE)
CombinedRearrangements$CD8 <- rowMeans(CombinedRearrangements[, cols_cd8_all, drop = FALSE], na.rm = TRUE)

# --- 3) resolve once on pooled columns ---
CombinedRearrangements <- resolveambiguous(
  CombinedRearrangements, "CD4", "CD8",
  ratio = 5, cutoff = 1e-5, do_normalize = FALSE
)

# indices for pooled assignment
CD4 <- which(CombinedRearrangements$CD4 > 0)
CD8 <- which(CombinedRearrangements$CD8 > 0)

# --- 4) zero *all* lineage columns (Stim + Unstim) according to pooled call ---
all_rows <- seq_len(nrow(CombinedRearrangements))
CombinedRearrangements[setdiff(all_rows, CD4), cols_cd4_all] <- 0
CombinedRearrangements[setdiff(all_rows, CD8), cols_cd8_all] <- 0

#### Identify XDRTCC's ------
## --- CONFIG ---
ref_col <- "Unstim-XKLT001-POD0-PBMC"  # PreTx Unstim reference
cutoff  <- 1e-5
fold    <- 2

## --- HELPERS ---
# Safely get sequence IDs from rownames; fallback to the "Amino Acid" column if needed
.get_seq_ids <- function(df) {
  rn <- rownames(df)
  if (!is.null(rn) && length(rn) && any(nzchar(rn))) return(rn)
  if ("Amino Acid" %in% colnames(df)) return(df[["Amino Acid"]])
  stop("No rownames and no 'Amino Acid' column found for sequence IDs.")
}

# Collect XDRTCCs from a set of Stim columns against the PreTx Unstim reference
get_xdrtcc <- function(df, ref_col, stim_cols, cutoff = 1e-5, fold = 2) {
  stopifnot(ref_col %in% colnames(df))
  stopifnot(all(stim_cols %in% colnames(df)))
  seq_ids <- .get_seq_ids(df)
  
  hits <- lapply(stim_cols, function(sc) {
    # condition: present in stim > cutoff AND stim >= fold * ref
    idx <- which(df[[sc]] > cutoff & df[[sc]] >= fold * df[[ref_col]])
    if (length(idx)) seq_ids[idx] else character(0)
  })
  unique(unlist(hits, use.names = FALSE))
}

## --- SELECT COLUMNS BY PATTERN ---
# All Stim lineage columns
stim_cd4_cols <- grep("^Stim-.*-CD4$", colnames(CombinedRearrangements), value = TRUE)
stim_cd8_cols <- grep("^Stim-.*-CD8$", colnames(CombinedRearrangements), value = TRUE)

## --- COMPUTE UNIQUE XDRTCC SETS ---
cd4_DRTCCs <- get_xdrtcc(CombinedRearrangements, ref_col, stim_cd4_cols, cutoff, fold)
cd8_DRTCCs <- get_xdrtcc(CombinedRearrangements, ref_col, stim_cd8_cols, cutoff, fold)

## --- BUILD FILTERED DATA FRAMES ---
seq_ids <- .get_seq_ids(CombinedRearrangements)
cd4_XDRTCC_df <- CombinedRearrangements[seq_ids %in% cd4_DRTCCs, , drop = FALSE]
cd8_XDRTCC_df <- CombinedRearrangements[seq_ids %in% cd8_DRTCCs, , drop = FALSE]

## (Optional) If you want the Amino Acid column to be the key and also as rownames:
if ("Amino Acid" %in% colnames(CombinedRearrangements)) {
  rownames(cd4_XDRTCC_df) <- cd4_XDRTCC_df[["Amino Acid"]]
  rownames(cd8_XDRTCC_df) <- cd8_XDRTCC_df[["Amino Acid"]]
}

### Now you can track the XDRTCCs using simple commands, e.g. 
### for CD8 XDRTCC cumulative frequency at POD33

sum(cd8_XDRTCC_df$`Unstim-XKLT001-POD33-CD8`)





















