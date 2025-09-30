### Export Unstim Samples from Adaptive Website using v2 ----

library(readr)

# Path to zip file from Adaptive
zip_path <- "sampleExport.zip"

# List all files inside the zip
tsv_files <- unzip(zip_path, list = TRUE)$Name
tsv_files <- tsv_files[grepl("\\.tsv$", tsv_files, ignore.case = TRUE)]

# Extract and read all TSVs into a list
data_list <- lapply(tsv_files, function(f) {
  read_tsv(unz(zip_path, f))
})

# Give the list elements the same names as the files (without extension)
names(data_list) <- tools::file_path_sans_ext(basename(tsv_files))

#### Remove non-productive and known contaminants sequences ----

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

clean_tcr_list <- function(df_list, 
                           status_col = "sequenceStatus", # column specifying if rearrangement is productive
                           aa_col = "aminoAcid", # cdr3
                           contaminants = NULL) {
  lapply(df_list, function(df) {
    # skip filtering if required columns are missing
    if (!(status_col %in% names(df)) || !(aa_col %in% names(df))) {
      warning("Missing expected columns in one dataframe, skipping filter.")
      return(df)
    }
    # keep only sequenceStatus == "In"
    df <- df[df[[status_col]] == "In", ]
    # drop rows with NA in amino acid column
    df <- df[!is.na(df[[aa_col]]), ]
    # keep only amino acids starting with C and ending with F
    df <- df[grepl("^C.*F$", df[[aa_col]]), ]
    # remove known contaminants 
    if (!is.null(contaminants)) {
      df <- df[!(df[[aa_col]] %in% contaminants), ]
    }
    df
  })
}

data_list <- clean_tcr_list(
  data_list,
  status_col = "sequenceStatus",
  aa_col = "aminoAcid",
  contaminants = contaminants.sequences
)

### Renormalize frequencies ----

renormalize_freq_list <- function(df_list, freq_col = "frequencyCount (%)") {
  lapply(df_list, function(df) {
    if (!(freq_col %in% names(df))) {
      warning("Frequency column not found in one dataframe, skipping renormalization.")
      return(df)
    }
    if (nrow(df) > 0) {
      s <- sum(df[[freq_col]], na.rm = TRUE)
      if (s > 0) {
        df[[freq_col]] <- df[[freq_col]] / s
      }
    }
    df
  })
}


data_list <- renormalize_freq_list(data_list, freq_col = "frequencyCount (%)")

######## Compute nucleotide-per-aa ----

### For PostTx CD4 ----
# Combine the dataframes
CombinedRearrangementDetailsPostTxCD4 <- rbind(data_list$`Unstim-XKLT001-POD28-CD4`, 
                                               data_list$`Unstim-XKLT001-POD33-CD4`, 
                                               data_list$`Unstim-XKLT001-POD49-CD4`)

# Remove duplicate rows based on the "nucleotide" column
CombinedRearrangementDetailsPostTxCD4 <- CombinedRearrangementDetailsPostTxCD4[!duplicated(CombinedRearrangementDetailsPostTxCD4$nucleotide),]

# List of XDRTCCs can be obtained from Identifying_XDRTCCs.R script
XDRTCCsinPostTxCD4 = CombinedRearrangementDetailsPostTxCD4[CombinedRearrangementDetailsPostTxCD4$aminoAcid %in% cd4_XDRTCC_df$`Amino Acid`, ]
# compute nucleotide-per-aa for xdrtccs postTx CD4
length(unique(XDRTCCsinPostTxCD4$nucleotide))/length(unique(XDRTCCsinPostTxCD4$aminoAcid))

# For non-XDRTCCs, need to download TCRb data from Adaptive Biotech immuneACCESS Database in Combined Rearrangements format
# Non-XDRTCCs is defined as present in PreTx Unstim, but not in any Stim population, so we need to filter CombinedRearrangementDetailsPostTxCD4
# configurable column names 
aa_col_comb   <- "Amino Acid"                    # in CombinedRearrangements
pretx_col     <- "Unstim-XKLT001-POD0-PBMC"      # PreTx Unstim column
aa_col_detail <- "aminoAcid"                       # in CombinedRearrangementDetailsPostTxCD4

#  identify all Stim columns 
stim_cols <- grep("^Stim", colnames(CombinedRearrangements), value = TRUE)

#  amino acids that are present at PreTx (>0) and absent from ALL Stim columns 
eligible_aa <- with(CombinedRearrangements, {
  # count how many Stim columns have >0 for each row
  stim_pos <- if (length(stim_cols)) {
    rowSums(CombinedRearrangements[, stim_cols, drop = FALSE] > 0, na.rm = TRUE)
  } else {
    rep(0L, nrow(CombinedRearrangements))
  }
  keep <- (!is.na(CombinedRearrangements[[aa_col_comb]])) &
    (CombinedRearrangements[[pretx_col]] > 0) &
    (stim_pos == 0)
  CombinedRearrangements[[aa_col_comb]][keep]
})

NonDRTCCsinPostTxCD4 <- CombinedRearrangementDetailsPostTxCD4[
  CombinedRearrangementDetailsPostTxCD4[[aa_col_detail]] %in% eligible_aa,
  ,
  drop = FALSE
]

# compute nucleotide-per-aa for non-xdrtccs post tx cd4
length(unique(NonDRTCCsinPostTxCD4$nucleotide))/length(unique(NonDRTCCsinPostTxCD4$aminoAcid))

# Now look at rest -- or clones that are in the Post-Tx but neither XDRTCCs or non-XDRTCCs

# define column names explicitly 
aa_col_posttx <- "aminoAcid"    # in CombinedRearrangementDetailsPostTxCD4
aa_col_xdrtcc <- "Amino Acid"   # in cd4_XDRTCC_df
aa_col_nondrtcc <- "aminoAcid"  # in NonDRTCCsinPostTxCD4

# --- collect sequences to remove ---
to_remove <- c(
  cd4_XDRTCC_df[[aa_col_xdrtcc]],
  NonDRTCCsinPostTxCD4[[aa_col_nondrtcc]]
)

# --- subset, keeping only rows NOT in to_remove ---
RestinPostTxCD4 <- CombinedRearrangementDetailsPostTxCD4[
  !(CombinedRearrangementDetailsPostTxCD4[[aa_col_posttx]] %in% to_remove),
  ,
  drop = FALSE
]

# compute nucleotide-per-aa for rest in PostTx CD4
length(unique(RestinPostTxCD4$nucleotide))/length(unique(RestinPostTxCD4$aminoAcid))

### For PostTx CD8 -----

# Combine the dataframes
CombinedRearrangementDetailsPostTxCD8 <- rbind(data_list$`Unstim-XKLT001-POD28-CD8`, 
                                               data_list$`Unstim-XKLT001-POD33-CD8`, 
                                               data_list$`Unstim-XKLT001-POD49-CD8`)

# Remove duplicate rows based on the "nucleotide" column
CombinedRearrangementDetailsPostTxCD8 <- CombinedRearrangementDetailsPostTxCD8[!duplicated(CombinedRearrangementDetailsPostTxCD8$nucleotide),]

# List of XDRTCCs can be obtained from Identifying_XDRTCCs.R script
XDRTCCsinPostTxCD8 = CombinedRearrangementDetailsPostTxCD8[CombinedRearrangementDetailsPostTxCD8$aminoAcid %in% cd8_XDRTCC_df$`Amino Acid`, ]
# compute nucleotide-per-aa for xdrtccs postTx CD8
length(unique(XDRTCCsinPostTxCD8$nucleotide))/length(unique(XDRTCCsinPostTxCD8$aminoAcid))

# For non-XDRTCCs, need to download TCRb data from Adaptive Biotech immuneACCESS Database in Combined Rearrangements format
# Non-XDRTCCs is defined as present in PreTx Unstim, but not in any Stim population, so we need to filter CombinedRearrangementDetailsPostTxCD8
# configurable column names 
aa_col_comb   <- "Amino Acid"                    # in CombinedRearrangements
pretx_col     <- "Unstim-XKLT001-POD0-PBMC"      # PreTx Unstim column
aa_col_detail <- "aminoAcid"                       # in CombinedRearrangementDetailsPostTxCD8

#  identify all Stim columns 
stim_cols <- grep("^Stim", colnames(CombinedRearrangements), value = TRUE)

#  amino acids that are present at PreTx (>0) and absent from ALL Stim columns 
eligible_aa <- with(CombinedRearrangements, {
  # count how many Stim columns have >0 for each row
  stim_pos <- if (length(stim_cols)) {
    rowSums(CombinedRearrangements[, stim_cols, drop = FALSE] > 0, na.rm = TRUE)
  } else {
    rep(0L, nrow(CombinedRearrangements))
  }
  keep <- (!is.na(CombinedRearrangements[[aa_col_comb]])) &
    (CombinedRearrangements[[pretx_col]] > 0) &
    (stim_pos == 0)
  CombinedRearrangements[[aa_col_comb]][keep]
})

NonDRTCCsinPostTxCD8 <- CombinedRearrangementDetailsPostTxCD8[
  CombinedRearrangementDetailsPostTxCD8[[aa_col_detail]] %in% eligible_aa,
  ,
  drop = FALSE
]

# compute nucleotide-per-aa for non-xdrtccs post tx CD8
length(unique(NonDRTCCsinPostTxCD8$nucleotide))/length(unique(NonDRTCCsinPostTxCD8$aminoAcid))

# Now look at rest -- or clones that are in the Post-Tx but neither XDRTCCs or non-XDRTCCs

# define column names explicitly 
aa_col_posttx <- "aminoAcid"    # in CombinedRearrangementDetailsPostTxCD8
aa_col_xdrtcc <- "Amino Acid"   # in cd8_XDRTCC_df
aa_col_nondrtcc <- "aminoAcid"  # in NonDRTCCsinPostTxCD8

# --- collect sequences to remove ---
to_remove <- c(
  cd8_XDRTCC_df[[aa_col_xdrtcc]],
  NonDRTCCsinPostTxCD8[[aa_col_nondrtcc]]
)

# --- subset, keeping only rows NOT in to_remove ---
RestinPostTxCD8 <- CombinedRearrangementDetailsPostTxCD8[
  !(CombinedRearrangementDetailsPostTxCD8[[aa_col_posttx]] %in% to_remove),
  ,
  drop = FALSE
]

# compute nucleotide-per-aa for rest in PostTx CD8
length(unique(RestinPostTxCD8$nucleotide))/length(unique(RestinPostTxCD8$aminoAcid))

### For PreTx ---------

# cd4 xdrtccs
XDRTCCsinPreTxCD4 = data_list$`Unstim-XKLT001-POD0-PBMC`[data_list$`Unstim-XKLT001-POD0-PBMC`$aminoAcid %in% cd4_XDRTCC_df$`Amino Acid`, ]
# compute nucleotide-per-aa for cd4 xdrtccs pretx
length(unique(XDRTCCsinPreTxCD4$nucleotide))/length(unique(XDRTCCsinPreTxCD4$aminoAcid))
# cd8 xdrtccs
XDRTCCsinPreTxCD8 = data_list$`Unstim-XKLT001-POD0-PBMC`[data_list$`Unstim-XKLT001-POD0-PBMC`$aminoAcid %in% cd8_XDRTCC_df$`Amino Acid`, ]
# compute nucleotide-per-aa for cd8 xdrtccs pretx
length(unique(XDRTCCsinPreTxCD8$nucleotide))/length(unique(XDRTCCsinPreTxCD8$aminoAcid))
# non-xdrtccs
# --- column names ---
aa_col_comb <- "Amino Acid"               # in CombinedRearrangements
pretx_col   <- "Unstim-XKLT001-POD0-PBMC" # PreTx Unstim column
aa_col_pre  <- "aminoAcid"                # in data_list$Unstim-XKLT001-POD0-PBMC

# --- find stim columns ---
stim_cols <- grep("^Stim", colnames(CombinedRearrangements), value = TRUE)

# --- eligible amino acids: >0 at PreTx and ==0 in all Stim ---
eligible_aa <- with(CombinedRearrangements, {
  stim_pos <- if (length(stim_cols)) {
    rowSums(CombinedRearrangements[, stim_cols, drop = FALSE] > 0, na.rm = TRUE)
  } else {
    rep(0L, nrow(CombinedRearrangements))
  }
  keep <- (!is.na(CombinedRearrangements[[aa_col_comb]])) &
    (CombinedRearrangements[[pretx_col]] > 0) &
    (stim_pos == 0)
  CombinedRearrangements[[aa_col_comb]][keep]
})

# --- subset the PreTx dataframe with those amino acids ---
NonXDRTCCsinPreTx <- data_list$`Unstim-XKLT001-POD0-PBMC`[
  data_list$`Unstim-XKLT001-POD0-PBMC`[[aa_col_pre]] %in% eligible_aa,
  ,
  drop = FALSE
]

# compute nucleotide-per-aa for non-xdrtccs preTx
length(unique(NonXDRTCCsinPreTx$nucleotide)) / length(unique(NonXDRTCCsinPreTx$aminoAcid))



