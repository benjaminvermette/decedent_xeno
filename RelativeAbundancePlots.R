library(immunarch)
library(tidyverse)
library(readr)
library(cdr3tools)
### Export Unstim Samples from Adaptive Website using v2 ----

library(readr)

# Path to zip file from Adaptive
zip_path <- "C:/Users/bv2269/Downloads/sampleExport.zip"

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
#### Adapt dataframes to immunarch format ---
# need to define format_immunarch_adaptive function using format_immunarch_adaptive.R script available in Github repository

data_list <- format_immunarch_adaptive(data_list)

# Can then use data_list on various immunarch functions
# Number of Unique Clones

Clon_num  <- repExplore(
  data_list,
  .method = c("clones"),
  .col = c("nt"),
  .coding = TRUE)

vis(Clon_num)

# Expansion

imm_rare <- repClonality(data_list, .method = "rare")

vis(imm_rare)

imm_hom <- repClonality(data_list,
                        .method = "homeo",
                        .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1)
)
vis(imm_hom)
