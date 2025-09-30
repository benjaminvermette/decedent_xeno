# From CombinedRearrangements file created using Identifying_XDRTCCs script, 
# use the below code to identify Non-XDRTCCs

# PreTx column
pretx_col <- "Unstim-XKLT001-POD0-PBMC"

# Stim columns
stim_cols <- grep("^Stim-", names(CombinedRearrangements), value = TRUE)

# Non-XDRTCCs: >0 in PreTx, and ==0 in *all* Stim columns
NonXDRTCC_CombinedRearrangements <- CombinedRearrangements[
  CombinedRearrangements[[pretx_col]] > 0 &
    rowSums(CombinedRearrangements[, stim_cols, drop = FALSE] > 0, na.rm = TRUE) == 0,
  ,
  drop = FALSE
]

# using the above dataframe, you can easily track non-Xdrtccs in any population
# e.g. cumulative frequency in POD33 CD8

sum(NonXDRTCC_CombinedRearrangements$`Unstim-XKLT001-POD33-CD8`)

# or number of unique clones in POD33 CD8

length(which((NonXDRTCC_CombinedRearrangements$`Unstim-XKLT001-POD33-CD8`>0)))
       