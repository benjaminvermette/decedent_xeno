library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)

## --- define unmapped and mapped (short) timepoint columns ---
col_pre <- "Unstim-XKLT001-POD0-PBMC"    
col_p28 <- "Unstim-XKLT001-POD28-CD8"
col_p33 <- "Unstim-XKLT001-POD33-CD8"
col_p49 <- "Unstim-XKLT001-POD49-CD8"

# named vector: names = short labels used in plot; values = actual column names in the data
time_cols <- c(PreTx = col_pre, POD28 = col_p28, POD33 = col_p33, POD49 = col_p49)

## --- helper: column-wise normalization to sum 1 ---
norm1 <- function(x) { s <- sum(x, na.rm = TRUE); if (s > 0) x/s else x }

## --- normalize each unstim CD8 timepoint column independently ---
CombinedRearrangements[, time_cols] <- lapply(
  as.data.frame(CombinedRearrangements[, time_cols, drop = FALSE]),
  norm1
)

## --- pick top-5 per timepoint, then union ---
top_n_per_col <- function(v, n = 5) {
  idx <- order(v, decreasing = TRUE, na.last = NA)
  idx <- idx[v[idx] > 0]   # ignore zeros
  head(idx, n)
}

top_idx_list   <- lapply(time_cols, function(cn) top_n_per_col(CombinedRearrangements[[cn]], 5))
top_union_idx  <- sort(unique(unlist(top_idx_list)))

## --- build wide table with Amino Acid as ID and short timepoint names ---
entities <- CombinedRearrangements[["Amino Acid"]][top_union_idx]

AlluvialPlotTop5ClonesTry <- CombinedRearrangements[top_union_idx, time_cols, drop = FALSE]
# rename columns to short labels
colnames(AlluvialPlotTop5ClonesTry) <- names(time_cols)
# add ID column first
AlluvialPlotTop5ClonesTry <- tibble::tibble(Entity = entities) %>%
  bind_cols(as.data.frame(AlluvialPlotTop5ClonesTry))

# Multiply unsorted PreTx values by 5 based on flow cytometry proportion of CD8s within CD3s
AlluvialPlotTop5ClonesTry$PreTx <- AlluvialPlotTop5ClonesTry$PreTx * 5

## --- to long format ---
df_long <- tidyr::pivot_longer(
  AlluvialPlotTop5ClonesTry,
  cols = c("PreTx", "POD28", "POD33", "POD49"),
  names_to = "TimePoint",
  values_to = "Frequency"
)

## --- colors ---
unique_entities <- unique(df_long$Entity)
base_colors <- c(
  "#F7DC6F", "#FF1493", "#3357FF", "#FF5733", "#8B008B", "#48C9B0",
  "#F39C12", "#C70039", "#1ABC9C", "#E74C3C", "#33FF57", "#2E86C1"
)
custom_colors <- setNames(rep(base_colors, length.out = length(unique_entities)), unique_entities)

## --- plot ---
ggplot(df_long, aes(x = TimePoint, stratum = Entity, alluvium = Entity, y = Frequency, fill = Entity)) +
  geom_alluvium(alpha = 0.8) +
  geom_stratum() +
  scale_x_discrete(limits = c("PreTx", "POD28", "POD33", "POD49")) +
  scale_y_continuous(
    sec.axis = sec_axis(~ ., name = "Frequency"),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = custom_colors) +
  labs(x = NULL, y = NULL) +
  theme_void() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.right = element_line(color = "black"),
    axis.text.y.right = element_text(color = "black", size = 20),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 20)
  )

## ---- add XDRTCC status to alluvial ----

AlluvialPlotTop5ClonesTry$DRTCC_status <- ifelse(
  AlluvialPlotTop5ClonesTry$Entity %in% cd8_XDRTCC_df$`Amino Acid`, 
  "CD8 XDRTCC",
  ifelse(AlluvialPlotTop5ClonesTry$PreTx > 0, "non-XDRTCC", NA)
)

# drop Entities that are neither XDRTCC nor present at PreTx:
AlluvialPlotTop5ClonesTry2 <- AlluvialPlotTop5ClonesTry[!is.na(AlluvialPlotTop5ClonesTry$DRTCC_status), ]

## ---- pivot and plot with legend by DRTCC status ----
df_long <- tidyr::pivot_longer(
  AlluvialPlotTop5ClonesTry2,
  cols = c("PreTx", "POD28", "POD33", "POD49"),
  names_to = "TimePoint",
  values_to = "Frequency"
)

# Build a per-Entity ordering key from the PreTx axis
# - XDRTCCs get a higher base so they draw on top
order_tbl <- AlluvialPlotTop5ClonesTry2 %>%
  transmute(
    Entity,
    DRTCC_status,
    PreTx = replace_na(PreTx, 0),
    order_key = if_else(DRTCC_status == "CD8 XDRTCC", 2e6, 1e6) + rank(PreTx, ties.method = "first")
  )

# Join order + (re)fix factor levels
df_long <- df_long %>%
  select(-DRTCC_status) %>%                       # avoid duplicate column on join
  left_join(order_tbl[, c("Entity","DRTCC_status","order_key")], by = "Entity") %>%
  mutate(
    DRTCC_status = factor(DRTCC_status, levels = c("non-XDRTCC", "CD8 XDRTCC"))
  ) %>%
  arrange(TimePoint, order_key)                   # draw non-X first, XDRTCC last (on top)

ggplot(df_long, aes(
  x = TimePoint,
  stratum = Entity,
  alluvium = Entity,
  y = Frequency,
  fill = DRTCC_status,
  order = order_key          # <- ensures XDRTCCs stack above at PreTx
)) +
  geom_alluvium(alpha = 0.9, color = "white", size = 0.5) +
  geom_stratum(size = 0.6, color = "grey") +
  scale_x_discrete(limits = c("PreTx", "POD28", "POD33", "POD49")) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Frequency"), expand = c(0, 0)) +
  scale_fill_manual(values = c("CD8 XDRTCC" = "#002A99", "non-XDRTCC" = "#0099FF")) +
  labs(x = NULL, y = NULL, fill = "CD8 DRTCC Status") +
  theme_void() +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.line.y.right = element_line(color = "black"),
    axis.text.y.right = element_text(color = "black", size = 20),
    axis.title.y.right = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 20)
  )

