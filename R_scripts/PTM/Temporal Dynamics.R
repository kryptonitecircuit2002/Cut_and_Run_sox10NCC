library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggalluvial)

load_file <- function(file, ptm, tp) {
  df <- read_tsv(
    file,
    col_names = FALSE,
    na = c(".")
  )
  colnames(df) <- c("chr","start","end","b_chr","b_start","b_end","b_name","b_score", "b_strand", "signalValue", "pvalue",
                    "qvalue", "peak")
  df <- df %>%
    mutate(
      peak_id = paste(chr, start, end, sep="_"),
      PTM = ptm,
      TP = tp,
      active = ifelse(is.na(b_chr), 0, 1)
    ) %>%
    dplyr::select(peak_id, TP, PTM, active)
  
  return(df)
}

files <- list(
  list("../h3k27a_16hpf.bed","H3K27ac","16hpf"),
  list("../h3k27a_24hpf.bed","H3K27ac","24hpf"),
  list("../h3k27a_48hpf.bed","H3K27ac","48hpf"),
  list("../h3k27m_16hpf.bed","H3K27me3","16hpf"),
  list("../h3k27m_24hpf.bed","H3K27me3","24hpf"),
  list("../h3k27m_48hpf.bed","H3K27me3","48hpf"),
  list("../h3k4m_16hpf.bed","H3K4me3","16hpf"),
  list("../h3k4m_24hpf.bed","H3K4me3","24hpf"),
  list("../h3k4m_48hpf.bed","H3K4me3","48hpf")
)

ptm_list <- lapply(files, function(x) load_file(x[[1]], x[[2]], x[[3]]))
ptm_all <- bind_rows(ptm_list)
ptm_all <- ptm_all %>%
  group_by(peak_id, TP, PTM) %>%
  summarise(active = max(active), .groups = "drop")


ptm_wide <- ptm_all %>%
  pivot_wider(
    names_from = PTM,
    values_from = active,
    values_fill = 0
  )

classify_ptm <- function(H3K27ac, H3K4me3, H3K27me3) {
  if (H3K4me3 == 1 & H3K27me3 == 1) return("Bivalent")
  if (H3K27me3 == 0 & (H3K27ac == 1 | H3K4me3 == 1)) return("Active")
  if (H3K27ac == 0 & H3K4me3 == 0 & H3K27me3 == 1) return("Repressed")
  
  return("No PTM")
}

ptm_wide$State <- mapply(classify_ptm,
                         ptm_wide$H3K27ac,
                         ptm_wide$H3K4me3,
                         ptm_wide$H3K27me3)

ptm_state_matrix <- ptm_wide %>%
  select(peak_id, TP, State) %>%
  pivot_wider(
    names_from = TP,
    values_from = State,
    values_fill = "Unassigned"
  )

PTM_LEVELS <- c(
  "Bivalent",
  "Active",
  "Repressed",
  "no PTM"
)

state_colors <- c(
  "Bivalent" = "#C77CFF",
  "Active" = "#F8766D",
  "Repressed" = "#00BFC4",
  "No PTM" = "grey70"
)


df_alluvial <- ptm_state_matrix %>%
  count(`16hpf`, `24hpf`, `48hpf`, name = "Freq") %>%
  complete(
    `16hpf` = PTM_LEVELS,
    `24hpf` = PTM_LEVELS,
    `48hpf` = PTM_LEVELS,
    fill = list(Freq = 0)
  )


state_colors <- c(
  "H3K27ac" = "#F8766D",
  "H3K4me3" = "#7CAE00",
  "H3K27me3" = "#00BFC4",
  "Bivalent" = "#C77CFF",
  "Ambivalent" = "#E68613",
  "Unassigned" = "grey70"
)

ggplot(
  df_alluvial,
  aes(
    axis1 = `16hpf`,
    axis2 = `24hpf`,
    axis3 = `48hpf`,
    y = Freq
  )
) +
  geom_alluvium(aes(fill = `16hpf`), width = 1/12) +
  
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/12) +
  
  scale_fill_manual(values = state_colors, drop = FALSE) +
  
  scale_x_discrete(
    limits = c("16hpf", "24hpf", "48hpf"),
    labels = c("16 hpf", "24 hpf", "48 hpf")
  ) +
  
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "PTM State Dynamics Over Time",
    y = "Peaks",
    x = "",
    fill = "PTM State"
  )


