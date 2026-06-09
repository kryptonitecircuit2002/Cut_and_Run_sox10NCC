library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggalluvial) 

load_state_file <- function(file, state, tp) {
  df <- read_tsv(
    file,
    col_names = FALSE,
    na = c(".")
  )
  colnames(df) <- c("chr","start","end", "b_chr","b_start","b_end", "b_name", "b_score", "b_strand")
  
  df <- df %>%
    mutate(
      peak_id = paste(chr, start, end, sep="_"),
      STATE = state,
      TP = tp,
      active = ifelse(is.na(b_chr), 0, 1)
    ) %>%
    select(peak_id, TP, STATE, active)
  
  return(df)
}

files <- list(
  list(".../h2_activeenh_18hpf.bed", "H2AZ_ActiveEnh", "16hpf"),
  list(".../h2_activeenh_24hpf.bed", "H2AZ_ActiveEnh", "24hpf"),
  list(".../h2_activeenh_48hpf.bed", "H2AZ_ActiveEnh", "48hpf"),
  list(".../h2_activeprom_18hpf.bed", "H2AZ_ActiveProm", "16hpf"),
  list(".../h2_activeprom_24hpf.bed", "H2AZ_ActiveProm", "24hpf"),
  list(".../h2_activeprom_48hpf.bed", "H2AZ_ActiveProm", "48hpf"),
  list(".../h2_bivalent_18hpf.bed", "H2AZ_Bivalent", "16hpf"),
  list(".../h2_bivalent_24hpf.bed", "H2AZ_Bivalent", "24hpf"),
  list(".../h2_bivalent_48hpf.bed", "H2AZ_Bivalent", "48hpf"),
  list(".../h2_repressed_18hpf.bed", "H2AZ_Repressed", "16hpf"),
  list(".../h2_repressed_24hpf.bed", "H2AZ_Repressed", "24hpf"),
  list("../h2_repressed_48hpf.bed", "H2AZ_Repressed", "48hpf"),
  
  list(".../h3_activeenh_18hpf.bed", "H33_ActiveEnh", "16hpf"),
  list(".../h3_activeenh_24hpf.bed", "H33_ActiveEnh", "24hpf"),
  list(".../h3_activeenh_48hpf.bed", "H33_ActiveEnh", "48hpf"),
  list(".../h3_activeprom_18hpf.bed", "H33_ActiveProm", "16hpf"),
  list(".../h3_activeprom_24hpf.bed", "H33_ActiveProm", "24hpf"),
  list(".../h3_activeprom_48hpf.bed", "H33_ActiveProm", "48hpf"),
  list(".../h3_bivalent_18hpf.bed", "H33_Bivalent", "16hpf"),
  list(".../h3_bivalent_24hpf.bed", "H33_Bivalent", "24hpf"),
  list(".../h3_bivalent_48hpf.bed", "H33_Bivalent", "48hpf"),
  list(".../h3_repressed_18hpf.bed", "H33_Repressed", "16hpf"),
  list(".../h3_repressed_24hpf.bed", "H33_Repressed", "24hpf"),
  list(".../h3_repressed_48hpf.bed", "H33_Repressed", "48hpf"),
  
  list(".../h2h3_activeenh_18hpf.bed", "H2H3_ActiveEnh", "16hpf"),
  list(".../h2h3_activeenh_24hpf.bed", "H2H3_ActiveEnh", "24hpf"),
  list(".../h2h3_activeenh_48hpf.bed", "H2H3_ActiveEnh", "48hpf"),
  list(".../h2h3_activeprom_18hpf.bed", "H2H3_ActiveProm", "16hpf"),
  list(".../h2h3_activeprom_24hpf.bed", "H2H3_ActiveProm", "24hpf"),
  list(".../h2h3_activeprom_48hpf.bed", "H2H3_ActiveProm", "48hpf"),
  list(".../h2h3_bivalent_18hpf.bed", "H2H3_Bivalent", "16hpf"),
  list(".../h2h3_bivalent_24hpf.bed", "H2H3_Bivalent", "24hpf"),
  list(".../h2h3_bivalent_48hpf.bed", "H2H3_Bivalent", "48hpf"),
  list(".../h2h3_repressed_18hpf.bed", "H2H3_Repressed", "16hpf"),
  list(".../h2h3_repressed_24hpf.bed", "H2H3_Repressed", "24hpf"),
  list(".../h2h3_repressed_48hpf.bed", "H2H3_Repressed", "48hpf")
)

ptm_variant_list <- lapply(files, function(x) load_state_file(x[[1]], x[[2]], x[[3]]))
ptm_variant_all <- bind_rows(ptm_variant_list)
ptm_variant_all <- ptm_variant_all %>%
  group_by(peak_id, TP, STATE) %>%
  summarise(active = max(active), .groups = "drop")

ptm_variant_wide <- ptm_variant_all %>%
  pivot_wider(
    names_from = STATE,
    values_from = active,
    values_fill = 0
  )

classify_state <- function(H2AZ_ActiveEnh, H2AZ_ActiveProm, H2AZ_Bivalent, H2AZ_Repressed, 
                           H33_ActiveEnh, H33_ActiveProm, H33_Bivalent, H33_Repressed,
                           H2H3_ActiveEnh, H2H3_ActiveProm, H2H3_Bivalent, H2H3_Repressed) {
  total_active <- sum(H2AZ_ActiveEnh, H2AZ_ActiveProm,H2AZ_Bivalent, H2AZ_Repressed, 
                      H33_ActiveEnh, H33_ActiveProm, H33_Bivalent, H33_Repressed, 
                      H2H3_ActiveEnh, H2H3_ActiveProm, H2H3_Bivalent, H2H3_Repressed)
  if (total_active == 1) {
    if (H2AZ_ActiveEnh == 1 | H2AZ_ActiveProm == 1) return("H2A.Z Active")
    if (H2AZ_Bivalent == 1) return("H2A.Z Poised")
    if (H2AZ_Repressed == 1) return("H2A.Z Repressed")
    if (H33_ActiveEnh == 1 | H33_ActiveProm == 1) return("H3.3 Active")
    if (H33_Bivalent == 1) return("H3.3 Poised")
    if (H33_Repressed == 1) return("H3.3 Repressed")
    if (H2H3_ActiveEnh == 1 | H2H3_ActiveProm == 1) return("Double Variant Active")
    if (H2H3_Bivalent == 1) return("Double Variant Poised")
    if (H2H3_Repressed == 1) return("Double Variant Repressed")
  }
  return("Unassigned")
}
ptm_variant_wide$State <- mapply(classify_state,
                                 ptm_variant_wide$H2AZ_ActiveEnh,
                                 ptm_variant_wide$H2AZ_ActiveProm,
                                 ptm_variant_wide$H2AZ_Bivalent,
                                 ptm_variant_wide$H2AZ_Repressed,
                                 ptm_variant_wide$H33_ActiveEnh,
                                 ptm_variant_wide$H33_ActiveProm,
                                 ptm_variant_wide$H33_Bivalent,
                                 ptm_variant_wide$H33_Repressed,
                                 ptm_variant_wide$H2H3_ActiveEnh,
                                 ptm_variant_wide$H2H3_ActiveProm,
                                 ptm_variant_wide$H2H3_Bivalent,
                                 ptm_variant_wide$H2H3_Repressed)
state_levels <- c("H2A.Z Active", "H2A.Z Poised", "H2A.Z Repressed", "H3.3 Active", "H3.3 Poised", "H3.3 Repressed", "Double Variant Active", "Double Variant Poised", "Double Variant Repressed", "Unassigned")

ptm_variant_matrix <- ptm_variant_wide %>%
  dplyr::select(peak_id, TP, State) %>%
  pivot_wider(
    names_from = TP,
    values_from = State,
    values_fill = "Unassigned"
  )
ptm_variant_transitions <- ptm_variant_matrix %>%
  count(`16hpf`, `24hpf`, `48hpf`, name = "Freq")

ptm_variant_transitions <- ptm_variant_transitions %>%
  mutate(
    `16hpf` = factor(`16hpf`, levels = state_levels),
    `24hpf` = factor(`24hpf`, levels = state_levels),
    `48hpf` = factor(`48hpf`, levels = state_levels)
  )

state_colors <- c(
  "H2A.Z Active" = "#4C72B0",  
  "H2A.Z Poised" = "#8FAAD9", 
  "H2A.Z Repressed" = "#2F4B7C",  
  "H3.3 Active" = "#DD8452",  
  "H3.3 Poised" = "#F2B48C",
  "H3.3 Repressed" = "#B45A2A", 
  "Double Variant Active" = "#55A868",
  "Double Variant Poised" = "#9AD1A3",
  "Double Variant Repressed" = "#2E7D4F", 
  "Unassigned" = "#BDBDBD"
)

ggplot(
  ptm_variant_transitions,
  aes(
    axis1 = `16hpf`,
    axis2 = `24hpf`,
    axis3 = `48hpf`,
    y = Freq
  )
) +
  geom_alluvium(aes(fill = `24hpf`), width = 1/12, alpha = 0.8) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/12, color = "black", size = 0.5) +
  scale_fill_manual(values = state_colors, drop = FALSE) +
  scale_x_discrete(
    limits = c("18hpf", "24hpf", "48hpf"),
    labels = c("16 hpf", "24 hpf", "48 hpf")
  ) +
  theme_void() +
  theme(
    text = element_text(family = "Times New Roman", size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 11, face = "italic"),
    legend.position = "right",
    axis.text.x = element_text(size = 12)
  ) +
  labs(
    title = "PTM Variant Chromatin State Transitions During Development",
    subtitle = paste0("Following ", sum(ptm_variant_transitions$Freq), " ptm-variant peaks from 16 hpf"),
    y = "Number of Peaks",
    x = "",
    fill = "Variant+PTM State"
  )


