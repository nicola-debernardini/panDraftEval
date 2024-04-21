######################################################  @
# Script to reproduce the analysis in figure 3:
# 1. Selection of optimal MRF threshold for the UHGG dataset
# 2. Selection of optimal MRF threshold for the OMD dataset

# Here is provided the code used to compute the statistics and plots
# Due to the storage space required by GEMs only binary matrix (i.e. rxnXmod) will be provided here.
# Contact the corresponding author to have access to the GEMs.

# --- pseudo-code ---
#  load rxnXmod table of reference
#  load rxnXmod table of gapfilled pan-GEM
#  for every GEM in the reference GEM list:
#   for every pan-GEM generated for that specific species with different MRF thresholds:
#       compute statistics to quantify its structural quality 

# Load R-packages
library(tidyr)
library(dplyr)
library(ggpubr)
library(data.table)
library(ggplot2)
library(RColorBrewer)

# Little helpers
source("./scr/function_collection.R")

# compare the pan-GEMs with the reference models and compute statistics
compute_structural_QC <- function(iso.binary_table.list, rxnXmod_dir_path) {
  all_dist_df <- data.frame(comp_lv = numeric(0), mod_id = character(), ref_mod_id = character(), dist_tp = character(), mod_category = character(), value = numeric(0))
  for (mod in names(iso.binary_table.list)) {
      cat(paste("\n", mod))

      # determine species-level gold standard reactome
      inIso_dt <- iso.binary_table.list[[mod]]
      reference_reactions <- inIso_dt$rxn
      Isolate_cols <- names(inIso_dt)[names(inIso_dt)!="rxn"] 
      freq_ref_rxn <- inIso_dt %>% # frequency of rxns
          select(all_of(Isolate_cols)) %>%
          mutate(b_i = rowSums(.)/length(Isolate_cols),
              rxn = inIso_dt$rxn) %>%
          select(b_i, rxn)
      
      # load data for panModel
      sp_rxnXmod_path <- file.path(rxnXmod_dir_path, mod, "rxnXmod.tsv")
      sp_rxnXmod <- read.table(sp_rxnXmod_path, sep="\t", header = TRUE)
      sp_rxnXmod <- as.data.table(sp_rxnXmod)

      for (th in seq(0, 1, 0.01)) {
          if (th==0.01) { next }
          # check if the file exist, which means if originally the num of MAG in the subset was higher than the 30 MAGs
          panModel_name <- paste0("X", th, "_panModel")
          col_to_extract <- c(panModel_name, "rxn")
          sp_rxnXmod <- as.data.table(sp_rxnXmod)
          pan.model_rxn <- sp_rxnXmod[, ..col_to_extract]

          # Compute statictics
          # Tsgb <- Tsgb_cal(pan.model_rxn, freq_ref_rxn, panModel_name)
          # dist_val  <- distance_cal(pan.model_rxn, freq_ref_rxn, dist_type)   
          # dist_val <- dist_val[1]
          predicted_reactions <- pan.model_rxn[["rxn"]][as.logical(pan.model_rxn[[panModel_name]])]
          metrics_table <- compute_contigency_table(reference_reactions, predicted_reactions) # compute contingency table for a pan.mod
          stat_val <- metrics_table[which.stat][[1]]
          
          # Populate the dataframe
          # new_row_Tsgb <- list(comp_lv = th,
          #                 mod_id = panModel_name, 
          #                 ref_mod_id = mod,
          #                 dist_tp = "Tsgb", 
          #                 mod_category = "pan.mod", 
          #                 value = Tsgb)
          # new_row_dist <- list(comp_lv = th,
          #                 mod_id = panModel_name, 
          #                 ref_mod_id = mod,
          #                 dist_tp = dist_type, 
          #                 mod_category = "pan.mod", 
          #                 value = dist_val)
          new_row_stat <- list(comp_lv = th,
                          mod_id = panModel_name, 
                          ref_mod_id = mod,
                          dist_tp = which.stat, 
                          mod_category = "pan.mod", 
                          value = stat_val)
          all_dist_df <- rbind(all_dist_df, new_row_stat) # new_row_Tsgb, new_row_dist
      }
  }
  return(all_dist_df)
}

# ---------------------------------------------------------------------------------------------------
# Parameters:
dist_type <- "binary"
which.stat <- "f1_score" # fB_score f1_score recall precision accuracy FDR

output.dir <- "dat/fig/"

# reference-GEMs file paths
iso.binary.rxn.table_omd      <- "./dat/omd/rxnXmod/iso/allISO_rxnXmod.txt"       # OMD
iso.binary.rxn.table_uhgg      <- "./dat/uhgg/rxnXmod/iso/allISO_rxnXmod.txt"    # UHGG
# pan-GEMs file paths
rxnXmod_dir_path_omd <- "./dat/omd/rxnXmod/panDraft"        # OMD
rxnXmod_dir_path_uhgg <- "./dat/uhgg/rxnXmod/panDraft"     # UHGG

# load in bulk the binary matrix (i.e. rxnXmod) of all the reference genomes
iso.binary_table.list_omd     <- load_files_from_paths(iso.binary.rxn.table_omd, ".tsv")
iso.binary_table.list_uhgg     <- load_files_from_paths(iso.binary.rxn.table_uhgg, ".tsv")

# Process the UHGG --------------------------------------------------------------------------------
all_dist_df <- compute_structural_QC(iso.binary_table.list_uhgg, rxnXmod_dir_path_uhgg)

# Read metadata
completeness.level.table  <- "./dat/uhgg/metadata/genomes-all_metadata.csv"
metadata <- read.table(completeness.level.table, sep = ",", header = TRUE, row.names = "Genome") # dec = ","
columns_to_convert_numeric <- c("Completeness")
metadata[columns_to_convert_numeric] <- lapply(metadata[columns_to_convert_numeric], as.numeric)
metadata$Species_name <- rownames(metadata)

# merge the dataset
merged_df <- merge(all_dist_df, metadata, by.x="ref_mod_id", by.y="Species_name")
plot_prob_dist_df <- parse_GTDBtaxonomy(merged_df, "Lineage")
plot_prob_dist_df <- refine_GTDBtaxonomy(plot_prob_dist_df)

# prepare dataset for plotting
desired_order <- plot_prob_dist_df %>% # get the desired order of phylum
  group_by(phylum) %>%
  summarise(count = length(unique(species))) %>%
  arrange(desc(count)) %>%
  select(phylum)
desired_order <- as.list(desired_order)$phylum
# assign the value other if more that eight phylum where present
plot_dt <- plot_prob_dist_df %>%
    mutate(phylum_color = ifelse(phylum %in% desired_order, phylum, "other"),
    phylum_color_noStat = ifelse(phylum %in% desired_order, phylum, "other")) 
plot_dt <- plot_dt %>% 
    group_by(phylum_color) %>% 
    mutate(count_SpXphylum = length(unique(species))) %>%
    arrange(desc(count_SpXphylum)) 
plot_dt$phylum_color <- paste0(plot_dt$phylum_color, " (s=", plot_dt$count_SpXphylum, ")")
plot_dt$phylum_color <- as.factor(plot_dt$phylum_color)
levels(plot_dt$phylum_color) <- unique(plot_dt$phylum_color)
# levels(plot_dt$phylum_color) <- c("Firmicutes_A (s=22)", "Firmicutes (s=9)", "Bacteroidota (s=4)", "Proteobacteria (s=2)", "Actinobacteriota (s=2)", "Verrucomicrobiota (s=1)") # UHGG

# Plot
set1_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#F781BF") # UHGG
# set1_palette <- brewer.pal(9, "Set1")

line_plot <- ggplot(plot_dt, aes(x = comp_lv, y = value, group = species, color = phylum_color)) +
  geom_line(alpha = 0.8, show.legend = FALSE, linewidth = .3) + 
  labs(x = "MRF threshold (%)", y = "F1 score") +
  theme_minimal()+ 
  scale_color_manual(values = set1_palette) + 
  guides(color = guide_legend(title = "Phylum")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
  axis.text.y = element_text(size = 8),
  legend.text = element_text(size = 6)
  )

# determine the MRF threshold maximizing the F1-score as average of all the species
n_mod <- length(unique(all_dist_df$ref_mod_id))
find_best_perf <- all_dist_df %>%
  group_by(comp_lv) %>%
  summarise(avgPerformance = sum(value)/n_mod,
    stdPerformance = sd(value)/n_mod) %>%
  mutate(species = NA,
    genera = NA)
max_index <- which.max(find_best_perf$avgPerformance)
find_best_perf$avgPerformance[max_index]
find_best_perf$stdPerformance[max_index]
max_index <- find_best_perf$comp_lv[max_index]
# determine the MRF threshold maximizing the F1-score of each species
find_best_perf_Xsp <- all_dist_df %>%
  group_by(ref_mod_id) %>%
  mutate(bestXsp = max(value),
  whichBest = value==bestXsp) %>%
  mutate(species = NA,
    genera = NA)
find_best_perf_Xsp <- find_best_perf_Xsp[find_best_perf_Xsp$whichBest,] %>%
  arrange(comp_lv) %>%
  distinct(ref_mod_id, .keep_all = TRUE)
# Add vertical lines for each MRF thresholds
comp_lv_values <- unique(find_best_perf_Xsp$comp_lv)
for (comp_lv_value in comp_lv_values) {
  line_plot <- line_plot +
    geom_vline(xintercept = comp_lv_value, linetype = "dashed", color = "grey", linewidth = .3)
}
line_plot <- line_plot +
  geom_line(data=find_best_perf, aes(x = comp_lv, y = avgPerformance), color = "#000000", linewidth = .6) + # red: "#ad1d00"
  geom_vline(xintercept = max_index, linetype = "dashed", color = "black", linewidth = .4) 

# save plot
lineplot_UHGG <- line_plot 
# cat(file.path(output.dir, paste0(which.stat, "_UHGG_wlegend.pdf")))
# ggsave(file.path(output.dir, paste0(which.stat, "_UHGG_wlegend.pdf")), lineplot_UHGG, width = 8, height = 8, units = "cm", dpi = 300)

# process the OMD -----------------------------------------------------------------------------------
all_dist_df <- compute_structural_QC(iso.binary_table.list_omd, rxnXmod_dir_path_omd)

# Read metadata
oceans_metadata_fn        <- "./dat/omd/metadata/genomes-summary.csv"
oceans_metadata <- read.table(oceans_metadata_fn, sep = ",", header = TRUE, dec = ".")
oceans_metadata <- oceans_metadata %>% # separate MAGs from reference genome and SAGs
  mutate(Genome_type = if_else(grepl("_METAG_", Genome), "MAG", if_else(grepl("_REFG_", Genome), "REF", if_else(grepl("_SAGS_", Genome), "SAG", ""))))
oceans_metadata <- oceans_metadata %>% # specify the sample name
  mutate(sample = sapply(strsplit(oceans_metadata[["Genome"]], "_"), function(x) paste(x[1], x[2], x[3], sep="_")))
columns_to_convert_numeric <- c("Completeness")
oceans_metadata[columns_to_convert_numeric] <- lapply(oceans_metadata[columns_to_convert_numeric], as.numeric)
oceans_metadata$Species_name <- oceans_metadata$Genome
metadata <- oceans_metadata

# merge the dataset
merged_df <- merge(all_dist_df, metadata, by.x="ref_mod_id", by.y="Species_name")
plot_prob_dist_df <- parse_GTDBtaxonomy(merged_df, "GTDB.Taxonomy")
plot_prob_dist_df <- refine_GTDBtaxonomy(plot_prob_dist_df)

# prepare dataset for plotting
desired_order <- plot_prob_dist_df %>% # get the desired order of phylum
  group_by(phylum) %>%
  summarise(count = length(unique(species))) %>%
  arrange(desc(count)) %>%
  select(phylum)
desired_order <- as.list(desired_order)$phylum
# assign the value other if more that eight phylum where present
plot_dt <- plot_prob_dist_df %>%
    mutate(phylum_color = ifelse(phylum %in% desired_order, phylum, "other"),
    phylum_color_noStat = ifelse(phylum %in% desired_order, phylum, "other")) 
plot_dt <- plot_dt %>% 
    group_by(phylum_color) %>% 
    mutate(count_SpXphylum = length(unique(species))) %>%
    arrange(desc(count_SpXphylum)) 
plot_dt$phylum_color <- paste0(plot_dt$phylum_color, " (s=", plot_dt$count_SpXphylum, ")")
plot_dt$phylum_color <- as.factor(plot_dt$phylum_color)
levels(plot_dt$phylum_color) <- unique(plot_dt$phylum_color)

# Plot
set1_palette <- c("#984EA3" ,"#999999") # OMD

line_plot <- ggplot(plot_dt, aes(x = comp_lv, y = value, group = species, color = phylum_color)) +
  geom_line(alpha = 0.8, show.legend = FALSE, linewidth = .3) + 
  labs(x = "MRF threshold (%)", y = "F1 score") +
  theme_minimal()+ 
  scale_color_manual(values = set1_palette) + 
  guides(color = guide_legend(title = "Phylum")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
  axis.text.y = element_text(size = 8),
  legend.text = element_text(size = 6)
  )

# determine the MRF threshold maximizing the F1-score as average of all the species
n_mod <- length(unique(all_dist_df$ref_mod_id))
find_best_perf <- all_dist_df %>%
  group_by(comp_lv) %>%
  summarise(avgPerformance = sum(value)/n_mod,
    stdPerformance = sd(value)/n_mod) %>%
  mutate(species = NA,
    genera = NA)
max_index <- which.max(find_best_perf$avgPerformance)
find_best_perf$avgPerformance[max_index]
find_best_perf$stdPerformance[max_index]
max_index <- find_best_perf$comp_lv[max_index]
# determine the MRF threshold maximizing the F1-score of each species
find_best_perf_Xsp <- all_dist_df %>%
  group_by(ref_mod_id) %>%
  mutate(bestXsp = max(value),
  whichBest = value==bestXsp) %>%
  mutate(species = NA,
    genera = NA)
find_best_perf_Xsp <- find_best_perf_Xsp[find_best_perf_Xsp$whichBest,] %>%
  arrange(comp_lv) %>%
  distinct(ref_mod_id, .keep_all = TRUE)
# Add vertical lines for each MRF thresholds
comp_lv_values <- unique(find_best_perf_Xsp$comp_lv)
for (comp_lv_value in comp_lv_values) {
  line_plot <- line_plot +
    geom_vline(xintercept = comp_lv_value, linetype = "dashed", color = "grey", linewidth = .3)
}
line_plot <- line_plot +
  geom_line(data=find_best_perf, aes(x = comp_lv, y = avgPerformance), color = "#000000", linewidth = .6) + # red: "#ad1d00"
  geom_vline(xintercept = max_index, linetype = "dashed", color = "black", linewidth = .4) 

# save plot
lineplot_omd <- line_plot
# cat(file.path(output.dir, paste0(which.stat, "_OMD_wlegend.pdf")))
# ggsave(file.path(output.dir, paste0(which.stat, "_OMD_wlegend.pdf")), line_plot, width = 8, height = 8, units = "cm", dpi = 300)

### fig. 3 -----------------------------------------------------------------------------------
# Combine the two plots
combined_plot <- ggarrange( lineplot_UHGG + labs(title = "a."), lineplot_omd + labs(title = "b."),
    ncol = 2, nrow = 1) # , widths = c(2, 1.2))
cat(file.path(output.dir, paste0(which.stat, "_UHGG_OMD.pdf")))
ggsave(file.path(output.dir, paste0(which.stat, "_UHGG_OMD.pdf")), combined_plot,  units = "cm", width = 13, height = 8, dpi = 300)