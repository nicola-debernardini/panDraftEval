######################################################  @
# Script to reproduce the analysis in figure 4:
# Structural comparison of reference models and pan-GEMs generated for MAGs with different completeness levels

# Load R-packages
suppressMessages(library(ggpubr))
suppressMessages(library(data.table)); setDTthreads(1)
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(sybil))

# Little helpers
source("./scr/function_collection.R")

mag.binary.rxn.table      <- "./dat/uhgg/rxnXmod/mag/allMAG_rxnXmod.txt"
iso.binary.rxn.table      <- "./dat/uhgg/rxnXmod/iso/allISO_rxnXmod.txt" 
panMAG.CompLV.gapfill.model.path  <- "./dat/uhgg/GEM/pan/varing_compLv/panFromMAG_Gapfilled_Path.txt" # file path list of all gapfilled pan-GEMs generated using MAG with varing completeness lv. 
completeness.level.table          <- "./dat/uhgg/metadata/genomes-all_metadata.csv"
output.dir                        <- "./dat/fig"

# Parameters:
dist_type <- "binary"
which.stat <- "f1_score" # fB_score f1_score recall precision accuracy FDR
n.iter <- 10

# Read metadata
metadata <- read.table(completeness.level.table, sep = ",", header = TRUE, row.names = "Genome") # dec = ","
columns_to_convert_numeric <- c("Completeness")
metadata[columns_to_convert_numeric] <- lapply(metadata[columns_to_convert_numeric], as.numeric)
metadata$Species_name <- rownames(metadata)

# create binary tables for RXNs of gapfilled pan-GEMs
panMAG.CompLV.gapfill.mods <- load_files_from_paths_for_RDS(panMAG.CompLV.gapfill.model.path, ".RDS")
for (mod_idx in names(panMAG.CompLV.gapfill.mods)) {
  mod_to_rename <- panMAG.CompLV.gapfill.mods[[mod_idx]]
  n_name <- paste0(mod_idx, "_", mod_to_rename@mod_id)
  cat(paste("\nModel ID:", mod_to_rename@mod_id, "is duplicated\n"))
  cat(paste("renames as:", n_name))
  mod_to_rename@mod_id <- n_name
  panMAG.CompLV.gapfill.mods[[mod_idx]] <- mod_to_rename
  # Add annotation column to model attributes if not already there
  mod_to_rename <- add_annotation_column_to_attributes(mod_to_rename)
}
res <- build_rxn2mod_dt(panMAG.CompLV.gapfill.mods)
panMAG.CompLV.gapfill.rxn2mod_dt <- res[[1]]
rm(list = c("panMAG.CompLV.gapfill.mods", "res")) # Remove some variables

# load data for MAG and ISO since I have alreacy calculated
mag.binary_table.list <- load_files_from_paths(mag.binary.rxn.table, ".tsv")
iso.binary_table.list <- load_files_from_paths(iso.binary.rxn.table, ".tsv")

# computer statictics 
all_dist_df <- data.frame(comp_lv = numeric(0), mod_id = character(), ref_mod_id = character(), dist_tp = character(), mod_category = character(), value = numeric(0), iter = numeric(0))
for (mod in names(iso.binary_table.list)) { 
  cat(paste("\n", mod, "\n"))   
  inIso_dt <- iso.binary_table.list[[mod]]
  inMAG_dt <- mag.binary_table.list[[mod]]

  reference_reactions <- inIso_dt$rxn # reactions in reference genomes
  MAG_cols <- colnames(inMAG_dt)[-1]

  # frequency of rxn in all the reference genomes
  Isolate_cols <- names(inIso_dt)[names(inIso_dt)!="rxn"]
  freq_ref_rxn <- inIso_dt %>%
    select(all_of(Isolate_cols)) %>%
    mutate(b_i = rowSums(.)/length(Isolate_cols),
        rxn = inIso_dt$rxn) %>%
    select(b_i, rxn)
  
  for (th in seq(50, 90, 10)){
    th_Species_name <- metadata %>% # subset MAGs based on completeness
        filter(Species_name %in% MAG_cols) %>%
        filter(Completeness >= th & Completeness < (th+10)) %>%
        select(Species_name, Completeness)
    cat(paste("num MAG in subset:", dim(th_Species_name)[1], " - Comp.th.: ", th,"\n"))
    
    # pan-GEMs ---------------------------------------------------------
    # check if the file exist, which depend if originally the num of MAGs in the subset was at least 30
    if (paste0(mod,"_iter_1_Compth_", th, "_panModel_panDraf_model") %in% names(panMAG.CompLV.gapfill.rxn2mod_dt)) {  
      for (i in seq_len(n.iter)){
        panModel_name <- paste0(mod,"_iter_", i, "_Compth_", th, "_panModel_panDraf_model")
        col_to_extract <- c(panModel_name, "rxn")
        pan.model_rxn <- panMAG.CompLV.gapfill.rxn2mod_dt[, ..col_to_extract]
        
        # Compute statictics
        Tsgb <- Tsgb_cal(pan.model_rxn, freq_ref_rxn, panModel_name)
        dist_val  <- distance_cal(pan.model_rxn, freq_ref_rxn, dist_type)   
        dist_val <- dist_val[1]
        predicted_reactions <- pan.model_rxn[["rxn"]][as.logical(pan.model_rxn[[panModel_name]])]
        metrics_table <- compute_contigency_table(reference_reactions, predicted_reactions) # compute contingency table for a pan.mod
        stat_val <- metrics_table[which.stat][[1]]
        
        # Populate the dataframe
        new_row_Tsgb <- list(comp_lv = th,
                        mod_id = panModel_name, 
                        ref_mod_id = mod,
                        dist_tp = "Tsgb", 
                        mod_category = "pan.mod", 
                        value = Tsgb,
                        iter = i)
        new_row_dist <- list(comp_lv = th,
                        mod_id = panModel_name, 
                        ref_mod_id = mod,
                        dist_tp = dist_type, 
                        mod_category = "pan.mod", 
                        value = dist_val,
                        iter = i)
        new_row_stat <- list(comp_lv = th,
                        mod_id = panModel_name, 
                        ref_mod_id = mod,
                        dist_tp = which.stat, 
                        mod_category = "pan.mod", 
                        value = stat_val,
                        iter = i)
        all_dist_df <- rbind(all_dist_df, new_row_Tsgb, new_row_dist, new_row_stat)
      }
    }

    # best MAG GEMs ---------------------------------------------------------
    for (mag in th_Species_name$Species_name){
      mag_data <- inMAG_dt %>%
        select(all_of(mag), rxn)

      # Compute statictics
      Tsgb <- Tsgb_cal(mag_data, freq_ref_rxn, mag)
      dist_val  <- distance_cal(mag_data, freq_ref_rxn, dist_type)   
      dist_val <- dist_val[1]
      predicted_reactions <- mag_data[["rxn"]][as.logical(mag_data[[mag]])]
      metrics_table <- compute_contigency_table(reference_reactions, predicted_reactions) # compute contingency table for a pan.mod
      stat_val <- metrics_table[which.stat][[1]]
      
      # Populate the dataframe
      new_row_Tsgb <- list(comp_lv = th,
                      mod_id = mag, 
                      ref_mod_id = mod,
                      dist_tp = "Tsgb", 
                      mod_category = "mag", 
                      value = Tsgb,
                      iter = NA)
      new_row_dist <- list(comp_lv = th,
                      mod_id = mag, 
                      ref_mod_id = mod,
                      dist_tp = dist_type, 
                      mod_category = "mag", 
                      value = dist_val,
                      iter = NA)
      new_row_stat <- list(comp_lv = th,
                      mod_id = mag, 
                      ref_mod_id = mod,
                      dist_tp = which.stat, 
                      mod_category = "mag", 
                      value = stat_val,
                      iter = NA)
      all_dist_df <- rbind(all_dist_df, new_row_Tsgb, new_row_dist, new_row_stat)
    }

    best_mag_id <- th_Species_name %>% # identify the best mag in the set
        filter(Completeness == max(Completeness)) %>%
        select(Species_name) %>%
        slice(1) %>%
        as.character(.)
    mag_data <- inMAG_dt %>%
        select(all_of(best_mag_id), rxn)
  
    Tsgb_bestmag <- Tsgb_cal(mag_data, freq_ref_rxn, best_mag_id)
    dist_bestmag  <- distance_cal(mag_data, freq_ref_rxn, dist_type)
    dist_bestmag  <- dist_bestmag[1]
    predicted_reactions <- mag_data[["rxn"]][as.logical(mag_data[[best_mag_id]])]
    metrics_table <- compute_contigency_table(reference_reactions, predicted_reactions) # compute contingency table for a pan.mod
    stat_bestmag <- metrics_table[which.stat][[1]]
          
    # populate the dataframe
    new_row_Tsgb <- list(comp_lv = th,
                    mod_id = best_mag_id, 
                    ref_mod_id = mod,
                    dist_tp = "Tsgb", 
                    mod_category = "best_mag", 
                    value = Tsgb_bestmag,
                    iter = NA)
    new_row_dist <- list(comp_lv = th,
                    mod_id = best_mag_id, 
                    ref_mod_id = mod,
                    dist_tp = dist_type, 
                    mod_category = "best_mag", 
                    value = dist_bestmag,
                    iter = NA)
    new_row_stat <- list(comp_lv = th,
                    mod_id = best_mag_id, 
                    ref_mod_id = mod,
                    dist_tp = which.stat, 
                    mod_category = "best_mag", 
                    value = stat_bestmag,
                    iter = NA)
    all_dist_df <- rbind(all_dist_df, new_row_Tsgb, new_row_dist, new_row_stat)
  }

  # reference GEMs --------------------------------------------------------- 
  for (iso in Isolate_cols){
    iso_data <- inIso_dt %>%
      select(all_of(iso), rxn)

    # Compute statictics
    Tsgb <- Tsgb_cal(iso_data, freq_ref_rxn, iso)
    dist_val  <- distance_cal(iso_data, freq_ref_rxn, dist_type)   
    dist_val <- dist_val[1]
    predicted_reactions <- iso_data[["rxn"]][as.logical(iso_data[[iso]])]
    metrics_table <- compute_contigency_table(reference_reactions, predicted_reactions) # compute contingency table for a pan.mod
    stat_val <- metrics_table[which.stat][[1]]
    
    # populate the dataframe
    new_row_Tsgb <- list(comp_lv = 100,
                    mod_id = iso, 
                    ref_mod_id = mod,
                    dist_tp = "Tsgb", 
                    mod_category = "iso", 
                    value = Tsgb,
                    iter = NA)
    new_row_dist <- list(comp_lv = 100,
                    mod_id = iso, 
                    ref_mod_id = mod,
                    dist_tp = dist_type, 
                    mod_category = "iso", 
                    value = dist_val,
                    iter = NA)
    new_row_stat <- list(comp_lv = 100,
                    mod_id = iso, 
                    ref_mod_id = mod,
                    dist_tp = which.stat, 
                    mod_category = "iso", 
                    value = stat_val,
                    iter = NA)
    all_dist_df <- rbind(all_dist_df, new_row_Tsgb, new_row_dist, new_row_stat)
  }
}

# Plots ------------------------------------------------------------------------------------------
what <- which.stat # which.stat dist_type "Tsgb"
best_mag_df <- all_dist_df %>% 
    filter(dist_tp == what) %>%
    filter(mod_category == "best_mag")
# summarize the results for the pan-GEMs generated using MAGs with different completeness lv.
pan.mod_df <- all_dist_df %>% 
    filter(dist_tp == what) %>%
    filter(mod_category == "pan.mod") %>%
    group_by(ref_mod_id, comp_lv) %>%
    summarise(avg_value = mean(value),
        std_value = sd(value),
        .groups = "drop")
# put back together the two dataset
merged_df <- merge(best_mag_df, pan.mod_df, by = c("ref_mod_id", "comp_lv"), all.x=TRUE)
merged_df <- merged_df %>%
    mutate(diff_pan.modVsBest.mag = avg_value - value)
merged_df <- merge(merged_df, metadata, by.x="ref_mod_id", by.y="Species_name")
merged_df <- parse_GTDBtaxonomy(merged_df, "Lineage")
merged_df <- refine_GTDBtaxonomy(merged_df)
# remove the species for which any pan-GEMs was reconstructed
merged_filtered_df <- merged_df %>%
  group_by(Species_rep) %>%
  filter(any(!is.na(diff_pan.modVsBest.mag)))

# Heatmap
cols <- RColorBrewer::brewer.pal(6,'OrRd')[c(1,6)]
heatmap_plot <- ggplot(merged_filtered_df, aes(x = species, y = factor(comp_lv), fill = diff_pan.modVsBest.mag)) +
  geom_tile() +
  labs(
    x = "Reference genome taxonomy", 
    y = "Competeness (%)", 
    fill = paste("F1 score difference"), 
    size = 8) +
  scale_fill_gradient(low=cols[1], high=cols[2], na.value = "gray") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 7),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 6)) 
# ggsave(file.path(output.dir, paste0("heatmap.pdf")), heatmap_plot,  units = "cm", width = 16, height = 10, dpi = 300)

# Boxplot
noBest_mag_df <- all_dist_df %>% 
    filter(dist_tp == what) %>%
    filter(mod_category != "best_mag") %>%
    group_by(comp_lv)

noBest_mag_df$comp_lv <- as.factor(noBest_mag_df$comp_lv)
box_plot <- ggplot(noBest_mag_df, aes(x=comp_lv, y=value, fill=mod_category)) +
  geom_boxplot() +
  labs(
    x = "Completeness (%)",
    y = "F1 score",
    size = 8) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 8),
  legend.text = element_text(size = 6),
  legend.title = element_text(size = 8))
# ggsave(file.path(output.dir, paste0("boxplot.pdf")), box_plot,  units = "cm", width = 16, height = 10, dpi = 300)

# Combine the two plots
combined_plot <- ggarrange(box_plot + labs(title = "a.") + theme(legend.position = "none"),
  heatmap_plot + labs(title = "b.") + theme(legend.position = "none"), 
  ncol = 2,  widths = c(1.2, 2)) 
cat(file.path(output.dir, paste0("boxPlot_all_MAG_CompLV_AND_", what, "_improvement_based_on_MAG_CompLV_gapfill_NOLEGEND_v2.pdf")))
ggsave(file.path(output.dir, paste0("boxPlot_all_MAG_CompLV_AND_", what, "_improvement_based_on_MAG_CompLV_gapfill_NOLEGEND_v2.pdf")), 
  combined_plot,  units = "cm", width = 13, height = 10, dpi = 300)
