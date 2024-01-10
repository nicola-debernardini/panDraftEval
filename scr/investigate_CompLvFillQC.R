# Load R-packages
suppressMessages(library(data.table)); setDTthreads(1)
suppressMessages(library(ggplot2)) # Graphical
suppressMessages(library(dplyr)) # Graphical
suppressMessages(library(sybil))
library(ggpubr)

# Little helpers
source("/home/bioinfo/users/niber/prj_panModel/scr/function_collection.R")
source("/home/bioinfo/users/niber/prj_panModel/scr/pan-draft_functions.R")

# Arguments:
iso.gapfill.model.pathList        <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft_18Oct/gapfill/ISO_Gapfilled_listsPath.txt"
mag.gapfill.model.pathList        <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft_18Oct/gapfill/MAG_Gapfilled_listsPath.txt"

panMAG.CompLV.gapfill.model.path  <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft_18Oct/CompLvFillQC/panFromMAG_Gapfilled_Path.txt"
completeness.level.table          <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/metadata/genomes-all_metadata.csv"
output.dir                        <- "/home/bioinfo/users/niber/prj_panModel/db/Paper/fig"

# Parameters:
dist_type <- "binary"
which.stat <- "f1_score" # fB_score f1_score recall precision accuracy FDR
n.iter <- 10

# Read metadata
metadata <- read.table(completeness.level.table, sep = ",", header = TRUE, row.names = "Genome") # dec = ","
columns_to_convert_numeric <- c("Completeness")
metadata[columns_to_convert_numeric] <- lapply(metadata[columns_to_convert_numeric], as.numeric)
metadata$Species_name <- rownames(metadata)

# Create binary tables for RXNs of gapfilled ISOLATE models and save results in a LIST
iso.binary_table.list <- list()
singleGeno_pathList <- readLines(iso.gapfill.model.pathList)
for (list in singleGeno_pathList) {
  dirty_spID <- tail(strsplit(list, "_")[[1]], n = 1)
  spID <- strsplit(dirty_spID, "\\.")[[1]][1]
  
  singleGeno_mods <- load_files_from_paths_for_RDS(list, ".RDS") # open all the models of the isolats of a single species
  res <- build_rxn2mod_dt(singleGeno_mods) # Build the data.table of 1|0 reaction in a list of models
  rxn2mod_dt <- res[[1]]
  iso.binary_table.list[[spID]] <- rxn2mod_dt
  rm(list = c("singleGeno_mods", "res")) # Remove some variables
}
length(iso.binary_table.list)
# SAVE
for (rxn2mod_dt_name in names(iso.binary_table.list)) {
  rxn2mod_dt <- iso.binary_table.list[[rxn2mod_dt_name]]
  fwrite(rxn2mod_dt, file = file.path(output.dir, paste0(rxn2mod_dt_name, "_ISO_rxnXmod.tsv")), sep = "\t", quote = FALSE)
}

# Create binary tables for RXNs of gapfilled MAG models and save results in a LIST
mag.binary_table.list <- list()
singleGeno_pathList <- readLines(mag.gapfill.model.pathList)
for (list in singleGeno_pathList) {
  dirty_spID <- tail(strsplit(list, "_")[[1]], n = 1)
  spID <- strsplit(dirty_spID, "\\.")[[1]][1]
  
  singleGeno_mods <- load_files_from_paths_for_RDS(list, ".RDS") # open all the models of the isolats of a single species
  res <- build_rxn2mod_dt(singleGeno_mods) # Build the data.table of 1|0 reaction in a list of models
  rxn2mod_dt <- res[[1]]
  mag.binary_table.list[[spID]] <- rxn2mod_dt
  rm(list = c("singleGeno_mods", "res")) # Remove some variables
}
length(mag.binary_table.list)
# SAVE
for (rxn2mod_dt_name in names(mag.binary_table.list)) {
  rxn2mod_dt <- mag.binary_table.list[[rxn2mod_dt_name]]
  fwrite(rxn2mod_dt, file = file.path(output.dir, paste0(rxn2mod_dt_name, "_MAG_rxnXmod.tsv")), sep = "\t", quote = FALSE)
}

# Create binary tables for RXNs of gapfilled panmodel from MAG and save results in a TABLE
panMAG.CompLV.gapfill.mods <- load_files_from_paths_for_RDS(panMAG.CompLV.gapfill.model.path, ".RDS") # list containing all the gapfilled models generated based on MAG completeness 
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
res <- build_rxn2mod_dt(panMAG.CompLV.gapfill.mods) # Build the data.table of 1|0 reaction in a list of models
panMAG.CompLV.gapfill.rxn2mod_dt <- res[[1]]
rm(list = c("panMAG.CompLV.gapfill.mods", "res")) # Remove some variables

# LOAD DATA for MAG and ISO since I have alreacy calculated
mag.binary.rxn.table      <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft_18Oct/gapfill/MAG/allMAG_rxnXmod.txt"
iso.binary.rxn.table      <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft_18Oct/gapfill/ISO/allISO_rxnXmod.txt" 
mag.binary_table.list <- load_files_from_paths(mag.binary.rxn.table, ".tsv")
iso.binary_table.list <- load_files_from_paths(iso.binary.rxn.table, ".tsv")


# Computer statictics for gapfilled panmodel based on MAG completeness
all_dist_df <- data.frame(comp_lv = numeric(0), mod_id = character(), ref_mod_id = character(), dist_tp = character(), mod_category = character(), value = numeric(0), iter = numeric(0))
for (mod in names(iso.binary_table.list)) { 
  cat(paste("\n", mod, "\n"))   
  inIso_dt <- iso.binary_table.list[[mod]] # iso rxn2mod_dt table
  inMAG_dt <- mag.binary_table.list[[mod]] # mag rxn2mod_dt table

  reference_reactions <- inIso_dt$rxn # reactions in reference genomes
  MAG_cols <- colnames(inMAG_dt)[-1] # MAG id in species

  # Calculate frequency of rxn in all the reference genomes 
  Isolate_cols <- names(inIso_dt)[names(inIso_dt)!="rxn"]
  freq_ref_rxn <- inIso_dt %>%
    select(all_of(Isolate_cols)) %>%
    mutate(b_i = rowSums(.)/length(Isolate_cols),
        rxn = inIso_dt$rxn) %>%
    select(b_i, rxn)
  
  for (th in seq(50, 90, 10)){
    th_Species_name <- metadata %>% # Subset MAGs based on completeness
        filter(Species_name %in% MAG_cols) %>%
        filter(Completeness >= th & Completeness < (th+10)) %>%
        select(Species_name, Completeness)
    cat(paste("num MAG in subset:", dim(th_Species_name)[1], " - Comp.th.: ", th,"\n"))
    # PAN MOD
    # check if the file exist, so if originally the num of MAG in the subset was higher than the TH: 30 mag
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

    # ALL MAG 
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

    # BEST MAG Tsgb of the most complete mag
    best_mag_id <- th_Species_name %>% # identify the best mag in the set
        filter(Completeness == max(Completeness)) %>%
        select(Species_name) %>%
        slice(1) %>%
        as.character(.)
    mag_data <- inMAG_dt %>%
        select(all_of(best_mag_id), rxn)
    # Computer statictics
    Tsgb_bestmag <- Tsgb_cal(mag_data, freq_ref_rxn, best_mag_id)
    dist_bestmag  <- distance_cal(mag_data, freq_ref_rxn, dist_type)
    dist_bestmag  <- dist_bestmag[1]
    predicted_reactions <- mag_data[["rxn"]][as.logical(mag_data[[best_mag_id]])]
    metrics_table <- compute_contigency_table(reference_reactions, predicted_reactions) # compute contingency table for a pan.mod
    stat_bestmag <- metrics_table[which.stat][[1]]
          
    # Populate the dataframe
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

  # ALL ISO
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
    
    # Populate the dataframe
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

dim(all_dist_df)
head(all_dist_df)

# fwrite(all_dist_df, file = file.path("/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft_18Oct/CompLvFillQC/fig/F1TsgbJacc_pan.mod_vs_bestMAG_gapfill.tsv"), sep = "\t")

# PLOT
what <- which.stat # which.stat dist_type "Tsgb"
best_mag_df <- all_dist_df %>% 
    filter(dist_tp == what) %>%
    filter(mod_category == "best_mag")
head(best_mag_df)
dim(best_mag_df)
# summarize the results for the different iterations of the pan models
pan.mod_df <- all_dist_df %>% 
    filter(dist_tp == what) %>%
    filter(mod_category == "pan.mod") %>%
    group_by(ref_mod_id, comp_lv) %>%
    summarise(avg_value = mean(value),
        std_value = sd(value),
        .groups = "drop")
head(pan.mod_df)
dim(pan.mod_df)
# put back together the two dataset
merged_df <- merge(best_mag_df, pan.mod_df, by = c("ref_mod_id", "comp_lv"), all.x=TRUE)
merged_df <- merged_df %>%
    mutate(diff_pan.modVsBest.mag = avg_value - value)
# add taxonomy
merged_df <- merge(merged_df, metadata, by.x="ref_mod_id", by.y="Species_name")
merged_df <- parse_GTDBtaxonomy(merged_df, "Lineage") # GTDB.Taxonomy Lineage
merged_df <- refine_GTDBtaxonomy(merged_df)
head(merged_df)
dim(merged_df)
# remove the species for which pangenome was never calculated (always less than 30 sp)
merged_filtered_df <- merged_df %>%
  group_by(Species_rep) %>%
  filter(any(!is.na(diff_pan.modVsBest.mag)))

# HEATMAP
cols <- RColorBrewer::brewer.pal(6,'OrRd')[c(1,6)] # OrRd PuBu
heatmap_plot <- ggplot(merged_filtered_df, aes(x = species, y = factor(comp_lv), fill = diff_pan.modVsBest.mag)) +
  geom_tile() +
  labs(x = "Reference Model ID", y = "Competeness level", 
  fill = paste("F1 score")) + # "Difference in", what
  scale_fill_gradient(low=cols[1], high=cols[2], na.value = "gray") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 80, hjust = 1),
  legend.title = element_text(size = 10)) # Rotate x-axis labels for readability
print(heatmap_plot)
ggsave(file.path(output.dir, paste0("only_heatmap.pdf")), 
  heatmap_plot,  units = "cm", width = 16, height = 10, dpi = 300)

# BOXPLOT
noBest_mag_df <- all_dist_df %>% 
    filter(dist_tp == what) %>%
    filter(mod_category != "best_mag") %>%
    group_by(comp_lv) 
dim(noBest_mag_df)

# noBest_mag_df$comp_lv <- noBest_mag_df$comp_lv
noBest_mag_df$comp_lv <- as.factor(noBest_mag_df$comp_lv)
box_plot <- ggplot(noBest_mag_df, aes(x=comp_lv, y=value, fill=mod_category)) + 
  geom_boxplot() + #position_nudge(x = 1)
  labs(
    x = "Completeness threshold",
    y = "F1 score"
    ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 80, hjust = 1),
  legend.text = element_text(size = 6),
  legend.title = element_text(size = 8)
  ) 
print(box_plot)

combined_plot <- ggarrange(box_plot + labs(title = "a."),
  heatmap_plot + labs(title = "b.") + theme(legend.position = "none"), 
  ncol = 2,  widths = c(1.2, 2)) # Combine the two plots
print(combined_plot)
ggsave(file.path(output.dir, paste0("boxPlot_all_MAG_CompLV_AND_", what, "_improvement_based_on_MAG_CompLV_gapfill_NOLEGEND.pdf")), 
  combined_plot,  units = "cm", width = 16, height = 10, dpi = 300)
