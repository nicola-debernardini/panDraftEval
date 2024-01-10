
############################################################################
### Improvement of pan.mod based on completeness subset compared to best MAG
# Parameters
min.rxn.freq.in.mods <- 0.08
num_MAGxiter <- 30
n.iter <- 10

# Helpers
suppressMessages(library(sybil))
suppressMessages(library(cplexAPI))
source("/home/bioinfo/users/niber/prj_panModel/scr/pan-draft_functions.R")
source("/home/bioinfo/users/niber/prj_panModel/scr/function_collection.R")

cat("\n\n### Estimating model improvement based on MAGs completeness\n")
cat(paste("Statistics calculated using", num_MAGxiter, "MAGs in", n.iter, "iteration and the", dist_type, "distance\n")) 

iso.binary.rxn.table      <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft_18Oct/gapfill/ISO/allISO_rxnXmod.txt" # "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft/allISO_rxnXmod.txt"
inIso_dt_list <- load_files_from_paths(iso.binary.rxn.table, ".tsv")

# Computer statictics for subcluster based on MAG completeness
all_dist_df <- data.frame(comp_lv = numeric(0), mod_id = character(), ref_mod_id = character(), dist_tp = character(), mod_category = character(), value = numeric(0), iter = numeric(0))
# list <- c("MGYG000002504")
for (mod in names(inIso_dt_list)) { #names(inIso_dt_list)
    cat(paste("\n", mod))   
    inIso_dt <- inIso_dt_list[[mod]]
    inMAG_dt <- inMAG_dt_list[[mod]]

    ### Process models to be able to reconstruct the pan-Draft 
    mod.path <- paste0("/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/input/MAG/input4modelMerging_", mod,".txt")
    rxn.weights.path <- paste0("/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/input/MAG/input4WeightsMerging_", mod,".txt")
    rxnXgene.table.path <- paste0("/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/input/MAG/input4XgenesMerging_", mod,".txt")
    pathways.table.path <- paste0("/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/input/MAG/Pathways_", mod,".txt")
    cat("\nLoading input data, be patient... \n")
    model_list <- load_files_from_paths_for_RDS(mod.path, "-draft.RDS")    
    weights_list <- load_files_from_paths_for_RDS(rxn.weights.path, "-rxnWeights.RDS")
    xgenes_list <- load_files_from_paths_for_RDS(rxnXgene.table.path, "-rxnXgenes.RDS")
    pathways_list <- load_files_from_paths_for_tbl(pathways.table.path, "-all-Pathways.tbl.gz")
    cat("\tcompleted\n") 

    # Build the data.table of presence/absence reaction in a list of models
    res <- build_rxn2mod_dt(model_list) 
    rxn2mod_dt <- res[[1]]
    mod_id2mod_dict <- res[[2]]
    first_mod_desc <- model_list[[1]]@mod_desc
    rm(list = c("model_list", "res")) # Remove some variables
    # Extract RXN from model, extract for all rxn the info from the first model having that rxn # SLOW PART
    info_all_rxns_mods <- list()
    for (rxn_id in rxn2mod_dt$rxn) {
        first_modID_with_rxn <- colnames(rxn2mod_dt)[rxn2mod_dt[rxn_id, ]==1][1] # first model ID having the reaction
        first_mod_with_rxn <- attr(mod_id2mod_dict, first_modID_with_rxn) # first model having the reaction
        rxn <- getReaction(first_mod_with_rxn, j = rxn_id) # extract the reaction from a model based on rxn_ID
        info_all_rxns_mods <- c(info_all_rxns_mods, rxn) # save rxn in a list
    }
    # find associations between MET_ID and MET_NAME/RXN_ID --- dictionary
    met_id2met_name_dict <- list()
    met_id2rxn_id_dict <- list()
    for (rxn in info_all_rxns_mods) {
        rxn_id = rxn@react_id
        for (idx in seq_along(rxn@met_id)) { # identify the metabolites associated to each rxn
            met_id = rxn@met_id[idx] 
            met_name = rxn@met_name[idx] 

            # MET: metabolite ID to metabolite name dictionary
            if (!(met_id %in% names(met_id2met_name_dict))) {
                met_id2met_name_dict[[met_id]] <- list(met_name)
            } else if (!(met_name %in% met_id2met_name_dict[[met_id]])) { # save all the met_name associated to a specific met_id 
                met_id2met_name_dict[[met_id]] <- c(met_id2met_name_dict[[met_id]], met_name)
            }
            # RXN: metabolite ID to rxn ID dictionary
            if (!(met_id %in% names(met_id2rxn_id_dict))) {
                met_id2rxn_id_dict[[met_id]] <- list(rxn_id)
            } else if (!(rxn_id %in% met_id2rxn_id_dict[[met_id]])) {
                met_id2rxn_id_dict[[met_id]] <- c(met_id2rxn_id_dict[[met_id]], rxn_id) # save all the rxn_id associated to a specific met_id 
            }
        }
    }
    met_id2rxn_id_df <- pad_dict_to_dataframe(met_id2rxn_id_dict)
    met_id2met_name_padded_df <- pad_dict_to_dataframe(met_id2met_name_dict)
    rm(list = c("met_id2rxn_id_dict", "met_id2met_name_dict")) # Remove some variables
    cat(paste("Let's standanrdize the name of the duplicated compounds:\n"))
    dupl.th <- dim(met_id2met_name_padded_df)[2]-2 # ids have duplicated names if at least 2 different names are associated with the same id
    met_id2duplicated_met_name_df <- met_id2met_name_padded_df[!(rowSums(is.na(met_id2met_name_padded_df)) > dupl.th),] # find metabolites that have a duplicated name 
    info_all_rxns_mods <- standardize_duplicated_met_name(met_id2duplicated_met_name_df, met_id2rxn_id_df, info_all_rxns_mods) # Standanrdize the duplicated MET_NAME 

    # Calculate frequency of rxn in all the reference genomes 
    Isolate_cols <- names(inIso_dt)[names(inIso_dt)!="rxn"]
    freq_ref_rxn <- inIso_dt %>%
    select(all_of(Isolate_cols)) %>%
    mutate(b_i = rowSums(.)/length(Isolate_cols),
        rxn = inIso_dt$rxn) %>%
    select(b_i, rxn)

    reference_reactions <- inIso_dt$rxn
    MAG_cols <- colnames(inMAG_dt)[-1]

    for (th in seq(50, 90, 10)){
        th_Species_name <- metadata %>% # Subset MAGs based on completeness
            filter(Species_name %in% MAG_cols) %>%
            filter(Completeness >= th & Completeness < (th+10)) %>%
            select(Species_name, Completeness)
        cat(paste("num MAG in subset:", dim(th_Species_name)[1], " - Comp.th.: ", th,"\n"))
        # check if the num of MAG in the subset is smaller that the number of MAGs desired to test the improvement of the pan-Draft
        if (dim(th_Species_name)[1] >= 30) { # if the number of MAGs per cluster is low recalibrate the number of num_MAGxiter 
            # Tsgb of pan.mod
            for (i in seq_len(n.iter)){
                shuffled_indices <- sample(rownames(th_Species_name))[1:num_MAGxiter] # shuffle the row indices
                subset_th_Species_name <- th_Species_name[shuffled_indices, ]
                # Define the pan.model based on a single th
                pan.model_rxn <- inMAG_dt %>% 
                    select(subset_th_Species_name$Species_name) %>%
                    mutate(freq_rxn = rowSums(.) / num_MAGxiter, 
                        pan.mod = ifelse(freq_rxn > min.rxn.freq.in.mods, 1, 0),
                        rxn = inMAG_dt$rxn) %>%
                    select(pan.mod, rxn)
                # Compute statictics
                Tsgb <- Tsgb_cal(pan.model_rxn, freq_ref_rxn, "pan.mod")
                dist_val  <- distance_cal(pan.model_rxn, freq_ref_rxn, dist_type)   
                dist_val <- dist_val[1]
                predicted_reactions <- pan.model_rxn[["rxn"]][as.logical(pan.model_rxn[["pan.mod"]])]
                metrics_table <- compute_contigency_table(reference_reactions, predicted_reactions) # compute contingency table for a pan.mod
                stat_val <- metrics_table[which.stat][[1]]
                
                # Populate the dataframe
                new_row_Tsgb <- list(comp_lv = th,
                                mod_id = paste("pan.mod", min.rxn.freq.in.mods), 
                                ref_mod_id = mod,
                                dist_tp = "Tsgb", 
                                mod_category = "pan.mod", 
                                value = Tsgb,
                                iter = i)
                new_row_dist <- list(comp_lv = th,
                                mod_id = paste("pan.mod", min.rxn.freq.in.mods), 
                                ref_mod_id = mod,
                                dist_tp = dist_type, 
                                mod_category = "pan.mod", 
                                value = dist_val,
                                iter = i)
                new_row_stat <- list(comp_lv = th,
                                mod_id = paste("pan.mod", min.rxn.freq.in.mods), 
                                ref_mod_id = mod,
                                dist_tp = which.stat, 
                                mod_category = "pan.mod", 
                                value = stat_val,
                                iter = i)
                all_dist_df <- rbind(all_dist_df, new_row_Tsgb, new_row_dist, new_row_stat)

                # Generate the pan-Draft
                # subset rxn2mod_dt based on min.rxn.freq.in.mods to reconstruct the pan-draft
                cat(paste("\niter:", i, "- Compth:", th,"- rxnth:", min.rxn.freq.in.mods))
                rxn_PresAbs_dt <- as.data.table(copy(pan.model_rxn))
                row_sub <- rxn_PresAbs_dt[,  "pan.mod"]==1   
                subSet_rxn_dt <- rxn_PresAbs_dt[row_sub[, 1], .SD, .SDcols = "rxn"] # select id of present reactions
                pan.mod <- build_panDraft(subSet_rxn_dt, info_all_rxns_mods, first_mod_desc) # build pan-Draft
                saveRDS(pan.mod, file.path(output.dir,paste0(mod, "_iter_", i, "_Compth_", th,"_panModel-draft.RDS")))

                ### Generate specific rxnWeights table scores 
                sub_weights_list <- weights_list[subset_th_Species_name$Species_name]
                weights_dt <- rbindlist(sub_weights_list, idcol = "model_id")
                weights_dt[, num.pan := .N, by = .(seed)] # Add lines for seed in order to obtain corrected median  
                # Calculate custom median of "weight" by grouping "seed"
                weights_dt[, weight.pan := custom_quartile_weight(weight, num.pan, num_MAGxiter, min.rxn.freq.in.mods), by = .(seed)]
                weights_dt[, num.pan := NULL] # drop the colum num.pan
                weights_dt <- weights_dt[order(seed, weight)] # alternative: "abs(weight - weigth.pan)" 
                weights_dt <- weights_dt[!duplicated(seed)] # keep the entry with the highest score.
                weights_dt[,c("weight", "weight.pan")] <- weights_dt[,c("weight.pan", "weight")] # Swap the values of "weight" and "weight.pan" columns
                colnames(weights_dt)[colnames(weights_dt)=="weight.pan"] <- "weight.old" # Rename column
                setkeyv(weights_dt, c("model_id", "seed")) # set data.table keys
                saveRDS(weights_dt, file.path(output.dir,paste0(mod, "_iter_", i, "_Compth_", th,"_panModel-rxnWeigths.RDS")))
                cat(paste("Save output in: ", file.path(output.dir,paste0(mod, "_iter_", i, "_Compth_", th,"_panModel-rxnWeigths.RDS\n"))))

                ### Update rxnXGenes table  
                sub_xgenes_list <- xgenes_list[subset_th_Species_name$Species_name]
                xgenes_dt <- rbindlist(sub_xgenes_list, idcol = "model_id")
                setkeyv(xgenes_dt, c("model_id", "seed"))
                xgenes_dt <- xgenes_dt[weights_dt] # extract only genes corresponding to the ref reactions
                saveRDS(xgenes_dt, file.path(output.dir,paste0(mod, "_iter_", i, "_Compth_", th,"_panModel-rxnXgenes.RDS")))
                cat(paste("Save output in: ", file.path(output.dir,paste0(mod, "_iter_", i, "_Compth_", th,"_panModel-rxnXgenes.RDS\n"))))

                ### Save the list of detected Pathways
                sub_pathways_list <- pathways_list[subset_th_Species_name$Species_name]
                pathways_dt <- do.call(rbind, pathways_list)
                pathways_dt <- pathways_dt[, .(Prediction = sum(Prediction)), by=c("ID", "Name")]
                pathways_dt$Prediction <- pathways_dt$Prediction > num_MAGxiter*min.rxn.freq.in.mods
                fwrite(pathways_dt, file = file.path(output.dir,paste0(mod, "_iter_", i, "_Compth_", th,"_panModel-tmp-Pathways.tbl")), sep = "\t", quote = FALSE)
                cat(paste("Save output in: ", file.path(output.dir,paste0(mod, "_iter_", i, "_Compth_", th,"_panModel-tmp-Pathways.tbl\n"))))
            }
        } 

        # Tsgb of the most complete mag
        # identify the best mag in the set
        best_mag_id <- th_Species_name %>%
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
}

# fwrite(all_dist_df, file = file.path("/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft/tbl/AccTsgbJacc_pan.mod_vs_bestMAG.tsv"), sep = "\t")

# PLOT
what <- which.stat # which.stat dist_type "Tsgb"
best_mag_df <- all_dist_df %>% 
    filter(dist_tp == what) %>%
    filter(mod_category == "best_mag")
head(best_mag_df)
dim(best_mag_df)

pan.mod_df <- all_dist_df %>% 
    filter(dist_tp == what) %>%
    filter(mod_category == "pan.mod") %>%
    group_by(ref_mod_id, comp_lv) %>%
    summarise(avg_value = mean(value),
        std_value = sd(value),
        .groups = "drop")
head(pan.mod_df)
dim(pan.mod_df)

merged_df <- merge(best_mag_df, pan.mod_df, by = c("ref_mod_id", "comp_lv"), all.x=TRUE)
merged_df <- merged_df %>%
    mutate(diff_pan.modVsBest.mag = avg_value - value)
# add taxonomy
merged_df <- merge(merged_df, metadata, by.x="ref_mod_id", by.y="Species_name")
merged_df <- parse_GTDBtaxonomy(merged_df, "Lineage") # GTDB.Taxonomy Lineage
merged_df <- refine_GTDBtaxonomy(merged_df)
head(merged_df)
dim(merged_df)

cols <- RColorBrewer::brewer.pal(6,'OrRd')[c(1,6)]
heatmap_plot <- ggplot(merged_df, aes(x = species, y = factor(comp_lv), fill = diff_pan.modVsBest.mag)) +
  geom_tile() +
  labs(x = "Reference Model ID", y = "Competeness level", fill = paste("Difference", what, "value")) +
  theme_minimal() +
  scale_fill_gradient(low=cols[1], high=cols[2], na.value = "gray") +
  theme(axis.text.x = element_text(angle = 80, hjust = 1)) # Rotate x-axis labels for readability
print(heatmap_plot)
        
out_fn <- paste0(what,"_improvement_based_on_MAG_CompLV.pdf")
ggsave(file.path(output.dir, out_fn), heatmap_plot, width = 10, height = 4, units = "in", dpi = 300)