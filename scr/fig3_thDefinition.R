######################################################  @
# Script to reproduce the analysis in figure 3:
# 1. Selection of optimal MRF threshold for UHGG
# 2. Selection of optimal MRF threshold for OMD

# Pseudo-code:
#  load rxnXmod table of reference
#  load rxnXmod table of gapfilled pan-GEM
#  for every GEM in the reference GEM list:
#   for every pan-GEM generated for that specific species with different threasholds:
#       calculate statistics and save them in a dataframe

library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)

# Little helpers
source("./scr/function_collection.R")

# Parameters:
dist_type <- "binary"
which.stat <- "f1_score" # fB_score f1_score recall precision accuracy FDR

# LOAD DATA since alreacy calculated
output.dir <- "/home/bioinfo/users/niber/prj_panModel/db/Paper/fig"
iso.binary.rxn.table      <- "/home/bioinfo/users/niber/prj_panModel/db/microbiomics_ocean/pan.draft/gapfill/ISO/allISO_rxnXmod.txt"
# iso.binary.rxn.table      <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft_18Oct/gapfill/ISO/allISO_rxnXmod.txt" 
iso.binary_table.list     <- load_files_from_paths(iso.binary.rxn.table, ".tsv")

all_dist_df <- data.frame(comp_lv = numeric(0), mod_id = character(), ref_mod_id = character(), dist_tp = character(), mod_category = character(), value = numeric(0))
for (mod in names(iso.binary_table.list)) { 
    cat(paste("\n", mod, "\n"))   
    # Prepare isolates data
    inIso_dt <- iso.binary_table.list[[mod]] # iso rxn2mod_dt table
    reference_reactions <- inIso_dt$rxn # reactions in reference genomes
    # Calculate frequency of rxn in all the reference genomes 
    Isolate_cols <- names(inIso_dt)[names(inIso_dt)!="rxn"]
    freq_ref_rxn <- inIso_dt %>%
        select(all_of(Isolate_cols)) %>%
        mutate(b_i = rowSums(.)/length(Isolate_cols),
            rxn = inIso_dt$rxn) %>%
        select(b_i, rxn)
    
    # Load data for panModel
    rxnXmod_dir_path <- "/home/bioinfo/users/niber/prj_panModel/db/microbiomics_ocean/pan.draft4minThValidation/rxnXmod/fill"
    # rxnXmod_dir_path <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft4minThValidation/rxnXmod"
    sp_rxnXmod_path <- file.path(rxnXmod_dir_path, mod, "rxnXmod.tsv")
    sp_rxnXmod <- read.table(sp_rxnXmod_path, sep="\t", header = TRUE)
    sp_rxnXmod <- as.data.table(sp_rxnXmod)


    for (th in seq(0, 1, 0.01)){
        if (th == 0.01){
            next 
        }
        # check if the file exist, so if originally the num of MAG in the subset was higher than the TH: 30 mag
        panModel_name <- paste0("X",th,"_panModel")
        col_to_extract <- c(panModel_name, "rxn")
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

# Read metadata
completeness.level.table  <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/metadata/genomes-all_metadata.csv"
metadata <- read.table(completeness.level.table, sep = ",", header = TRUE, row.names = "Genome") # dec = ","
columns_to_convert_numeric <- c("Completeness")
metadata[columns_to_convert_numeric] <- lapply(metadata[columns_to_convert_numeric], as.numeric)
metadata$Species_name <- rownames(metadata)
length(unique(metadata$Species_rep))

oceans_metadata_fn        <- "/home/bioinfo/users/niber/prj_panModel/db/microbiomics_ocean/metadata/genomes-summary.csv"  
oceans_metadata <- read.table(oceans_metadata_fn, sep = ",", header = TRUE, dec = ".")
oceans_metadata <- oceans_metadata %>% # separate MAGs from reference genome and SAGs
  mutate(Genome_type = if_else(grepl("_METAG_", Genome), "MAG", if_else(grepl("_REFG_", Genome), "REF", if_else(grepl("_SAGS_", Genome), "SAG", ""))))
oceans_metadata <- oceans_metadata %>% # specify the sample name
  mutate(sample = sapply(strsplit(oceans_metadata[["Genome"]], "_"), function(x) paste(x[1], x[2], x[3], sep="_")))
columns_to_convert_numeric <- c("Completeness")
oceans_metadata[columns_to_convert_numeric] <- lapply(oceans_metadata[columns_to_convert_numeric], as.numeric)
oceans_metadata$Species_name <- oceans_metadata$Genome
metadata <- oceans_metadata
head(metadata$dRep.Species.Cluster)
length(unique(metadata$dRep.Species.Cluster))

# merge the dataset
merged_df <- merge(all_dist_df, metadata, by.x="ref_mod_id", by.y="Species_name")
plot_prob_dist_df <- parse_GTDBtaxonomy(merged_df, "GTDB.Taxonomy") # GTDB.Taxonomy Lineage
plot_prob_dist_df <- refine_GTDBtaxonomy(plot_prob_dist_df)
colnames(plot_prob_dist_df)
head(plot_prob_dist_df)

#UHGG
# column_to_save <- c("species_ID", "Power_law_k_estimate", "Power_law_alpha_estimate",
#   "Heaps_law_k_estimate", "Heaps_law_gamma_estimate",
#   "nIsolates", "nRxn_in_Isolates", "nMAG", "nRxn_in_MAG", "Overlap_IsoVsMAG_size",
#   "strict_core", "core", "shell", "cloud",
#   "domain", "phylum", "class", "order", "family", "genera", "species", "Lineage",
#   "Genome_accession")

# OCEAN
# column_to_save <- c("species_ID", "Power_law_k_estimate", "Power_law_alpha_estimate", 
  # "Heaps_law_k_estimate", "Heaps_law_gamma_estimate", 
  # "nIsolates", "nRxn_in_Isolates", "nMAG", "nRxn_in_MAG", "Overlap_IsoVsMAG_size",
  # "strict_core", "core", "shell", "cloud",
  # "domain", "phylum", "class", "order", "family", "genera", "species", "GTDB.Taxonomy",                     
#   "dRep.Species.Cluster", "dRep.Representative.Genome", "SpecI.Species.Cluster", "mOTUs.Species.Cluster")       
# to_save_df <- plot_prob_dist_df[, ..column_to_save]
# fwrite(to_save_df, file = file.path(output.dir,"all_curated_pan-reactome_stat.tsv"), sep = "\t")

# Probability dist --- PLOT
library(RColorBrewer)
unique(plot_prob_dist_df$phylum)
set1_palette <- brewer.pal(2, "Set1")  # You can specify the number of colors you want

line_plot <- ggplot(plot_prob_dist_df, aes(x = comp_lv, y = value, group = species, color = phylum)) +
  geom_line(alpha = 0.8) +
  labs(x = "MRF threshold (%)", y = "F1 score") +
  theme_minimal()+ 
  scale_color_manual(values = set1_palette) + 
  guides(color = guide_legend(title = "Phylum")) 
#   theme(legend.position = "none")
  # scale_color_discrete()

# determine the maximum value
n_mod <- length(unique(all_dist_df$ref_mod_id))
find_best_perf <- all_dist_df %>%
  group_by(comp_lv) %>%
  summarise(avgPerformance = sum(value)/n_mod,
    stdPerformance = sd(value)/n_mod) %>%
  mutate(species = NA,
    genera = NA)
max_index <- which.max(find_best_perf$avgPerformance)
max_index <- find_best_perf$comp_lv[max_index]

line_plot <- line_plot +
  geom_line(data=find_best_perf, aes(x = comp_lv, y = avgPerformance), color = "#000000", linewidth = .2) + # red: "#ad1d00"
  geom_vline(xintercept = max_index, linetype = "dashed", color = "black", linewidth = .1) 

cat(file.path(output.dir, paste0(which.stat, "_maxLine_prob_dens_ocean.pdf")))
ggsave(file.path(output.dir, paste0(which.stat, "_maxLine_prob_dens_ocean.pdf")), line_plot, width = 16, height = 8, units = "cm", dpi = 300)










########################################################################################################
# VERSION ONLY DRAFT
########################################################################################################
# LOAD DATA since alreacy calculated
output.dir <- "/home/bioinfo/users/niber/prj_panModel/db/Paper/fig"
# gapfill
iso.binary.rxn.table      <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft_18Oct/gapfill/ISO/allISO_rxnXmod.txt" 
iso.binary.rxn.table      <- "/home/bioinfo/users/niber/prj_panModel/db/microbiomics_ocean/pan.draft/gapfill/ISO/allISO_rxnXmod.txt"
# draft
iso.binary.rxn.table      <- "/home/bioinfo/users/niber/prj_panModel/db/microbiomics_ocean/pan.draft/draft/ISO/allISO_rxnXmod.txt"
iso.binary.rxn.table      <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft_18Oct/draft/ISO/allISO_rxnXmod.txt" 
#
iso.binary_table.list     <- load_files_from_paths(iso.binary.rxn.table, ".tsv")

all_dist_df <- data.frame(comp_lv = numeric(0), mod_id = character(), ref_mod_id = character(), dist_tp = character(), mod_category = character(), value = numeric(0))
for (mod in names(iso.binary_table.list)) { 
    cat(paste("\n", mod, "\n"))   
    # Prepare isolates data
    inIso_dt <- iso.binary_table.list[[mod]] # iso rxn2mod_dt table
    reference_reactions <- inIso_dt$rxn # reactions in reference genomes
    # Calculate frequency of rxn in all the reference genomes 
    Isolate_cols <- names(inIso_dt)[names(inIso_dt)!="rxn"]
    freq_ref_rxn <- inIso_dt %>%
        select(all_of(Isolate_cols)) %>%
        mutate(b_i = rowSums(.)/length(Isolate_cols),
            rxn = inIso_dt$rxn) %>%
        select(b_i, rxn)
    
    # Load data for panModel
    # 7 percent
    rxnXmod_dir_path <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft4minThValidation_7perc/rxnXmod"
    rxnXmod_dir_path <- "/home/bioinfo/users/niber/prj_panModel/db/microbiomics_ocean/pan.draft4minThValidation_7perc/rxnXmod"
    # draft
    rxnXmod_dir_path <- "/home/bioinfo/users/niber/prj_panModel/db/microbiomics_ocean/pan.draft4minThValidation/rxnXmod/draft"
    rxnXmod_dir_path <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft4minThValidation/rxnXmod/draft"
    sp_rxnXmod_path <- file.path(rxnXmod_dir_path, mod, "rxnXmod.tsv")
    sp_rxnXmod <- read.table(sp_rxnXmod_path, sep="\t", header = TRUE)
    sp_rxnXmod <- as.data.table(sp_rxnXmod)


    for (th in seq(0, 1, 0.01)){
        # check if the file exist, so if originally the num of MAG in the subset was higher than the TH: 30 mag
        panModel_name <- paste0("X",th,"_panModel")
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

# Read metadata
# for UHGG
completeness.level.table  <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/metadata/genomes-all_metadata.csv"
metadata <- read.table(completeness.level.table, sep = ",", header = TRUE, row.names = "Genome") # dec = ","
columns_to_convert_numeric <- c("Completeness")
metadata[columns_to_convert_numeric] <- lapply(metadata[columns_to_convert_numeric], as.numeric)
metadata$Species_name <- rownames(metadata)
dim(metadata)

# pangenome info 
pangenome_metadata_fn <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/metadata/pangenome_metadata.csv"
pangenome_metadata <- read.table(pangenome_metadata_fn, sep = ";", header = TRUE, dec = ",")
metadata <- merge(metadata, pangenome_metadata, by.x = "Genome_accession", by.y = "Species.representative", all.x = TRUE)
rownames(metadata) <- metadata$Species_name

# for OCEAN
oceans_metadata_fn        <- "/home/bioinfo/users/niber/prj_panModel/db/microbiomics_ocean/metadata/genomes-summary.csv"  
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
plot_prob_dist_df <- parse_GTDBtaxonomy(merged_df, "Lineage") # GTDB.Taxonomy Lineage
plot_prob_dist_df <- refine_GTDBtaxonomy(plot_prob_dist_df)
colnames(plot_prob_dist_df)
head(plot_prob_dist_df)

#UHGG
# column_to_save <- c("species_ID", "Power_law_k_estimate", "Power_law_alpha_estimate",
#   "Heaps_law_k_estimate", "Heaps_law_gamma_estimate",
#   "nIsolates", "nRxn_in_Isolates", "nMAG", "nRxn_in_MAG", "Overlap_IsoVsMAG_size",
#   "strict_core", "core", "shell", "cloud",
#   "domain", "phylum", "class", "order", "family", "genera", "species", "Lineage",
#   "Genome_accession")

# OCEAN
# column_to_save <- c("species_ID", "Power_law_k_estimate", "Power_law_alpha_estimate", 
  # "Heaps_law_k_estimate", "Heaps_law_gamma_estimate", 
  # "nIsolates", "nRxn_in_Isolates", "nMAG", "nRxn_in_MAG", "Overlap_IsoVsMAG_size",
  # "strict_core", "core", "shell", "cloud",
  # "domain", "phylum", "class", "order", "family", "genera", "species", "GTDB.Taxonomy",                     
#   "dRep.Species.Cluster", "dRep.Representative.Genome", "SpecI.Species.Cluster", "mOTUs.Species.Cluster")       
# to_save_df <- plot_prob_dist_df[, ..column_to_save]
# fwrite(to_save_df, file = file.path(output.dir,"all_curated_pan-reactome_stat.tsv"), sep = "\t")

                                
# "Firmicutes_A" "#E41A1C"    "Bacteroidota (m. 11307; s. 75)" "#377EB8"               
# "Firmicutes" "#4DAF4A"  "Proteobacteria (m. 5531; s. 31)" "#984EA3" 
# "Actinobacteriota" "#FF7F00"     "Firmicutes_C (m. 4389; s. 25)" "#FFFF33" 
# "other" "#A65628" "Verrucomicrobiota (m. 515; s. 10)" "#F781BF" 
# "Cyanobacteria" "#999999" 

# Probability dist --- PLOT
# get the desired order of phylum 
desired_order <- plot_prob_dist_df %>%
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

set1_palette <- brewer.pal(9, "Set1")  # You can specify the number of colors you want
set1_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#F781BF") # UHGG
# set1_palette <- c("#984EA3" ,"#999999") # OMD

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

# determine the maximum value of the threshold as average of all the species
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


# determine the maximum value of the threshold of each species
find_best_perf_Xsp <- all_dist_df %>%
  group_by(ref_mod_id) %>%
  mutate(bestXsp = max(value),
  whichBest = value==bestXsp) %>%
  mutate(species = NA,
    genera = NA)
find_best_perf_Xsp <- find_best_perf_Xsp[find_best_perf_Xsp$whichBest,] %>%
  arrange(comp_lv) %>%
  distinct(ref_mod_id, .keep_all = TRUE)

# Add vertical lines for each comp_lv value
comp_lv_values <- unique(find_best_perf_Xsp$comp_lv)
for (comp_lv_value in comp_lv_values) {
  line_plot <- line_plot +
    geom_vline(xintercept = comp_lv_value, linetype = "dashed", color = "grey", linewidth = .3)
}
line_plot <- line_plot +
  geom_line(data=find_best_perf, aes(x = comp_lv, y = avgPerformance), color = "#000000", linewidth = .6) + # red: "#ad1d00"
  geom_vline(xintercept = max_index, linetype = "dashed", color = "black", linewidth = .4) 
# SAVE
print(line_plot)
cat(file.path(output.dir, paste0(which.stat, "_OMD_draft_legend.pdf"))) # "_UHGG_draft_legend.pdf"
# ggsave(file.path(output.dir, paste0(which.stat, "_OMD_draft_legend.pdf")), line_plot, width = 8, height = 8, units = "cm", dpi = 300) # "_UHGG_draft_legend.pdf"
# _Ocean_7percFillvsISOFill.pdf

##################################
### FIGURE 3 in paper
lineplot_ocean <- line_plot 
lineplot_UHGG <- line_plot 

print(lineplot_ocean)
print(lineplot_UHGG)

# Combine the two plots
library(ggpubr)
combined_plot <- ggarrange( lineplot_UHGG + labs(title = "a."), lineplot_ocean + labs(title = "b."),
    ncol = 2, nrow = 1) # , widths = c(2, 1.2))
print(combined_plot)
cat(file.path(output.dir, paste0(which.stat, "_UHGG_Ocean_draft.pdf")))
ggsave(file.path(output.dir, paste0(which.stat, "_UHGG_Ocean_draft.pdf")), combined_plot,  units = "cm", width = 13, height = 8, dpi = 300)
##################################




## CORRELATIONs 
# How variable are the TH             # UHGG    # Ocean
mean(find_best_perf_Xsp$comp_lv)*100  # 6.625  # 12.11111
sd(find_best_perf_Xsp$comp_lv)*100    # 4.823487 # 6.153138

# Generate the dataset with the average of information and single pangenome values

# ISOLATES code
names(iso.binary_table.list)
path_metadata_40sp <- "/home/bioinfo/users/niber/prj_panModel/db/Paper/fig/metadata_40sp.txt"
metadata_40sp_df <- fread(path_metadata_40sp)
colnames(metadata_40sp_df)

metadata_40sp_df <- metadata_40sp_df  %>%
  filter(Genome_type == "MAG") %>%
  group_by(Species_rep) %>%
  mutate(mean_Length = mean(Length),
    mean_N_contigs = mean(N_contigs), 
    mean_N50 = mean(N50),
    mean_GC_content = mean(GC_content),
    mean_Completeness = mean(Completeness), 
    mean_Contamination = mean(Contamination)
  ) %>%
  distinct(., Species_rep, .keep_all = TRUE) %>% 
  select(Species_rep, n_MAG, n_REF, Number_of_samples, mean_Length, mean_N_contigs, mean_N50, mean_GC_content, mean_Completeness, mean_Contamination)

head(metadata_40sp_df)
colnames(metadata_40sp_df)

# pangenome info 
pangenome_metadata_fn <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/metadata/pangenome_metadata.csv"
pangenome_metadata <- read.table(pangenome_metadata_fn, sep = ";", header = TRUE, dec = ",")
col_of_interest <- c("Species.representative", "Core.genes", "Pan.genome.size", "Number.of.SNVs", "Geographical.diversity")
pangenome_metadata <- pangenome_metadata[, col_of_interest]
colnames(pangenome_metadata)
metadata <- merge(metadata, pangenome_metadata, by.x = "Genome_accession", by.y = "Species.representative", all.x = TRUE)


meanCont_best_perf_Xsp <- merge(metadata_40sp_df, find_best_perf_Xsp[, c("comp_lv", "ref_mod_id")], by.x="Species_rep", by.y="ref_mod_id")
dim(meanCont_best_perf_Xsp)
head(meanCont_best_perf_Xsp)

meanCont_best_perf_Xsp <- as.data.table(meanCont_best_perf_Xsp)
df <- meanCont_best_perf_Xsp[,-c("Species_rep")]
df_standardized <- scale(df)

# Initialize a data frame to store results
columns_of_interest <- c("n_MAG", "n_REF", "Number_of_samples",
"mean_Length",        "mean_N_contigs",     "mean_N50",          
"mean_GC_content",    "mean_Completeness",  "mean_Contamination",
"comp_lv")
correlation_results <- data.frame(column = character(), correlation = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Loop through columns of interest
for (col in columns_of_interest) {
  # Perform correlation test
  correlation_test <- cor.test(df_standardized[, "comp_lv"], df_standardized[, col])
  
  # Store results in the data frame
  correlation_results <- rbind(correlation_results, 
                                data.frame(column = col, 
                                           correlation = correlation_test$estimate, 
                                           p_value = correlation_test$p.value,
                                           stringsAsFactors = FALSE))
}
# Correct p-values for multiple testing using Benjamini-Hochberg correction
correlation_results$p_value_adjusted <- p.adjust(correlation_results$p_value, method = "BH")
print(correlation_results)

meanCont_best_perf_Xsp_plot <- meanCont_best_perf_Xsp[, c("Species_rep", "countMAG", "comp_lv", "countISO")] %>%
  arrange(comp_lv)

# Create a scatter plot
scatter_plot <- ggplot(meanCont_best_perf_Xsp_plot, aes(x = comp_lv, y = countISO)) +
  geom_point() +
  geom_smooth(se = FALSE, color = "blue") +  # Add the regression line
  labs(x = "Threashold", y = "countISO") +
  annotate("text", x = 0.16, y = 120, label = paste0("Corr. ", round(estimate, 2), "; p-value: ", round(p_value, 3))) +
  theme_minimal()

print(scatter_plot)
ggsave(file.path(output.dir, paste0("scatter_plot_CORR_nISOvsTh.pdf")), scatter_plot, width = 16, height = 8, units = "cm", dpi = 300)


