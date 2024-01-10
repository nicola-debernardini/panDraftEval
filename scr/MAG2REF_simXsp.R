# Load R-packages
suppressMessages(library(data.table)); setDTthreads(1)
suppressMessages(library(ggplot2)) # Graphical
suppressMessages(library(dplyr)) # Graphical
library(ggpubr)
library(RColorBrewer)


# Little helpers
source("./scr/function_collection.R")

# Arguments:
mag.binary.rxn.table      <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft_18Oct/gapfill/MAG/allMAG_rxnXmod.txt" # "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft/allMAG_rxnXmod.txt"
iso.binary.rxn.table      <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft_18Oct/gapfill/ISO/allISO_rxnXmod.txt" # "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft/allISO_rxnXmod.txt"
completeness.level.table  <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/metadata/genomes-all_metadata.csv"
# output.dir                <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft_18Oct/fig"

# Parameters:
dist_type <- "binary"
which.stat <- "f1_score" # fB_score f1_score recall precision accuracy FDR

# Read metadata
metadata <- read.table(completeness.level.table, sep = ",", header = TRUE, row.names = "Genome") # dec = ","
columns_to_convert_numeric <- c("Completeness")
metadata[columns_to_convert_numeric] <- lapply(metadata[columns_to_convert_numeric], as.numeric)
metadata$Species_name <- rownames(metadata)

############################################################################
# Accuracy of GEMs compared to reference genome based on MAG completeness - NO PAN MODELs
### Read Input
inIso_dt_list <- load_files_from_paths(iso.binary.rxn.table, ".tsv")
inMAG_dt_list <- load_files_from_paths(mag.binary.rxn.table, ".tsv")

all_dist_df <- data.frame(mod_id = character(), ref_mod_id = character(),  dist_tp = character(), mod_category = character(), value = numeric(0))
for (mod in names(inIso_dt_list)) {
    cat(paste(mod,"\n"))   
    inIso_dt <- inIso_dt_list[[mod]]
    inMAG_dt <- inMAG_dt_list[[mod]]

    # Calculate frequency of rxn in all the reference genomes 
    Isolate_cols <- names(inIso_dt)[names(inIso_dt)!="rxn"]
    freq_ref_rxn <- inIso_dt %>%
    select(all_of(Isolate_cols)) %>%
    mutate(b_i = rowSums(.)/length(Isolate_cols),
        rxn = inIso_dt$rxn) %>%
    select(b_i, rxn)

    reference_reactions <- inIso_dt$rxn

    # Average Tsgb of the MAGs
    # dist_val  <- distance_cal(inMAG_dt, freq_ref_rxn, dist_type)   
    # dist_val_matr <- as.matrix(dist_val)
    
    for (mag in colnames(inMAG_dt)[-1]){ # exclude "rxn"
        # Tsgb <- Tsgb_cal(inMAG_dt, freq_ref_rxn, mag)
        # dist_val <- dist_val_matr["ref", mag]

        predicted_reactions <- inMAG_dt[["rxn"]][as.logical(inMAG_dt[[mag]])]
        metrics_table <- compute_contigency_table(reference_reactions, predicted_reactions) # compute contingency table for a pan.mod
        stat_val <- metrics_table[which.stat][[1]]
        
        # Populate the dataframe
        # new_row_Tsgb <- list(mod_id = mag, 
        #                 ref_mod_id = mod,
        #                 dist_tp = "Tsgb", 
        #                 mod_category = "mag", 
        #                 value = Tsgb)
        # new_row_dist <- list(mod_id = mag,
        #                 ref_mod_id = mod, 
        #                 dist_tp = dist_type, 
        #                 mod_category = "mag", 
        #                 value = dist_val)
        new_row_stat <- list(mod_id = mag,
                        ref_mod_id = mod, 
                        dist_tp = which.stat, 
                        mod_category = "mag", 
                        value = stat_val)
        all_dist_df <- rbind(all_dist_df, new_row_stat) # new_row_Tsgb, new_row_dist)
    }
}

merged_df <- merge(all_dist_df, metadata, by.x="mod_id", by.y="Species_name")
plot_df <- parse_GTDBtaxonomy(merged_df, "Lineage") # GTDB.Taxonomy Lineage
plot_df <- refine_GTDBtaxonomy(plot_df)
plot_df$dist_tp <-as.factor(plot_df$dist_tp)
# ALREADY CALCULATED
# fwrite(plot_df, file = file.path("/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft/tbl/F1TsgbJacc_allMAGs_vs_ISOs.tsv"), sep = "\t")
# plot_df <- read.table("/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft/tbl/F1TsgbJacc_allMAGs_vs_ISOs.tsv", sep = "\t", header = TRUE)
head(plot_df)
dim(plot_df)


# Correlation
cor_test_result <- cor.test(plot_df$Completeness, plot_df$value, method = "spearman",  exact = FALSE)
cor_test_result <- cor.test(plot_df$Completeness, plot_df$value, method = "pearson")
cor_coefficient <- cor_test_result$estimate
p_value <- cor_test_result$p.value
cat("Spearman's Rank Correlation Coefficient:", cor_coefficient, "\n")
cat("P-value:", p_value, "\n")

# Assess what is the average improvement level of MAGs with a 90% completeness level
mean(plot_df[plot_df$Completeness<71 & plot_df$Completeness>69, "value"]) 
sd(plot_df[plot_df$Completeness<71 & plot_df$Completeness>69, "value"]) 
mean(plot_df[plot_df$Completeness<81 & plot_df$Completeness>79, "value"]) 
sd(plot_df[plot_df$Completeness<81 & plot_df$Completeness>79, "value"]) 
mean(plot_df[plot_df$Completeness<91 & plot_df$Completeness>89, "value"]) 
sd(plot_df[plot_df$Completeness<91 & plot_df$Completeness>89, "value"]) 
mean(plot_df[plot_df$Completeness>=99, "value"]) 
sd(plot_df[plot_df$Completeness>=99, "value"]) 

# -------------
### EVALUATE MAGs with few or no Ref and high number of MAGs
metadata_REF_taxon <- metadata %>%
    filter(Species_name==Species_rep) %>%
    select("Species_name", "Lineage")
metadata_REF_taxon <- parse_GTDBtaxonomy(metadata_REF_taxon, "Lineage") # GTDB.Taxonomy Lineage
metadata_REF_taxon <- refine_GTDBtaxonomy(metadata_REF_taxon)

# calculate stat on MAG number and completeness level per species
plot_dt <- metadata %>%
    group_by(Species_rep) %>%
    mutate(n_Iso = sum(Genome_type=="Isolate"), 
    n_MAG = sum(Genome_type=="MAG"), 
    MAGComp_mean = mean(Completeness),
    MAGComp_std = sd(Completeness),
    bestMAG_Comp = max(Completeness)) %>%
    filter(n_Iso < 1) %>%
    filter(n_MAG >= 30) %>%
    arrange(n_MAG) %>%
    filter(Species_name==Species_rep) %>% # select only representative
    ungroup()

cat(paste0("Number of species with at least 30 MAG and no reference genome: ",sum(plot_dt$bestMAG_Comp < 90),"/",dim(plot_dt)[1], "\n"))

# sort phylum based on the number of MAGs
plot_dt <- merge(plot_dt, metadata_REF_taxon, by="Species_name")
desired_order <- plot_dt %>% 
    group_by(phylum) %>%
    summarise(count = n()) %>%
    arrange(desc(count)) 

# assign the value other if more that eight phylum where present
desired_order <- desired_order$phylum[1:8]
plot_dt <- plot_dt %>%
    mutate(phylum_color = ifelse(phylum %in% desired_order, phylum, "other"),
    phylum_color_noStat = ifelse(phylum %in% desired_order, phylum, "other")) 
plot_dt <- plot_dt %>% 
    group_by(phylum_color) %>% 
    mutate(count_SpXphylum = n(),
    count_MAGxphylum_count = sum(n_MAG)) %>%
    arrange(desc(count_SpXphylum)) 

plot_dt$phylum_color <- paste0(plot_dt$phylum_color, " (m. ", plot_dt$count_MAGxphylum_count, "; s. ", plot_dt$count_SpXphylum, ")")
plot_dt$phylum_color <- as.factor(plot_dt$phylum_color)
levels(plot_dt$phylum_color) <- unique(plot_dt$phylum_color)

# desired_order <- c("Actinobacteriota", "Bacteroidota", "Firmicutes", "Firmicutes_A", "Proteobacteria", "Verrucomicrobiota")
# plot_dt <- plot_dt %>%
#     mutate(phylum_color = ifelse(phylum %in% desired_order, phylum, "other")) 
# plot_dt$phylum_color <- as.factor(plot_dt$phylum_color)
# levels(plot_dt$phylum_color) <- c(desired_order, "other") # setdiff(all_levels, desired_order) 

### Relationship between completeness and accuracy or other metrics
set1_palette <- brewer.pal(9, "Set1")  # You can specify the number of colors you want
plot_dt_UHGG <- plot_dt

# PLOT BestCompleteness MAG x Species noIsolates
p <- ggplot(plot_dt, aes(x = n_MAG, y = bestMAG_Comp, color = phylum_color)) + # factor(n_Iso)
  geom_point(alpha = 0.4, size = 2) + #shape = 1
  scale_x_log10() +
  labs(x = "Number of MAG", y = "Completeness (%)") +
  labs(color = "Phylum") + # "#Isolates"
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
  legend.position = "left",
  legend.text = element_text(size = 8)
  ) + 
  scale_color_manual(values = set1_palette)  # Use the Set1 palette 
print(p)
ggsave(file.path(output.dir, "BestCompletenessMAG_x_Sp_noIsolates.pdf"), p, width = 6, height = 6, units = "in", dpi = 300)


### PLOT: Relationship between completeness and accuracy or other metrics
what <- which.stat # which.stat dist_type "Tsgb"
single_metric_plot_df <- plot_df[plot_df$dist_tp==what,]
common_values <- intersect(single_metric_plot_df$ref_mod_id, single_metric_plot_df$mod_id)
print(common_values)

single_metric_plot_df[single_metric_plot_df$mod_id==single_metric_plot_df$ref_mod_id]

phylum_in_plot_p <- unique(plot_dt$phylum_color_noStat)
palette <- unique(set1_palette[match(single_metric_plot_df$phylum, phylum_in_plot_p)])

sum_plot <- ggplot(single_metric_plot_df, aes(x=Completeness, y=value, color=phylum)) + 
    geom_point(alpha = 0.3, show.legend = FALSE, size = 1) +
    geom_smooth(method = "auto", se = TRUE, show.legend = FALSE) +  # Add trend lines "lm"
    labs(
        x = "Completeness (%)",
        y = "F1 score") + # what
        # title = "Fit a generalized additive model 'mgcv::gam' with formula = 'y ~ s(x, bs = 'cs')") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+             
    scale_color_manual(values = palette) # Use the Set1 palette
print(sum_plot)

# out_fn <- paste0(what,"_MAGs_dotplot.pdf")
# ggsave(file.path(output.dir, out_fn), sum_plot, width = 10, height = 8, units = "in", dpi = 300)


#####################
### Ocean
#####################
# Input rxnXmod.txt
iso.binary.rxn.table      <- "/home/bioinfo/users/niber/prj_panModel/db/microbiomics_ocean/pan.draft/gapfill/ISO/allISO_rxnXmod.txt"
mag.binary.rxn.table      <- "/home/bioinfo/users/niber/prj_panModel/db/microbiomics_ocean/pan.draft/gapfill/MAG/allMAG_rxnXmod.txt"
oceans_metadata_fn        <- "/home/bioinfo/users/niber/prj_panModel/db/microbiomics_ocean/metadata/genomes-summary.csv"  
output.dir                <- "/home/bioinfo/users/niber/prj_panModel/db/microbiomics_ocean/pan.draft/fig"

############################################################################
# Accuracy of GEMs compared to reference genome based on MAG completeness - NO PAN MODELs
### Read Input
inIso_dt_list <- load_files_from_paths(iso.binary.rxn.table, ".tsv")
inMAG_dt_list <- load_files_from_paths(mag.binary.rxn.table, ".tsv")

# Input metadata
oceans_metadata <- read.table(oceans_metadata_fn, sep = ",", header = TRUE, dec = ".")
oceans_metadata <- oceans_metadata %>% # separate MAGs from reference genome and SAGs
  mutate(Genome_type = if_else(grepl("_METAG_", Genome), "MAG", if_else(grepl("_REFG_", Genome), "REF", if_else(grepl("_SAGS_", Genome), "SAG", ""))))
oceans_metadata <- oceans_metadata %>% # specify the sample name
  mutate(sample = sapply(strsplit(oceans_metadata[["Genome"]], "_"), function(x) paste(x[1], x[2], x[3], sep="_")))
columns_to_convert_numeric <- c("Completeness")
oceans_metadata[columns_to_convert_numeric] <- lapply(oceans_metadata[columns_to_convert_numeric], as.numeric)
oceans_metadata$Species_name <- oceans_metadata$Genome
dim(oceans_metadata)
head(oceans_metadata)

head(oceans_metadata[oceans_metadata$detected.in.the.water.column=="true",])
dim(oceans_metadata[oceans_metadata$detected.in.the.water.column=="true" & oceans_metadata$depth.layer=="EPI",])
dim(oceans_metadata[oceans_metadata$detected.in.the.water.column=="false" & oceans_metadata$depth.layer=="EPI",])
# 33647/34815 have been detected in the water columns
# 20482/34815 have been detected in the water columns
sum(oceans_metadata$oxygen..µmol.kg. > 200, na.rm = TRUE)


all_dist_df <- data.frame(mod_id = character(), ref_mod_id = character(),  dist_tp = character(), mod_category = character(), value = numeric(0))
for (mod in names(inIso_dt_list)) {
    cat(paste(mod,"\n"))   
    inIso_dt <- inIso_dt_list[[mod]]
    inMAG_dt <- inMAG_dt_list[[mod]]

    # Calculate frequency of rxn in all the reference genomes 
    Isolate_cols <- names(inIso_dt)[names(inIso_dt)!="rxn"]
    freq_ref_rxn <- inIso_dt %>%
    select(all_of(Isolate_cols)) %>%
    mutate(b_i = rowSums(.)/length(Isolate_cols),
        rxn = inIso_dt$rxn) %>%
    select(b_i, rxn)

    reference_reactions <- inIso_dt$rxn

    # Average Tsgb of the MAGs
    # dist_val  <- distance_cal(inMAG_dt, freq_ref_rxn, dist_type)   
    # dist_val_matr <- as.matrix(dist_val)
    for (mag in colnames(inMAG_dt)[-1]){ # exclude "rxn"
        # Tsgb <- Tsgb_cal(inMAG_dt, freq_ref_rxn, mag)
        # dist_val <- dist_val_matr["ref", mag]

        predicted_reactions <- inMAG_dt[["rxn"]][as.logical(inMAG_dt[[mag]])]
        metrics_table <- compute_contigency_table(reference_reactions, predicted_reactions) # compute contingency table for a pan.mod
        stat_val <- metrics_table[which.stat][[1]]
        
        # Populate the dataframe
        # new_row_Tsgb <- list(mod_id = mag, 
        #                 ref_mod_id = mod,
        #                 dist_tp = "Tsgb", 
        #                 mod_category = "mag", 
        #                 value = Tsgb)
        # new_row_dist <- list(mod_id = mag,
        #                 ref_mod_id = mod, 
        #                 dist_tp = dist_type, 
        #                 mod_category = "mag", 
        #                 value = dist_val)
        new_row_stat <- list(mod_id = mag,
                        ref_mod_id = mod, 
                        dist_tp = which.stat, 
                        mod_category = "mag", 
                        value = stat_val)
        all_dist_df <- rbind(all_dist_df, new_row_stat) # new_row_Tsgb, new_row_dist)
    }
}
# 
merged_df <- merge(all_dist_df, oceans_metadata, by.x="mod_id", by.y="Species_name")
plot_df <- parse_GTDBtaxonomy(merged_df, "GTDB.Taxonomy") # GTDB.Taxonomy Lineage
plot_df <- refine_GTDBtaxonomy(plot_df)
plot_df$dist_tp <-as.factor(plot_df$dist_tp)
# ALREADY CALCULATED
# fwrite(plot_df, file = file.path("/home/bioinfo/users/niber/prj_panModel/db/microbiomics_ocean/pan.draft/tbl/F1TsgbJacc_allMAGs_vs_ISOs.tsv"), sep = "\t")
# plot_df <- read.table("/home/bioinfo/users/niber/prj_panModel/db/microbiomics_ocean/pan.draft/tbl/F1TsgbJacc_allMAGs_vs_ISOs.tsv", sep = "\t", header = TRUE)
head(plot_df)
dim(plot_df)

# Assess what is the average improvement level of MAGs with a 90% completeness level
dim(plot_df[plot_df$Completeness<71 & plot_df$Completeness>69, ]) 
mean(plot_df[plot_df$Completeness<71 & plot_df$Completeness>69, "value"]) 
sd(plot_df[plot_df$Completeness<71 & plot_df$Completeness>69, "value"]) 
dim(plot_df[plot_df$Completeness<81 & plot_df$Completeness>79, ]) 
mean(plot_df[plot_df$Completeness<81 & plot_df$Completeness>79, "value"]) 
sd(plot_df[plot_df$Completeness<81 & plot_df$Completeness>79, "value"]) 
dim(plot_df[plot_df$Completeness<91 & plot_df$Completeness>89, ]) 
mean(plot_df[plot_df$Completeness<91 & plot_df$Completeness>89, "value"]) 
sd(plot_df[plot_df$Completeness<91 & plot_df$Completeness>89, "value"]) 
dim(plot_df[plot_df$Completeness>=99, ]) 
mean(plot_df[plot_df$Completeness>=99, "value"]) 
sd(plot_df[plot_df$Completeness>=99, "value"]) 

### EVALUATE MAGs with few or no Ref and high number of MAGs
metadata_REF_taxon <- oceans_metadata %>%
    filter(dRep.Representative.Genome) %>%
    select("Species_name", "GTDB.Taxonomy")
metadata_REF_taxon <- parse_GTDBtaxonomy(metadata_REF_taxon, "GTDB.Taxonomy") # GTDB.Taxonomy Lineage
metadata_REF_taxon <- refine_GTDBtaxonomy(metadata_REF_taxon)
head(metadata_REF_taxon)

plot_dt <- oceans_metadata %>%
    group_by(dRep.Species.Cluster) %>%
    mutate(n_Iso = sum(Genome_type=="REF"),
    n_MAG = sum(Genome_type=="MAG"), 
    MAGComp_mean = mean(Completeness),
    MAGComp_std = sd(Completeness),
    bestMAG_Comp = max(Completeness)) %>%
    filter(n_Iso < 1) %>%
    filter(n_MAG >= 30) %>%
    arrange(n_MAG) %>%
    filter(dRep.Representative.Genome) %>%
    ungroup()
head(plot_dt)

cat(paste0("Number of species with at least 30 MAG and no reference genome: ",sum(plot_dt$bestMAG_Comp < 90),"/",dim(plot_dt)[1], "\n"))

plot_dt <- merge(plot_dt, metadata_REF_taxon, by="Species_name")
desired_order <- plot_dt %>% 
    group_by(phylum) %>%
    summarise(count = n()) %>%
    arrange(desc(count)) 

desired_order <- desired_order$phylum[1:8]
plot_dt <- plot_dt %>%
    mutate(phylum_color = ifelse(phylum %in% desired_order, phylum, "other"),
    phylum_color_noStat = ifelse(phylum %in% desired_order, phylum, "other")) 
plot_dt <- plot_dt %>% 
    group_by(phylum_color) %>% 
    mutate(count_SpXphylum = n(),
    count_MAGxphylum_count = sum(n_MAG)) %>%
    arrange(desc(count_SpXphylum)) 

plot_dt$phylum_color <- paste0(plot_dt$phylum_color, " (m. ", plot_dt$count_MAGxphylum_count, "; s. ", plot_dt$count_SpXphylum, ")")
plot_dt$phylum_color <- as.factor(plot_dt$phylum_color)
levels(plot_dt$phylum_color) <- unique(plot_dt$phylum_color)

### Relationship between completeness and accuracy or other metrics
set1_palette <- brewer.pal(9, "Set1")  # You can specify the number of colors you want

# PLOT BestCompleteness MAG x Species noIsolates
p_ocean <- ggplot(plot_dt, aes(x = n_MAG, y = bestMAG_Comp, color = phylum_color)) + # factor(n_Iso)
  geom_point(alpha = 0.4, size = 2) + #shape = 1
  scale_x_log10() +
  labs(x = "Number of MAG", y = "Completeness (%)") +
  labs(color = "Phylum") + # "#Isolates"
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
  legend.position = "left",
  legend.text = element_text(size = 8)
  ) + 
  scale_color_manual(values = set1_palette)  # Use the Set1 palette 
print(p_ocean)
# ggsave(file.path(output.dir, "BestCompletenessMAG_x_Sp_noIsolates.pdf"), p_ocean, width = 6, height = 6, units = "in", dpi = 300)


### Relationship between completeness and accuracy or other metrics
single_metric_plot_df <- plot_df[plot_df$dist_tp==what,]
head(single_metric_plot_df)
dim(single_metric_plot_df)

# calculate number of species and MAG per phylum
result <- single_metric_plot_df %>%
  group_by(ref_mod_id) %>%
  summarise(count_mod_id = n_distinct(mod_id))
result_merged_df <- merge(result, oceans_metadata, by.x="ref_mod_id", by.y="Species_name")
result_merged_df <- parse_GTDBtaxonomy(result_merged_df, "GTDB.Taxonomy") # GTDB.Taxonomy Lineage
result_merged_df <- refine_GTDBtaxonomy(result_merged_df)

dim(result_merged_df)
result_mag <- result_merged_df %>%
  group_by(phylum) %>%
  summarise(count_mod_id = sum(count_mod_id)) %>%
  arrange(count_mod_id)
result_sp <- result_merged_df %>%
  group_by(phylum) %>%
  summarise(count_mod_id = n_distinct(ref_mod_id)) %>%
  arrange(count_mod_id)
# which Species has enough MAGs in the ocean in order to calculate the panDraft? 
single_metric_plot_df %>%
  group_by(ref_mod_id) %>%
  summarize(count_mod_id = n_distinct(mod_id),
  Comp50_60 = sum(Completeness>=50 & Completeness<60),
  Comp60_70 = sum(Completeness>=60 & Completeness<70),
  Comp70_80 = sum(Completeness>=70 & Completeness<80),
  Comp80_90 = sum(Completeness>=80 & Completeness<90),
  Comp90_100 = sum(Completeness>=90 & Completeness<100))
# OUTPUT --none-- 
# ref_mod_id     count_mod_id Comp50_60 Comp60_70 Comp70_80 Comp80_90 Comp90_100
# MALA_SAMN0542…           64         7         5         7        19         26
# MARD_SAMEA470…           40         7         9        14        10          0
# MARD_SAMN0449…           52         9         4         9        18         12
# MARD_SAMN0521…           37         4         5        11         8          9
# MARD_SAMN0532…           50         6         7         6        13         18
# MARD_SAMN0876…           33         9         5         9         5          4
# TARA_SAMEA262…           32         4         7         4         6         11
# TARA_SAMEA262…          115        14        20        17        44         19
# TARA_SAMEA439…           31         5         5         9         3          5

phylum_in_plot_p <- unique(plot_dt$phylum_color_noStat)
palette <- unique(set1_palette[match(single_metric_plot_df$phylum, phylum_in_plot_p)])

sum_plot_ocean <- ggplot(single_metric_plot_df, aes(x=Completeness, y=value, color=phylum)) + 
    geom_point(alpha = 0.3, show.legend = FALSE, size = 1) +
    geom_smooth(method = "auto", se = TRUE, show.legend = FALSE) +  # Add trend lines "lm"
    labs(
        x = "Completeness (%)",
        y = "F1 score") + # what
        # title = "Fit a generalized additive model 'mgcv::gam' with formula = 'y ~ s(x, bs = 'cs')") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +             
    scale_color_manual(values = palette) 
    
print(sum_plot_ocean)
# out_fn <- paste0(what,"_MAGs_dotplot.pdf")
# ggsave(file.path(output.dir, out_fn), sum_plot_ocean, width = 10, height = 8, units = "in", dpi = 300)

##################################
### FIGURE 2 in paper
output.dir <- "/home/bioinfo/users/niber/prj_panModel/db/Paper/fig"
# Combine the two plots
combined_plot <- ggarrange( p + labs(title = "a."), sum_plot + labs(title = "b."),
    p_ocean + labs(title = "c."), sum_plot_ocean + labs(title = "d."),
    ncol = 2, nrow = 2, widths = c(2, 1.2))
print(combined_plot)
ggsave(file.path(output.dir, paste0("BestCompletenessMAG_x_Sp_noIsolates_AND_", what, "_MAGs_dotplot.pdf")), combined_plot,  units = "cm", width = 16, height = 16, dpi = 300)
##################################




##############################
# ADDITIONAL PLOTs GUT.... devi ricaricare i dati se hai lanciato quelli del Ocean
# PLOT Avg. MAG x Species noIsolates
p <- ggplot(plot_dt, aes(x = n_MAG, y = MAGComp_mean, color = phylum_color, size = MAGComp_std)) + # color = factor(n_Iso)
  geom_point(shape = 1) +
  scale_size(range = c(0.5, 4)) +  # Adjust the range as needed
  scale_x_log10() +
  labs(x = "number of MAG", y = "Avg. Completeness") +
  labs(color = "#Isolates", size = "Comp. std") +
  theme_minimal() +
  scale_color_manual(values = set1_palette)  # Use the Set1 palette 
print(p)
ggsave(file.path(output.dir, "Avg.CompletenessMAG_x_Sp_noIsolates.pdf"), p, width = 6, height = 6, units = "in", dpi = 300)

### Comparison between metrics 
what <- c(which.stat, dist_type, "Tsgb") # which.stat dist_type "Tsgb"
single_metric_plot_df <- plot_df[plot_df$dist_tp %in% what,]
head(single_metric_plot_df)

sum_plot <- ggplot(single_metric_plot_df, aes(x=species, y=value, fill=dist_tp)) + 
    geom_boxplot() +
    labs(
        x = "Species",
        y = "Similarities level") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text = element_text(size = 14), 
            axis.title = element_text(size = 14))
print(sum_plot)
out_fn <- paste0("F1AndJaccAndTsgb_boxplot.pdf")
ggsave(file.path(output.dir, out_fn), sum_plot, width = 10, height = 8, units = "in", dpi = 300)