######################################################  @
# Script to reproduce the analysis in figure 2:
# 1. Completeness level of the "best MAG" within frequently detected SGBs
# 2. GEMs quality quantification with respect to reference genome models

# Load R-packages
suppressMessages(library(data.table)); setDTthreads(1)
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr)) 
library(RColorBrewer)
library(ggpubr)

# Little helpers
source("./scr/function_collection.R")

#####################
### UHGG Catalog
#####################

# Arguments:
mag.binary.rxn.table      <- "./dat/uhgg/rxnXmod/mag/allMAG_rxnXmod.txt"
iso.binary.rxn.table      <- "./dat/uhgg/rxnXmod/iso/allISO_rxnXmod.txt"
completeness.level.table  <- "./dat/uhgg/metadata/genomes-all_metadata.csv"
output.dir                <- "./dat/uhgg/out"

# Parameters - type of distance and score to compute between the reference-genome model and the individual GEMs:
dist_type <- "binary"
which.stat <- "f1_score" # other options are: fB_score, f1_score, recall, precision, accuracy, FDR

# Read metadata
metadata <- read.table(completeness.level.table, sep = ",", header = TRUE, row.names = "Genome")
columns_to_convert_numeric <- c("Completeness")
metadata[columns_to_convert_numeric] <- lapply(metadata[columns_to_convert_numeric], as.numeric)
metadata$Species_name <- rownames(metadata)

# Read Input
inIso_dt_list <- load_files_from_paths(iso.binary.rxn.table, ".tsv")
inMAG_dt_list <- load_files_from_paths(mag.binary.rxn.table, ".tsv")

# -------------
# Compute statistics: distance between GEMs to reference-genome models 
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

# combine metadata and model statistics
merged_df <- merge(all_dist_df, metadata, by.x="mod_id", by.y="Species_name")
plot_df <- parse_GTDBtaxonomy(merged_df, "Lineage") # GTDB.Taxonomy Lineage
plot_df <- refine_GTDBtaxonomy(plot_df)
plot_df$dist_tp <-as.factor(plot_df$dist_tp)

# -------------
# Correlation between completeness level and GEM quality
cor_test_result <- cor.test(plot_df$Completeness, plot_df$value, method = "spearman",  exact = FALSE)
cor_test_result <- cor.test(plot_df$Completeness, plot_df$value, method = "pearson")
cor_coefficient <- cor_test_result$estimate
p_value <- cor_test_result$p.value
cat("Spearman's Rank Correlation Coefficient:", cor_coefficient, "\n")
cat("P-value:", p_value, "\n")

# -------------
# Calculate how many SGB have few or no Ref and a high number of MAGs
metadata_REF_taxon <- metadata %>%
    filter(Species_name==Species_rep) %>%
    select("Species_name", "Lineage")
metadata_REF_taxon <- parse_GTDBtaxonomy(metadata_REF_taxon, "Lineage") # parse the GTDB taxonomy
metadata_REF_taxon <- refine_GTDBtaxonomy(metadata_REF_taxon) # refine the GTDB taxonomy names

# process the metadata to calculate stat on MAG number and completeness level for each SGB
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

cat(paste0("The number of species with at least 30 MAG, no reference genome and a best-MAG completeness level < 90 is: ",sum(plot_dt$bestMAG_Comp < 90),"/",dim(plot_dt)[1], "\n"))

# Read non-redundant genome sequence metadata


# Calculate how many SGB have few or no Ref and a high number of MAGs
metadata_REF_taxon <- metadata %>%
    filter(Species_name==Species_rep) %>%
    select("Species_name", "Lineage")
metadata_REF_taxon <- parse_GTDBtaxonomy(metadata_REF_taxon, "Lineage") # parse the GTDB taxonomy
metadata_REF_taxon <- refine_GTDBtaxonomy(metadata_REF_taxon) # refine the GTDB taxonomy names

# process the metadata to calculate stat on MAG number and completeness level for each SGB
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

# -------------
# prepare dataset for the first plot:
plot_dt <- merge(plot_dt, metadata_REF_taxon, by="Species_name")
# sort phylum based on the number of MAGs
desired_order <- plot_dt %>% 
    group_by(phylum) %>%
    summarise(count = n()) %>%
    arrange(desc(count)) 

# for graphical representation purposes: set the phylum order for the first 8 based on the MAG number and assigne the remaining as "other"
desired_order <- desired_order$phylum[1:8]
plot_dt <- plot_dt %>%
    mutate(phylum_color = ifelse(phylum %in% desired_order, phylum, "other"),
    phylum_color_noStat = ifelse(phylum %in% desired_order, phylum, "other")) 
plot_dt <- plot_dt %>% 
    group_by(phylum_color) %>% 
    mutate(count_SpXphylum = n(),
    count_MAGxphylum_count = sum(n_MAG)) %>%
    arrange(desc(count_SpXphylum)) 
plot_dt$phylum_color <- paste0(plot_dt$phylum_color, " (m=", plot_dt$count_MAGxphylum_count, "; s=", plot_dt$count_SpXphylum, ")")
plot_dt$phylum_color <- as.factor(plot_dt$phylum_color)
levels(plot_dt$phylum_color) <- unique(plot_dt$phylum_color)

plot_dt_UHGG <- plot_dt
plot_df_UHGG <- plot_df
phylum_in_plot_uhgg <- unique(plot_dt_UHGG$phylum_color_noStat)

# -------------
# prepare dataset for the second plot:
# relationship between completeness level of genome and quality of GEMs
single_metric_plot_df <- plot_df_UHGG[plot_df_UHGG$dist_tp==which.stat,] # extract the values of the desired metric
# single_metric_plot_df$phylum_color <- paste0(single_metric_plot_df$phylum_color, " (m=", single_metric_plot_df$count_MAGxphylum_count, "; s=", single_metric_plot_df$count_SpXphylum, ")")
# plot_dt$phylum_color <- as.factor(plot_dt$phylum_color)
levels(single_metric_plot_df$phylum)
single_metric_plot_df

#####################
### OMD dataset
#####################
# Input rxnXmod.txt
iso.binary.rxn.table      <- "./dat/omd/rxnXmod/iso/allISO_rxnXmod.txt"
mag.binary.rxn.table      <- "./dat/omd/rxnXmod/mag/allMAG_rxnXmod.txt"
oceans_metadata_fn        <- "./dat/omd/metadata/genomes-summary.csv"
output.dir                <- "./dat/omd/out"

# Read input
inIso_dt_list <- load_files_from_paths(iso.binary.rxn.table, ".tsv")
inMAG_dt_list <- load_files_from_paths(mag.binary.rxn.table, ".tsv")

# Input metadata
oceans_metadata <- read.table(oceans_metadata_fn, sep = ",", header = TRUE, dec = ".")
oceans_metadata <- oceans_metadata %>%  # distinguish genomes base on annotation (MAGs, reference genome and SAGs)
  mutate(Genome_type = if_else(grepl("_METAG_", Genome), "MAG", 
                        if_else(grepl("_REFG_", Genome), "REF", 
                        if_else(grepl("_SAGS_", Genome), "SAG", "")))
        )
oceans_metadata <- oceans_metadata %>% # specify the sample name
  mutate(sample = sapply(strsplit(oceans_metadata[["Genome"]], "_"), function(x) paste(x[1], x[2], x[3], sep="_")))
columns_to_convert_numeric <- c("Completeness")
oceans_metadata[columns_to_convert_numeric] <- lapply(oceans_metadata[columns_to_convert_numeric], as.numeric)
oceans_metadata$Species_name <- oceans_metadata$Genome
dim(oceans_metadata)

# -------------
# Compute statistics: distance between GEMs to reference-genome models
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

# combine metadata and model statistics
merged_df <- merge(all_dist_df, oceans_metadata, by.x="mod_id", by.y="Species_name")
plot_df <- parse_GTDBtaxonomy(merged_df, "GTDB.Taxonomy")
plot_df <- refine_GTDBtaxonomy(plot_df)
plot_df$dist_tp <-as.factor(plot_df$dist_tp)

# -------------
# Calculate how many SGB have few or no Ref and a high number of MAGs
metadata_REF_taxon <- oceans_metadata %>%
    filter(dRep.Representative.Genome) %>%
    select("Species_name", "GTDB.Taxonomy")
metadata_REF_taxon <- parse_GTDBtaxonomy(metadata_REF_taxon, "GTDB.Taxonomy") # GTDB.Taxonomy Lineage
metadata_REF_taxon <- refine_GTDBtaxonomy(metadata_REF_taxon)

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

cat(paste0("Number of species with at least 30 MAG and no reference genome: ",sum(plot_dt$bestMAG_Comp < 90),"/",dim(plot_dt)[1], "\n"))

# -------------
# prepare dataset for the first plot:
plot_dt <- merge(plot_dt, metadata_REF_taxon, by="Species_name")
desired_order <- plot_dt %>%
    group_by(phylum) %>%
    summarise(count = n()) %>%
    arrange(desc(count))

# for graphical representation purposes: set the phylum order for the first 8 based on the MAG number and assigne the remaining as "other"
desired_order <- desired_order$phylum[1:8]
plot_dt <- plot_dt %>%
    mutate(phylum_color = ifelse(phylum %in% desired_order, phylum, "other"),
    phylum_color_noStat = ifelse(phylum %in% desired_order, phylum, "other")) 
plot_dt <- plot_dt %>% 
    group_by(phylum_color) %>% 
    mutate(count_SpXphylum = n(),
    count_MAGxphylum_count = sum(n_MAG)) %>%
    arrange(desc(count_SpXphylum)) 
plot_dt$phylum_color <- paste0(plot_dt$phylum_color, " (m=", plot_dt$count_MAGxphylum_count, "; s=", plot_dt$count_SpXphylum, ")")
plot_dt$phylum_color <- as.factor(plot_dt$phylum_color)
levels(plot_dt$phylum_color) <- unique(plot_dt$phylum_color)

plot_dt_OMD <- plot_dt
plot_df_OMD <- plot_df
phylum_in_plot_OMD <- unique(plot_dt_OMD$phylum_color_noStat)

# prepare dataset for the second plot:
single_metric_plot_OMD <- plot_df_OMD[plot_df_OMD$dist_tp==which.stat,]
dim(single_metric_plot_OMD)

head(single_metric_plot_df)

# -------------
# PLOTs
palette <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Dark2"))  # You can specify the number of colors you want
phylum_in_plot <- unique(c(phylum_in_plot_uhgg, phylum_in_plot_OMD)) # you have to compute this before

# fig. 2a best MAG vs. number of MAGs per SGB (no isolates)
set1_palette <- unique(palette[match(plot_dt_UHGG$phylum_color_noStat, phylum_in_plot)])
new_order <- unique(plot_dt_UHGG$phylum_color_noStat)
plot_dt_UHGG$phylum_color_noStat <- factor(plot_dt_UHGG$phylum_color_noStat, levels = new_order)
levels(plot_dt_UHGG$phylum_color_noStat)

p <- ggplot(plot_dt_UHGG, aes(x = n_MAG, y = bestMAG_Comp, color = phylum_color)) + # phylum_color_noStat phylum_color # factor(n_Iso)
    geom_point(alpha = 0.4, size = 2) +
    scale_x_log10(limits = c(30, 6000)) +
    labs(
        x = "Number of MAGs",
        y = "Maximum completeness (%)",
        size = 8) +
    labs(color = "UHGG") + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    legend.position =  c(.8, 0.4),
    legend.text = element_text(size = 6)
    ) +
    scale_color_manual(values = set1_palette) +
    ylim(75, 100)
p
# ggsave(file.path(output.dir, "BestCompletenessMAG_x_Sp_noIsolates.pdf"), p, width = 6, height = 6, units = "in", dpi = 300)


# fig. 2b MAG's completeness level vs. GEMs quality
set1_palette <- unique(palette[match(single_metric_plot_df$phylum, phylum_in_plot)])
new_order <- unique(single_metric_plot_df$phylum)
single_metric_plot_df$phylum <- factor(single_metric_plot_df$phylum, levels=new_order)
sum_plot <- ggplot(single_metric_plot_df, aes(x=Completeness, y=value, color=phylum)) + 
    geom_point(alpha = 0.3, show.legend = TRUE, size = 1) +
    geom_smooth(method = "auto", se = TRUE, show.legend = FALSE) +  # Add trend lines "lm"
    labs(
        x = "Completeness (%)",
        y = "F1 score",
        size = 8) + # which.stat
        # title = "Fit a generalized additive model 'mgcv::gam' with formula = 'y ~ s(x, bs = 'cs')") +
    labs(color = "UHGG") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    legend.position =  c(.8, 0.3),
    legend.text = element_text(size = 6)) +
    scale_color_manual(values = set1_palette)
sum_plot
# ggsave(file.path(output.dir, paste0(which.stat,"_MAGsCLv.vs.GEMsQ.pdf")), sum_plot, width = 10, height = 8, units = "in", dpi = 300)


# fig. 2c best MAG vs. number of MAGs per SGB (no isolates)
set1_palette <- unique(palette[match(plot_dt_OMD$phylum_color_noStat, phylum_in_plot)])
new_order <- unique(plot_dt_OMD$phylum_color_noStat)
plot_dt_OMD$phylum_color_noStat <- factor(plot_dt_OMD$phylum_color_noStat, levels = new_order)

p_ocean <- ggplot(plot_dt_OMD, aes(x = n_MAG, y = bestMAG_Comp, color = phylum_color)) + # factor(n_Iso)
    geom_point(alpha = 0.4, size = 2) +
    scale_x_log10(limits = c(30, 400)) +
    labs(
        x = "Number of MAGs", 
        y = "Maximum completeness (%)",
        size = 8) +
    labs(color = "OMD") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    legend.position =  c(.8, 0.3),
    legend.text = element_text(size = 6)
    ) + 
    scale_color_manual(values = set1_palette) +
    ylim(60, 100)
p_ocean
# ggsave(file.path(output.dir, "BestCompletenessMAG_x_Sp_noIsolates.pdf"), p_ocean, width = 6, height = 6, units = "in", dpi = 300)


# fig. 2d MAG's completeness level vs. GEMs quality
set1_palette <- unique(palette[match(single_metric_plot_OMD$phylum, phylum_in_plot)])
new_order <- unique(single_metric_plot_OMD$phylum)
single_metric_plot_OMD$phylum <- factor(single_metric_plot_OMD$phylum, levels = new_order)

sum_plot_ocean <- ggplot(single_metric_plot_OMD, aes(x=Completeness, y=value, color=phylum)) + 
    geom_point(alpha = 0.3, show.legend = TRUE, size = 1) +
    geom_smooth(method = "auto", se = TRUE, show.legend = FALSE) +  # Add trend lines "lm"
    labs(
        x = "Completeness (%)",
        y = "F1 score",
        size = 8) + # which.stat
        # title = "Fit a generalized additive model 'mgcv::gam' with formula = 'y ~ s(x, bs = 'cs')") +
    labs(color = "OMD") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    legend.position =  c(.8, 0.3),
    legend.text = element_text(size = 6)) +
    scale_color_manual(values = set1_palette) 
# out_fn <- paste0(which.stat,"_MAGsCLv.vs.GEMsQ.pdf")
# ggsave(file.path(output.dir, out_fn), sum_plot_ocean, width = 10, height = 8, units = "in", dpi = 300)


# -------------
### FIGURE 2 in paper
output.dir <- "./dat/fig"
# Combine the two plots
combined_plot <- ggarrange( p + labs(title = "a."), p_ocean + labs(title = "b."),
    sum_plot + labs(title = "c."), sum_plot_ocean + labs(title = "d."),
    ncol = 2, nrow = 2, widths = c(1, 1))
ggsave(file.path(output.dir, "fig2.pdf"), combined_plot,  units = "cm", width = 16, height = 16, dpi = 300)
