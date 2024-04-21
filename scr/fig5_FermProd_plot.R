######################################################
# Script to reproduce the fermentation test (figure 5):

# Load R-packages
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(viridis)
library(tidyr)

# Little helpers
source("./scr/function_collection.R")

calculate_confusion_matrix <- function(true_values, predicted_values) {
  # Ensure both vectors are of the same length
  if (length(true_values) != length(predicted_values)) {
    stop("Vectors must have the same length.")
  }
  
  # Convert the vectors to logical
  true_values <- as.logical(true_values)
  predicted_values <- as.logical(predicted_values)
  
  # Calculate TP, TN, FP, FN
  TP <- sum(true_values & predicted_values)
  TN <- sum(!true_values & !predicted_values)
  FP <- sum(!true_values & predicted_values)
  FN <- sum(true_values & !predicted_values)
  
  # Return the confusion matrix
  return(c(TP = TP, TN = TN, FP = FP, FN = FN))
}

calculate_normalized_confusion_matrix <- function(true_values, predicted_values) {

  # Ensure both vectors are of the same length
  if (length(true_values) != length(predicted_values)) {
    stop("Vectors must have the same length.")
  }
  
  # Convert the vectors to logical
  true_values <- as.logical(true_values)
  TP <- sum(predicted_values[true_values])
  FP <- sum(predicted_values[!true_values])
  realP <- sum(true_values)
  realN <- sum(!true_values)

  # Return the confusion matrix
  return(c(TP = TP, TN = realN-FP, FP = FP, FN = realP-TP))
}

calculate_performance_array <- function(true_values, predicted_values) {

  # Ensure both vectors are of the same length
  if (length(true_values) != length(predicted_values)) {
    stop("Vectors must have the same length.")
  }

  # Convert the vectors to logical
  true_values <- as.logical(true_values)
  predicted_values <- as.logical(predicted_values)
  
  # Calculate TP, TN, FP, FN
  TP <- as.numeric(true_values & predicted_values)
  TN <- as.numeric(!true_values & !predicted_values)
  FP <- as.numeric(!true_values & predicted_values)
  FN <- as.numeric(true_values & !predicted_values)
  
  res <- data.frame(TP=TP, TN=TN, FP=FP, FN=FN)
  # Return the confusion matrix
  return(res)
}

# Arguments:
output.dir                <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/fermProdTest/out_7Feb_2024"
gs.fermprod_pan <- read.table(file.path(output.dir, "panDraft_fermProd_31Gen.tsv"), header = TRUE, sep = "\t", quote = "")
gs.fermprod <- read.table(file.path(output.dir, "mag_fermProd_6Feb.tsv"), header = TRUE, sep = "\t", quote = "")
dim(gs.fermprod_pan)
dim(gs.fermprod)

# pan-Draft
panMAG.CompLV.gapfill.model.path  <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/fermProdTest/panFromMAG_GapfilledFT_Path.txt"
panMAG.CompLV.gapfill.mods.fn <- readLines(panMAG.CompLV.gapfill.model.path)

mag.gapfill.model.pathList        <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/fermProdTest/MAG_GapfilledFT_listsPath.txt"
singleGeno_pathList <- readLines(mag.gapfill.model.pathList)

# fermentation product
fermentation_product_fn <- "/home/bioinfo/projects/panGenome/panDraftEval/data/phenotypeEval/fermentation_product_info.csv"
fermentation_product_binary <- fread(fermentation_product_fn, header = TRUE, sep = ",", skip=2)
# true fermentation product
cpd_colName <- fermentation_product_binary[, grepl("cpd", colnames(fermentation_product_binary)), with = FALSE]
uhgg_id_4_fermentation_test <- fermentation_product_binary$UHGG_id
cpd_colName$uhgg_id <- uhgg_id_4_fermentation_test
cpd_colName_heatmap <- cpd_colName %>%
  melt(id.vars = "uhgg_id", variable.name = "cpd", value.name = "fermProd")

# Read metadata MAG
completeness.level.table          <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/metadata/genomes-all_metadata.csv"
metadata <- read.table(completeness.level.table, sep = ",", header = TRUE, row.names = "Genome")
columns_to_convert_numeric <- c("Completeness")
metadata[columns_to_convert_numeric] <- lapply(metadata[columns_to_convert_numeric], as.numeric)
metadata$Species_name <- rownames(metadata)
metadata <- parse_GTDBtaxonomy(metadata, "Lineage") 
metadata <- refine_GTDBtaxonomy(metadata)
head(metadata)

# EXISTING PAN DRAFT
suf <- ".RDS"
exiting_pan_list <- list()
for (fn in panMAG.CompLV.gapfill.mods.fn) {
    filename = basename(fn)
    filename_without_suffix <- sub(suf, "", filename)
    complLv <- strsplit(filename_without_suffix, "_")[[1]][5]
    uhgg_id <- strsplit(filename_without_suffix, "_")[[1]][1]
    list_id <- paste0(uhgg_id,"_",complLv)
    exiting_pan_list[[list_id]] <- c(uhgg_id, complLv)
}
exiting_pan_df <- data.frame(
  uhgg_id = character(),
  compLv = character(),
  stringsAsFactors = FALSE
)
# Fill the data frame with values from exiting_pan_list
for (list_id in names(exiting_pan_list)) {
  values <- unlist(exiting_pan_list[[list_id]])
  exiting_pan_df[nrow(exiting_pan_df) + 1, ] <- values
}
exiting_pan_df$compLv <- as.numeric(exiting_pan_df$compLv)
exiting_pan_dt <- as.data.table(exiting_pan_df)
exiting_pan_dt$uhgg_compLv_id <- paste0(exiting_pan_dt$uhgg_id, "_", exiting_pan_dt$compLv) 
head(exiting_pan_df)
dim(exiting_pan_df)

# EXIST MAG
exiting_mag_list <- list()
for (fn in singleGeno_pathList) {
    mod_pathList <- readLines(fn)
    ref_filename = basename(fn)
    filename_without_suffix <- sub(".txt", "", ref_filename)
    ref_id <- strsplit(filename_without_suffix, "_")[[1]][3]
    for (mod_fn in mod_pathList){
        mag_filename = basename(mod_fn)
        mag_id <- sub(".RDS", "", mag_filename)
        exiting_mag_list[[mag_id]] <- c(mag_id, ref_id)
    }
}
exiting_mag_df <- data.frame(
  uhgg_id = character(),
  ref_id = character(),
  stringsAsFactors = FALSE
)
# Fill the data frame with values from exiting_pan_list
for (list_id in names(exiting_mag_list)) {
  values <- unlist(exiting_mag_list[[list_id]])
  exiting_mag_df[nrow(exiting_mag_df) + 1, ] <- values
}
head(exiting_mag_df)
dim(exiting_mag_df)

# PARAMETERS: percentage for TP/TN/... across species on percentage of prediction per single MAG
# predictions are computed per species based on the predictions computed for subset of single MAGs

# extract from the exchange name the cpd identifier
fermprod_pan_stat <- gs.fermprod_pan %>%
    mutate(cpd = sapply(strsplit(ex, "_"), "[", 2))
colnames(fermprod_pan_stat)[1] <- "pan_id"
colnames(fermprod_pan_stat)[7] <- "uhgg_id"
# filter: only exchange reaction of cpd of interest are keeped
# EX_rxn that were not predicted have been added to the table 
compLv_df <- data.frame(compLv = c(50, 60, 70, 80, 90), dummy = 1)
iteration_df <- data.frame(iter_num = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), dummy = 1)
cpd_colName_heatmap$dummy = 1 
cartesian_cpd_colName_heatmap <- merge(cpd_colName_heatmap, compLv_df, by = "dummy", allow.cartesian=TRUE)
cartesian_cpd_colName_heatmap <- merge(cartesian_cpd_colName_heatmap, iteration_df, by = "dummy", allow.cartesian=TRUE)
cartesian_cpd_colName_heatmap <- cartesian_cpd_colName_heatmap[,-1]
fermprod_pan_stat_merg <- merge(fermprod_pan_stat, cartesian_cpd_colName_heatmap, by = c("uhgg_id", "cpd", "compLv", "iter_num"), all.y = TRUE)
# filter: keep only info for the model that have been generated
fermprod_pan_stat_merg$uhgg_id_compLv <- paste0(fermprod_pan_stat_merg$uhgg_id,"_",fermprod_pan_stat_merg$compLv)
exiting_pan_df$uhgg_id_compLv <- paste0(exiting_pan_df$uhgg_id,"_",exiting_pan_df$compLv)
fermprod_pan_stat_merg <- fermprod_pan_stat_merg[fermprod_pan_stat_merg$uhgg_id_compLv %in% exiting_pan_df$uhgg_id_compLv, ] 

# set to zero the missing values
fermprod_pan_stat_merg$mtf.flux[is.na(fermprod_pan_stat_merg$mtf.flux)] <- 0
fermprod_pan_stat_merg$u[is.na(fermprod_pan_stat_merg$u)] <- 0
fermprod_pan_stat_merg$gr.rate[is.na(fermprod_pan_stat_merg$gr.rate)] <- 1
# set to zero ref. export flux annotated as NA 
fermprod_pan_stat_merg$fermProd[is.na(fermprod_pan_stat_merg$fermProd)] <- 0
# filter: only flux larger than 10^-4 are considered as export
fermprod_pan_stat_merg$u_binary <- fermprod_pan_stat_merg$u/fermprod_pan_stat_merg$gr.rate > 1e-4
fermprod_pan_stat_merg$mtf.flux_binary <- fermprod_pan_stat_merg$mtf.flux/fermprod_pan_stat_merg$gr.rate > 1e-4

# compute: statistics based on ref. fermentation product per EX_rxn
r_tmp <- calculate_performance_array(fermprod_pan_stat_merg$fermProd, fermprod_pan_stat_merg$mtf.flux_binary)
colnames(r_tmp) <- c("mtf.flux_TP", "mtf.flux_TN", "mtf.flux_FP", "mtf.flux_FN")
fermprod_pan_stat_merg <- cbind(fermprod_pan_stat_merg, r_tmp)
r_tmp <- calculate_performance_array(fermprod_pan_stat_merg$fermProd, fermprod_pan_stat_merg$u_binary)
colnames(r_tmp) <- c("u_TP", "u_TN", "u_FP", "u_FN")
fermprod_pan_stat_merg <- cbind(fermprod_pan_stat_merg, r_tmp)
head(fermprod_pan_stat_merg)
dim(fermprod_pan_stat_merg)

# compute: summarize statistics per genome
num_cpd_to_be_test <- length(unique(cpd_colName_heatmap$cpd)) # number of tested ref. cpd

# compute perc
fermprod_pan_stat_merg_pl_perc <- fermprod_pan_stat_merg %>%
    group_by(uhgg_id, compLv) %>%
    reframe(
        mtf.flux_tot_TP = sum(mtf.flux_TP),
        mtf.flux_tot_TN = sum(mtf.flux_TN),
        mtf.flux_tot_FP = sum(mtf.flux_FP),
        mtf.flux_tot_FN = sum(mtf.flux_FN),
        mtf.flux_perc_TP = sum(mtf.flux_TP)/(10*num_cpd_to_be_test),
        mtf.flux_perc_TN = sum(mtf.flux_TN)/(10*num_cpd_to_be_test),
        mtf.flux_perc_FP = sum(mtf.flux_FP)/(10*num_cpd_to_be_test),
        mtf.flux_perc_FN = sum(mtf.flux_FN)/(10*num_cpd_to_be_test),
        mtf.flux_TPR = mtf.flux_tot_TP / (mtf.flux_tot_TP + mtf.flux_tot_FN),
        mtf.flux_precision = mtf.flux_tot_TP / (mtf.flux_tot_TP + mtf.flux_tot_FP),
        mtf.flux_f1 = (2*mtf.flux_tot_TP) / (2*mtf.flux_tot_TP + mtf.flux_tot_FP + mtf.flux_tot_FN),
        u_tot_TP = sum(u_TP),
        u_tot_TN = sum(u_TN),
        u_tot_FP = sum(u_FP),
        u_tot_FN = sum(u_FN),
        u_perc_TP = sum(u_TP)/(10*num_cpd_to_be_test),
        u_perc_TN = sum(u_TN)/(10*num_cpd_to_be_test),
        u_perc_FP = sum(u_FP)/(10*num_cpd_to_be_test),
        u_perc_FN = sum(u_FN)/(10*num_cpd_to_be_test),
        u_TPR = u_tot_TP / (u_tot_TP + u_tot_FN),
        u_precision = u_tot_TP / (u_tot_TP + u_tot_FP),
        u_f1 = (2*u_tot_TP) / (2*u_tot_TP + u_tot_FP + u_tot_FN)) 
dim(fermprod_pan_stat_merg_pl_perc)
head(data.frame(fermprod_pan_stat_merg_pl_perc))

#################
# MAG
# add ref. genome id
gs.fermprod_wCmLv <- merge(gs.fermprod, metadata[, c("Species_name", "Species_rep")], by.x=".id", by.y="Species_name")
# extract from the exchange name the cpd identifier
fermprod_mag_stat <- gs.fermprod_wCmLv %>%
    mutate(cpd = sapply(strsplit(ex, "_"), "[", 2))
colnames(fermprod_mag_stat)[9] <- "uhgg_id"
# filter: only exchange reaction of cpd of interest are keeped
tested_mag_withatleast1posEX <- unique(fermprod_mag_stat[, c(".id", "uhgg_id")])
tested_mag_withatleast1posEX$dummy <- 1 
cpd_colName_heatmap$dummy <- 1 
cartesian_cpd_colName_heatmap <- merge(cpd_colName_heatmap, tested_mag_withatleast1posEX, by = c("uhgg_id", "dummy"), allow.cartesian=TRUE)
cartesian_cpd_colName_heatmap <- cartesian_cpd_colName_heatmap[,-2]
fermprod_mag_stat_merg <- merge(fermprod_mag_stat, cartesian_cpd_colName_heatmap, by = c("uhgg_id", "cpd", ".id"), all.y = TRUE)
# add: completeness Lv
fermprod_mag_stat_merg <- merge(fermprod_mag_stat_merg, metadata[, c("Species_name", "Completeness")], by.x=".id", by.y="Species_name")
fermprod_mag_stat_merg <- fermprod_mag_stat_merg %>%
                mutate(compLv = ifelse(Completeness < 60, 50,                         # assign completeness level
                ifelse(Completeness >= 60 & Completeness < 70, 60,
                ifelse(Completeness >= 70 & Completeness < 80, 70,
                ifelse(Completeness >= 80 & Completeness < 90, 80,
                ifelse(Completeness >= 90 & Completeness <= 100, 90, NA))))))
# filter: keep only info for the pan-model that have been generated
fermprod_mag_stat_merg$uhgg_id_compLv <- paste0(fermprod_mag_stat_merg$uhgg_id,"_",fermprod_mag_stat_merg$compLv)
exiting_pan_df$uhgg_id_compLv <- paste0(exiting_pan_df$uhgg_id,"_",exiting_pan_df$compLv)
fermprod_mag_stat_merg <- fermprod_mag_stat_merg[fermprod_mag_stat_merg$uhgg_id_compLv %in% exiting_pan_df$uhgg_id_compLv, ] 
# set to zero the missing values
fermprod_mag_stat_merg$mtf.flux[is.na(fermprod_mag_stat_merg$mtf.flux)] <- 0
fermprod_mag_stat_merg$u[is.na(fermprod_mag_stat_merg$u)] <- 0
fermprod_mag_stat_merg$gr.rate[is.na(fermprod_mag_stat_merg$gr.rate)] <- 1
# set to zero ref. export flux annotated as NA 
fermprod_mag_stat_merg$fermProd[is.na(fermprod_mag_stat_merg$fermProd)] <- 0
# filter: only flux larger than 10^-4 are considered as export
fermprod_mag_stat_merg$u_binary <- fermprod_mag_stat_merg$u/fermprod_mag_stat_merg$gr.rate > 1e-4
fermprod_mag_stat_merg$mtf.flux_binary <- fermprod_mag_stat_merg$mtf.flux/fermprod_mag_stat_merg$gr.rate > 1e-4

# compute: statistics based on ref. fermentation product per EX_rxn
r_tmp <- calculate_performance_array(fermprod_mag_stat_merg$fermProd, fermprod_mag_stat_merg$mtf.flux_binary)
colnames(r_tmp) <- c("mtf.flux_TP", "mtf.flux_TN", "mtf.flux_FP", "mtf.flux_FN")
fermprod_mag_stat_merg <- cbind(fermprod_mag_stat_merg, r_tmp)
r_tmp <- calculate_performance_array(fermprod_mag_stat_merg$fermProd, fermprod_mag_stat_merg$u_binary)
colnames(r_tmp) <- c("u_TP", "u_TN", "u_FP", "u_FN")
fermprod_mag_stat_merg <- cbind(fermprod_mag_stat_merg, r_tmp)
head(fermprod_mag_stat_merg)
dim(fermprod_mag_stat_merg)

# compute: summarize statistics per genome
num_cpd_to_be_test <- length(unique(cpd_colName_heatmap$cpd)) # number of tested ref. cpd
# compute perc
fermprod_mag_stat_merg_pl_perc <- fermprod_mag_stat_merg %>%
    group_by(uhgg_id, compLv) %>%
    reframe(
        mtf.flux_tot_TP = sum(mtf.flux_TP),
        mtf.flux_tot_TN = sum(mtf.flux_TN),
        mtf.flux_tot_FP = sum(mtf.flux_FP),
        mtf.flux_tot_FN = sum(mtf.flux_FN),
        mtf.flux_perc_TP = sum(mtf.flux_TP)/(length(unique(.id))*num_cpd_to_be_test),
        mtf.flux_perc_TN = sum(mtf.flux_TN)/(length(unique(.id))*num_cpd_to_be_test),
        mtf.flux_perc_FP = sum(mtf.flux_FP)/(length(unique(.id))*num_cpd_to_be_test),
        mtf.flux_perc_FN = sum(mtf.flux_FN)/(length(unique(.id))*num_cpd_to_be_test),
        mtf.flux_TPR = mtf.flux_tot_TP / (mtf.flux_tot_TP + mtf.flux_tot_FN),
        mtf.flux_precision = mtf.flux_tot_TP / (mtf.flux_tot_TP + mtf.flux_tot_FP),
        mtf.flux_f1 = (2*mtf.flux_tot_TP) / (2*mtf.flux_tot_TP + mtf.flux_tot_FP + mtf.flux_tot_FN),
        u_tot_TP = sum(u_TP),
        u_tot_TN = sum(u_TN),
        u_tot_FP = sum(u_FP),
        u_tot_FN = sum(u_FN),
        u_perc_TP = sum(u_TP)/(length(unique(.id))*num_cpd_to_be_test),
        u_perc_TN = sum(u_TN)/(length(unique(.id))*num_cpd_to_be_test),
        u_perc_FP = sum(u_FP)/(length(unique(.id))*num_cpd_to_be_test),
        u_perc_FN = sum(u_FN)/(length(unique(.id))*num_cpd_to_be_test),
        u_TPR = u_tot_TP / (u_tot_TP + u_tot_FN),
        u_precision = u_tot_TP / (u_tot_TP + u_tot_FP),
        u_f1 = (2*u_tot_TP) / (2*u_tot_TP + u_tot_FP + u_tot_FN)) 
dim(fermprod_mag_stat_merg_pl_perc)
head(data.frame(fermprod_mag_stat_merg_pl_perc))

# PAN & MAG plot
# PLOT ON PERCENTAGE
fermprod_mag_stat_merg_pl_perc$dummy_grup <- "mag"
fermprod_pan_stat_merg_pl_perc$dummy_grup <- "pan"
mag_pan_stat_perc <- rbind(fermprod_mag_stat_merg_pl_perc, fermprod_pan_stat_merg_pl_perc)
head(data.frame(mag_pan_stat_perc))
dim(mag_pan_stat_perc)

mag_pan_stat_perc_pl <- mag_pan_stat_perc[,c("uhgg_id", "compLv",
                    "mtf.flux_perc_TP", "mtf.flux_perc_TN", 
                    "mtf.flux_perc_FP", "mtf.flux_perc_FN", 
                    "u_perc_TP", "u_perc_TN",
                    "u_perc_FP", "u_perc_FN", "dummy_grup")]

melted_data <- gather(mag_pan_stat_perc_pl, key = "variable", value = "value", -uhgg_id, -compLv, -dummy_grup)
melted_data$stat <- sapply(strsplit(melted_data$variable, "_"), "[", 3)
melted_data$predtype_modtype <- sapply(strsplit(melted_data$variable, "_"), "[", 1)
melted_data <- subset(melted_data, select = -variable)
expanded_data <- spread(melted_data, key = stat, value = value)
expanded_data[expanded_data$dummy_grup == "mag", "compLv"] <- expanded_data[expanded_data$dummy_grup == "mag", "compLv"] + 3
expanded_data$dummy_grup_predtype_modtype <- paste0(expanded_data$dummy_grup,"_",expanded_data$predtype_modtype)


pl_TP <- ggplot(expanded_data, aes(x = compLv)) +
    geom_point(aes(y = TP, color = "gray", shape = predtype_modtype, size = 1.)) +
    stat_summary(aes(y = TP, group = dummy_grup_predtype_modtype, color = dummy_grup, linetype = dummy_grup), 
        fun = "mean", geom = "line") +
    stat_summary(aes(y = TP, group = dummy_grup_predtype_modtype, color = dummy_grup), 
        fun.data = "mean_sd", geom = "errorbar", # errorbar crossbar
        width = 3, position=position_dodge(2.)) +
    stat_summary(aes(y = TP, group = dummy_grup_predtype_modtype, color = dummy_grup, shape = predtype_modtype, size = 2.), 
        fun = "mean", geom = "point", position=position_dodge(2.)) +
    # stat_summary(aes(y = TP, group = predtype_modtype), fun.data = "mean_se", geom = "errorbar", width = 0.2) + 
    scale_shape_manual(values = c("u" = 1, "mtf.flux" = 2)) +
    scale_color_manual(values = c("mag" = viridis::plasma(5)[1], "pan" = viridis::plasma(5)[3])) +
    scale_linetype_manual(values = c("mag" = 1, "pan" = 2)) +
    scale_size(range = c(1, 2)) +
    labs(x = "Completeness level", y = "Percentage on total prediction", color = "Prediction", shape = "Model type") +
    # ylim(-0.03, 1.01) +
    theme_minimal() 

pl_TN <- ggplot(expanded_data, aes(x = compLv)) +
    geom_point(aes(y = TN, color = "gray", shape = predtype_modtype, size = 1.)) +
    stat_summary(aes(y = TN, group = dummy_grup_predtype_modtype, color = dummy_grup, linetype = dummy_grup), 
        fun = "mean", geom = "line") +
    stat_summary(aes(y = TN, group = dummy_grup_predtype_modtype, color = dummy_grup), 
        fun.data = "mean_sd", geom = "errorbar", # errorbar crossbar
        width = 3, position=position_dodge(2.)) +
    stat_summary(aes(y = TN, group = dummy_grup_predtype_modtype, color = dummy_grup, shape = predtype_modtype, size = 2.), 
        fun = "mean", geom = "point", position=position_dodge(2.)) +
    # stat_summary(aes(y = TN, group = predtype_modtype), fun.data = "mean_se", geom = "errorbar", width = 0.2) + 
    scale_shape_manual(values = c("u" = 1, "mtf.flux" = 2)) +
    scale_color_manual(values = c("mag" = viridis::plasma(5)[1], "pan" = viridis::plasma(5)[3])) +
    scale_linetype_manual(values = c("mag" = 1, "pan" = 2)) +
    scale_size(range = c(1, 2)) +
    labs(x = "Completeness level", y = "Percentage on total prediction", color = "Prediction", shape = "Model type") +
    # ylim(-0.03, 1.01) +
    theme_minimal() 

pl_FP <- ggplot(expanded_data, aes(x = compLv)) +
    geom_point(aes(y = FP, color = "gray", shape = predtype_modtype, size = 1.)) +
    stat_summary(aes(y = FP, group = dummy_grup_predtype_modtype, color = dummy_grup, linetype = dummy_grup), 
        fun = "mean", geom = "line") +
    stat_summary(aes(y = FP, group = dummy_grup_predtype_modtype, color = dummy_grup), 
        fun.data = "mean_sd", geom = "errorbar", # errorbar crossbar
        width = 3, position=position_dodge(2.)) +
    stat_summary(aes(y = FP, group = dummy_grup_predtype_modtype, color = dummy_grup, shape = predtype_modtype, size = 2.), 
        fun = "mean", geom = "point", position=position_dodge(2.)) +
    # stat_summary(aes(y = FP, group = predtype_modtype), fun.data = "mean_se", geom = "errorbar", width = 0.2) + 
    scale_shape_manual(values = c("u" = 1, "mtf.flux" = 2)) +
    scale_color_manual(values = c("mag" = viridis::plasma(5)[1], "pan" = viridis::plasma(5)[3])) +
    scale_linetype_manual(values = c("mag" = 1, "pan" = 2)) +
    scale_size(range = c(1, 2)) +
    labs(x = "Completeness level", y = "Percentage on total prediction", color = "Prediction", shape = "Model type") +
    # ylim(-0.03, 1.01) +
    theme_minimal() 

pl_FN <- ggplot(expanded_data, aes(x = compLv)) +
    geom_point(aes(y = FN, color = "gray", shape = predtype_modtype, size = 1.)) +
    stat_summary(aes(y = FN, group = dummy_grup_predtype_modtype, color = dummy_grup, linetype = dummy_grup), 
        fun = "mean", geom = "line") +
    stat_summary(aes(y = FN, group = dummy_grup_predtype_modtype, color = dummy_grup), 
        fun.data = "mean_sd", geom = "errorbar", # errorbar crossbar
        width = 3, position=position_dodge(2.)) +
    stat_summary(aes(y = FN, group = dummy_grup_predtype_modtype, color = dummy_grup, shape = predtype_modtype, size = 2.), 
        fun = "mean", geom = "point", position=position_dodge(2.)) +
    # stat_summary(aes(y = FN, group = predtype_modtype), fun.data = "mean_se", geom = "errorbar", width = 0.2) + 
    scale_shape_manual(values = c("u" = 1, "mtf.flux" = 2)) +
    scale_color_manual(values = c("mag" = viridis::plasma(5)[1], "pan" = viridis::plasma(5)[3])) +
    scale_linetype_manual(values = c("mag" = 1, "pan" = 2)) +
    scale_size(range = c(1, 2)) +
    labs(x = "Completeness level", y = "Percentage on total prediction", color = "Prediction", shape = "Model type") +
    # ylim(-0.03, 1.01) +
    theme_minimal() 

combined_plot <- ggarrange(
    pl_TP + labs(title = "a.") + 
        theme(legend.position = "none",
        axis.title.x = element_blank()) +
        labs(y = "Percentage TP on total prediction"),
    pl_TN + labs(title = "b.") +
        theme(legend.position = "none",
        axis.title.x = element_blank()) +
        labs(y = "Percentage TN on total prediction"),
    pl_FP + labs(title = "c.") +
        theme(legend.position = "none") +
        labs(y = "Percentage FP on total prediction"),
    pl_FN + labs(title = "d.") +
        theme(legend.position = "none") +
        labs(y = "Percentage FN on total prediction"),
    ncol = 2, nrow = 2, widths = c(1, 1))

ggsave(file.path(output.dir, "confusion_lineplot_allStat_percOnSpPred.pdf"),
    plot = combined_plot,
    units = "cm", width = 16, height = 16, dpi = 300)

# F1 
mag_pan_sumStat_perc_pl <- mag_pan_stat_perc[, c("uhgg_id", "compLv", "dummy_grup", "mtf.flux_f1", "u_f1")] # "mtf.flux_TPR", "mtf.flux_precision", "u_TPR", "u_precision",
melted_data <- gather(mag_pan_sumStat_perc_pl, key = "variable", value = "value", -uhgg_id, -compLv, -dummy_grup)
melted_data$stat <- sapply(strsplit(melted_data$variable, "_"), "[", 2)
melted_data$predtype_modtype <- sapply(strsplit(melted_data$variable, "_"), "[", 1)
melted_data <- subset(melted_data, select = -variable)
expanded_data <- spread(melted_data, key = stat, value = value)
expanded_data[expanded_data$dummy_grup == "mag", "compLv"] <- expanded_data[expanded_data$dummy_grup == "mag", "compLv"] + 3
expanded_data$dummy_grup_predtype_modtype <- paste0(expanded_data$dummy_grup,"_",expanded_data$predtype_modtype)

pl_f1 <- ggplot(expanded_data, aes(x = compLv)) +
    geom_point(aes(y = f1, color = "gray", shape = predtype_modtype, size = 1.)) +
    stat_summary(aes(y = f1, group = dummy_grup_predtype_modtype, color = dummy_grup, linetype = dummy_grup), 
        fun = "mean", geom = "line") +
    stat_summary(aes(y = f1, group = dummy_grup_predtype_modtype, color = dummy_grup), 
        fun.data = "mean_sd", geom = "errorbar", # errorbar crossbar
        width = 3, position=position_dodge(1.3)) +
    stat_summary(aes(y = f1, group = dummy_grup_predtype_modtype, color = dummy_grup, shape = predtype_modtype, size = 2.), 
        fun = "mean", geom = "point", position=position_dodge(1.3)) +
    # stat_summary(aes(y = f1, group = predtype_modtype), fun.data = "mean_se", geom = "errorbar", width = 0.2) + 
    scale_shape_manual(values = c("u" = 1, "mtf.flux" = 2)) +
    scale_color_manual(values = c("mag" = viridis::plasma(5)[1], "pan" = viridis::plasma(5)[3])) +
    scale_linetype_manual(values = c("mag" = 1, "pan" = 2)) +
    scale_size(range = c(1, 2)) +
    labs(x = "Completeness level", y = "Percentage on total prediction", color = "Prediction", shape = "Model type") +
    # ylim(0, 0.5) +
    guides(size = FALSE) +  
    theme_minimal() 

ggsave(file.path(output.dir, "f1_mtf.flux_percOnSpPred_withLegend.pdf"),
    plot = pl_f1,
    units = "cm", width = 16, height = 10, dpi = 300)

expanded_data$uhgg_id_dummy_grup <- paste0(expanded_data$uhgg_id,"_",expanded_data$dummy_grup)
u <- expanded_data[expanded_data$predtype_modtype == "u",]
pl_f1 <- ggplot(u, aes(x = compLv)) +
    geom_point(aes(y = f1, color = uhgg_id, shape = predtype_modtype, size = 1.)) +
    geom_line(aes(y = f1, color = uhgg_id, group = uhgg_id_dummy_grup, linetype = dummy_grup)) +
    scale_color_viridis_d(option = "plasma") +
    scale_linetype_manual(values = c("mag" = 1, "pan" = 2)) +
    scale_size(range = c(1, 2)) +
    labs(x = "Completeness level", y = "Percentage on total prediction", color = "Prediction", shape = "Model type") +
    guides(size = FALSE) +  
    theme_minimal() 

ggsave(file.path(output.dir, "f1_U_percOnSpPred_withLegend_spLine.pdf"),
    plot = pl_f1,
    units = "cm", width = 16, height = 10, dpi = 300)

combined_plot_f1andStat <- ggarrange(
    combined_plot + 
        theme(legend.position = "none"),
    pl_f1 + labs(title = "e.") +
        labs(y = "F1 score"),
    ncol = 1, nrow = 2, widths = c(1, 1))

ggsave(file.path(output.dir, "confusion_lineplot_allStat_percOnSpPred_andF1_noyAxisLim.pdf"),
    plot = combined_plot_f1andStat,
    units = "cm", width = 16, height = 24, dpi = 300)
    

# STATISTICAL TEST
x <- mag_pan_sumStat_perc_pl[mag_pan_sumStat_perc_pl$dummy_grup == "mag", "mtf.flux_f1"][["mtf.flux_f1"]]
y <- mag_pan_sumStat_perc_pl[mag_pan_sumStat_perc_pl$dummy_grup == "pan", "mtf.flux_f1"][["mtf.flux_f1"]]
result_f1 <- wilcox.test(x, y, alternative = "two.sided")
print(result_f1)

x <- mag_pan_sumStat_perc_pl[mag_pan_sumStat_perc_pl$dummy_grup == "mag", "u_f1"][["u_f1"]]
y <- mag_pan_sumStat_perc_pl[mag_pan_sumStat_perc_pl$dummy_grup == "pan", "u_f1"][["u_f1"]]
result_f1 <- wilcox.test(x, y, alternative = "two.sided")
print(result_f1)


comp_th <- c(50, 60, 70, 80, 90)
res_t.test <- data.frame()
for (th in comp_th){ 
    x <- mag_pan_sumStat_perc_pl[mag_pan_sumStat_perc_pl$dummy_grup == "mag" & mag_pan_sumStat_perc_pl$compLv == th, c("uhgg_id", "mtf.flux_f1")][["mtf.flux_f1"]]
    y <- mag_pan_sumStat_perc_pl[mag_pan_sumStat_perc_pl$dummy_grup == "pan" & mag_pan_sumStat_perc_pl$compLv == th, c("uhgg_id", "mtf.flux_f1")][["mtf.flux_f1"]]
    res_t_mtf.flux <- t.test(x, y, paired = TRUE) # alternative = "less"
    res_t.test <- rbind(res_t.test, data.frame(compLv = th, p_val = res_t_mtf.flux$p.value, predtype = "mtf.flux"))

    x <- mag_pan_sumStat_perc_pl[mag_pan_sumStat_perc_pl$dummy_grup == "mag" & mag_pan_sumStat_perc_pl$compLv == th, "u_f1"][["u_f1"]]
    y <- mag_pan_sumStat_perc_pl[mag_pan_sumStat_perc_pl$dummy_grup == "pan" & mag_pan_sumStat_perc_pl$compLv == th, "u_f1"][["u_f1"]]
    res_t_u <- t.test(x, y, paired = TRUE) # alternative = "less"
    res_t.test <- rbind(res_t.test, data.frame(compLv = th, p_val = res_t_u$p.value, predtype = "u"))
}

res_t.test$sig <- ifelse(res_t.test$p_val<0.05, "*", "")
"""
> res_t.test
   compLv       p_val predtype sig
1      50 0.049090165 mtf.flux   *
2      50 0.059862780        u    
3      60 0.026116386 mtf.flux   *
4      60 0.066046252        u    
5      70 0.002766801 mtf.flux   *
6      70 0.012932514        u   *
7      80 0.008421519 mtf.flux   *
8      80 0.074464329        u    
9      90 0.011402641 mtf.flux   *
10     90 0.054537524        u   
"""

# HEATMAP
fermprod_pan_stat_percxCPD <- fermprod_pan_stat_merg %>%
    group_by(uhgg_id, compLv, cpd) %>%
    reframe(
        mtf.flux_TPxCPD = sum(mtf.flux_TP)/length(uhgg_id_compLv),
        mtf.flux_TNxCPD = sum(mtf.flux_TN)/length(uhgg_id_compLv),
        mtf.flux_FPxCPD = sum(mtf.flux_FP)/length(uhgg_id_compLv),
        mtf.flux_FNxCPD = sum(mtf.flux_FN)/length(uhgg_id_compLv),
        u_TPxCPD = sum(u_TP)/length(uhgg_id_compLv),
        u_TNxCPD = sum(u_TN)/length(uhgg_id_compLv),
        u_FPxCPD = sum(u_FP)/length(uhgg_id_compLv),
        u_FNxCPD = sum(u_FN)/length(uhgg_id_compLv))

fermprod_mag_stat_percxCPD <- fermprod_mag_stat_merg %>%
    group_by(uhgg_id, compLv, cpd) %>%
    reframe(
        mtf.flux_TPxCPD = sum(mtf.flux_TP)/length(uhgg_id_compLv),
        mtf.flux_TNxCPD = sum(mtf.flux_TN)/length(uhgg_id_compLv),
        mtf.flux_FPxCPD = sum(mtf.flux_FP)/length(uhgg_id_compLv),
        mtf.flux_FNxCPD = sum(mtf.flux_FN)/length(uhgg_id_compLv),
        u_TPxCPD = sum(u_TP)/length(uhgg_id_compLv),
        u_TNxCPD = sum(u_TN)/length(uhgg_id_compLv),
        u_FPxCPD = sum(u_FP)/length(uhgg_id_compLv),
        u_FNxCPD = sum(u_FN)/length(uhgg_id_compLv))

fermprod_mag_stat_percxCPD <- merge(fermprod_mag_stat_percxCPD, unique(fermprod_mag_stat_merg[,c("fermProd", "uhgg_id", "cpd")]), by = c("uhgg_id", "cpd"))
fermprod_mag_stat_percxCPD$mod_type <- "mag"
dim(fermprod_mag_stat_percxCPD)
fermprod_pan_stat_percxCPD <- merge(fermprod_pan_stat_percxCPD, unique(fermprod_pan_stat_merg[,c("fermProd", "uhgg_id", "cpd")]), by = c("uhgg_id", "cpd"))
fermprod_pan_stat_percxCPD$mod_type <- "pan"
dim(fermprod_pan_stat_percxCPD)

# add annotations: taxa and cpd name
merged_plot <- rbind(fermprod_pan_stat_percxCPD, fermprod_mag_stat_percxCPD)
dim(merged_plot)
mapping_cpd_df <- data.frame(
                cpd = c("cpd00211", "cpd00029", "cpd11640", "cpd00363", "cpd00141", "cpd00221", "cpd00159", "cpd00047", "cpd00036"), 
                compound = c("Butyrate", "Acetate", "H2", "Ethanol", "Propionate", "D-Lactate", "L-Lactate", "Formate", "Succinate"))
merged_plot <- merge(merged_plot, mapping_cpd_df, by = "cpd")
metadata_uhggname <- metadata[metadata$Species_name %in% merged_plot$uhgg_id, c("Species_name", "species")]
merged_plot <- merge(merged_plot, metadata_uhggname, by.x = "uhgg_id", by.y = "Species_name")

# heatmap
for (compth in c(50, 60, 70, 80, 90)){
    for (what_stat in list(c("u_TPxCPD", "u_FPxCPD"), c("mtf.flux_TPxCPD", "mtf.flux_FPxCPD"))){
        for (mt in c("mag", "pan")){

            what_stat_TP <- what_stat[1]
            what_stat_FP <- what_stat[2]
            print(what_stat_FP)
            print(what_stat_TP)

            pl_data <- merged_plot[merged_plot$mod_type == mt & merged_plot$compLv == compth, ]

            pl1 <-  ggplot(pl_data, 
                    aes(y = factor(species),  x = factor(compound))) +
                    geom_tile(aes(fill = factor(fermProd))) +
                    scale_fill_manual(values = c("white", "#386ab361"), 
                        na.value = "#bbbbbb",
                        name = "prod") +
                    labs(x = "Cpd", y = "UHGG id")
            
            pl2 <- pl1 + geom_point(aes(size = get(what_stat_TP)))  +
                    geom_point(aes(size = get(what_stat_FP))) +
                    scale_size(range = c(-1, 3))+
                    theme_bw() +
                    # theme_minimal() +
                    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                    labs(size = "perc") 
            
            ggsave(file.path(output.dir, paste0(what_stat,"_",compth,"compLv_",mt,"_fermProd.pdf")), 
                pl2,  units = "cm", width = 10, height = 10, dpi = 300)
        }
    }
}

















# OTHER
####################################
# boxplot
stat <- "mean_mtf.flux_TPR" # mean_mtf.flux_precision mean_mtf.flux_f1 mean_mtf.flux_TPR
pl_pan <- ggplot(fermprod_pan_stat_merg_pl_summary, aes(x = factor(compLv), y = get(stat))) +
  geom_boxplot() +
  labs(x = "Completeness level", y = "stat") +
  ylim(0, 1) +
  theme_minimal()
ggsave(file.path(output.dir, paste0("statxSp_",stat,"_PAN_box.pdf")),
    pl_pan,  units = "cm", width = 8, height = 8, dpi = 300)

# lineplot
stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data = fun, colour = "red", geom = geom, width = 0.2, ...)
}


# boxplot
stat <- "mean_u_f1" # mean_mtf.flux_precision mean_mtf.flux_f1 mean_mtf.flux_TPR mean_u_TPR mean_u_f1 mean_u_precision
pl_pan <- ggplot(fermprod_mag_stat_merg_pl_summary, aes(x = factor(compLv), y = get(stat))) +
  geom_boxplot() +
  labs(x = "Completeness level", y = "stat") +
  ylim(0, 1) +
  theme_minimal()
ggsave(file.path(output.dir, paste0("statxSp_",stat,"_MAG_box.pdf")),
    pl_pan,  units = "cm", width = 8, height = 8, dpi = 300)

# lineplot
stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data = fun, colour = "red", geom = geom, width = 0.2, ...)
}

#################################
# mtf.flux_tot_FN can be predicted having stat for the other three parameters. 
# Reason: rows for non prediced exchanges are missing in the data.frame 
fermprod_mag_stat_merg_pl <- fermprod_mag_stat_merg %>%
    group_by(uhgg_id, compLv, .id) %>%
    reframe(
        mtf.flux_tot_TP = sum(mtf.flux_TP),
        mtf.flux_tot_TN = sum(mtf.flux_TN),
        mtf.flux_tot_FP = sum(mtf.flux_FP),
        mtf.flux_tot_FN = sum(mtf.flux_FN),
        mtf.flux_TPR = mtf.flux_tot_TP / (mtf.flux_tot_TP + mtf.flux_tot_FN),
        mtf.flux_precision = mtf.flux_tot_TP / (mtf.flux_tot_TP + mtf.flux_tot_FP),
        mtf.flux_f1 = (2*mtf.flux_tot_TP) / (2*mtf.flux_tot_TP + mtf.flux_tot_FP + mtf.flux_tot_FN),
        u_tot_TP = sum(u_TP),
        u_tot_TN = sum(u_TN),
        u_tot_FP = sum(u_FP),
        u_tot_FN = sum(u_FN),
        u_TPR = u_tot_TP / (u_tot_TP + u_tot_FN),
        u_precision = u_tot_TP / (u_tot_TP + u_tot_FP),
        u_f1 = (2*u_tot_TP) / (2*u_tot_TP + u_tot_FP + u_tot_FN)) 
dim(fermprod_mag_stat_merg_pl)
# some genome didn't have any export reaction for the ref. ferm. products
table(fermprod_mag_stat_merg_pl[, c("uhgg_id", "compLv")])
# check total number of prediction per genome (should be 7)
# rowSums(fermprod_mag_stat_merg_pl[,c("mtf.flux_tot_TP", "mtf.flux_tot_TN", "mtf.flux_tot_FP", "mtf.flux_tot_FN")])
# pr <- merge(exiting_mag_df[exiting_mag_df$ref_id == "MGYG000000074",], metadata[, c("Species_name", "Completeness")], by.x="uhgg_id", by.y="Species_name")
# pr <- merge(unique(gs.fermprod_wCmLv[gs.fermprod_wCmLv$Species_rep == "MGYG000000074", c(".id", "Species_rep")]),
#     metadata[, c("Species_name", "Completeness")], by.x=".id", by.y="Species_name")
# dim(pr[pr$Completeness<60,])

fermprod_mag_stat_merg_pl_summary <- fermprod_mag_stat_merg_pl %>%
    group_by(uhgg_id, compLv) %>%
    reframe(
        mean_mtf.flux_TPR = mean(mtf.flux_TPR),
        mean_mtf.flux_precision = mean(mtf.flux_precision),
        mean_mtf.flux_f1 = mean(mtf.flux_f1),
        std_mtf.flux_TPR = sd(mtf.flux_TPR),
        std_mtf.flux_precision = sd(mtf.flux_precision),
        std_mtf.flux_f1 = sd(mtf.flux_f1),
        mean_u_TPR = mean(u_TPR),
        mean_u_precision = mean(u_precision),
        mean_u_f1 = mean(u_f1),
        std_u_TPR = sd(u_TPR),
        std_u_precision = sd(u_precision),
        std_u_f1 = sd(u_f1))
fermprod_mag_stat_merg_pl_summary$dummy_grup <- "mag"
# some have any FP or TP so the precision cannot be computed
head(fermprod_mag_stat_merg_pl_summary)

####################################
# mtf.flux_tot_FN can be predicted having stat for the other three parameters. 
# Reason: rows for non prediced exchanges are missing in the data.frame 
fermprod_pan_stat_merg_pl <- fermprod_pan_stat_merg %>%
    group_by(uhgg_id, compLv, iter_num) %>%
    reframe(
        mtf.flux_tot_TP = sum(mtf.flux_TP),
        mtf.flux_tot_TN = sum(mtf.flux_TN),
        mtf.flux_tot_FP = sum(mtf.flux_FP),
        mtf.flux_tot_FN = sum(mtf.flux_FN),
        mtf.flux_TPR = mtf.flux_tot_TP / (mtf.flux_tot_TP + mtf.flux_tot_FN),
        mtf.flux_precision = mtf.flux_tot_TP / (mtf.flux_tot_TP + mtf.flux_tot_FP),
        mtf.flux_f1 = (2*mtf.flux_tot_TP) / (2*mtf.flux_tot_TP + mtf.flux_tot_FP + mtf.flux_tot_FN),
        u_tot_TP = sum(u_TP),
        u_tot_TN = sum(u_TN),
        u_tot_FP = sum(u_FP),
        u_tot_FN = sum(u_FN),
        u_TPR = u_tot_TP / (u_tot_TP + u_tot_FN),
        u_precision = u_tot_TP / (u_tot_TP + u_tot_FP),
        u_f1 = (2*u_tot_TP) / (2*u_tot_TP + u_tot_FP + u_tot_FN)) 
dim(fermprod_pan_stat_merg_pl)
# two pan-draft (MGYG000001300) do not have any export reaction for the ref. ferm. products but were included anyway
table(fermprod_pan_stat_merg_pl[, c("uhgg_id", "compLv")])
# rowSums(fermprod_pan_stat_merg_pl[,c("mtf.flux_tot_TP", "mtf.flux_tot_TN", "mtf.flux_tot_FP", "mtf.flux_tot_FN")])

# average the performance through the species to give the same weight (irrelevant for pangenome because they are 10)
fermprod_pan_stat_merg_pl_summary <- fermprod_pan_stat_merg_pl %>%
    group_by(uhgg_id, compLv) %>%
    reframe(
        mean_mtf.flux_TPR = mean(mtf.flux_TPR),
        mean_mtf.flux_precision = mean(mtf.flux_precision),
        mean_mtf.flux_f1 = mean(mtf.flux_f1),
        std_mtf.flux_TPR = sd(mtf.flux_TPR),
        std_mtf.flux_precision = sd(mtf.flux_precision),
        std_mtf.flux_f1 = sd(mtf.flux_f1),
        mean_u_TPR = mean(u_TPR),
        mean_u_precision = mean(u_precision),
        mean_u_f1 = mean(u_f1),
        std_u_TPR = sd(u_TPR),
        std_u_precision = sd(u_precision),
        std_u_f1 = sd(u_f1))
fermprod_pan_stat_merg_pl_summary$dummy_grup <- "pan"
head(fermprod_pan_stat_merg_pl_summary)

################################################
# PAN & MAG plot
colnames(fermprod_mag_stat_merg_pl_summary)
colnames(fermprod_pan_stat_merg_pl_summary)

mag_pan_stat <- rbind(fermprod_mag_stat_merg_pl_summary, fermprod_pan_stat_merg_pl_summary)
head(data.frame(mag_pan_stat))

f1_big <- data.frame(mag_pan_stat[mag_pan_stat$mean_u_f1 > 0.9,])
f1_big_sp <- merge(f1_big, metadata[metadata$Species_name == metadata$Species_rep, c("Species_rep", "species")], by.x="uhgg_id", by.y="Species_rep")
f1_big_sp[,c("uhgg_id", "compLv","species", "mean_mtf.flux_f1" , "std_mtf.flux_f1", "mean_u_f1", "std_u_f1", "dummy_grup")]
colnames(f1_big_sp)

stat <- "mean_mtf.flux_f1" # mean_u_f1
# mean_mtf.flux_precision mean_mtf.flux_f1 mean_mtf.flux_TPR mean_u_TPR mean_u_f1 mean_u_precision
pl_pan <- ggplot(mag_pan_stat, aes(x = factor(compLv), y = get(stat), color = dummy_grup)) + 
        geom_point() +
        stat_summary(aes(y = get(stat), group = dummy_grup), fun = "mean", geom = "line") +
        stat_summary(fun.data = "mean_sd", geom = "crossbar", width = 0.2, mapping = aes(group = dummy_grup)) +
        # stat_sum_df("mean_sd", mapping = aes(group = compLv)) + # fun.args = list(mult = 1)
        labs(x = "Completeness level", y = "stat") +
        ylim(0, 1) +
        scale_color_manual(values = c("mag" = "red", "pan" = "blue")) + 
        theme_minimal()

ggsave(file.path(output.dir, paste0("statxSp_",stat,"_MAGandPAN_line.pdf")),
    pl_pan,  units = "cm", width = 8, height = 8, dpi = 300)

# STATISTICAL TEST
x <- fermprod_mag_stat_merg_pl_summary[["mean_u_f1"]]
y <- fermprod_pan_stat_merg_pl_summary[["mean_u_f1"]]
result_mean_u_f1 <- wilcox.test(x, y, alternative = "two.sided")
print(result_mean_u_f1)

x <- fermprod_mag_stat_merg_pl_summary[["mean_mtf.flux_f1"]]
y <- fermprod_pan_stat_merg_pl_summary[["mean_mtf.flux_f1"]]
result_mean_mtf.flux_f1 <- wilcox.test(x, y, alternative = "two.sided")
print(result_mean_mtf.flux_f1)

#####################################
# PARAMETERS: percentage for TP/TN/... across species on percentage of prediction per species
# predictions are computed per species based on the predictions computed for all the MAGs
results <- data.frame()
for (what_stat in c("u", "mtf.flux")){
    for (compth in c(50, 60, 70, 80, 90)){
        flux_in_more_than <- .0

        # parse PAN results
        gs.fermprod_pan_heatmap <- gs.fermprod_pan %>%
        group_by(ref_uhgg_id, compLv, ex) %>%
        summarise(perc_above_zero = sum(get(what_stat)/gr.rate > 0) / 10) %>% #, na.rm = TRUE
        mutate(cpd = sapply(strsplit(ex, "_"), "[", 2))
        colnames(gs.fermprod_pan_heatmap)[1] <- "uhgg_id"
        dim(gs.fermprod_pan_heatmap)

        # parse MAG results
        # link Species_rep to the Species_name
        gs.fermprod_heatmap <- merge(gs.fermprod, metadata[, c("Species_name", "Species_rep", "Completeness")], by.x=".id", by.y="Species_name")
        gs.fermprod_heatmap <- gs.fermprod_heatmap %>%
        mutate(compLv = ifelse(Completeness < 60, 50,                         # assign completeness level
                        ifelse(Completeness >= 60 & Completeness < 70, 60,
                        ifelse(Completeness >= 70 & Completeness < 80, 70,
                        ifelse(Completeness >= 80 & Completeness < 90, 80,
                        ifelse(Completeness >= 90 & Completeness <= 100, 90, NA)))))) %>%
        group_by(Species_rep, compLv) %>%
        mutate(mag_in_subset = length(unique(.id))) %>%                      # number of genomes for that species
        ungroup() %>%
        group_by(Species_rep, compLv, ex) %>%
        reframe(perc_above_zero = sum(get(what_stat)/gr.rate > 0) / mag_in_subset) %>% # na.rm = TRUE
        distinct() %>%
        mutate(cpd = sapply(strsplit(ex, "_"), "[", 2))
        colnames(gs.fermprod_heatmap)[4] <- "perc_above_zero_mag"
        colnames(gs.fermprod_heatmap)[1] <- "uhgg_id"
        dim(gs.fermprod_heatmap)

        # MERGE RESULTS MAG, PAN and original prediction 
        compLv_df <- data.frame(compLv = c(50, 60, 70, 80, 90), dummy = 1)
        cpd_colName_heatmap$dummy = 1 
        merged_data <- merge(cpd_colName_heatmap, compLv_df, by = "dummy", allow.cartesian=TRUE)
        merged_data <- merged_data[, -1]
        merged_data <- merge(merged_data, gs.fermprod_pan_heatmap, by = c("uhgg_id", "cpd", "compLv"), all.x=TRUE)
        merged_data <- merge(merged_data, gs.fermprod_heatmap, by = c("uhgg_id", "cpd", "compLv"), all.x=TRUE)
        dim(merged_data)

        # filter fluxes less frequent than "flux_in_more_than"
        merged_data$perc_above_zero[merged_data$perc_above_zero < flux_in_more_than] <- 0
        merged_data$perc_above_zero_mag[merged_data$perc_above_zero_mag < flux_in_more_than] <- 0

        # SET NA TO THE SP WITHOUT A PAN GENOME
        merged_plot <- merged_data[merged_data$compLv == compth]
        merged_plot$fermProd[is.na(merged_plot$fermProd)] <- 0
        merged_plot$uhgg_compLv_id <- paste0(merged_plot$uhgg_id, "_", merged_plot$compLv) 
        merged_plot$fermProd[!(merged_plot$uhgg_compLv_id %in% exiting_pan_dt$uhgg_compLv_id)] <- NA
        merged_plot$perc_above_zero_mag[!(merged_plot$uhgg_compLv_id %in% exiting_pan_dt$uhgg_compLv_id)] <- NA

        # Drop missing pan models
        merged_plot <- merged_plot[merged_plot$uhgg_compLv_id %in% exiting_pan_dt$uhgg_compLv_id, ]

        mapping_df <- data.frame(
                        cpd = c("cpd00211", "cpd00029", "cpd11640", "cpd00363", "cpd00141", "cpd00221", "cpd00159", "cpd00047"), 
                        compound = c("Butyrate", "Acetate", "H2", "Ethanol", "Propionate", "D-Lactate", "L-Lactate", "Formate")
                    )
        merged_plot <- merge(merged_plot, mapping_df, by = "cpd")
        metadata_uhggname <- metadata[metadata$Species_name %in% merged_plot$uhgg_id, c("Species_name", "species")]
        merged_plot <- merge(merged_plot, metadata_uhggname, by.x = "uhgg_id", by.y = "Species_name")

        pl1 <-  ggplot(merged_plot, 
                aes(y = factor(species),  x = factor(compound))) +
                geom_tile(aes(fill = factor(fermProd))) +
                scale_fill_manual(values = c("white", "#386ab361"), 
                    na.value = "#bbbbbb",
                    name = "prod") +
                labs(x = "Cpd", y = "UHGG id")
        # pan 
        pl2 <- pl1 + geom_point(aes(size = perc_above_zero))  +
                scale_size(range = c(-1, 3))+
                # theme_bw() +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                labs(size = "perc") 
        ggsave(file.path(output.dir, paste0(what_stat,"_",compth,"compLv_pan_fermProd_",flux_in_more_than,"minFreq.pdf")), 
        pl2,  units = "cm", width = 10, height = 10, dpi = 300)
        cat(paste0(what_stat,"_",compth,"compLv_pan_fermProd_",flux_in_more_than,"minFreq.pdf\n"))

        # mag
        pl3 <- pl1 + geom_point(aes(size = perc_above_zero_mag))  +
                # scale_fill_manual("black") +
                scale_size(range = c(-1, 3))+
                theme_bw() +
                # theme_classic() +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                labs(size = "perc") 
        ggsave(file.path(output.dir, paste0(what_stat,"_",compth,"compLv_mag_fermProd_",flux_in_more_than,"minFreq.pdf")), 
        pl3,  units = "cm", width = 10, height = 10, dpi = 300)
        cat(paste0(what_stat,"_",compth,"compLv_mag_fermProd_",flux_in_more_than,"minFreq.pdf\n"))


        # PERFORMANCE
        merged_plot_perf <- merged_plot[(merged_plot$uhgg_compLv_id %in% exiting_pan_dt$uhgg_compLv_id)] # keep only species for which there is a pan-Draft
        merged_plot_perf$perc_above_zero_mag[is.na(merged_plot_perf$perc_above_zero_mag)] <- 0 # set to 0 exchange for which no flux was predicted
        merged_plot_perf$perc_above_zero[is.na(merged_plot_perf$perc_above_zero)] <- 0
        # compute normalized confusion matrix based on frequency of prediction
        norm_performance_mag <- calculate_normalized_confusion_matrix(merged_plot_perf$fermProd, merged_plot_perf$perc_above_zero_mag) 
        norm_performance_pan <- calculate_normalized_confusion_matrix(merged_plot_perf$fermProd, merged_plot_perf$perc_above_zero)
        norm_performance_mag["compLv"] <- compth
        norm_performance_mag["type"] <- "mag"
        norm_performance_mag["tot_pred"] <- length(merged_plot_perf$fermProd)
        norm_performance_mag["sp_num"] <- length(unique(merged_plot_perf$uhgg_id))
        norm_performance_mag["pred_type"] <- what_stat
        norm_performance_pan["compLv"] <- compth
        norm_performance_pan["type"] <- "pan"
        norm_performance_pan["tot_pred"] <- length(merged_plot_perf$fermProd)
        norm_performance_pan["sp_num"] <- length(unique(merged_plot_perf$uhgg_id))
        norm_performance_pan["pred_type"] <- what_stat
        # append results
        results <- rbind(results, t(data.frame(value = norm_performance_mag)))
        results <- rbind(results, t(data.frame(value = norm_performance_pan)))
    }
}

head(results)
results$TP <- as.numeric(results$TP)
results$TN <- as.numeric(results$TN)
results$FP <- as.numeric(results$FP)
results$FN <- as.numeric(results$FN)
results$compLv <- as.numeric(results$compLv)
results$tot_pred <- as.numeric(results$tot_pred)
results$sp_num <- as.numeric(results$sp_num)
results$pred_spNum <- paste0("(p.", results$tot_pred, ", s.", results$sp_num, ")")
results$predtype_modtype <- paste0(results$pred_type, " ", results$type)

# PLOT
pl_TP <- ggplot(results, aes(x = compLv)) +
        geom_line(aes(y = TP, color = "TP", linetype = predtype_modtype)) +
        scale_linetype_manual(values = c("u mag" = 1, "u pan" = 2, "mtf.flux mag" = 1, "mtf.flux pan" = 2)) +
        geom_point(aes(y = TP, color = "TP", shape = predtype_modtype)) +
        geom_text(aes(y = TP, label = pred_spNum), nudge_x = 0.5, nudge_y = -20., size = 2) + # Add labels based on sp_num
        scale_shape_manual(values = c("u mag" = 1, "u pan" = 1, "mtf.flux mag" = 2, "mtf.flux pan" = 2)) +
        labs(x = "Completeness level", y = "Percentage on total prediction", color = "Prediction", shape = "Model type") +
        # scale_color_viridis_d(option = "plasma") +
        scale_color_manual(values = c("TP" = viridis::plasma(5)[1])) +
        ylim(0, 50) +
        theme_minimal() 
        # theme(legend.position = "none")

pl_FP <- ggplot(results, aes(x = compLv)) +
        geom_line(aes(y = FP, color = "FP", linetype = predtype_modtype)) +
        scale_linetype_manual(values = c("u mag" = 1, "u pan" = 2, "mtf.flux mag" = 1, "mtf.flux pan" = 2)) +
        geom_point(aes(y = FP, color = "FP", shape = predtype_modtype)) +
        # geom_text(aes(y = FP, label = pred_spNum), nudge_x = 0.5, nudge_y = -1., size = 2) + # Add labels based on sp_num
        scale_shape_manual(values = c("u mag" = 1, "u pan" = 1, "mtf.flux mag" = 2, "mtf.flux pan" = 2)) +
        labs(x = "Completeness level", y = "Percentage on total prediction", color = "Prediction", shape = "Model type") +
        # scale_color_viridis_d(option = "plasma") +
        scale_color_manual(values = c("FP" = viridis::plasma(5)[2])) +
        ylim(0, 50) +
        theme_minimal() 
        # theme(legend.position = "none")

pl_FN <- ggplot(results, aes(x = compLv)) +
    geom_line(aes(y = FN, color = "FN", linetype = predtype_modtype)) +
    scale_linetype_manual(values = c("u mag" = 1, "u pan" = 2, "mtf.flux mag" = 1, "mtf.flux pan" = 2)) +
    geom_point(aes(y = FN, color = "FN", shape = predtype_modtype)) +
    # geom_text(aes(y = FN, label = pred_spNum), nudge_x = 0.5, nudge_y = -1., size = 2) + # Add labels based on sp_num
    scale_shape_manual(values = c("u mag" = 1, "u pan" = 1, "mtf.flux mag" = 2, "mtf.flux pan" = 2)) +
    labs(x = "Completeness level", y = "Percentage on total prediction", color = "Prediction", shape = "Model type") +
    # scale_color_viridis_d(option = "plasma") +
    scale_color_manual(values = c("FN" = viridis::plasma(5)[3])) +
    ylim(0, 50) +
    theme_minimal() 
    # theme(legend.position = "none")

pl_TN <- ggplot(results, aes(x = compLv)) +
    geom_line(aes(y = TN, color = "TN", linetype = predtype_modtype)) +
    scale_linetype_manual(values = c("u mag" = 1, "u pan" = 2, "mtf.flux mag" = 1, "mtf.flux pan" = 2)) +
    geom_point(aes(y = TN, color = "TN", shape = predtype_modtype)) +
    geom_text(aes(y = TN, label = pred_spNum), nudge_x = 0.5, nudge_y = -20., size = 2) + # Add labels based on sp_num
    scale_shape_manual(values = c("u mag" = 1, "u pan" = 1, "mtf.flux mag" = 2, "mtf.flux pan" = 2)) +
    labs(x = "Completeness level", y = "Percentage on total prediction", color = "Prediction", shape = "Model type") +
    # scale_color_viridis_d(option = "plasma") +
    scale_color_manual(values = c("TN" = viridis::plasma(5)[4])) +
    ylim(0, 50) +
    theme_minimal() 
    # theme(legend.position = "none")
        
combined_plot <- ggarrange(
    pl_TP + labs(title = "a.") + 
        theme(legend.position = "none",
        axis.title.x = element_blank()) +
        labs(y = "Percentage TP on total prediction"),
    pl_TN + labs(title = "b.") +
        theme(legend.position = "none",
        axis.title.x = element_blank()) +
        labs(y = "Percentage TN on total prediction"),
    pl_FP + labs(title = "c.") +
        theme(legend.position = "none") +
        labs(y = "Percentage FP on total prediction"),
    pl_FN + labs(title = "d.") +
        theme(legend.position = "none") +
        labs(y = "Percentage FN on total prediction"),
    ncol = 2, nrow = 2, widths = c(1, 1))

ggsave(file.path(output.dir, "confusion_lineplot_allStat.pdf"),
    plot = combined_plot,
    units = "cm", width = 14, height = 16, dpi = 300)

# ggsave(file.path(output.dir, "confusion_lineplot_TP.pdf"),
#     pl_TP,  units = "cm", width = 8, height = 8, dpi = 300)
# ggsave(file.path(output.dir, "confusion_lineplot_FP.pdf"),
#     pl_FP,  units = "cm", width = 8, height = 8, dpi = 300)
# ggsave(file.path(output.dir, "confusion_lineplot_FN.pdf"),
#     pl_FN,  units = "cm", width = 8, height = 8, dpi = 300)
# ggsave(file.path(output.dir, "confusion_lineplot_TN.pdf"),
#     pl_TN,  units = "cm", width = 8, height = 8, dpi = 300)
