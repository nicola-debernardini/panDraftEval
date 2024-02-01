# HEATMAP
library(dplyr)
library(data.table)
library(ggplot2) # Graphical

# Little helpers
source("/home/bioinfo/projects/panGenome/panDraftEval/scr/function_collection.R")

# HELPER
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

# HELPER
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

# INPUT
output.dir                <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/fermProdTest/out"
gs.fermprod_pan <- read.table(file.path(output.dir, "panDraft_fermProd_31Gen.tsv"), header = TRUE, sep = "\t", quote = "")
gs.fermprod <- read.table(file.path(output.dir, "mag_fermProd_31Gen.tsv"), header = TRUE, sep = "\t", quote = "")
dim(gs.fermprod_pan)
dim(gs.fermprod)

# pan-Draft
panMAG.CompLV.gapfill.model.path  <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/fermProdTest/panFromMAG_GapfilledFT_Path.txt"
panMAG.CompLV.gapfill.mods.fn <- readLines(panMAG.CompLV.gapfill.model.path)

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
dim(exiting_pan_df)


# PARAMETERS
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
library(viridis)
pl <- ggplot(results, aes(x = compLv)) +
        geom_line(aes(y = TP, color = "TP", linetype = predtype_modtype)) +
        geom_line(aes(y = TN, color = "TN", linetype = predtype_modtype)) +
        geom_line(aes(y = FN, color = "FN", linetype = predtype_modtype)) +
        geom_line(aes(y = FP, color = "FP", linetype = predtype_modtype)) +
        scale_linetype_manual(values = c("u mag" = 1, "u pan" = 2, "mtf.flux mag" = 1, "mtf.flux pan" = 2)) +
        geom_point(aes(y = TP, color = "TP", shape = predtype_modtype)) +
        geom_point(aes(y = TN, color = "TN", shape = predtype_modtype)) +
        geom_point(aes(y = FN, color = "FN", shape = predtype_modtype)) +
        geom_point(aes(y = FP, color = "FP", shape = predtype_modtype)) +
        geom_text(aes(y = TN, label = pred_spNum), nudge_x = 0.5, nudge_y = -1., size = 2) + # Add labels based on sp_num
        scale_shape_manual(values = c("u mag" = 1, "u pan" = 1, "mtf.flux mag" = 2, "mtf.flux pan" = 2)) +
        labs(x = "Completeness level", y = "Percentage on total prediction", color = "Prediction", shape = "Model type") +
        # scale_color_viridis_d(option = "plasma") +
        scale_color_manual(values = c("TP" = viridis::plasma(5)[1], 
                                    "TN" = viridis::plasma(5)[4],
                                    "FN" = viridis::plasma(5)[3],
                                    "FP" = viridis::plasma(5)[2])) +
        theme_minimal()
        
ggsave(file.path(output.dir, "confusion_lineplot.pdf"),
    pl,  units = "cm", width = 16, height = 10, dpi = 300)




    # plot
    # confusion_df$value_perc <- confusion_df$value / sum(confusion_df$value)
    # confusion_plot <- ggplot(confusion_df, aes(x = "Metrics", y = value_perc, fill = metric)) +
    # geom_bar(stat = "identity") +
    # labs(title = "Confusion Matrix (Percentage)",
    #     x = "",
    #     y = "Percentage") +
    # scale_y_continuous(labels = scales::percent_format()) +
    # theme_minimal()
    # ggsave(file.path(output.dir, paste0(what_stat,"_",compth,"_confusion_plot_mag_",flux_in_more_than,"minFreq.pdf")),
    # confusion_plot,  units = "cm", width = 6, height = 10, dpi = 300)

    # confusion_df <- data.frame(metric = names(norm_performance_pan), value = norm_performance_pan)
    # confusion_df$value_perc <- confusion_df$value / sum(confusion_df$value)
    # confusion_plot <- ggplot(confusion_df, aes(x = "Metrics", y = value_perc, fill = metric)) +
    # geom_bar(stat = "identity") +
    # labs(title = "Confusion Matrix (Percentage)",
    #     x = "",
    #     y = "Percentage") +
    # scale_y_continuous(labels = scales::percent_format()) +
    # theme_minimal()
    # ggsave(file.path(output.dir, paste0(what_stat,"_",compth,"_confusion_plot_pan_",flux_in_more_than,"minFreq.pdf")),
    # confusion_plot,  units = "cm", width = 6, height = 10, dpi = 300)


# PERFORMANCE
head(data.frame(merged_data))
merged_data$fermProd[is.na(merged_data$fermProd)] <- 0
merged_data$predmag <- as.numeric(merged_data$perc_above_zero_mag > 0)
merged_data$predpan <- as.numeric(merged_data$perc_above_zero > 0)

confusion_matrix_mag <- calculate_confusion_matrix(merged_data$fermProd, merged_data$predmag)
confusion_matrix_pan <- calculate_confusion_matrix(merged_data$fermProd, merged_data$predpan)

confusion_df <- data.frame(metric = names(confusion_matrix_mag), value = confusion_matrix_mag)
confusion_df$value_perc <- confusion_df$value / sum(confusion_df$value) 
confusion_plot <- ggplot(confusion_df, aes(x = "Metrics", y = value_perc, fill = metric)) +
  geom_bar(stat = "identity") +
  labs(title = "Confusion Matrix (Percentage)",
       x = "",
       y = "Percentage") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal()
ggsave(file.path(output.dir, paste0(flux_in_more_than,"less_",what_stat,compLv,"_confusion_plot_mag.pdf")), 
  confusion_plot,  units = "cm", width = 6, height = 10, dpi = 300)

confusion_df <- data.frame(metric = names(confusion_matrix_pan), value = confusion_matrix_pan)
confusion_df$value_perc <- confusion_df$value / sum(confusion_df$value) 
confusion_plot <- ggplot(confusion_df, aes(x = "Metrics", y = value_perc, fill = metric)) +
  geom_bar(stat = "identity") +
  labs(title = "Confusion Matrix (Percentage)",
       x = "",
       y = "Percentage") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal()
ggsave(file.path(output.dir, paste0(flux_in_more_than,"less_",what_stat,compLv,"_confusion_plot_pan.pdf")), 
  confusion_plot,  units = "cm", width = 6, height = 10, dpi = 300)













# COMPARE RESULTS
library(tidyr)
library(dplyr)
str(t.all_panDraft_pred_res)
str(t.all_pred_res)

numeric_columns <- c("TP", "FP", "TN", "FN", "CompLv_th")
t.all_pred_res[numeric_columns] <- lapply(t.all_pred_res[numeric_columns], as.numeric)
t.all_panDraft_pred_res[numeric_columns] <- lapply(t.all_panDraft_pred_res[numeric_columns], as.numeric)

numeric_columns <- c("TP", "FP", "TN", "FN")
plot_df <- data.frame()
for (uhgg_id in unique(t.all_panDraft_pred_res$ref_Id)) {
  cat(paste(uhgg_id, "\n"))
  CompLv_th <- 80
  if (dim(t.all_panDraft_pred_res[t.all_panDraft_pred_res$ref_Id == uhgg_id & t.all_panDraft_pred_res$compLv == CompLv_th,])[1] > 0) {
    panPheno_betweenTh <- t.all_panDraft_pred_res[t.all_panDraft_pred_res$ref_Id == uhgg_id & t.all_panDraft_pred_res$compLv == CompLv_th, c("TP", "FP", "TN", "FN")]

    MAG_id_betweenTh <-  metadata[metadata$Species_rep == uhgg_id & metadata$Completeness < (CompLv_th+10), "Species_name"]
    length(MAG_id_betweenTh)
    magPheno_betweenTh <- t.all_pred_res[t.all_pred_res$mag_Id %in% MAG_id_betweenTh, c("TP", "FP", "TN", "FN")]
    dim(magPheno_betweenTh)
    # something is missing, WHY ???????

    colnames(panPheno_betweenTh)

    mean_panPheno_betweenTh <- colMeans(panPheno_betweenTh[numeric_columns])
    mean_magPheno_betweenTh <- colMeans(magPheno_betweenTh[numeric_columns])
    mean_panPheno_betweenTh["type"] = "pan"
    mean_magPheno_betweenTh["type"] = "mag"

    mean_Pheno_sp <- as.data.frame(rbind(mean_panPheno_betweenTh, mean_magPheno_betweenTh))
    mean_Pheno_sp[numeric_columns] <- lapply(mean_Pheno_sp[numeric_columns], as.numeric)
    plot_df <- rbind(plot_df, mean_Pheno_sp)
  } else {
    cat(paste(uhgg_id, "do not have a pan-model in the completeness th", CompLv_th, CompLv_th+10, "\n"))
    next
  }
}

num_sp <- dim(plot_df)[1]/2 
cat(paste("The number of species consider for this statistics is", num_sp, "\n"))

plot_df <- gather(plot_df, key = "statistic", value = "value", -type)
summarized_df <- plot_df %>%
  group_by(type, statistic) %>%
  summarize(sum_value = sum(value, na.rm = TRUE))


# Plot
plot <- ggplot(summarized_df, aes(x = statistic, y = sum_value, fill = type)) +
  geom_col(position = "dodge") +  # Use geom_col() instead of geom_bar()
  labs(x = "Type",
      y = "Count") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_minimal()
print(plot)

ggsave(file.path(output.dir, paste0("prova_bar_", CompLv_th, "_totalSum.pdf")), 
  plot,  units = "cm", width = 16, height = 10, dpi = 300)



uhgg_id <- "MGYG000000113"
CompLv_th <- 50
panPheno_betweenTh <- t.all_panDraft_pred_res[t.all_panDraft_pred_res$ref_Id == uhgg_id & t.all_panDraft_pred_res$compLv == as.character(CompLv_th), c("TP", "FP", "TN", "FN")]

MAG_id_betweenTh <-  metadata[metadata$Species_rep == uhgg_id & metadata$Completeness < (CompLv_th+10), "Species_name"]
length(MAG_id_betweenTh)
magPheno_betweenTh <- t.all_pred_res[t.all_pred_res$mag_Id %in% MAG_id_betweenTh, c("TP", "FP", "TN", "FN")]
dim(magPheno_betweenTh)
# something is missing, WHY ???????

magPheno_betweenTh["type"] = "mag"
panPheno_betweenTh["type"] = "pan"
numeric_columns <- c("TP", "FP", "TN", "FN")
plot_df <- rbind(panPheno_betweenTh, magPheno_betweenTh)
plot_df[numeric_columns] <- lapply(plot_df[numeric_columns], as.numeric)
plot_df <- gather(plot_df, key = "statistic", value = "value", -type)
dim(plot_df)

# Plot
plot <- ggplot(plot_df, aes(x = statistic, y = value, fill = type)) +
  geom_boxplot() +
  labs(x = "Type",
       y = "Count") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_minimal()
print(plot)

ggsave(file.path(output.dir, paste0("prova.pdf")), 
  plot,  units = "cm", width = 16, height = 10, dpi = 300)



### CORRELATION
metadata_compLv <- metadata[,c("Species_name", "Completeness")]
confusion_matrix_res_df$Species_name <- rownames(confusion_matrix_res_df)
merged_data <- merge(metadata_compLv, confusion_matrix_res_df, by = "Species_name")

# Scatter plot for Completeness vs. TP
ggplot(merged_data, aes(x = Completeness, y = TP)) +
  geom_point() +
  labs(title = "Completeness vs. True Positives", x = "Completeness", y = "True Positives")

# Scatter plot for Completeness vs. TN
ggplot(merged_data, aes(x = Completeness, y = TN)) +
  geom_point() +
  labs(title = "Completeness vs. True Negatives", x = "Completeness", y = "True Negatives")

# Scatter plot for Completeness vs. FP
ggplot(merged_data, aes(x = Completeness, y = FP)) +
  geom_point() +
  labs(title = "Completeness vs. False Positives", x = "Completeness", y = "False Positives")

# Scatter plot for Completeness vs. FN
ggplot(merged_data, aes(x = Completeness, y = FN)) +
  geom_point() +
  labs(title = "Completeness vs. False Negatives", x = "Completeness", y = "False Negatives")

library(corrplot)
cor_data <- merged_data[, c("Completeness", "TP", "TN", "FP", "FN")]
cor_matrix <- cor(cor_data)
p_values <- matrix(NA, nrow = ncol(cor_matrix), ncol = ncol(cor_matrix))
# Loop through each pair of variables and calculate the correlation and p-value
for (i in 1:(ncol(cor_matrix) - 1)) {
  for (j in (i + 1):ncol(cor_matrix)) {
    cor_test_result <- cor.test(cor_data[, i], cor_data[, j])
    cor_matrix[i, j] <- cor_test_result$estimate
    cor_matrix[j, i] <- cor_test_result$estimate
    p_values[i, j] <- cor_test_result$p.value
    p_values[j, i] <- cor_test_result$p.value
  }
}

print(p_values)
# Plot the correlation matrix using a heatmap
corrplot(cor_matrix, method = "color", type = "upper", addCoef.col = "black", sig.level = 0.05, insig = "blank", p.mat = p_values)


