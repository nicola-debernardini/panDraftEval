#
#
############################################################################

# Load R-packages
library(data.table)
library(ggplot2) # Graphical
library(sybil)

# select solver
if( "cplexAPI" %in% rownames(installed.packages()) ){
  sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1
}else{
  warning("glpkAPI is used but cplexAPI is recommended because it is much faster")
  sybil::SYBIL_SETTINGS("SOLVER","glpkAPI"); ok <- 5
}

# Little helpers
# source("./scr/function_collection.R")
source("/home/bioinfo/projects/panGenome/panDraftEval/scr/function_collection.R")
source("/home/bioinfo/bin/gapseq_panDraft/src/constrain.model.R")

# Function to compute confusion matrix
compute_confusion_matrix <- function(predicted, actual) {
  tp <- sum(predicted == 1 & actual == 1)
  fp <- sum(predicted == 1 & actual == 0)
  tn <- sum(predicted == 0 & actual == 0)
  fn <- sum(predicted == 0 & actual == 1)
  
  return(c(TP = tp, FP = fp, TN = tn, FN = fn))
}

setObjFun_altEnergy_sources <- function(mod) {
  mod@obj_coef <- rep(0, mod@react_num) # zero the objective
  # add biolog like test
  mql  <- "cpd15499[c0]"; mqn   <- "cpd15500[c0]" # menaquinone
  uql  <- "cpd15561[c0]"; uqn   <- "cpd15560[c0]" # ubiquinone
  h    <- "cpd00067[c0]"
  nad  <- "cpd00003[c0]"; nadh  <- "cpd00004[c0]"
  fdox <- "cpd11621[c0]"; fdred <- "cpd11620[c0]" # ferredoxin
  pql  <- "cpd27796[c0]"; pqn   <- "cpd27797[c0]" # plastoquinone
  mod <- addReact(mod, "ESP1", met=c(mql,h,mqn), Scoef=c(-1,2,1), lb=0, ub=1000, metComp = rep(1,3))
  mod <- addReact(mod, "ESP2", met=c(uql,h,uqn), Scoef=c(-1,2,1), lb=0, ub=1000, metComp = rep(1,3))
  mod <- addReact(mod, "ESP3", met=c(nadh,h,nad), Scoef=c(-1,1,1), lb=0, ub=1000, metComp = rep(1,3))
  mod <- addReact(mod, "ESP4", met=c(fdred,fdox), Scoef=c(-1,1), lb=0, ub=1000, metComp = rep(1,2))
  mod <- addReact(mod, "ESP5", met=c(pql,h,pqn), Scoef=c(-1,2,1), lb=0, ub=1000, metComp = rep(1,3))
  mod <- changeObjFunc(mod, react=c("ESP1", "ESP2", "ESP3", "ESP4", "ESP5"), obj_coef=c(1,1,1,1,1))

  return(mod)
}

# Arguments:
output.dir                <- "/home/bioinfo/projects/panGenome/panDraftEval/data/fig"

UHGG2proTraitId_fn <- "/home/bioinfo/projects/panGenome/panDraftEval/data/phenotypeEval/40IDs2spName_plus_NCBItaxa.tsv"
UHGG2proTraitId <- read.table(UHGG2proTraitId_fn, header = TRUE, sep = "\t")

proTrait_binary_fn <- "/home/bioinfo/projects/panGenome/panDraftEval/data/phenotypeEval/ProTraits_binaryIntegratedPr0.95.tsv"
proTrait_binary <- fread(proTrait_binary_fn, header = TRUE, sep = "\t")

mag.gapfill.model.pathList        <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft_18Oct/gapfill/MAG_Gapfilled_listsPath.txt"
singleGeno_pathList <- readLines(mag.gapfill.model.pathList)

panMAG.CompLV.gapfill.model.path  <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft_18Oct/CompLvFillQC/panFromMAG_Gapfilled_Path.txt"
# panMAG_pathList <- readLines(panMAG.CompLV.gapfill.model.path)

# metadata
completeness.level.table          <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/metadata/genomes-all_metadata.csv"
metadata <- read.table(completeness.level.table, sep = ",", header = TRUE, row.names = "Genome")

cpd4test.pathList <- "/home/bioinfo/projects/panGenome/panDraftEval/data/phenotypeEval/ProTraitcpd2ModelSEEDcpd.csv"
cpd4test <- read.csv(cpd4test.pathList, sep = ",", header = TRUE, row.names = "ProTrait_cpd")

MM_medium_fn <- "/home/bioinfo/bin/gapseq_panDraft/dat/media/MM_glu.csv"
MM_medium <- fread(MM_medium_fn)
MM_medium <- MM_medium[!(name %in% c("D-Glucose", "O2")), ] # Subsets the data table to include only rows where the value in the 'name' column is not equal to "D-Glucose"

# Parameters:

# Read metadata
columns_to_convert_numeric <- c("Completeness")
metadata[columns_to_convert_numeric] <- lapply(metadata[columns_to_convert_numeric], as.numeric)
metadata$Species_name <- rownames(metadata)

# Process input table
# subset the species of interest: 20 sp.
proTrait_binary <- proTrait_binary[proTrait_binary$Tax_ID %in% UHGG2proTraitId$ProTrait_ID,]
# remove empty columns
proTrait_binary[proTrait_binary == "?"] <- NA
col_atLeast_one_pred <- sapply(proTrait_binary, function(col) any(!is.na(col)))
proTrait_binary <- proTrait_binary[, ..col_atLeast_one_pred]

# select carbon sources 
carbon_sources <- names(proTrait_binary)[1:53]
elements_to_drop <- c("alkaline_phosphatase", "esterase_lipase__c8", "esterase__c4","gelatin_hydrolysis", "β-galactosidase") # also drop
carbon_sources <- carbon_sources[!carbon_sources %in% elements_to_drop]
carbon_sources_df <- proTrait_binary[, ..carbon_sources]
# exclude species without pred. phenotypes
row_atLeast_one_pred <- apply(carbon_sources_df[,c(-1,-2)], 1,function(row) any(!is.na(row)))
carbon_sources_df <- carbon_sources_df[row_atLeast_one_pred,]

# How many test can I perform
sum(!is.na(carbon_sources_df))
sum(is.na(carbon_sources_df))
result_0 <- rowSums(sapply(carbon_sources_df[, c(-1, -2)], function(x) as.numeric(!is.na(x) & x == "0")))
sum(result_0) # negative
result_1 <- rowSums(sapply(carbon_sources_df[, c(-1, -2)], function(x) as.numeric(!is.na(x) & x == "1")))
sum(result_1) # positive
# Which cpd are have many data
result_1 <- colSums(sapply(carbon_sources_df[, c(-1, -2)], function(x) as.numeric(!is.na(x) & x == "1")))
result_1[result_1 > 0] # positive in more than 1 sp.
result_0 <- colSums(sapply(carbon_sources_df[, c(-1, -2)], function(x) as.numeric(!is.na(x) & x == "0")))
result_0[result_0 > 1] # negative in more than 1 sp.

### Load models
uhgg_id_2_protrait_id <- UHGG2proTraitId[UHGG2proTraitId$ProTrait_ID %in% carbon_sources_df$Tax_ID, c("Species_name", "ProTrait_ID")]
rownames(uhgg_id_2_protrait_id) <- uhgg_id_2_protrait_id$Species_name

# GENOME OF MAG
singleGeno_pathList <- singleGeno_pathList[c(6)]
all_pred_res <- data.frame()
for (list in singleGeno_pathList) {
  dirty_spID <- tail(strsplit(list, "_")[[1]], n = 1)
  spID <- strsplit(dirty_spID, "\\.")[[1]][1]
  proTraitID <- uhgg_id_2_protrait_id[spID,"ProTrait_ID"]

  # tailor the medium to the species: select only cpd meaningfull to be tested
  cpds_to_be_pred_df <- !is.na(carbon_sources_df[carbon_sources_df$Tax_ID == proTraitID,])
  true_indices <- which(cpds_to_be_pred_df, arr.ind = TRUE)
  cpds_to_be_pred_list <- colnames(cpds_to_be_pred_df)[true_indices[, 2]][c(-1,-2)]
  gold_std_carbon_usage <- carbon_sources_df[carbon_sources_df$Tax_ID == proTraitID, ..cpds_to_be_pred_list]
  # cat("\n", "Testing row names:", paste(cpds_to_be_pred_list, collapse = ", "), "\n")
  cpd4test4sp <- cpd4test[cpds_to_be_pred_list,] # all cpd to be tested
  media2test <- data.table(compounds=cpd4test4sp$ModelSEED_cpd, name=cpd4test4sp$ModelSEED_name, maxFlux=100)

  # load only the model of genome of interest
  if (spID %in% uhgg_id_2_protrait_id$Species_name) {
    cat("\n", "Processing:", spID, "\n")
    singleGeno_mods <- load_files_from_paths_for_RDS(list, ".RDS") # open all the models of a single species
    sp_carbon_source_usage_res <- data.frame()
    counter <- 0
  }else{
    next # procede to the next reference (set of genomes)  
  }

  for (geno in names(singleGeno_mods)) {
    counter <- counter + 1
    # cat("\n", "Processing:", geno)
    mod <- singleGeno_mods[names(singleGeno_mods) == geno][[1]]
    mod <- setObjFun_altEnergy_sources(mod) # change ObjFun

    media2test.res <- data.table(compounds=cpd4test4sp$ModelSEED_cpd, sustainReducingComp=0)
    # add one carbon source at the time and perform the simulation
    for (i in 1:nrow(media2test)) {
      # update the medium
      tmp_MM_medium <- MM_medium 
      tmp_MM_medium <- rbind(tmp_MM_medium, media2test[i, ])
      mod <- constrain.model(mod, media = tmp_MM_medium) # constrain model
      sol <- optimizeProb(mod, retOptSol=F)

      if(sol$stat == ok & sol$obj >= 1e-7){
        media2test.res[match(media2test[i, ]$compounds, media2test.res$compounds),"sustainReducingComp"] <- 1
        # cat(paste(media2test[i, "name"], ":", sol$obj), "\n")
      }
      # else
        # cat(paste(media2test[i, "name"], ":", NA), "\n")
    }
    
    # store results single genome 
    names(media2test.res) <- c("compounds", geno)
    t_media2test.res <- as.data.frame(t(media2test.res[,-1]))
    names(t_media2test.res) <- media2test.res$compounds
    confusion_matrix_res_df <- compute_confusion_matrix(t_media2test.res, gold_std_carbon_usage)
    
    t_media2test.res$proTrait_Id <- proTraitID
    t_media2test.res$ref_Id <- spID
    t_media2test.res$mag_Id <- geno
    sp_carbon_source_usage_res <- rbind(sp_carbon_source_usage_res, t_media2test.res)  
  }

  merge_df <- merge(sp_carbon_source_usage_res, confusion_matrix_res_df, by="mag_Id")
    mod <- 
  t.merge_df <- as.data.frame(t(merge_df))
  t.merge_df$Rownames <- rownames(t.merge_df)
  if (all(dim(all_pred_res) == 0)) { # if the final dataset is still empty append the restults
    all_pred_res <- rbind(all_pred_res, t.merge_df)
  }else{ # otherwise merge them
    all_pred_res <- merge(all_pred_res, t.merge_df, by = "Rownames", all = TRUE)
  }
}

# Transpose final results 
rownames(all_pred_res) <- all_pred_res$Rownames
all_pred_res$Rownames <- NULL
colnames(all_pred_res) <- all_pred_res["mag_Id",]
t.all_pred_res <- as.data.frame(t(all_pred_res))
head(t.all_pred_res)
dim(t.all_pred_res)



### PAN-DRAFT
all_pred_res <- data.frame()
panMAG.CompLV.gapfill.mods <- load_files_from_paths_for_RDS(panMAG.CompLV.gapfill.model.path, ".RDS") 
for (mod_idx in names(panMAG.CompLV.gapfill.mods)) {
  mod_idx <- "MGYG000002492_iter_7_Compth_80_panModel"
  spID <- strsplit(mod_idx, "_")[[1]][1]
  iter_num <- strsplit(mod_idx, "_")[[1]][3]
  compLv <- strsplit(mod_idx, "_")[[1]][5]
  proTraitID <- uhgg_id_2_protrait_id[spID, "ProTrait_ID"]

  # load only the model of genome of interest
  if (spID %in% uhgg_id_2_protrait_id$Species_name) {
    cat("\n", "Processing:", mod_idx)
    mod <- panMAG.CompLV.gapfill.mods[mod_idx]
  }else{
    next # procede to the next reference (set of genomes)  
  }

  # tailor the medium to the species: select only cpd meaningfull to be tested
  cpds_to_be_pred_df <- !is.na(carbon_sources_df[carbon_sources_df$Tax_ID == proTraitID,])
  true_indices <- which(cpds_to_be_pred_df, arr.ind = TRUE)
  cpds_to_be_pred_list <- colnames(cpds_to_be_pred_df)[true_indices[, 2]][c(-1,-2)]
  gold_std_carbon_usage <- carbon_sources_df[carbon_sources_df$Tax_ID == proTraitID, ..cpds_to_be_pred_list]
  # cat("\n", "Testing row names:", paste(cpds_to_be_pred_list, collapse = ", "), "\n")
  cpd4test4sp <- cpd4test[cpds_to_be_pred_list,] # all cpd to be tested
  media2test <- data.table(compounds=cpd4test4sp$ModelSEED_cpd, name=cpd4test4sp$ModelSEED_name, maxFlux=100)

  # add one carbon source at the time and perform the simulation
  mod <- setObjFun_altEnergy_sources(mod) # change ObjFun
  media2test.res <- data.table(compounds=cpd4test4sp$ModelSEED_cpd, sustainReducingComp=0)
  for (i in 1:nrow(media2test)) {
    # update the medium
    tmp_MM_medium <- MM_medium 
    tmp_MM_medium <- rbind(tmp_MM_medium, media2test[i, ])
    mod <- constrain.model(mod, media = tmp_MM_medium) # constrain model
    sol <- optimizeProb(mod, retOptSol = F)

    if(sol$stat == ok & sol$obj >= 1e-7){
      media2test.res[match(media2test[i, ]$compounds, media2test.res$compounds),"sustainReducingComp"] <- 1
      # cat(paste(media2test[i, "name"], ":", sol$obj), "\n")
    }
    # else
      # cat(paste(media2test[i, "name"], ":", NA), "\n")
  }
    
  # store results single genome 
  names(media2test.res) <- c("compounds", mod_idx)
  t_media2test.res <- as.data.frame(t(media2test.res[,-1]))
  names(t_media2test.res) <- media2test.res$compounds
  confusion_matrix_res_df <- compute_confusion_matrix(t_media2test.res, gold_std_carbon_usage)

  t_media2test.res$proTrait_Id <- proTraitID
  t_media2test.res$ref_Id <- spID
  t_media2test.res$compLv <- compLv
  t_media2test.res$iter_num <- iter_num
  t_media2test.res$mag_Id <- mod_idx
  merge_df <- cbind(t_media2test.res, t(confusion_matrix_res_df))
   
  t.merge_df <- as.data.frame(t(merge_df))
  t.merge_df$Rownames <- rownames(t.merge_df)
  if (all(dim(all_pred_res) == 0)) { # if the final dataset is still empty append the restults
    all_pred_res <- rbind(all_pred_res, t.merge_df)
  }else{ # otherwise merge them
    all_pred_res <- merge(all_pred_res, t.merge_df, by = "Rownames", all = TRUE)
  }
}

# Transpose final results 
rownames(all_pred_res) <- all_pred_res$Rownames
all_pred_res$Rownames <- NULL
colnames(all_pred_res) <- all_pred_res["mag_Id",]
t.all_pred_res <- as.data.frame(t(all_pred_res))
head(t.all_pred_res)
dim(t.all_pred_res)








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


