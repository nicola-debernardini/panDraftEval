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
source("/home/bioinfo/projects/panGenome/panDraftEval/scr/sybil_toolkit.R")

# Function to compute confusion matrix
compute_confusion_matrix <- function(predicted, actual) {
  tp <- sum(predicted == 1 & actual == 1)
  fp <- sum(predicted == 1 & actual == 0)
  tn <- sum(predicted == 0 & actual == 0)
  fn <- sum(predicted == 0 & actual == 1)
  
  return(c(TP = tp, FP = fp, TN = tn, FN = fn))
}

# Arguments:
output.dir                <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/fermProdTest/out"

fermentation_product_fn <- "/home/bioinfo/projects/panGenome/panDraftEval/data/phenotypeEval/fermentation_product_info.csv"
fermentation_product_binary <- fread(fermentation_product_fn, header = TRUE, sep = ",", skip=2)

mag.gapfill.model.pathList        <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/fermProdTest/MAG_GapfilledFT_listsPath.txt"
singleGeno_pathList <- readLines(mag.gapfill.model.pathList)

panMAG.CompLV.gapfill.model.path  <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/fermProdTest/panFromMAG_GapfilledFT_Path.txt"

# metadata
completeness.level.table          <- "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/metadata/genomes-all_metadata.csv"
metadata <- read.table(completeness.level.table, sep = ",", header = TRUE, row.names = "Genome")

# The simulations of anaerobic growth were all performed assuming the same growth medium ("FT.csv"). 
# It comprised several organic compounds (i.e. carbohydrates, polyols, nucleotides, amino acids, organic acids) as potential energy sources
medium_fn <- "/home/bioinfo/bin/gapseq_panDraft/dat/media/FT.csv"
medium <- fread(medium_fn)
medium <- medium[!(name %in% c("O2")), ] # Subsets the data table to include only rows where the value in the 'name' column is not equal to "D-Glucose"

# Read metadata
columns_to_convert_numeric <- c("Completeness")
metadata[columns_to_convert_numeric] <- lapply(metadata[columns_to_convert_numeric], as.numeric)
metadata$Species_name <- rownames(metadata)

### Load models
uhgg_id_4_fermentation_test <- fermentation_product_binary$UHGG_id
# uhgg_id_4_fermentation_test <- uhgg_id_4_fermentation_test[1:6]

# GENOME OF MAG
limit_mags_to <- 1000
prod.mets <- list()
growth    <- list()
for (list in singleGeno_pathList) {
  counter <- 0
  # list = "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/input/MAG/input4modelMerging_Gapfilled_MGYG000001292.txt"
  dirty_spID <- tail(strsplit(list, "_")[[1]], n = 1)
  spID <- strsplit(dirty_spID, "\\.")[[1]][1]

  # load only the model of genome of interest
  if (spID %in% uhgg_id_4_fermentation_test) {
    cat("\n", "Processing:", spID, "\n")
    singleGeno_mods <- load_files_from_paths_for_RDS(list, ".RDS") # open all the models of a single species
  }else{
    next # procede to the next reference (set of genomes)  
  }

  for (geno in names(singleGeno_mods)) {
    # limit the analysis to 100 GEMs per species
    if (counter >= limit_mags_to) {
      next
    }

    # geno = "MGYG000284276"
    cat("\n", geno, "\n")
    mod <- singleGeno_mods[names(singleGeno_mods) == geno][[1]]
    mod <- constrain.model(mod, media = medium) # constrain model
    
    # Growth rate
    growth[[geno]] <- optimizeProb(mod)@lp_obj
    if (growth[[geno]] == 0){
      cat(paste(geno, "is not able to grow\n"))
      dt <- data.table(matrix(NA, nrow = 1, ncol = 5))
      colnames(dt) <- c("ex", "rxn.name", "l", "u", "mtf.flux")
      prod.mets[[geno]] <- dt
      next
    }
    counter <- counter + 1

    # Flux-Variability-Analysis (FVA) AND Flux Balance Analysis (FBA) coupled with a minimisation of total flux
    prod.mets[[geno]] <- get.produced.metabolites(mod)
  }
}

# Metabolites with a positive exchange flux (i.e. outflow) are considered as fermentation products
gs.fermprod <- rbindlist(prod.mets, idcol = T)

for(i in 1:length(growth)) {
  gs.fermprod[.id == names(growth[i]), gr.rate := growth[[i]]]
}
gs.fermprod[, recon.method := "gapseq"]

# We normalised the outflow of the individual fermentation products by the predicted growth rate of the respective organism
# SAVE
fwrite(gs.fermprod, file = file.path(output.dir, "mag_fermProd_31Gen.tsv"), sep = "\t", quote = FALSE)
# gs.fermprod <- read.table(file.path(output.dir, "mag_fermProd_31Gen.tsv"), header = TRUE, sep = "\t", quote = "")
head(gs.fermprod)


### PAN-DRAFT
panMAG.CompLV.gapfill.mods <- load_files_from_paths_for_RDS(panMAG.CompLV.gapfill.model.path, ".RDS") 

prod.mets <- list()
growth    <- list()
for (mod_idx in names(panMAG.CompLV.gapfill.mods)) {
  spID <- strsplit(mod_idx, "_")[[1]][1]
  iter_num <- strsplit(mod_idx, "_")[[1]][3]
  compLv <- strsplit(mod_idx, "_")[[1]][5]

  # load only the model for genome of interest
  if (spID %in% uhgg_id_4_fermentation_test) {
    cat("\n", "Processing:", mod_idx)
  }else{
    next # procede to the next reference (set of genomes)  
  }

  mod <- panMAG.CompLV.gapfill.mods[mod_idx][[1]]
  mod <- constrain.model(mod, media = medium) # constrain model
  
  # Growth rate
  growth[[mod_idx]] <- optimizeProb(mod)@lp_obj
  if (growth[[mod_idx]] == 0) {
    cat(paste("\n", mod_idx, "is not able to grow\n"))
    dt <- data.table(matrix(NA, nrow = 1, ncol = 8))
    colnames(dt) <- c("ex", "rxn.name",
                      "l", "u", "mtf.flux",
                      "ref_uhgg_id", "iter_num",
                      "compLv")
    prod.mets[[mod_idx]] <- dt

    next
  }

  # Flux-Variability-Analysis (FVA) AND Flux Balance Analysis (FBA) coupled with a minimisation of total flux
  prod.mets[[mod_idx]] <- get.produced.metabolites(mod)
  prod.mets[[mod_idx]][, ref_uhgg_id := spID]
  prod.mets[[mod_idx]][, iter_num := iter_num]
  prod.mets[[mod_idx]][, compLv := compLv]
}

# Metabolites with a positive exchange flux (i.e. outflow) are considered as fermentation products
gs.fermprod_pan <- rbindlist(prod.mets, idcol = T)
for(i in 1:length(growth)) {
  gs.fermprod_pan[.id == names(growth[i]), gr.rate := growth[[i]]]
}
gs.fermprod_pan[, recon.method := "pan"]

# We normalised the outflow of the individual fermentation products by the predicted growth rate of the respective organism
# SAVE
fwrite(gs.fermprod_pan, file = file.path(output.dir, "panDraft_fermProd_31Gen.tsv"), sep = "\t", quote = FALSE)
# gs.fermprod_pan <- read.table(file.path(output.dir, "panDraft_fermProd_31Gen.tsv"), header = TRUE, sep = "\t", quote = "")
head(gs.fermprod_pan)