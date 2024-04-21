###################################################### @
# Script to reproduce the fermentation test (figure 5):
# !!! due to space limitation: only a subset of all the GEMs used in the paper are provided in the github directory !!!

# Load R-packages
library(data.table)
library(ggplot2) 
library(sybil)

# select solver
if( "cplexAPI" %in% rownames(installed.packages()) ){
  sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1
}else{
  warning("glpkAPI is used but cplexAPI is recommended because it is much faster")
  sybil::SYBIL_SETTINGS("SOLVER","glpkAPI"); ok <- 5
}

# Little helpers
source("./scr/function_collection.R")
source("~/bin/gapseq_panDraft/src/constrain.model.R") # path to the constrain.model.R script in your gapseq directory 
source("./scr/sybil_toolkit.R")

# Function to compute confusion matrix
compute_confusion_matrix <- function(predicted, actual) {
  tp <- sum(predicted == 1 & actual == 1)
  fp <- sum(predicted == 1 & actual == 0)
  tn <- sum(predicted == 0 & actual == 0)
  fn <- sum(predicted == 0 & actual == 1)
  return(c(TP = tp, FP = fp, TN = tn, FN = fn))
}

# Arguments:
output.dir              <- "./out"
fermentation_product_fn <- "./dat/uhgg/ferm_prod_test/fermentation_product_info.csv" # known fermentation products
mag.gapfill.model.pathList        <- "./dat/uhgg/ferm_prod_test/MAG_GapfilledFT_listsPath.txt" # location of GEMs for the MAGs divided by SGB
panMAG.CompLV.gapfill.model.path  <- "./dat/uhgg/ferm_prod_test/panFromMAG_GapfilledFT_Path.txt" # location of pan-GEMs
completeness.level.table          <- "./dat/uhgg/metadata/genomes-all_metadata.csv" # genomes metadat
medium_fn <- "./dat/uhgg/ferm_prod_test/media/FT.csv" # growth medium

# read reference set of fermentation products
fermentation_product_binary <- fread(fermentation_product_fn, header = TRUE, sep = ",", skip=1)
# read GEM location
singleGeno_pathList <- readLines(mag.gapfill.model.pathList)

# The simulations of anaerobic growth were all performed assuming the same growth medium ("FT.csv"). 
# It comprised several organic compounds (i.e. carbohydrates, polyols, nucleotides, amino acids, organic acids) as potential energy sources
medium <- fread(medium_fn)
medium <- medium[!(name %in% c("O2")), ] # exclude the oxigen from the medium

# read metadata
metadata <- read.table(completeness.level.table, sep = ",", header = TRUE, row.names = "Genome")
metadata$Species_name <- rownames(metadata)

### Load models
uhgg_id_4_fermentation_test <- fermentation_product_binary$UHGG_id
### MAGs ---------------
# limit_mags_to <- 2 # compute the simulation for a limited number of GEMs per SGB
prod.mets <- list()
growth    <- list()
# Flux-Variability-Analysis (FVA) AND Flux Balance Analysis (FBA) coupled with a minimisation of total flux
for (list in singleGeno_pathList) {
  counter <- 0
  dirty_spID <- tail(strsplit(list, "_")[[1]], n = 1)
  spID <- strsplit(dirty_spID, "\\.")[[1]][1]

  # load the GEMs of an SGB only if the SGB has any known associated fermentation product
  if (spID %in% uhgg_id_4_fermentation_test) {
    cat("\n", "Processing:", spID, "\n")
    singleGeno_mods <- load_files_from_paths_for_RDS(list, ".RDS") # load all the models of the SGB
  }else{
    next # procede to the next reference (set of genomes)  
  }

  for (geno in names(singleGeno_mods)) {
    # limit the analysis to n GEMs per SGB
    # if (counter >= limit_mags_to) { 
    #   next
    # }
    cat("\n", geno, "\n")
    mod <- singleGeno_mods[names(singleGeno_mods) == geno][[1]]
    mod <- constrain.model(mod, media = medium) # constrain the model
    
    # Growth rate
    growth[[geno]] <- optimizeProb(mod)@lp_obj
    if (growth[[geno]] == 0){
      cat(paste(geno, "is not able to grow\n"))
      dt <- data.table(matrix(NA, nrow = 1, ncol = 5))
      colnames(dt) <- c("ex", "rxn.name", "l", "u", "mtf.flux")
      prod.mets[[geno]] <- dt
      next
    }
    prod.mets[[geno]] <- get.produced.metabolites(mod)
    counter <- counter + 1
  }
}

# We normalised the outflow of the individual fermentation products by the predicted growth rate of the respective organism
# Metabolites with a positive exchange flux (i.e. outflow) are considered as fermentation products
gs.fermprod <- rbindlist(prod.mets, idcol = T)
for(i in 1:length(growth)) {
  gs.fermprod[.id == names(growth[i]), gr.rate := growth[[i]]]
}
gs.fermprod[, recon.method := "gapseq"]

# save
fwrite(gs.fermprod, file = file.path(output.dir, "mag_fermProd.tsv"), sep = "\t", quote = FALSE)

### PAN-DRAFT ---------------
panMAG.CompLV.gapfill.mods <- load_files_from_paths_for_RDS(panMAG.CompLV.gapfill.model.path, ".RDS") 

prod.mets <- list()
growth    <- list()

# Flux-Variability-Analysis (FVA) AND Flux Balance Analysis (FBA) coupled with a minimisation of total flux
for (mod_idx in names(panMAG.CompLV.gapfill.mods)) {
  spID <- strsplit(mod_idx, "_")[[1]][1]
  iter_num <- strsplit(mod_idx, "_")[[1]][3]
  compLv <- strsplit(mod_idx, "_")[[1]][5]

  # load the GEMs of an SGB only if the SGB has any known associated fermentation product
  if (spID %in% uhgg_id_4_fermentation_test) {
    cat("\n", "Processing:", mod_idx)
  }else{
    next 
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

  prod.mets[[mod_idx]] <- get.produced.metabolites(mod)
  prod.mets[[mod_idx]][, ref_uhgg_id := spID]
  prod.mets[[mod_idx]][, iter_num := iter_num]
  prod.mets[[mod_idx]][, compLv := compLv]
}

gs.fermprod_pan <- rbindlist(prod.mets, idcol = T)
for(i in 1:length(growth)) {
  gs.fermprod_pan[.id == names(growth[i]), gr.rate := growth[[i]]]
}
gs.fermprod_pan[, recon.method := "pan"]

# save
fwrite(gs.fermprod_pan, file = file.path(output.dir, "panDraft_fermProd.tsv"), sep = "\t", quote = FALSE)