###################################################### @
# Functions collection

compute_contigency_table <- function(ref_rxn, pred_rxn) {
  all_reactions <- unique(c(ref_rxn, pred_rxn)) # Create a vector of all unique reactions
  common_reactions <- intersect(ref_rxn, pred_rxn)
  contingency_table <- table(
    Reference = factor(all_reactions %in% ref_rxn, levels = c(FALSE, TRUE)),
    Predicted = factor(all_reactions %in% pred_rxn, levels = c(FALSE, TRUE))
  )

  tp <- contingency_table[2, 2]
  fn <- contingency_table[2, 1]
  fp <- contingency_table[1, 2]
  tn <- contingency_table[1, 1]
  
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  f1_score <- 2 * (precision * recall) / (precision + recall)
  FDR <- fp / (tp + fp)
  
  metrics.list <- list(
                    accuracy=accuracy, 
                    precision=precision, 
                    recall=recall, 
                    f1_score=f1_score,
                    FDR=FDR)
  return(metrics.list)
}


compute_stat4diffMAGsetSize <- function(data, MAG_cols, th, reference_reactions, n.iter, n.MAG.cols, which.stat) {
  prob_dist_df <- data.frame(step = numeric(0), stat = numeric(0), group = numeric(0)) # df to save prob
  # compute statistics for n iteration
  for (n_iter in seq_len(n.iter)){
    cat(paste("\niteration",n_iter))
    rand_MAG_cols <- sample(MAG_cols, replace = FALSE) # randomize input
    for (n_MAG in seq_len(n.MAG.cols)){
      # given a randomized list of MAG select a subset, the predicted rxn are all the reaction in the subset
      sum_sub_data <- data %>%
        select(rand_MAG_cols[1:n_MAG]) %>%
        mutate(freq_rxn = rowSums(.),
          rxn = data$rxn) %>%
        filter(freq_rxn >= 1) %>%
        filter((freq_rxn/n_MAG) >= th)

      predicted_reactions <- sum_sub_data$rxn
      metrics_table <- compute_contigency_table(reference_reactions, predicted_reactions, fb.score.par) # computer contingency table for a pan.mod
      prob_dist_df <- rbind(prob_dist_df, c(n_MAG, metrics_table[which.stat][[1]], n_iter))
    }
  }
  colnames(prob_dist_df) <- c("num_MAGs", "stat", "group") 
  return(prob_dist_df)
}

compute_stat4diffMAGsetSize_parallel <- function(data, MAG_cols, th, reference_reactions, n.iter, n.MAG.cols, which.stat) {
  cl <- makeCluster(num.cores) # Initialize a parallel backend with n cores
  registerDoParallel(cl)
  clusterEvalQ(cl, { library(dplyr) })  # Load the necessary packages in the parallel workers
  clusterExport(cl, c("compute_contigency_table", "fb.score.par")) # Define your custom function compute_contigency_table inside the parallel workers

  prob_dist_df <- data.frame(step = numeric(0), stat = numeric(0), group = numeric(0))
  iter_list <- seq_len(n.iter)
  # Parallelize the loop for each iteration 
  result_list <- foreach(n_iter = iter_list, .combine = rbind) %dopar% { 
    rand_MAG_cols <- sample(MAG_cols, replace = FALSE) # randomize
    result_df <- data.frame()
    for (n_MAG in seq_len(n.MAG.cols)){
      # given a randomized list of MAG select a subset, the predicted rxn are all the reaction in the subset
      sum_sub_data <- data %>%
        select(rand_MAG_cols[1:n_MAG]) %>%
        mutate(freq_rxn = rowSums(.),
          rxn = data$rxn) %>%
        filter(freq_rxn >= 1) %>%
        filter((freq_rxn/n_MAG) >= th)
    
      predicted_reactions <- rownames(sum_sub_data)
      metrics_table <- compute_contigency_table(reference_reactions, predicted_reactions, fb.score.par) # compute contingency table for a pan.mod
      result_df <- rbind(result_df, c(n_MAG, metrics_table[which.stat][[1]], n_iter))
    }
    colnames(result_df) <- c("num_MAGs", "stat", "group")
    return(result_df)
  }
  prob_dist_df <- do.call(rbind, result_list)  # Combine the results from all iterations
  stopCluster(cl) # Stop the parallel backend
  return(prob_dist_df)
}

# FUNCTION - compute the Tsgb score.
# similarity score between a vector containing the frequency of reaction in reference genome and a vector of present/absent reaction in a MAG   
Tsgb_cal <- function(mag_data, freq_ref_rxn, mag) {
    Tsgb_df <- merge(mag_data, freq_ref_rxn, by.x = "rxn", by.y = "rxn", all.y = TRUE)
    Tsgb_df[is.na(Tsgb_df)] <- 0 # Replace NA values with 0
    Tsgb_df$numerator_i <- Tsgb_df[[mag]] * Tsgb_df[,"b_i"]
    Tsgb <- sum(Tsgb_df$numerator_i) / sum(Tsgb_df$b_i)
    return(Tsgb)
}

# FUNCTION - compute the distance score
distance_cal <- function(mag_data, freq_ref_rxn, dist_type){
    ref_mtrx <- data.frame(matrix(1, nrow = nrow(freq_ref_rxn), ncol = 1))
    names(ref_mtrx) <- "ref"
    ref_mtrx$rxn <- freq_ref_rxn$rxn

    # Merge mag_data and ref_mtrx, including missing rows and columns
    data_mtrx <- merge(mag_data, ref_mtrx, by = "rxn", all = TRUE)
    data_mtrx[is.na(data_mtrx)] <- 0 # Replace NA values with 0
    col_name <- data_mtrx$rxn 
    data_mtrx$rxn <- NULL # Remove the "Row.names" column
    t_data_mtrx <- t(data_mtrx) 
    colnames(t_data_mtrx) <- col_name # Set column names to match the original row names
    dist_val <- dist(t_data_mtrx, method = dist_type)
    return(dist_val)
}

### load files (.RDS) from paths
# specify suffix to save the prefix as id in a list
load_files_from_paths <- function(list_fn, suf) {
  cat("\tloading file from:", list_fn, "\n")
  file_paths <- readLines(list_fn)
  file_list <- list() # Create an empty list to store the loaded files
  # Iterate over the file paths and read the corresponding files
  for (path in file_paths) {
    file_data <- read.table(path, sep="\t", header = TRUE)
    filename <- head(tail(strsplit(path, "/")[[1]], 2), 1)
    file_list[[filename]] <- file_data
  }
  return(file_list)
}

### load files (.RDS) from paths
# specify suffix to save the prefix as id in a list
load_files_from_paths_for_RDS <- function(list_fn, suf) {
  cat("\tloading file from:", list_fn, "\n")
  file_paths <- readLines(list_fn)
  file_list <- list() # Create an empty list to store the loaded files
  # Iterate over the file paths and read the corresponding files
  for (path in file_paths) {
    file_data <- readRDS(path)
    filename <- basename(path)
    filename_without_suffix <- sub(suf, "", filename)
    file_list[[filename_without_suffix]] <- file_data
  }
  return(file_list)
}

### load files (.tbl) from paths
load_files_from_paths_for_tbl <- function(list_fn, suf) {
  cat("\tloading file from:", list_fn, "\n")
  file_paths <- readLines(list_fn)
  file_list <- list() # Create an empty list to store the loaded files
  # Iterate over the file paths and read the corresponding files
  for (path in file_paths) {
    file_data <- fread(path)
    filename <- basename(path)
    filename_without_suffix <- sub(suf, "", filename)
    file_list[[filename_without_suffix]] <- file_data
  }
  return(file_list)
}

# extract taxonomy form GTDB annotation
parse_GTDBtaxonomy <- function(df, colGTDBtaxonomy) {
    df %>%
        mutate(domain = sapply(strsplit(df[[colGTDBtaxonomy]], ";"), function(x) x[1])) %>%
        mutate(phylum = sapply(strsplit(df[[colGTDBtaxonomy]], ";"), function(x) x[2])) %>%
        mutate(class = sapply(strsplit(df[[colGTDBtaxonomy]], ";"), function(x) x[3])) %>%
        mutate(order = sapply(strsplit(df[[colGTDBtaxonomy]], ";"), function(x) x[4])) %>%
        mutate(family = sapply(strsplit(df[[colGTDBtaxonomy]], ";"), function(x) x[5])) %>%
        mutate(genera = sapply(strsplit(df[[colGTDBtaxonomy]], ";"), function(x) x[6])) %>%
        mutate(species = sapply(strsplit(df[[colGTDBtaxonomy]], ";"), function(x) x[7])) %>%   
        return(df)
}
# refine GTDB annotation by removing the "s__"
refine_GTDBtaxonomy <- function(df) {
    df %>%
        mutate(domain = sapply(strsplit(df[["domain"]], "__"), function(x) x[2])) %>%
        mutate(phylum = sapply(strsplit(df[["phylum"]], "__"), function(x) x[2])) %>%
        mutate(class = sapply(strsplit(df[["class"]], "__"), function(x) x[2])) %>%
        mutate(order = sapply(strsplit(df[["order"]], "__"), function(x) x[2])) %>%
        mutate(family = sapply(strsplit(df[["family"]], "__"), function(x) x[2])) %>%
        mutate(genera = sapply(strsplit(df[["genera"]], "__"), function(x) x[2])) %>%
        mutate(species = sapply(strsplit(df[["species"]], "__"), function(x) x[2])) %>%
        return(df)
}