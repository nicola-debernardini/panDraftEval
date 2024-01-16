# conda activate prj_anne2 (ion)
import pandas as pd
import numpy as np
import cobra as cb
import seaborn as sns
import os

config = cb.Configuration()
config.solver = 'cplex'

# Load metadata
completeness_level_table = "/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/metadata/genomes-all_metadata.csv"
metadata = pd.read_csv(completeness_level_table, sep=",", index_col="Genome")
metadata["Completeness"] = pd.to_numeric(metadata["Completeness"], errors="coerce")
metadata["Species_name"] = metadata.index

# List of compound to be tested
cpd4test_path_list = "/home/bioinfo/projects/panGenome/panDraftEval/data/phenotypeEval/ProTraitcpd2ModelSEEDcpd.csv"
cpd4test = pd.read_csv(cpd4test_path_list, sep=",", index_col="ProTrait_cpd")

# Read input
UHGG2proTraitId_fn = "/home/bioinfo/projects/panGenome/panDraftEval/data/phenotypeEval/40IDs2spName_plus_NCBItaxa.tsv"
UHGG2proTraitId = pd.read_table(UHGG2proTraitId_fn, sep="\t", index_col="Species_name")

proTrait_binary_fn = "/home/bioinfo/projects/panGenome/panDraftEval/data/phenotypeEval/ProTraits_binaryIntegratedPr0.95.tsv"
proTrait_binary = pd.read_csv(proTrait_binary_fn, sep="\t")

gut_medium_fn = "/home/bioinfo/bin/gapseq_panDraft/dat/media/gut.csv"
gut_medium = pd.read_csv(gut_medium_fn)
gut_medium['exchange_reaction'] = 'EX_' + gut_medium['compounds'] + '_e0'
gut_medium_dict = dict(zip(gut_medium['exchange_reaction'], gut_medium['maxFlux']))

# Process input table
proTrait_binary = proTrait_binary[proTrait_binary["Tax_ID"].isin(UHGG2proTraitId["ProTrait_ID"])]
proTrait_binary.replace("?", pd.NA, inplace=True)
col_at_least_one_pred = proTrait_binary.columns[proTrait_binary.notna().any()]
proTrait_binary = proTrait_binary[col_at_least_one_pred]

# Select carbon sources
carbon_sources = proTrait_binary.columns[1:54]
elements_to_drop = ["alkaline_phosphatase", "esterase_lipase__c8", "esterase__c4", "gelatin_hydrolysis"]
carbon_sources = carbon_sources.difference(elements_to_drop)
carbon_sources_df = proTrait_binary[list(carbon_sources)]

# Exclude species without pred. phenotypes
row_at_least_one_pred = carbon_sources_df.iloc[:, 1:].notna().any(axis=1)
carbon_sources_df = carbon_sources_df[row_at_least_one_pred]

# How many tests can be performed
print(f"Null tests: {pd.isna(carbon_sources_df.iloc[:, 2:]).sum().sum()}")
print(f"Total tests: {pd.notna(carbon_sources_df.iloc[:, 2:]).sum().sum()}")
num_positive_tests = (carbon_sources_df.iloc[:, 2:] == "1").sum(axis=1)
num_negative_tests = (carbon_sources_df.iloc[:, 2:] == "0").sum(axis=1)
print(f"Positive tests: {num_positive_tests.sum()}")  # Positive
print(f"Negative tests: {num_negative_tests.sum()}")  # Negative

# Which cpds have many data
positive_cpds = (carbon_sources_df.iloc[:, 2:] == "1").sum(axis=0)
negative_cpds = (carbon_sources_df.iloc[:, 2:] == "0").sum(axis=0)
print("Positive CPDs:")
print(positive_cpds[positive_cpds > 1])  # Positive in more than 1 sp.
print("Negative CPDs:")
print(negative_cpds[negative_cpds > 1])  # Negative in more than 1 sp.

# Load models
def load_file_from_path_for_xml(file_path, f_format):
  new_mod_fn = os.path.splitext(file_path)[0] + f"{f_format}"
  model = cb.io.read_sbml_model(new_mod_fn)
  dirty_sp_id = file_path.split("/")[-1]
  sp_id = dirty_sp_id.split(".")[0]
    
  return (sp_id, model)

uhgg_id_sp_of_interest = UHGG2proTraitId.loc[UHGG2proTraitId["ProTrait_ID"].isin(carbon_sources_df["Tax_ID"]),].index
# single_geno_path_list = [line.strip() for line in open("/home/bioinfo/users/niber/prj_panModel/db/UHGG_v2.1_db/pan.draft_18Oct/gapfill/MAG_Gapfilled_listsPath.txt")]
single_geno_path_list = [line.strip() for line in open("/home/bioinfo/projects/panGenome/panDraftEval/data/phenotypeEval/MAG_listPath_prova.txt")]
single_geno_path_list = single_geno_path_list[0:2] # tmp
for file_path in single_geno_path_list:
    dirty_sp_id = file_path.split("_")[-1]
    rep_id = dirty_sp_id.split(".")[0]
    print(f"Working on referrance: {rep_id}")
    
    for mod_fn in [line.strip() for line in open(file_path)]:  
      sp_id, single_geno_mod = load_file_from_path_for_xml(mod_fn, ".xml.gz") # Assuming models are saved in MATLAB format
      print(f"Processing model: {sp_id}")
      
      # # load the GUT medium ("gut.csv")
      # medium_rxn_present_dict = {}
      # for rxn in gut_medium_dict.keys():
      #     if single_geno_mod.reactions.has_id(rxn):
      #         medium_rxn_present_dict[rxn] = gut_medium_dict[rxn]
      #     else: 
      #         pass
      # single_geno_mod.medium = medium_rxn_present_dict

      single_geno_medium = single_geno_mod.medium 
      # add one at the time the compounds of interest
      for cpd in cpd4test["ModelSEED_cpd"]:
        test4cpd = 'EX_' + cpd + '_e0'
        if single_geno_mod.reactions.has_id(test4cpd):
          single_geno_medium_cpdTest = dict(single_geno_medium)         
          single_geno_medium_cpdTest[test4cpd] = 100
          single_geno_mod.medium = single_geno_medium_cpdTest
          # print(f"added: {cpd}")        
          # select the objective function: ESP1 + ESP2 + ESP3
          
          # run the simulation
          single_geno_sol = single_geno_mod.optimize()
          print(single_geno_sol.objective_value)

        else:
          pass

# rxn17049  
single_geno_mod.reactions
single_geno_mod.reactions.has_id("rxn17049_c0")


single_geno_mod.medium["EX_cpd00054_e0"]
single_geno_medium["EX_cpd00054_e0"]

# Display reactions
print("Reactions:")
for reaction in single_geno_mod.reactions:
    print(reaction.id)

# Display metabolites
print("\nMetabolites:")
for metabolite in model.metabolites:
    print(metabolite.id)