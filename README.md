# QSAR_predictions_database

This repository contains the script(s) and output data described in the "Dataset on ecotoxicity predictions of 2676 chemicals, using three ecotoxicological quantitative structure-activity relationship platforms" journal article.

These scripts may be used under the Creative Commons Attribution-NonCommercial (CC-BY-NC) license, if you need to cite it please use the CITATION.cff. 

## Main output

The main output can be found in the "Output" folder.
They can also be accessed through these links:

|Output file|Description|
|-----------|-----------|
[QSAR_predictions.tsv](Output/QSAR_predictions.tsv)|QSAR predictions database, wide format|
[experimental_dataset.tsv](Output/experimental_dataset.tsv)|Curated empirical data, long format|
[identifiers.tsv](Output/identifiers.tsv)|The identifiers and physico-chemical data of compounds in the QSAR database and experimental data|
[QSAR_predictions_database_content_descriptions.xlsx](Output/QSAR_predictions_database_content_descriptions.xlsx)|Excel sheets describing the contents of the database, empirical data and identifiers|


## Below follows a short guide to using the script

QSAR_data_collection_script.R is the main script used to produce the dataset. It calls a number of helper functions, located under the subfolder “Functions”. The main script handles directory management and packages, standardization and merging of empirical data, extraction of identifiers from empirical data, gathering of physicochemical properties, output data formatting and output data export.

Specific functions are called for import of empirical data, cleaning and filtering of empirical data, collection of SMILES and InChIKeys, running QSAR platforms, importing and curating QSAR outputs and calculating compounded QSAR predictions. A full list of function files and a short description of their uses can be found below.

A handful of intermediate files are created by the script to reduce the load on third party APIs and to avoid rerunning time consuming steps of the script, unless instructed to. These files can be found under the “Intermediate files” folder. A full list of intermediate files and a short description of their purpose can be found below.

To run the script, first you need to install and load the required packages, and specify certain paths (working directories, QSAR directories, empirical data locations and QSAR tool executables), finally the functions should be imported and a few key parameters should be set. All of this is done in the main script under "1. Packages, working directories and housekeeping".

The script should then be run section by section, with additional user interaction needed for running the QSAR tools.

The QSAR tools are launched by the script, but the tools needs to be run by the user, with instructions available through the script. If the user need additional help running the QSAR platforms, please find the documentation for each tool from their corresponding developers.

The script was written in R version 4.1.3 using the RStudio editor version 2022.12.0 Build 353

## The main script and all helper functions rely on the following R-packages:

|Package|Version
|-----------|------------|
|data.table      |1.14.2|
|stringr         |1.4.0|
|dplyr           |1.0.8|
|tidyr           |1.2.0|
|readr           |2.1.2|
|webchem         |1.1.3|
|readxl          |1.4.0|

## The script uses data from the following sources:

The empirical dataset contains data from the following sources:
US EPA ECOTOX [1] ASCII-file direct link:
https://gaftp.epa.gov/ecotox/ecotox_ascii_09_15_2022.zip
EFSA pesticide report [2]:
Pierobon, E., Neri, M. C., Marroncelli, S., & Croce, V. (2012). Completion of data entry of pesticide ecotoxicology Tier 1 study endpoints in a XML schema–database. EFSA Supporting Publications, 9(11), 326E

In addition the following file is needed:
"ECOTOX-Term-Appendix-C.csv"" available at:
https://cfpub.epa.gov/ecotox/help.cfm?sub=term-appendix#
Navigate to Appendix C, and "Export as..." to download the CSV version of the file


## The script uses the following QSAR tools:

### ECOSAR
ECOSAR 2.2 [5] application available at:
https://www.epa.gov/tsca-screening-tools/ecological-structure-activity-relationships-ecosar-predictive-model

### VEGA
VEGA 1.1.5 [6] application available at:
https://www.vegahub.eu/portfolio-item/vega-qsar/


### Toxicity Estimation Software Tool (T.E.S.T.)
T.E.S.T. 5.1.1.0 [7] application available at:
https://www.epa.gov/chemical-research/toxicity-estimation-software-tool-test




## Function files:
|Function|Description|
|-----------|------------|
|"ECOTOX_build_function.R"            |A function to build a dataframe from the ECOTOX ascii-download|
|"ECOTOX_cleanup_function.R"          |A function to clean up the ECOTOX dataframe (remove special characters, fixing CAS errors, translating concentrations to mg/l, translating species to correct names and converting durations to hours)|
|"ECOTOX_filter_function.R"           |A function to filter the ECOTOX database based on filter settings|
|"ECOTOX_import_function.R"           |A function that wraps the other three ECOTOX-functions, and serves to prepare the ECOTOX database for the current project|
|"EFSA_cleanup_function.R"            |A function to clean up the EFSA dataframe (fixing CAS errors, translating concentrations to mg/l, translating species to correct names and converting durations to hours)|
|"EFSA_filter_function.R"             |A function to filter the ECOTOX database based on filter settings|
|"EFSA_import_function.R"             |A function that wraps the other three EFSA-functions, and serves to prepare the EFSA data for the current project|
|"QSAR_add_inchikey_function.R"       |A function that wraps the chemical identifier resolver (CIR) query functon for InChIKeys in webchem for the QSAR predictions database project|
|"QSAR_add_smiles_function.R"         |A function that wraps the chemical identifier resolver (CIR) query functon for SMILES in webchem for the QSAR predictions database project, in addition manually curates some SMILES|
|"QSAR_processing_function.R"         |A function that starts, assists with and collects the results of the three QSAR platforms used in the QSAR predictions database project|
|"QSAR_subset_reduction_function.R"   |A function that reduces QSAR predictions to single predictions per subset (or keeps both the raw and calculated predictions, as per arguments)|

## Intermediate files:
|File|Description|
|-----------|------------|
|"ECOTOX_filtered_12.Rda"                       |the ECOTOX database post building, cleanup and filtering|
|"ECOTOX_identifiers_12.Rda"                    |compound metadata from the ECOTOX database|
|"ECOTOX_identifiers_cir_lookup.Rda"              |compound metadata from the ECOTOX database post cir_query|
|"EFSA_CIR_lookup.Rda"                            |compound metadata from the EFSA database post cir_query|
|"EFSA_filtered_12.Rda"                         |the EFSA database post cleanup and filtering|
|"EFSA_identifiers_12.Rda"                      |compound metadata from the EFSA database|
|"experimental_dataset_post_merge_v12.Rda"      |The merged empirical data|
|"identifiers_12.Rda"                           |The chemical identifiers and physicochemical data used in the QSAR predictions datbase project|
|"identifiers_cid_lookup.Rda"                     |The chemical identifiers used in the QSAR predictions datbase project, post collection of CID|
|"identifiers_logkow_pka_lookup.Rda"              |The chemical identifiers used in the QSAR predictions datbase project, post collection of pka and logp (logkow)|
|"identifiers_post_merge_v12.Rda"               |The chemical identifiers used in the QSAR predictions datbase project, post merge|
|"inchikey_lookupfile.Rda"                        |The InChIKeys used in the QSAR predictions datbase project|
|"logp_lookup.Rda"                                |A lookup table for collected logp (logkow)|
|"pka_lookup.Rda"                                 |A lookup table for collected logp (pka)|
|"QSAR_all_processec_v12.Rda"                   |A backup of the processed QSAR predictions, long format|
|"qsar_data_12.Rda"                             |A backup of the "raw" QSAR predictions, long format|







