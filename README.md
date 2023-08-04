# QSAR_predictions_database

This repository contains the script(s) and output data described in the "Dataset on ecotoxicity predictions of 2676 chemicals, using three ecotoxicological quantitative structure-activity relationship platforms" journal article.

These scripts may be used under the Creative Commons Attribution-NonCommercial (CC-BY-NC) license, if you need to cite it please use the CITATION.cff. 


## Below follows a short guide to usaing the script

QSAR_data_collection_script_v11.R is the main script used to produce the dataset. It calls a number of helper functions, located under the subfolder “Functions”. The main script handles directory management and packages, standardization and merging of empirical data, extraction of identifiers from empirical data, gathering of physicochemical properties, output data formatting and output data export.

Specific functions are called for import of empirical data, cleaning and filtering of empirical data, collection of SMILES and InChIKeys, running QSAR platforms, importing and curating QSAR outputs and calculating compounded QSAR predictions. For a full list of function files and a short description of their uses can be found below.

A handful of intermediate files are created by the script to reduce the load on third party APIs and to avoid rerunning time consuming steps of the script, unless instructed to. These files can be found under the “Intermediate files” folder. For a full list of intermediate files and a short description of their purpose can be found below.

The script was written (and is functional) in R version 4.1.3 using the RStudio editor version 2022.12.0 Build 353

###The main script and all helper functions rely on the following R-packages:

|Package|Version
|-----------|------------|
|data.table      |1.14.2|
|stringr         |1.4.0|
|dplyr           |1.0.8|
|tidyr           |1.2.0|
|readr           |2.1.2|
|webchem         |1.1.3|
|readxl          |1.4.0|


The empirical dataset contains data from the following sources:
US EPA ECOTOX [1] ASCII-file direct link:
https://gaftp.epa.gov/ecotox/ecotox_ascii_09_15_2022.zip
EFSA pesticide report [2]:
Pierobon, E., Neri, M. C., Marroncelli, S., & Croce, V. (2012). Completion of data entry of pesticide ecotoxicology Tier 1 study endpoints in a XML schema–database. EFSA Supporting Publications, 9(11), 326E

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
|"ECOTOX_filtered_11.Rda"                       |the ECOTOX database post building, cleanup and filtering|
|"ECOTOX_identifiers_11.Rda"                    |compound metadata from the ECOTOX database|
|"ECOTOX_identifiers_cir_dump.Rda"              |compound metadata from the ECOTOX database post cir_query|
|"EFSA_CIR_dump.Rda"                            |compound metadata from the EFSA database post cir_query|
|"EFSA_filtered_11.Rda"                         |the EFSA database post cleanup and filtering|
|"EFSA_identifiers_11.Rda"                      |compound metadata from the EFSA database|
|"experimental_dataset_post_merge_v11.Rda"      |The merged empirical data|
|"identifiers_11.Rda"                           |The chemical identifiers and physicochemical data used in the QSAR predictions datbase project|
|"identifiers_cid_dump.Rda"                     |The chemical identifiers used in the QSAR predictions datbase project, post collection of CID|
|"identifiers_logkow_pka_dump.Rda"              |The chemical identifiers used in the QSAR predictions datbase project, post collection of pka and logp (logkow)|
|"identifiers_post_merge_v11.Rda"               |The chemical identifiers used in the QSAR predictions datbase project, post merge|
|"inchikey_dumpfile.Rda"                        |The InChIKeys used in the QSAR predictions datbase project|
|"logp_dump.Rda"                                |A lookup table for collected logp (logkow)|
|"pka_dump.Rda"                                 |A lookup table for collected logp (pka)|
|"QSAR_all_processec_v11.Rda"                   |A backup of the processed QSAR predictions, long format|
|"qsar_data_11.Rda"                             |A backup of the "raw" QSAR predictions, long format|







