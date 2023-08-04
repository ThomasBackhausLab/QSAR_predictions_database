# QSAR_predictions_database
A repository for the QSAR prediction database project

These scripts may be used under the Creative Commons Attribution-NonCommercial (CC-BY-NC) license. Below follows a short guide to their usage.

QSAR_data_collection_script_v11.R is the main script used to produce the dataset. It calls a number of helper functions, located under the subfolder “Functions”. The main script handles directory management and packages, standardization and merging of empirical data, extraction of identifiers from empirical data, gathering of physicochemical properties, output data formatting and output data export.

Specific functions are called for import of empirical data, cleaning and filtering of empirical data, collection of SMILES and InChIKeys, running QSAR platforms, importing and curating QSAR outputs and calculating compounded QSAR predictions.

A handful of intermediate files are created by the script to reduce the load on third party APIs and to avoid rerunning time consuming steps of the script, unless instructed to. These files can be found under the “Intermediate files” folder.

The script was written (and is functional) in R version 4.1.3 using the RStudio editor version 2022.12.0 Build 353

The main script and all helper functions rely on the following R-packages:

data.table      1.14.2

stringr         1.4.0

dplyr           1.0.8

tidyr           1.2.0

readr           2.1.2

webchem         1.1.3

readxl          1.4.0


TEMP readme considerations:

Short script guide
The scripts used to produce the dataset are available from the git repository https://github.com/ThomasBackhausLab/QSAR_predictions_database/tree/main
The scripts may be used under the Creative Commons Attribution-NonCommercial (CC-BY-NC) license. Below follows a short guide to their usage.
QSAR_data_collection_script_v11.R is the main script used to produce the dataset. It calls a number of helper functions, located under the subfolder “Functions”. The main script handles directory management and packages, standardization and merging of empirical data, extraction of identifiers from empirical data, gathering of physicochemical properties, output data formatting and output data export.
Specific functions are called for import of empirical data, cleaning and filtering of empirical data, collection of SMILES and InChIKeys, running QSAR platforms, importing and curating QSAR outputs and calculating compounded QSAR predictions.
A handful of intermediate files are created by the script to reduce the load on third party APIs and to avoid rerunning time consuming steps of the script, unless instructed to. These files can be found under the “Intermediate files” folder.
The main script and all helper functions rely on packages described in table 2.






