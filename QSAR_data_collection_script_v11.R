################################################################################
#    A script for collecting and standardizing experimental data and           #
#                           QSAR predictions                                   #
#                                                                              #
################################################################################
#
# Original author: Patrik Svedberg
#
# Based on the output of QSAR models:
# 
# Vega 1.1.5 LC50/EC50/NOEC models
# ECOSAR 2.2
# T.E.S.T 5.1.1.0
# 
# And the empirical databases ECOTOX from the US EPA and the EFSA pesticides
# dataset
#
################################################################################
#                  Table of contents                                           #
#                                                                              #
# 0. Descriptions, disclaimers and datafiles
#   0.1. Versions and changelog
#   0.2. Planned changes and WIP
#   0.3. Input and output files scheme
# 1. Packages, working directories and housekeeping   <- Adapt to your own system here
# 2. Loading empirical data
#   2.1. Importing and filtering EFSA
#   2.2. Importing and filtering USEPA ECOTOX
#   2.3. Merge data and identifiers
#   2.4. Add chemical properties
# 3. Running QSARs
#   3.1. Calculate metadata on QSAR predictions
#   3.2. Fix ECOSAR ChV -> NOEC
# 4. calculating mean/low of "quality" predictions
#   4.1. Reformatting data to wide format "big QSAR database"
#   4.2. Exporting output data
#
################################################################################
#         0. Descriptions, disclaimers and datafiles                           #
################################################################################
#
#
################################################################################
#         0.1. Versions and changelog                                          #
################################################################################

# Version handling in text in this script was implemented for version 4.0
# Any previous versions or changes are not logged here, but can be tracked in
# GIT repository

## Version 4.0 changes

# - Changed some ToC and administrative things (Added changelog, added 
# input/output file descriptions)
# - Changed how pesticide_class works with the identifier merger 
# - Added exp database specific accuracy graphs
# - improved and changed a few graphs

## Version 5.0 changes

# - moved plotting to separate script
# - changed the input version of EFSA to "original" data instead of data processed by Mikael Gustavsson

## Version 6.0 changes

# - incorporated the most recent ecotox and did some minor changes in variables
# - Changed duration handling in EFSA (using a single column: exposureDuration instead of also using expDuration)
# - Added multiple settings for vega calculations

## Version 7.0 changes

# - Rearranged and improved experimental data filtering and processing
# - introduced "forced_rerun" to enable running all of the script even though saved dump files exists
# - defining OECD species at start of script

## Version 8.0 changes

# - Adding scenarios, rearranging to facilitate
# - externalized tons of functions (not chem props)

## Version 9.0 changes

# - Name change from QSARmerger_vX to QSAR_main_script_vX
# - New system for the base frames


## Version 10.0 changes

# - logP filter set to 5


## Version 11.0 changes

# Split off data collection (and QSAR prediction) from analysis
# Scenarios, summaries and analysis now part of a different script (and article)


################################################################################
#         0.2. Planned changes and WIP                                         #
################################################################################



################################################################################
#         0.3. Input and output files scheme                                   #
################################################################################

####### Input files:
#
# ------'Full table.xlsx'------------------
# 
# The EFSA pesticides database
#
#
# ------'ECOTOX Summary_v11_Rev12_09_15_2022.csv'-------------------------------
#
# The US EPA ECOTOX database
#
#
# ------'ECOTOX-Term-Appendix-C.csv'--------------------------------------------
#
# Appendix C from US EPA ECOTOX, containing descriptsions of some terms 
# for the ECOTOX database
#
#
# ------------------------------------------------------------------------------
#
# QSAR outputs are handled within the QSAR_processing_function and will not be
# specified here
#
#
################################################################################


####### Output files:
#
# --------'EFSA_filtered_', version, '.Rda'-------------------------------------
#
# The EFSA data filtered to fit the usage within this script (endpoints matching
# the QSAR models, removal of low quality data points etc)
#
# --------'ECOTOX_combined_filtered_', version, '.Rda'--------------------------
#
# The ECOTOX data filtered to fit the usage within this script (endpoints 
# matching the QSAR models etc)
#
# --------'identifiers_', version, '.Rda'---------------------------------------
#
# A list of all identifiers collected from the original data sources
#
#
# --------'qsar_data_', version, '.Rda'-----------------------------------------
#
# The output from QSAR_processing saved to allow work from systems without
# possibility to run the QSAR_processing function (i.e. no QSARs installed)
#
#
# --------'/identifiers/CAS.txt'------------------------------------------------
#
# The CAS numbers field from the identifiers frame for input into QSAR models
#
#
# --------'/identifiers/SMILES.txt'---------------------------------------------
#
# The SMILES field from the identifiers frame for input into QSAR models
#
#
# --------'/identifiers/CASSMILES.txt'------------------------------------------
#
# The SMILES and CAS numbers fields from the identifiers frame for input into 
# QSAR models
#
#
################################################################################

####### Dump files:
#
# Dump files are produced and/or loaded within the script to reduce script
# times when rerunning the script. In general they contain collected data
# from various APIs which are generally query-capped. They can take anywhere
# from a few minutes to half a day to reproduce
#
#
################################################################################


################################################################################
#         1. Packages, working directories and housekeeping                    #
################################################################################

# A small piece of code to stop "run all" to function here (mainly because sometimes I misclick)
stop('Please do not "run all" in this script - run section by section')

### Clear environment, memory and console
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.



#### Prepare loading packages if this is the first time ####
#library(BiocManager)
#BiocManager::install(c('ChemmineOB', 'ChemmineR'))

#### Load packages ####
packages <- c('readxl', 
              'dplyr',
              'data.table',
              'readr',
              'stringr',
              'webchem',
              "ChemmineR",
              "ChemmineOB",
              'tidyr',
              'data.table',
              'ggplot2',
              'forcats',
              'gridExtra',
              'gridtext',
              'grid',
              'hms',
              'MASS',
              'R.utils')

lapply(packages, require, character.only = TRUE)

#Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jdk1.7.0_51\\jre') # Needed to get Java running for rcdk
library(rcdk)
library(rJava)

#######################

## Set working directory (sub directories will be created by function):

# Git directory
git_directory <- 'C:/Git/QSAR_output_merger'

# Input directory (for input files other than very large databases)
input_directory <- 'C:/Git/QSAR_output_merger/Input'

# output directory (Main outputs)
output_directory <- 'C:/Git/QSAR_output_merger/Output'

# intermediate files directory (Dumpfiles, working files etc)
intermediate_directory <- 'C:/Git/QSAR_output_merger/Intermediate files'


## Paths to raw data (These can be very large files and may not be suited for storage in git directory)
EFSA_filepath = 'E:/Datasets/EFSA/Full table.xlsx'
ECOTOX_filepath = 'E:/Datasets/US_EPA_ECOTOX/ECOTOX 09_15_2022 v8.Rda'


# QSAR directory (function will make subfolders here)
QSAR_directory <- 'C:/Projects/QSARmerger/QSAR'

## Set paths to your QSAR models
vegapath <- 'C:/Program Files/Vega 1.1.5/vega-1.1.5-b48/VEGA.jar'
ecosarpath <- 'C:/Program Files/Ecosar 2.2/ecosarapplication/bin/ecosarapplication64.exe'
testpath <- 'C:/Program Files/TEST 5.1.1.0/TEST.exe' 



# # Load functions in order of first use (if any of these fails, pull a new version of the script)
EFSA_import_function = dget('Functions/EFSA_import_function.R')
EFSA_cleanup_function = dget('Functions/EFSA_cleanup_function.R')
EFSA_filter_function = dget('Functions/EFSA_filter_function.R')
QSAR_add_smiles_function = dget('Functions/QSAR_add_smiles_function.R')
ECOTOX_import_function = dget('Functions/ECOTOX_import_function.R')
ECOTOX_cleanup_function = dget('Functions/ECOTOX_cleanup_function.R')
ECOTOX_filter_function = dget('Functions/ECOTOX_filter_function.R')
QSAR_add_inchikey_function = dget('Functions/QSAR_add_inchikey_function.R')
QSAR_processing_function = dget('Functions/QSAR_processing_function.R')
QSAR_subset_reduction_function = dget('Functions/QSAR_subset_reduction_function.R')


## Set script version
version = 11

## Set forced rerun
# If True   will regenerate all dump-files    (Slow option, has to be done the first time)
# If False  will use dumpfiles if they exist (Faster option)
forced_rerun = F


## Set standard species for fish and algae (synonyms on same rows, unique species on new row)
oecd_fish_species = c('Danio rerio', 'Cyprinus rerio', 'Brachydanio rerio',
                      'Pimephales promelas',
                      'Cyprinus carpio',
                      'Oryzias latipes', 'Poecilia latipes',
                      'Poecilia reticulata', 'Acanthophacelus reticulata', 'Poecilia (Acanthophacelus) reticulata', 'Poecilia latipinna reticulata',
                      'Lepomis macrochirus',
                      'Oncorhynchus mykiss', 'Salmo mykiss', 'Oncorhynchus nerka mykiss', 'Parasalmo mykiss')

oecd_algae_species = c('Pseudokirchneriella subcapitata', 'Selenastrum capricornutum', 'Raphidocelis subcapitata',
                       'Desmodesmus subspicatus', 'Scenedesmus subspicatus')

daphnia_species = c('Daphnia magna', 
                    'Daphnia pulex')


# A small piece of code to stop "run all" to function here, run section by section
stop('Please do not "run all" in this script - run section by section')


################################################################################
#           2. Loading empirical data                                          #
################################################################################

################################################################################
#           2.1. Importing and filtering EFSA                                  #
################################################################################

# Make sure there is a filepath
if(!is.null(EFSA_filepath) & !is.na(EFSA_filepath)){
  
  # Set filters of EFSA data
  EFSA_filters = list('operator' = '=',                          # Set wich operators to include (EFSA_cleanup sets "=" and "range" as options)
                      'reliability' = '1|2',                     # Set which reliability scores to include (regexp, options 1 to 4 separated by |)
                      'active substance' = 'active substance',   # Set how to handle active substances vs formulation (regexp, check EFSA database for alternatives)
                      'species' = c(daphnia_species,             # Set which species to include (vector of species names)
                                    oecd_fish_species, 
                                    oecd_algae_species),
                      'salts' = '\\.',                           # Set regexp for removal of salts (regexp, will remove all data matching the regexp in the SMILES field)
                      'medium' = 'freshwater')                   # Set which mediums to use (regexp, check EFSA database for alternatives)
                                                                 # Additionally, duration filters, effect filters and endpoint filters can be set for the EFSA_filter_function
  
  ## Functions required for EFSA import:
  # EFSA_cleanup_function = dget('EFSA_cleanup_function.R')
  # EFSA_filter_function = dget('EFSA_filter_function.R')
  # QSAR_add_smiles_function = dget('QSAR_add_smiles_function.R')
  
  # Import EFSA
  EFSA_handle = EFSA_import_function(EFSA_filepath,
                                     version = version,
                                     rerun = T,
                                     filters = EFSA_filters,
                                     settings = NULL)
  
  # Separate data from unique identifiers frame
  EFSA_filtered = EFSA_handle[['EFSA_data']]
  EFSA_identifiers = EFSA_handle[['EFSA_identifiers']]
  rm('EFSA_handle')
  
  # Save identifiers separately
  save(EFSA_identifiers, file = paste0(intermediate_directory, '/EFSA_identifiers.Rda'))
  
} else {
  
  print('No EFSA filepath provided, the rest of the script may not run properly without it')
  
}



################################################################################
#           2.2. Importing and filtering USEPA ECOTOX                          #
################################################################################

if(!is.null(ECOTOX_filepath) & !is.na(ECOTOX_filepath)){
  
  # Set ECOTOX dat filters
  ECOTOX_filters = list('operator' = 'None',                            # Set wich operators to include (ECOTOX has 'None' in place of '=')
                        'active substance' = 'active substance',        # Set how to handle active substances vs formulation (regexp, check ECOTOX database for alternatives)
                        'medium' = '(FW)|(NONE)',                       # Set which mediums to use (regexp, check ECOTOX database for alternatives)
                        'species' = c(daphnia_species,                  # Set which species to include (vector of species names)
                                      oecd_fish_species, 
                                      oecd_algae_species),
                        'endpoint' = c('NOEC',                          # Set which endpoints to include (vector of endpoint codes)
                                       'EC50', 
                                       'LC50', 
                                       'IC50', 
                                       'LD50'),
                        settings = NULL)                                # Settings to send to cleanup function
  
  # Set molweights file (if any) for the cleanup function
  ECOTOX_molweights_file = 'Additional data/molweight_dump.Rda'
  
  ## Functions required for ECOTOX import:
  # ECOTOX_cleanup_function = dget('ECOTOX_cleanup_function.R')
  # ECOTOX_filter_function = dget('ECOTOX_filter_function.R')
  # QSAR_add_smiles_function = dget('QSAR_add_smiles_function.R')
  
  # Import ECOTOX
  ECOTOX_handle = ECOTOX_import_function(ECOTOX_filepath,
                                         version = version,
                                         rerun = F,
                                         filters = ECOTOX_filters,
                                         molweight_dump = ECOTOX_molweights_file,
                                         settings = NULL)
  
  # Separate data from unique identifiers frame
  ECOTOX_filtered = ECOTOX_handle[['ECOTOX_data']]
  ECOTOX_identifiers = ECOTOX_handle[['ECOTOX_identifiers']]
  rm('ECOTOX_handle')
  
  
  save(ECOTOX_identifiers, file = 'ECOTOX_identifiers.Rda')
  #load(file = 'ECOTOX_identifiers.Rda')

} else {
  
  print('No ECOTOX filepath provided, the rest of the script may not run properly without it')
  
}
################################################################################
#           2.3. Merge data and identifiers                                    #
################################################################################


### Merge the datasets

# Standardize column names
# colnames(ECOTOX_filtered)
colnames(EFSA_filtered)

ECOTOX_filtered = ECOTOX_filtered %>%
  rename(Species_group = species_group) %>%
  rename(original_CAS = cas_number) %>%
  rename(Endpoint = endpoint) %>%
  rename(Conc_sign = conc1_mean_op) %>%
  rename(Conc_mgL = mgperL) %>%
  rename(Media_type = media_type) %>%
  rename(Species_name = latin_name) %>%
  rename(year = publication_year)

EFSA_filtered = EFSA_filtered %>%
  rename(Endpoint = endpointType) %>%
  rename(Conc_mgL = expValue) %>%
  rename(Conc_sign = expValueOp) %>%
  rename(Media_type = medium) %>%
  rename(Species_name = organism)



# Add empty columns so we can rbind dataframes, for columns not standardizable
ECOTOX_filtered[,colnames(EFSA_filtered)[!colnames(EFSA_filtered) %in% colnames(ECOTOX_filtered)]] <- NA
EFSA_filtered[,colnames(ECOTOX_filtered)[!colnames(ECOTOX_filtered) %in% colnames(EFSA_filtered)]] <- NA


# rbind the datasets
experimental_dataset = rbind(EFSA_filtered, ECOTOX_filtered)



### Merge the identifiers


# Match columns
# colnames(ECOTOX_identifiers)
# colnames(EFSA_identifiers)

ECOTOX_identifiers[,colnames(EFSA_identifiers)[!colnames(EFSA_identifiers) %in% colnames(ECOTOX_identifiers)]] <- NA
EFSA_identifiers[,colnames(ECOTOX_identifiers)[!colnames(ECOTOX_identifiers) %in% colnames(EFSA_identifiers)]] <- NA




# rbind the identifiers
identifiers = rbind(EFSA_identifiers, ECOTOX_identifiers)

# Check total
nrow(identifiers)

# Remove duplicates
identifiers = identifiers %>% distinct(original_CAS, .keep_all = T)

# Check total
nrow(identifiers)

# Remove entries with dots (salts) or empty fields (no SMILES found) or special characters
identifiers = identifiers[!grepl(identifiers$original_SMILES, pattern = '\\.') & !identifiers$original_SMILES == '' & !grepl(identifiers$original_SMILES, pattern = '\\%') & !grepl(identifiers$original_SMILES, pattern = '\\|'),]

# Check total
nrow(identifiers)

# Remove entries without SMILES
identifiers = identifiers[!is.na(identifiers$original_SMILES),]

# Check total
nrow(identifiers)


# Get inchikeys for all compounds (if it isnt there already)
identifiers = QSAR_add_inchikey_function(
  identifiers, 
  settings = NULL,
  local_dumpfile = 'inchikey_dumpfile.Rda')



# Add molecule_number to identifiers
identifiers$molecule_number = seq(1, nrow(identifiers), 1)


## Save the identifiers to files (CAS and SMILES)
# CAS
write.table(identifiers[,"original_CAS"],
            file = paste0(directory, '/identifiers/CAS.txt'),
            col.names = F,
            row.names = F,
            quote = F)
# SMILES
write.table(identifiers[,"original_SMILES"],
            file = paste0(directory, '/identifiers/SMILES.txt'),
            col.names = F,
            row.names = F,
            quote = F)
# CAS and SMILES
write.table(identifiers[,c("original_CAS", "original_SMILES")],
            file = paste0(directory, '/identifiers/CASSMILES.txt'),
            sep = '\t',
            col.names = F,
            row.names = F,
            quote = F)

# Special case for ECOSAR if the identifiers list is more than 1000 entries
if(nrow(identifiers) > 1000){
  
  # Save the first 1k entries
  write.table(identifiers[1:1000,"original_SMILES"],
              file = paste0(directory, '/identifiers/SMILES_1to1000.txt'),
              col.names = F,
              row.names = F,
              quote = F)
  
  # If there are more than 2k entries, save the first 1k to 2k and move on
  if(nrow(identifiers) > 2000){
    
    write.table(identifiers[1001:2000,"original_SMILES"],
                file = paste0(directory, '/identifiers/SMILES_1001to2000.txt'),
                col.names = F,
                row.names = F,
                quote = F)
    
    # If there are more than 3k entries, save the first 2k to 3k and move on
    if(nrow(identifiers) > 3000){
      
      write.table(identifiers[2001:3000,"original_SMILES"],
                  file = paste0(directory, '/identifiers/SMILES_2001to3000.txt'),
                  col.names = F,
                  row.names = F,
                  quote = F)
      
      # Save the 3k to end entries
      # Elongate this if more than 4k entries
      write.table(identifiers[3001:nrow(identifiers),"original_SMILES"],
                  file = paste0(directory, '/identifiers/SMILES_3001toEnd.txt'),
                  col.names = F,
                  row.names = F,
                  quote = F)
      
      
    } else {
      # Not more than 3k entries, save 2k to end
      write.table(identifiers[2001:nrow(identifiers),"original_SMILES"],
                  file = paste0(directory, '/identifiers/SMILES_2001toEnd.txt'),
                  col.names = F,
                  row.names = F,
                  quote = F)
      
    }  
    
  } else {
    # Not more than 2k entries, save 1k to end
    write.table(identifiers[1001:nrow(identifiers),"original_SMILES"],
                file = paste0(directory, '/identifiers/SMILES_1001toEnd.txt'),
                col.names = F,
                row.names = F,
                quote = F)
    
    
  }
  
}

save(experimental_dataset, file = paste0('experimental_dataset_post_merge_v', version, '.Rda'))
save(identifiers, file = paste0('identifiers_post_merge_v', version, '.Rda'))

#load(file = paste0('experimental_dataset_post_merge_v', version, '.Rda'))
#load(file = paste0('identifiers_post_merge_v', version, '.Rda'))




################################################################################
#           2.4. Add chemical properties                                       #
################################################################################

# Get CID from webchem (to query pubchem through webchem for pKa, logkow etc)

# Set save counter, to save every 10 queries (in case webchem throws us out)
j = 0
identifiers_backup = identifiers

# This part was rewritten to be a for loop instead of apply since the limiting factor isnt R's speed issues with for loops, but the API
if(!file.exists(file = 'identifiers_cid_dump.Rda')){
  
  # Add CID as all NA to identifiers
  identifiers$CID = NA
  
  # If we don't have a dupfile/backup we do a fresh round of queries
  for(i in 1:nrow(identifiers)){
    
    # If we have a CID, move on
    if(!is.na(identifiers[i, 'CID'])){
      next
    }
    
    # Get current SMILES
    current_smiles = identifiers[i, "original_SMILES"]
    
    # Collect the query results
    temp = get_cid(current_smiles, from = 'smiles', domain = 'compound', match = 'first')
    
    # Extract CID
    identifiers[i, 'CID'] = temp$cid
    
    # In this loop we save every 10 times we add a new CID due to the volatility of the API
    if(j %% 10 == 0){
      print(paste0('Saving progress. Row number ', i))
      save(identifiers, file = 'identifiers_cid_dump.Rda')
    }
    
    # Save a final time upon completion
    if(i == nrow(identifiers)){
      print('Collection finished, saving...')
      save(identifiers, file = 'identifiers_cid_dump.Rda')
    }
    
    j = j + 1
    
    
  }
  
  
} else {
  # If we have a dumpfile/backup we do things differently
  
  # Do a little jig to avoid overwriting
  temp_identifiers =  identifiers
  
  # Get the dumpfile/backup
  load(file = 'identifiers_cid_dump.Rda')
  
  # Rename dumpfile identifiers
  old_identifiers =  identifiers
  
  # take our current identifiers back
  identifiers =  temp_identifiers
  
  # Remove current identifiers CID
  identifiers$CID = NULL
  
  # Add CID from dumpfile identifiers
  identifiers = merge(identifiers, old_identifiers[,c("original_CAS", 'CID')], by = 'original_CAS', all.x = T)
  
  for(i in 1:nrow(identifiers)){
    
    # For any identifiers not covered in our dupfile, we do the query as in the previous part

    if(!is.na(identifiers[i, 'CID'])){
      next
    }
    
    current_smiles = identifiers[i, "original_SMILES"]
    
    temp = get_cid(current_smiles, from = 'smiles', domain = 'compound', match = 'first')
    
    identifiers[i, 'CID'] = temp$cid
    
    # In this loop we save every 10 times we add a new CID due to the volatility of the API
    if(j %% 10 == 0){
      print(paste0('Saving progress. Row number ', i))
      save(identifiers, file = 'identifiers_cid_dump.Rda')
    }
    
    if(i == nrow(identifiers)){
      print('Collection finished, saving...')
      save(identifiers, file = 'identifiers_cid_dump.Rda')
    }
    
    j = j + 1
    
    
  }
  
}

# Remove any all NA rows
identifiers = identifiers[rowSums(is.na(identifiers)) != ncol(identifiers) - 1,]

# Check the number of NA among CID
sum(is.na(identifiers$CID))

# Fix them by hand if possible using pubchem (searching for CAS instead of SMILES)
print(identifiers[is.na(identifiers$CID),"original_CAS"])


# # This part needs to be run once per new iteration of the dump file, add more as needed
# identifiers[identifiers$original_CAS == "110488-70-5", 'CID'] = '5889665'
# identifiers[identifiers$original_CAS == "16752-77-5", 'CID'] = '5353758'
# identifiers[identifiers$original_CAS == "11141-17-6", 'CID'] = '5281303'
# identifiers[identifiers$original_CAS == "82657-04-3", 'CID'] = '6442842'
# identifiers[identifiers$original_CAS == "131860-33-8", 'CID'] = '3034285'
# identifiers[identifiers$original_CAS == "101007-06-1", 'CID'] = '6436606'
# identifiers[identifiers$original_CAS == "141517-21-7", 'CID'] = '11664966'
# identifiers[identifiers$original_CAS == "51596-11-3", 'CID'] = '9959038'
# identifiers[identifiers$original_CAS == "143390-89-0", 'CID'] = '6112114'
# identifiers[identifiers$original_CAS == "79538-32-2", 'CID'] = '11534837'
# identifiers[identifiers$original_CAS == "131983-72-7", 'CID'] = '6537961'
# identifiers[identifiers$original_CAS == "542-75-6", 'CID'] = '24883'
# identifiers[identifiers$original_CAS == "60-54-8" , 'CID'] = '54675776'
# identifiers[identifiers$original_CAS == "8025-81-8", 'CID'] = '6419898'
# identifiers[identifiers$original_CAS == "34010-21-4" , 'CID'] = '5364711'
# identifiers[identifiers$original_CAS == "4170-30-3"  , 'CID'] = '447466'
# identifiers[identifiers$original_CAS == "80214-83-1" , 'CID'] = '6915744'
# identifiers[identifiers$original_CAS == "153719-23-4", 'CID'] = '5821911'
# identifiers[identifiers$original_CAS == "91465-08-6", 'CID'] = '6440557'
# identifiers[identifiers$original_CAS == "73851-70-4", 'CID'] = '3033889'
# identifiers[identifiers$original_CAS == "17924-92-4" , 'CID'] = '5281576'
# identifiers[identifiers$original_CAS == "117704-25-3" , 'CID'] = '9832750'
# identifiers[identifiers$original_CAS == "141-66-2"    , 'CID'] = '5371560'
# identifiers[identifiers$original_CAS == "25875-51-8"  , 'CID'] = '9570438'
# identifiers[identifiers$original_CAS == "54739-18-3"  , 'CID'] = '5324346'
# identifiers[identifiers$original_CAS == "68085-85-8" , 'CID'] = '5281873'
# identifiers[identifiers$original_CAS == "59-87-0" , 'CID'] = '5447130'
# identifiers[identifiers$original_CAS == "7786-34-7", 'CID'] = '5355863'
# identifiers[identifiers$original_CAS == "470-90-6" , 'CID'] = '10107'
# identifiers[identifiers$original_CAS == "6923-22-4"  , 'CID'] = '5371562'
# identifiers[identifiers$original_CAS == "76703-62-3" , 'CID'] = '6440554'
# identifiers[identifiers$original_CAS == "31218-83-4" , 'CID'] = '5372405'
# identifiers[identifiers$original_CAS == "76703-65-6" , 'CID'] = '6440557'
# identifiers[identifiers$original_CAS == "13171-21-6" , 'CID'] = '25750'
# identifiers[identifiers$original_CAS == "7166-19-0"  , 'CID'] = '6508331'
# identifiers[identifiers$original_CAS == "112-80-1", 'CID'] = '445639'
# identifiers[identifiers$original_CAS == "19902-04-6"  , 'CID'] = '6439694'
# identifiers[identifiers$original_CAS == "361377-29-9"  , 'CID'] = '11048796'
# identifiers[identifiers$original_CAS == "13411-16-0" , 'CID'] = '6436061'
# identifiers[identifiers$original_CAS == "20056-92-2"  , 'CID'] = '5362794'
# identifiers[identifiers$original_CAS == "928-97-2"  , 'CID'] = '5284503'
# identifiers[identifiers$original_CAS == "592-46-1" , 'CID'] = '638071'
# identifiers[identifiers$original_CAS == "928-96-1" , 'CID'] = '5281167'
# identifiers[identifiers$original_CAS == "79-77-6"  , 'CID'] = '638014'
# identifiers[identifiers$original_CAS == "15271-41-7" , 'CID'] = '76970394'
# identifiers[identifiers$original_CAS == "141-05-9" , 'CID'] = '5271566'
# identifiers[identifiers$original_CAS == "94-62-2", 'CID'] = '638024'
# identifiers[identifiers$original_CAS == "94-67-7"     , 'CID'] = '135408751'
# identifiers[identifiers$original_CAS == "107-29-9" , 'CID'] = '5324279'
# identifiers[identifiers$original_CAS == "22910-86-7", 'CID'] = '5281852'
# identifiers[identifiers$original_CAS == "62037-80-3", 'CID'] = '51342034'
# identifiers[identifiers$original_CAS == "101043-37-2", 'CID'] = '445434'
# identifiers[identifiers$original_CAS == "10028-15-6" , 'CID'] = '24823'

# Double check that we are done
sum(is.na(identifiers$CID))
print(identifiers[is.na(identifiers$CID),"original_CAS"])

save(identifiers, file = 'identifiers_cid_dump.Rda')
load(file = 'identifiers_cid_dump.Rda')


### Get logp and pka

# Make empty frames first time
logp_frame = data.frame('CID' = NA, 'Name' = NA, 'Result' = NA, 'SourceName' = NA, 'SourceID' = NA)
pka_frame = data.frame('CID' = NA, 'Name' = NA, 'Result' = NA, 'SourceName' = NA, 'SourceID' = NA)

# If we have run this script before we load the output
if(file.exists('logp_dump.Rda') & file.exists('pka_dump.Rda')){
  
  load(file = 'logp_dump.Rda')
  load(file = 'pka_dump.Rda')
  
}

for(i in 1:nrow(identifiers)){  
  
  current_CID = identifiers[i, "CID"]
  
  # If we have processed this CID already go to next
  if(current_CID %in% logp_frame$CID){
    
    next
    
  }
  
  temp_logp = pc_sect(current_CID, section = c('logP'), domain = 'compound')
  
  if(ncol(temp_logp) > 0){
    
    logp_frame = rbind(logp_frame, temp_logp)
    
  } else {
    
    temp_logp2 = pc_prop(current_CID, properties = c('XlogP'))
    
    if(ncol(temp_logp2) > 1){
      
      temp_logp_row = data.frame('CID' = temp_logp2$CID, 'Name' = NA, 'Result' = temp_logp2$XLogP, 'SourceName' = NA, 'SourceID' = NA)
      
      logp_frame = rbind(logp_frame, temp_logp_row)
      
    } else {
      
      temp_logp_row = data.frame('CID' = current_CID, 'Name' = NA, 'Result' = NA, 'SourceName' = NA, 'SourceID' = NA)
      
      logp_frame = rbind(logp_frame, temp_logp_row)
      
    }
    
    
  }
  
  temp_pKa = pc_sect(current_CID, section = c('Dissociation Constants'), domain = 'compound')
  
  
  if(ncol(temp_pKa) > 0){
    
    pka_frame = rbind(pka_frame, temp_pKa)
    
  }
  
  # Every one hundred compounds we save the results
  if(i %% 100 == 0){
    
    save(logp_frame, file = 'logp_dump.Rda')
    save(pka_frame, file = 'pka_dump.Rda')
    
  }
  
}

# Drop the NA rows in the beginning and do a last save
logp_frame = logp_frame[rowSums(is.na(logp_frame)) != ncol(logp_frame),]
pka_frame = pka_frame[rowSums(is.na(pka_frame)) != ncol(pka_frame),]
save(logp_frame, file = 'logp_dump.Rda')
save(pka_frame, file = 'pka_dump.Rda')



# Unpack the contents of the pka and logp frames into the identifiers

identifiers$logkow = NA
identifiers$logkow_source = NA
identifiers$pka = NA

# Create a list to "fix by hand" for entries with multiple values for different pH etc
temp_logp_fix_by_hand = list()


if(!file.exists(file = 'identifiers_logkow_pka_dump.Rda') | forced_rerun){
  
  for(i in 1:nrow(identifiers)){  
    
    current_CID = identifiers[i, "CID"]
    
    # Get any entries matching the CID from logp and pka frames
    
    temp_logp = logp_frame[logp_frame$CID == current_CID,]
    temp_pka = pka_frame[pka_frame$CID == current_CID,]
    
    
    ###### Process logP results
    
    if(nrow(temp_logp) > 1){
      # If mutluiple logP values
      
      # Gather the values
      temp_values = tryCatch(
        {as.numeric(sub(unlist(str_extract_all(temp_logp$Result, pattern = '(?<=([Ll]og Kow ?[=:] ?)|^|([;,] )|(and ))(?<!pH )-?[0-9]{1,2}([\\.,][0-9]{1,3})?(?=( )|($))(?!( °C)|( C)|( at pH [^67])|( \\(pH [^67]))')), pattern = ',', replacement = '.'))
        },
        error=function(cond){
          message(paste("Error logp in the following CID: ", current_CID))
          message(temp_logp)
          message(cond)
          return(NA)
        },
        warning=function(cond){
          message(paste("Warning logp in the following CID: ", current_CID))
          message(temp_logp)
          message(cond)
          return(NA)
        }
      )
      
      # Take median value (use two decimals)
      temp_median_logp = format(round(median(temp_values), digits = 2), nsmall = 2)
      
      # Save median
      identifiers[i, "logkow"] = temp_median_logp
      
      # If all are the same, and there is a sourceID (which isnt there for computed)
      if((var(temp_values, na.rm = T) == 0)|(length(temp_values) == 1) & !all(is.na(temp_logp$SourceID))){
        
        identifiers[i, "logkow_source"] = 'experimental'
        
      } else if(!all(is.na(temp_logp$SourceID))){
        
        identifiers[i, "logkow_source"] = 'median experimental'  
        
      } else {
        
        identifiers[i, "logkow_source"] = 'computed'  
        
      }
      
      
    } else if(nrow(temp_logp) == 0){
      
      identifiers[i, "logkow"] = NA
      
    } else if(nrow(temp_logp) == 1){
      # If only one result, save it
      
      temp_values = tryCatch(
        {as.numeric(sub(unlist(str_extract_all(temp_logp$Result, pattern = '(?<=([Ll]og Kow ?[=:] ?)|^|([;,] )|(and ))(?<!pH )-?[0-9]{1,2}([\\.,][0-9]{1,3})?(?=( )|($))(?!( °C)|( C)|( at pH [^67])|( \\(pH [^67]))')), pattern = ',', replacement = '.'))
        },
        error=function(cond){
          message(paste("Error logp in the following CID: ", current_CID))
          message(temp_pka)
          message(cond)
          return(NA)
        },
        warning=function(cond){
          message(paste("Warning logp in the following CID: ", current_CID))
          message(temp_pka)
          message(cond)
          return(NA)
        }
      )
      
      
      if(length(temp_values) == 0){
        
        identifiers[i, "logkow"] = NA
        
      } else {
        
        identifiers[i, "logkow"] = format(round(median(temp_values), digits = 2), nsmall = 2)
        
        if(!is.na(temp_logp$SourceID)){
          
          identifiers[i, "logkow_source"] = 'experimental'  
          
        } else {
          
          identifiers[i, "logkow_source"] = 'computed'  
          
        } 
      }
    } 
    
    
    # Process pKa results (here we can have multiple values per line)
    if(nrow(temp_pka) > 1){
      # If mutluiple lines
      
      # Gather the values
      temp_values = tryCatch(
        {unlist(str_extract_all(temp_pka$Result, pattern = '(?<=(^ ?)|([,:;=]) |(and )|(p[Kk]a[1234]? i?s? ?))(?<![^p][Kk]a?[1234]= ?)-?[0-9]{1,2}([\\.,][0-9]{1,3})?(?!( °C)|( C)|(X10)|[:\\.])'))
        },
        error=function(cond){
          message(paste("Error pka in the following CID: ", current_CID))
          message(temp_pka)
          message(cond)
          return(NA)
        },
        warning=function(cond){
          message(paste("Warning pka in the following CID: ", current_CID))
          message(temp_pka)
          message(cond)
          return(NA)
        }
      )
      
      if(length(temp_values) > 0){
        
        # Remove mutliples of the same value and sort by value
        
        temp_values = tryCatch(
          {sort(as.numeric(unique(temp_values)), decreasing = F)
          },
          error=function(cond){
            message(paste("Error numeric in the following CID: ", current_CID))
            message(temp_values)
            message(cond)
            return(NA)
          },
          warning=function(cond){
            message(paste("Warning numeric in the following CID: ", current_CID))
            message(temp_values)
            message(cond)
            return(NA)
          }
        )
        
        # Format into two decimals
        temp_values = unlist(sapply(temp_values, FUN = function(x){format(round(median(x), digits = 2), nsmall = 2)}))
        
        # Remove duplicates again after rounding
        temp_values = unique(temp_values)
        
        # paste them pipe separated for later processing
        temp_multiple_pka = paste(temp_values, collapse = '|')
        
        # Save median
        identifiers[i, "pka"] = temp_multiple_pka
        
        
      }
      
    } else if(ncol(temp_pka) == 0){
      # If no result save NA
      
      identifiers[i, "pka"] = NA
      
    } else if(nrow(temp_pka) == 1){
      # Sometimes multiple pKa values come in a single string, we split that up
      temp_values = tryCatch(
        {unlist(str_extract_all(temp_pka$Result, pattern = '(?<=(^ ?)|([,:;=]) |(and )|(p[Kk]a[1234]? i?s? ?))-?[0-9]{1,2}([\\.,][0-9]{1,3})?(?!( °C)|( C)|(X10)|[:\\.])'))
        },
        error=function(cond){
          message(paste("Error pka in the following CID: ", current_CID))
          message(temp_pka)
          message(cond)
          break
          return(NA)
        },
        warning=function(cond){
          message(paste("Warning pka in the following CID: ", current_CID))
          message(temp_pka)
          message(cond)
          break
          return(NA)
        }
      )
      
      if(length(temp_values) > 0){
        
        # Remove mutliples of the same value
        temp_values = tryCatch(
          {sort(as.numeric(unique(temp_values)), decreasing = F)
          },
          error=function(cond){
            message(paste("Error numeric in the following CID: ", current_CID))
            message(temp_values)
            message(cond)
            return(NA)
          },
          warning=function(cond){
            message(paste("Warning numeric in the following CID: ", current_CID))
            message(temp_values)
            message(cond)
            return(NA)
          }
        )
        
        # Format into two decimals
        temp_values = unlist(sapply(temp_values, FUN = function(x){format(round(median(x), digits = 2), nsmall = 2)}))
        
        # Remove duplicates again after rounding
        temp_values = unique(temp_values)
        
        if(length(temp_values) > 1){
          
          identifiers[i, "pka"] = paste(temp_values, collapse = '|')
          
        } else {
          
          # If only one result, save it
          identifiers[i, "pka"] = temp_values
          
        }
      }
    }
    
    if(is.na(identifiers[i, "pka"]) & nrow(temp_pka) > 0){
      
      print('pka results:')
      print(temp_pka)
      print('pKa saved:')
      print(identifiers[i, "pka"])
      
      
    }
    
  }
  
  
  ## Fix some logP troublemakers by hand
  
  # Spiromesifen log Kow = 4.55 (pH 2 and 7.5) Hazardous Substances Data Bank (HSDB)
  identifiers[identifiers$CID == '9907412', "logkow"] = 4.55
  identifiers[identifiers$CID == '9907412', "logkow_source"] = 'experimental'  
  
  # Mesotrione log Kow (at 20 °C): 0.11 (unbuffered water), 0.90 (pH 5), <-1.0 (pH 7 and pH 9) Hazardous Substances Data Bank (HSDB)
  identifiers[identifiers$CID == '175967', "logkow"] = 0.7        # Picking computed value
  identifiers[identifiers$CID == '175967', "logkow_source"] = 'computed'  
  
  # Kathon 930 log Kow = 2.8, also reported as 4.5 Hazardous Substances Data Bank (HSDB)
  identifiers[identifiers$CID == '91688', "logkow"] = 2.8
  identifiers[identifiers$CID == '91688', "logkow_source"] = 'experimental'  
  
  # Cybutryne log Kow = 3.95; also reported as 2.8 Hazardous Substances Data Bank (HSDB)
  identifiers[identifiers$CID == '91590', "logkow"] = 3.95
  identifiers[identifiers$CID == '91590', "logkow_source"] = 'experimental'  
  
  # Tralomethrin log Kow = approximately 5 at 25 °C Hazardous Substances Data Bank (HSDB)
  identifiers[identifiers$CID == '48132', "logkow"] = 5.00
  identifiers[identifiers$CID == '48132', "logkow_source"] = 'experimental'  
  
  # Sulfluramid log Kow >6.80 (un-ionized) Hazardous Substances Data Bank (HSDB)
  identifiers[identifiers$CID == '77797', "logkow"] = 6.80
  identifiers[identifiers$CID == '77797', "logkow_source"] = 'experimental'  
  
  # Toxaphene The median log Kow for toxaphene is 5.90. Hazardous Substances Data Bank (HSDB)
  identifiers[identifiers$CID == '5284469', "logkow"] = 5.90
  identifiers[identifiers$CID == '5284469', "logkow_source"] = 'experimental'  
  
  # 2-Ethylhexanal log Kow = 3.07, measured OECD Method 107 Hazardous Substances Data Bank (HSDB)
  identifiers[identifiers$CID == '31241', "logkow"] = 3.07
  identifiers[identifiers$CID == '31241', "logkow_source"] = 'experimental'  
  
  # 2,4,5-T butoxyethyl ester log Kow of 5.32 (est) Hazardous Substances Data Bank (HSDB)
  identifiers[identifiers$CID == '17349', "logkow"] = 5.32
  identifiers[identifiers$CID == '17349', "logkow_source"] = 'experimental'  
  
  
  
  ## Finally some I was not able to fix
  
  # Triethylenetetramine -1.4/-1.66 ILO International Chemical Safety Cards (ICSC)
  identifiers[identifiers$CID == '5565', "logkow"] = NA         # Unsure how to handle this
  identifiers[identifiers$CID == '5565', "logkow_source"] = 'experimental'  
  
  # Naphthenic acids 5->6 (calculated) ILO International Chemical Safety Cards (ICSC)
  identifiers[identifiers$CID == '20849290', "logkow"] = NA         # Unsure how to handle this
  identifiers[identifiers$CID == '20849290', "logkow_source"] = 'experimental'  
  
  # Ethylene glycol diacetate Log Kow= 0.10/0.38 (calc) Hazardous Substances Data Bank (HSDB)
  identifiers[identifiers$CID == '8121', "logkow"] = NA         # Unsure how to handle this
  identifiers[identifiers$CID == '8121', "logkow_source"] = 'experimental'  
  
  
  
  
  
  ## Fix some pka troublemakers by hand
  
  # Yellow prisms from water. Becomes anhydrous at 140 °C, decomposes at 313-314 °C. UV max (0.1N NaOH): 230,310 nm (epsilon 1400, 19600); (0.1N HCl): 222, 237 nm (epsilon 9240, 21300); (methanol): 216, 239 nm (epsilon 8940, 19300). pKa1 7.77, pKa2 11.7. Insoluble in water, acetone, ether. Soluble in hot ethanol, in alkaline solutions with slow decomposition. /6-Mercaptopurine monohydrate/
  identifiers[identifiers$CID == '667490', "pka"] = '7.77|11.7'
  
  # The pKa of the conjugate acid of propiconazole is 1.09. Hazardous Substances Data Bank (HSDB)
  identifiers[identifiers$CID == '43234', "pka"] = '1.09'
  
  # PKA 4.62 Hazardous Substances Data Bank (HSDB)
  identifiers[identifiers$CID == '6950', "pka"] = '4.62'
  
  # 2,3,4,5-Tetrachlorophenol pKa=6.35 Hazardous Substances Data Bank (HSDB)     6765
  identifiers[identifiers$CID == '21013', "pka"] = '6.35'
  
  save(identifiers, file = 'identifiers_logkow_pka_dump.Rda')
  
} else {
  
  load(file = 'identifiers_logkow_pka_dump.Rda')
  
}

identifiers$logkow = ifelse(identifiers$logkow == 'NA', NA, identifiers$logkow)

# Some substances we fail to get lokkow for
nrow(identifiers[is.na(identifiers$logkow),])

# Save all identifiers, for work on other systems
save(identifiers, file = paste0('identifiers_', version, '.Rda'))
load(file = paste0('identifiers_', version, '.Rda'))


################################################################################
#             3. Running QSAR_processing                                       #
################################################################################

# Option to turn off running the QSARS (the output files are already there, and you just want to recompile them)
run_QSARs = TRUE

QSAR_processing_function = dget('Functions/QSAR_processing_function.R')

# Run the QSARs using the processing script. More infor in the script.
QSAR_output = QSAR_processing_function(identifiers = identifiers, 
                                       QSAR_paths = c(vegapath, ecosarpath, testpath),
                                       working_directory = directory,
                                       run_QSARs = run_QSARs,
                                       format = 'long',
                                       ecosar_fix = T)

# Save the merged dataframe to avoid having to reprocess it
save(QSAR_output, file = paste0('qsar_data_', version, '.Rda'))

load(file = paste0('qsar_data_', version, '.Rda'))



################################################################################
#               4. calculating aggregated prediction for QSAR platforms        #
################################################################################

# Load QSAR_output dataframe to avoid extra work
load(file = paste0('qsar_data_', version, '.Rda'))

QSAR_daphnia_acute = QSAR_subset_reduction_function(QSAR_output, 
                                                    species = 'daphnia', 
                                                    in_format = "long",
                                                    out_format = "long",
                                                    endpoint = 'acute',
                                                    vega_filter = c('good', 'moderate'),
                                                    vega_filter2 = c('low'),
                                                    vega_method = c('geo_mean', 'low'), 
                                                    vega_method2 = c('geo_mean', 'low'), 
                                                    ecosar_method = c('geo_mean', 'drop baseline'),
                                                    ecosar_method2 = c('low', 'drop baseline'),
                                                    data_filter = 'all')

QSAR_daphnia_chronic = QSAR_subset_reduction_function(QSAR_output, 
                                                      species = 'daphnia', 
                                                      in_format = "long",
                                                      out_format = "long",
                                                      endpoint = 'chronic',
                                                      vega_filter = c('good', 'moderate'),
                                                      vega_filter2 = c('low'),
                                                      vega_method = c('geo_mean', 'low'), 
                                                      vega_method2 = c('geo_mean', 'low'), 
                                                      ecosar_method = c('geo_mean', 'drop baseline'),
                                                      ecosar_method2 = c('low', 'drop baseline'),
                                                      data_filter = 'all')



QSAR_fish_acute = QSAR_subset_reduction_function(QSAR_output, 
                                                 species = 'fish', 
                                                 in_format = "long",
                                                 out_format = "long",
                                                 endpoint = 'acute',
                                                 vega_filter = c('good', 'moderate'),
                                                 vega_filter2 = c('low'),
                                                 vega_method = c('geo_mean', 'low'), 
                                                 vega_method2 = c('geo_mean', 'low'), 
                                                 ecosar_method = c('geo_mean', 'drop baseline'),
                                                 ecosar_method2 = c('low', 'drop baseline'),
                                                 data_filter = 'all')


QSAR_fish_chronic = QSAR_subset_reduction_function(QSAR_output, 
                                                   species = 'fish', 
                                                   in_format = "long",
                                                   out_format = "long",
                                                   endpoint = 'chronic',
                                                   vega_filter = c('good', 'moderate'),
                                                   vega_filter2 = c('low'),
                                                   vega_method = c('geo_mean', 'low'), 
                                                   vega_method2 = c('geo_mean', 'low'), 
                                                   ecosar_method = c('geo_mean', 'drop baseline'),
                                                   ecosar_method2 = c('low', 'drop baseline'),
                                                   data_filter = 'all')


QSAR_algae_acute = QSAR_subset_reduction_function(QSAR_output, 
                                                  species = 'algae', 
                                                  in_format = "long",
                                                  out_format = "long",
                                                  endpoint = 'acute',
                                                  vega_filter = c('good', 'moderate'),
                                                  vega_filter2 = c('low'),
                                                  vega_method = c('geo_mean', 'low'), 
                                                  vega_method2 = c('geo_mean', 'low'), 
                                                  ecosar_method = c('geo_mean', 'drop baseline'),
                                                  ecosar_method2 = c('low', 'drop baseline'),
                                                  data_filter = 'all')

QSAR_algae_chronic = QSAR_subset_reduction_function(QSAR_output, 
                                                    species = 'algae', 
                                                    in_format = "long",
                                                    out_format = "long",
                                                    endpoint = 'chronic',
                                                    vega_filter = c('good', 'moderate'),
                                                    vega_filter2 = c('low'),
                                                    vega_method = c('geo_mean', 'low'), 
                                                    vega_method2 = c('geo_mean', 'low'), 
                                                    ecosar_method = c('geo_mean', 'drop baseline'),
                                                    ecosar_method2 = c('low', 'drop baseline'),
                                                    data_filter = 'all')

QSAR_all = do.call("rbind", list(QSAR_daphnia_acute, QSAR_daphnia_chronic, QSAR_fish_acute, QSAR_fish_chronic, QSAR_algae_acute, QSAR_algae_chronic))



# Sometimes some ChV stays, but we fixed that, lets just rename them to NOEC here
QSAR_all$model_endpoint = ifelse(QSAR_all$model_endpoint == 'NOEC ChV', 'NOEC', QSAR_all$model_endpoint)

save(QSAR_all, file = paste0('QSAR_all_processec_v', version, '.Rda'))
#load(file = paste0('QSAR_all_processec_v', version, '.Rda'))


################################################################################
#               4.1. Reformatting data to wide format "big QSAR database"      #
################################################################################

# ECOSAR outputs some alternative salts for some substances, which makes pivot
# not work correctly, we remove these

# Save an "ld" version to see how much was removed in the process
ECOSAR_subset_old = QSAR_all[QSAR_all$QSAR_tool == 'ECOSAR' & is.na(QSAR_all$reliability),]

ECOSAR_subset = distinct_at(QSAR_all[QSAR_all$QSAR_tool == 'ECOSAR' & is.na(QSAR_all$reliability),], .vars = c('original_SMILES', 'model', 'value'), .keep_all = T)

ECOSAR_CAS = unique(ECOSAR_subset$original_CAS)

# Loop over CAS to fix them
for(i in 1:length(ECOSAR_CAS)){
  
  # Extract CAS
  current_CAS = ECOSAR_CAS[i]
  
  # Subset non-calculated values for aid CAS
  current_subset = ECOSAR_subset[ECOSAR_subset$original_CAS == current_CAS,]
  
  # Distinct at input SMILES, model and value fields
  temp_frame = distinct_at(current_subset, .vars = c('original_SMILES', 'model', 'value'), .keep_all = T)
  
  # Chechk if we have duplicates
  if(any(duplicated(temp_frame$model))){
    
    #print(paste0('CAS with an issue: ', current_CAS))
    
    # If duplicates fixed by CAS identification (most of the time)
    if(current_CAS %in% as.cas(unique(temp_frame$ecosar_CAS))){
      
      #print(paste0('CAS matched, options: ', paste(unique(as.cas(unique(temp_frame$ecosar_CAS))), collapse = ', ')))
      
      # Remove old rows
      ECOSAR_subset = ECOSAR_subset[!(ECOSAR_subset$original_CAS == current_CAS & is.na(ECOSAR_subset$reliability)),]
      
      # Add new rows
      ECOSAR_subset = rbind(ECOSAR_subset, temp_frame[temp_frame$original_CAS == as.cas(temp_frame$ecosar_CAS),])
      
    } else{
      
      print(paste0('CAS mismatch, current_SMILES: ', unique(unique(temp_frame$original_SMILES))))
      print(paste0('Matching metals, SMILES options: ', paste(unique(unique(temp_frame$ecosar_SMILES)), collapse = ', ')))
      
      # Special cases where depreciated CAS are used
      
      # Find which metal (if any) is causing the issue and match it with ecosar SMILES
      current_metal = str_extract(unique(temp_frame$original_SMILES), pattern = '\\[[A-Za-z]{1,2}\\]')
      
      if(!is.na(current_metal)){
        
        print(paste0('Chosen SMILES: ', unique(temp_frame[grepl(temp_frame$ecosar_SMILES, pattern = current_metal, fixed = T), 'ecosar_SMILES'])))
        
        # Use this metal to get the right rows from temp_frame
        
        # Remove old rows
        ECOSAR_subset = ECOSAR_subset[!(ECOSAR_subset$original_CAS == current_CAS & is.na(ECOSAR_subset$reliability)),]
        
        # Add new rows
        ECOSAR_subset = rbind(ECOSAR_subset, temp_frame[grepl(temp_frame$ecosar_SMILES, pattern = current_metal, fixed = T),])
        
        
      } else if(unique(temp_frame$original_SMILES) %in% unique(temp_frame$ecosar_SMILES)){
        
        # Then we check if we have matching SMILES instead
        print(paste0('Chosen SMILES: ', unique(temp_frame[temp_frame$ecosar_SMILES == unique(temp_frame$original_SMILES), 'ecosar_SMILES'])))
        
        # Remove old rows
        ECOSAR_subset = ECOSAR_subset[!(ECOSAR_subset$original_CAS == current_CAS & is.na(ECOSAR_subset$reliability)),]
        
        # Add new rows
        ECOSAR_subset = rbind(ECOSAR_subset, temp_frame[temp_frame$ecosar_SMILES == unique(temp_frame$original_SMILES),])
        
      } else if(current_CAS == '72-20-8'){
        
        # Then we check if we have matching SMILES instead
        print(paste0('Manually fixing HEXADRIN/ENDRIN mismatch'))
        
        # Remove old rows
        ECOSAR_subset = ECOSAR_subset[!(ECOSAR_subset$original_CAS == current_CAS & is.na(ECOSAR_subset$reliability)),]
        
        # Add new rows
        ECOSAR_subset = rbind(ECOSAR_subset, temp_frame[temp_frame$ecosar_CAS == '128109',])
        
      } else {
        print('Something else')
        break
        
      }
      
      
    }
    
    
  }
  
}

# Remove all ECOSAR predictions
QSAR_all_deduped = QSAR_all[!(QSAR_all$QSAR_tool == 'ECOSAR' & is.na(QSAR_all$reliability)),]

# Add the new filtered down predictions
QSAR_all_deduped = rbind(QSAR_all_deduped, ECOSAR_subset)

# Remove EXPERIMENTAL T.E.S.T. predictions
QSAR_all_deduped = QSAR_all_deduped[!(QSAR_all_deduped$QSAR_tool == 'T.E.S.T.' & QSAR_all_deduped$reliability == 'EXPERIMENTAL'),]

# Final dedupe
QSAR_all_deduped = distinct_at(QSAR_all_deduped, .vars = c('original_SMILES', 'model', 'value'), .keep_all = T)

# We also need to fix VEGA reliabilities, so we remove VEGA initially
QSAR_all_deduped_no_vega = QSAR_all_deduped[!(QSAR_all_deduped$QSAR_tool == 'VEGA' & QSAR_all_deduped$reliability != 'calculated'),]
QSAR_all_deduped_vega_raw_only = QSAR_all_deduped[QSAR_all_deduped$QSAR_tool == 'VEGA' & QSAR_all_deduped$reliability != 'calculated',]

# We pivot VEGA
QSAR_vega_wide = as.data.frame(pivot_wider(QSAR_all_deduped_vega_raw_only, names_from = c(model), values_from = c(value, reliability), id_cols = original_CAS, names_glue = '{model}_{.value}'))

# Rearrange the columns
QSAR_vega_wide = QSAR_vega_wide[,order(colnames(QSAR_vega_wide))]

# We do this based on CAS due to HEXADRIN and DIELDRIN having the same SMILES :O :O but ECOSAR manages to provide predictions for the correct substances...
QSAR_all_wide = as.data.frame(pivot_wider(QSAR_all_deduped_no_vega, names_from = c(model), values_from = c(value), id_cols = original_CAS, ))

# Add VEGA
QSAR_all_wide = merge(QSAR_all_wide, QSAR_vega_wide, by = 'original_CAS')



## Fix column names
QSAR_all_wide = QSAR_all_wide[, order(colnames(QSAR_all_wide))]

QSAR_all_wide = QSAR_all_wide %>%
  relocate(original_CAS)

# Add identifier info
QSAR_all_wide = merge(QSAR_all_wide, identifiers, by = 'original_CAS')

QSAR_all_wide <- QSAR_all_wide %>%
  dplyr::select(colnames(identifiers), everything())

# Remove column number (if its in there)
QSAR_all_wide = QSAR_all_wide[,!colnames(QSAR_all_wide) %in% c('molecule_number')]

# Save colnames if something goes wrong
colnames_old = colnames(QSAR_all_wide)

# Change column names
colnames(QSAR_all_wide) = str_replace_all(colnames(QSAR_all_wide), pattern = '^(?![VTE])', replacement = 'META_')
colnames(QSAR_all_wide) = str_replace_all(colnames(QSAR_all_wide), pattern = '(?<=ECOSAR_)(?=[afd])', replacement = 'calculated_')
colnames(QSAR_all_wide) = str_replace_all(colnames(QSAR_all_wide), pattern = '(?<=ECOSAR_)(?=[A-Z])', replacement = 'raw_') 
colnames(QSAR_all_wide) = str_replace_all(colnames(QSAR_all_wide), pattern = '(?<=VEGA_)(?=[afd])', replacement = 'calculated_')
colnames(QSAR_all_wide) = str_replace_all(colnames(QSAR_all_wide), pattern = '(?<=VEGA_)(?=[A-Z])', replacement = 'raw_') 
colnames(QSAR_all_wide) = str_replace_all(colnames(QSAR_all_wide), pattern = '(?<=TEST_)', replacement = 'raw_') 

# Replace all NA with -7777
QSAR_all_wide[is.na(QSAR_all_wide)] = -7777

# Add internal database identifier for each compound
QSAR_all_wide$META_QSARn = paste0('QSARn', as.numeric(rownames(QSAR_all_wide)))

# Move this col first
QSAR_all_wide = QSAR_all_wide %>%
  relocate(META_InChIKey, .after = META_original_SMILES) %>%
  relocate(META_QSARn)


################################################################################
#               4.2. Exporting output data                                     #
################################################################################


# Save the wide format dataframe in tab separated csv
write.table(QSAR_all_wide, file = paste0('QSAR_all_wide_v', version, '.csv'), sep = '\t', col.names = T, row.names = F, quote = F)

# Double check save
# QSAR_all_wide_new = fread(file = paste0('QSAR_all_wide_v', version, '.csv'), sep = '\t')















################################################################################
#           2.4. Calculate overlap between empirical sets                      #
################################################################################

species_vector = c('daphnia',
                   'fish',
                   'algae')

substance_overview_frame = data.frame(species = species_vector,
                                      total_chemicals_unfiltered = NA,
                                      total_chemicals_filtered = NA,
                                      EFSA_chemicals_unfiltered = NA,
                                      EFSA_chemicals_filtered = NA,
                                      ECOTOX_chemicals_unfiltered = NA,
                                      ECOTOX_chemicals_filtered = NA,
                                      overlap_unfiltered = NA,
                                      overlap_filtered = NA)

# For this part we need the uncurated version of the data
EFSA = as.data.frame(read_xlsx(EFSA_filepath))
load(ECOTOX_filepath)

for(i in 1:length(species_vector)){
  
  current_species = species_vector[i]
  
  if(current_species == 'daphnia'){
    
    current_species_vector = daphnia_species
    
  } else if(current_species == 'fish'){
    
    current_species_vector = oecd_fish_species
    
  } else if(current_species == 'algae'){
    
    current_species_vector = oecd_algae_species
    
  } 
  
  current_total_chemicals_unfiltered = length(unique(append(EFSA[EFSA$organism %in% current_species_vector, "CASNO"],
                                                            ECOTOX[ECOTOX$latin_name %in% current_species_vector, "cas_number"])))
  
  current_total_chemicals_filtered = length(unique(append(EFSA_filtered[EFSA_filtered$Species_name %in% current_species_vector, "original_CAS"],
                                                          ECOTOX_filtered[ECOTOX_filtered$Species_name %in% current_species_vector, "original_CAS"])))
  
  current_EFSA_chemicals_unfiltered = length(unique(EFSA[grepl(EFSA$organism, pattern = paste0('(',paste(current_species_vector, collapse = ')|('), ')')), "CASNO"]))
  
  current_EFSA_chemicals_filtered = length(unique(EFSA_filtered[EFSA_filtered$Species_name %in% current_species_vector, "original_CAS"]))
  
  current_ECOTOX_chemicals_unfiltered = length(unique(ECOTOX[grepl(ECOTOX$latin_name, pattern = paste0('(',paste(current_species_vector, collapse = ')|('), ')')), "cas_number"]))
  
  current_ECOTOX_chemicals_filtered = length(unique(ECOTOX_filtered[ECOTOX_filtered$Species_name %in% current_species_vector, "original_CAS"]))
  
  
  
  current_overlap_unfiltered = length(unique(EFSA[grepl(EFSA$organism, pattern = paste0('(',paste(current_species_vector, collapse = ')|('), ')')) & EFSA$CASNO %in% ECOTOX[ECOTOX$latin_name %in% current_species_vector, "cas_number"], 'CASNO']))
  
  current_overlap_filtered = length(unique(EFSA_filtered[EFSA_filtered$Species_name %in% current_species_vector & EFSA_filtered$original_CAS %in% ECOTOX_filtered[ECOTOX_filtered$Species_name %in% current_species_vector, "original_CAS"], 'original_CAS']))
  
  
  substance_overview_frame[substance_overview_frame$species == current_species, "total_chemicals_unfiltered"] = current_total_chemicals_unfiltered
  substance_overview_frame[substance_overview_frame$species == current_species, "total_chemicals_filtered"] = current_total_chemicals_filtered
  substance_overview_frame[substance_overview_frame$species == current_species, "EFSA_chemicals_unfiltered"] = current_EFSA_chemicals_unfiltered
  substance_overview_frame[substance_overview_frame$species == current_species, "EFSA_chemicals_filtered"] = current_EFSA_chemicals_filtered
  substance_overview_frame[substance_overview_frame$species == current_species, "ECOTOX_chemicals_unfiltered"] = current_ECOTOX_chemicals_unfiltered
  substance_overview_frame[substance_overview_frame$species == current_species, "ECOTOX_chemicals_filtered"] = current_ECOTOX_chemicals_filtered
  substance_overview_frame[substance_overview_frame$species == current_species, "overlap_unfiltered"] = current_overlap_unfiltered
  substance_overview_frame[substance_overview_frame$species == current_species, "overlap_filtered"] = current_overlap_filtered
  
}

# Finally add total row
substance_overview_frame = rbind(substance_overview_frame, c("species" = 'all',
                                                             "total_chemicals_unfiltered" = length(unique(c(EFSA$CASNO, ECOTOX$cas_number))),
                                                             "total_chemicals_filtered" = nrow(identifiers), 
                                                             "EFSA_chemicals_unfiltered" = length(unique(c(EFSA$CASNO))), 
                                                             "EFSA_chemicals_filtered" = nrow(EFSA_identifiers),
                                                             "ECOTOX_chemicals_unfiltered" = length(unique(c(ECOTOX$cas_number))), 
                                                             "ECOTOX_chemicals_filtered" = nrow(ECOTOX_identifiers), 
                                                             "overlap_unfiltered" = length(unique(EFSA[EFSA$CASNO %in% ECOTOX[, "cas_number"], 'CASNO'])), 
                                                             "overlap_filtered" = sum(EFSA_identifiers$original_CAS %in% ECOTOX_identifiers$original_CAS)))

# Save the overlap as a CSV-file
write_delim(substance_overview_frame, file = 'compounds_per_species_table.csv', delim = '\t', col_names = T, quote = 'none')




################################################################################
#             3.1. Calculate metadata on QSAR predictions                      #
################################################################################

vega_unique_substances_daphnia_acute = length(unique(QSAR_output[QSAR_output$QSAR_tool == 'VEGA' & grepl(QSAR_output$model_organism, pattern = 'daphni', ignore.case = T) & QSAR_output$duration == '48h', 'original_CAS']))
vega_unique_substances_daphnia_chronic = length(unique(QSAR_output[QSAR_output$QSAR_tool == 'VEGA' & grepl(QSAR_output$model_organism, pattern = 'daphni', ignore.case = T) & QSAR_output$duration == '21d', 'original_CAS']))
vega_unique_substances_fish_acute = length(unique(QSAR_output[QSAR_output$QSAR_tool == 'VEGA' & grepl(QSAR_output$model_organism, pattern = 'fish', ignore.case = T) & QSAR_output$duration == '96h', 'original_CAS']))
vega_unique_substances_fish_chronic = length(unique(QSAR_output[QSAR_output$QSAR_tool == 'VEGA' & grepl(QSAR_output$model_organism, pattern = 'fish', ignore.case = T) & QSAR_output$duration == 'ELS', 'original_CAS']))
vega_unique_substances_algae_acute = length(unique(QSAR_output[QSAR_output$QSAR_tool == 'VEGA' & grepl(QSAR_output$model_organism, pattern = 'algae', ignore.case = T) & QSAR_output$duration == '72h' & QSAR_output$model_endpoint == 'EC50', 'original_CAS']))
vega_unique_substances_algae_chronic = length(unique(QSAR_output[QSAR_output$QSAR_tool == 'VEGA' & grepl(QSAR_output$model_organism, pattern = 'algae', ignore.case = T) & QSAR_output$duration == '72h' & QSAR_output$model_endpoint == 'NOEC', 'original_CAS']))

ecosar_unique_substances_daphnia_acute = length(unique(QSAR_output[QSAR_output$QSAR_tool == 'ECOSAR' & grepl(QSAR_output$model_organism, pattern = 'daphni', ignore.case = T) & QSAR_output$duration == '48h', 'original_CAS']))
ecosar_unique_substances_daphnia_chronic = length(unique(QSAR_output[QSAR_output$QSAR_tool == 'ECOSAR' & grepl(QSAR_output$model_organism, pattern = 'daphni', ignore.case = T) & QSAR_output$model_endpoint == 'ChV', 'original_CAS']))
ecosar_unique_substances_fish_acute = length(unique(QSAR_output[QSAR_output$QSAR_tool == 'ECOSAR' & grepl(QSAR_output$model_organism, pattern = 'fish', ignore.case = T) & QSAR_output$duration == '96h', 'original_CAS']))
ecosar_unique_substances_fish_chronic = length(unique(QSAR_output[QSAR_output$QSAR_tool == 'ECOSAR' & grepl(QSAR_output$model_organism, pattern = 'fish', ignore.case = T) & QSAR_output$model_endpoint == 'ChV', 'original_CAS']))
ecosar_unique_substances_algae_acute = length(unique(QSAR_output[QSAR_output$QSAR_tool == 'ECOSAR' & grepl(QSAR_output$model_organism, pattern = 'algae', ignore.case = T) & QSAR_output$duration == '96h' & QSAR_output$model_endpoint == 'EC50', 'original_CAS']))
ecosar_unique_substances_algae_chronic = length(unique(QSAR_output[QSAR_output$QSAR_tool == 'ECOSAR' & grepl(QSAR_output$model_organism, pattern = 'algae', ignore.case = T) & QSAR_output$model_endpoint == 'ChV', 'original_CAS']))

test_unique_substances_daphnia_acute = length(unique(QSAR_output[QSAR_output$QSAR_tool == 'T.E.S.T.' & grepl(QSAR_output$model_organism, pattern = 'daphni', ignore.case = T) & QSAR_output$duration == '48h', 'original_CAS']))
test_unique_substances_fish_acute = length(unique(QSAR_output[QSAR_output$QSAR_tool == 'T.E.S.T.' & grepl(QSAR_output$model_organism, pattern = 'fish', ignore.case = T) & QSAR_output$duration == '96h', 'original_CAS']))

