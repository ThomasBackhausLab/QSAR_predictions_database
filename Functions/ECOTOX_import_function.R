################################################################################
#    A function for importing and cleaning up ECOTOX for QSAR_main             #
################################################################################
#
# Original author: Patrik Svedberg
#
# Based on the script QSARmerger_vX.R:
# 
#
################################################################################
#                  Table of contents                                           #
# 
# 0. Planned changes and WIP
#   0.1. Versions and changelog
# 1. Initiating function and parsing arguments
# 2. Packages and housekeeping


################################################################################
# 0. Planned changes and WIP
################################################################################
#

################################################################################
#         0.1. Versions and changelog                                          #
################################################################################

## Version 1.0

# Externalized the function

## Version 2.0

# Minor feedback changes

################################################################################
# 1. Initiating function and parsing arguments
################################################################################

function(ECOTOX_filepath,
         version = NULL,
         rerun = FALSE,
         filters = list('operator' = 'None',
                        'active substance' = 'active substance',
                        'medium' = '(FW)|(NONE)',
                        'species' = c(daphnia_species, oecd_fish_species, oecd_algae_species),
                        'endpoint' = c('NOEC', 'EC50', 'LC50', 'IC50', 'LD50')),
         molweight_dump = NULL,
         settings = NULL){
  
  ECOTOX_cleanup_function = dget('Functions/ECOTOX_cleanup_function.R')
  ECOTOX_filter_function = dget('Functions/ECOTOX_filter_function.R')
  QSAR_add_smiles_function = dget('Functions/QSAR_add_smiles_function.R')
  
  # Check that the input makes sense
  if(is.null(ECOTOX_filepath)){
    print('ECOTOX_filepath not provided')
    print('Rerun function with correct arguments')
    return(FALSE)
  }
  
  print('Parsing settings')
  if(is.null(settings)){
    print('No additional settings provided')
  }
  
  load_cache_file = F
  if(rerun){
    
    print('Rerunning function, this might take some time due to CIR_query')
    
  } else if(!is.null(version)) {
    
    load_cache_file = T
    print('Loading cache_file')
    
  } else {
    
    print('No rerun, but no cache file chosen. Exiting function.')
    return(FALSE)
    
  }
  
  
  ################################################################################
  # 2. function body
  ################################################################################
  
  if(!load_cache_file){
    
    print('Loading ECOTOX file')
    
    # Load usepa ecotox
    if(file.exists(ECOTOX_filepath)){
      load(ECOTOX_filepath)
    } else {
      print('Unknown ECOTOX file or filepath')
    }
    
    # Initial filtering using filter function (to reduce computational load of further steps)
    print('Running ECOTOX_filter_function for the first time')
    ECOTOX_filtered_initial = ECOTOX_filter_function(ECOTOX, filters = filters)
    
    
    # Clean up ECOTOX with function
    print('Running ECOTOX_cleanup_function')
    ECOTOX_clean = ECOTOX_cleanup_function(ECOTOX_filtered_initial, 
                                           molweights = molweight_dump,
                                           settings = c('fix species', 'fix molweights'))
    
    
    # Extract identifiers from ECOTOX
    
    ECOTOX_unique_cas = unique(ECOTOX_clean$cas_number)
    ECOTOX_unique_cas_fixed = as.cas(ECOTOX_unique_cas)
    ECOTOX_identifiers = data.frame(original_CAS = ECOTOX_unique_cas_fixed)
    
    # Add SMILES from cir/pubchem using add_SMILES function
    ECOTOX_identifiers = QSAR_add_smiles_function(ECOTOX_identifiers, local_dumpfile = 'ECOTOX_identifiers_cir_dump.Rda', )
    
    
    # Add smiles from the identifiers frame to the original ECOTOX frame
    ECOTOX_clean = merge(ECOTOX_clean, ECOTOX_identifiers, by.x = 'cas_number', by.y = 'original_CAS', all.x = T)
    
    # Remove entries without SMILES
    ECOTOX_clean = ECOTOX_clean[!ECOTOX_clean$original_SMILES == '' & !is.na(ECOTOX_clean$original_SMILES),]
    
    # Filter ECOTOX with function again, this time removing some additional
    print('Running ECOTOX_filter_function for the second time')
    ECOTOX_filtered = ECOTOX_filter_function(ECOTOX_clean, filters = list(salts = '\\.'))
    
    # Add pesticide class as NA to match with EFSA (Not needed anymore)
    # ECOTOX_filtered$pesticide_class = NA
    
    # Add Database = 'ECOTOX'
    ECOTOX_filtered$Database = 'ECOTOX'
    
    # Save ECOTOX_filteredfiltered, for work on other systems
    save(ECOTOX_filtered, file = paste0('ECOTOX_filtered_', version, '.Rda'))
    save(ECOTOX_identifiers, file = paste0('ECOTOX_identifiers_', version, '.Rda'))
    
  } else {
    
    load(file = paste0('ECOTOX_filtered_', version, '.Rda'))
    load(file = paste0('ECOTOX_identifiers_', version, '.Rda'))
    
  }
  
  return(list('ECOTOX_data' = ECOTOX_filtered,
              'ECOTOX_identifiers' = ECOTOX_identifiers))
  
}