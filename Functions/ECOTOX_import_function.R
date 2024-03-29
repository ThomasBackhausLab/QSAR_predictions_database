################################################################################
#    A function for importing and cleaning up ECOTOX for QSAR_main             #
################################################################################
#
# Original author: Patrik Svedberg
#
# Contact email: patrik.svedberg@bioenv.gu.se
# (if the above doesnt work (i.e. years down the line) try p.a.svedberg@gmail.com)
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
         ecotox_cir_lookupfile = paste0(intermediate_directory, '/ECOTOX_identifiers_cir_lookup.Rda'),
         molweight_lookup = NULL,
         ECx_to_NOEC = NULL,
         settings = NULL){
  
  ECOTOX_build_function = dget('Functions/ECOTOX_build_function.R')
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
  
  ECx_to_NOEC_bool = F
  if(!is.null(ECx_to_NOEC)){
    ECx_to_NOEC_bool = T
    print(paste0('ECx_to_NOEC set, translating ', ECx_to_NOEC, ' into NOECs'))
  }
  
  ################################################################################
  # 2. function body
  ################################################################################
  
  if(!load_cache_file){
    
    print('Building ECOTOX from ascii')
    
    ECOTOX = ECOTOX_build_function(ECOTOX_ascii_path = ECOTOX_filepath,
                                   settings = NULL)
    
    print('Finished building ECOTOX from ascii')
    
    # reporting loss of data due to filtering etc, this is the first report
    print('Number of data points in ECOTOX raw:')
    print(nrow(ECOTOX))
    print('Number of unique CAS in ECOTOX raw:')
    print(length(unique(ECOTOX$cas_number)))
    
    
    # If some ECx_to_NOEC is set, we run this here before initial filter
    if(ECx_to_NOEC_bool){
      
      # Check if single value, or vector and handle as such
      if(length(ECx_to_NOEC) == 1){
        
        # Single value
        current_ECx = ECx_to_NOEC
        
        # Make a regexp out of the ECx (for LCx, ICx etc)
        current_ECx = str_replace(current_ECx, pattern = '^[EIL]{1}[DC]{1}', replacement = '[EIL]{1}[DC]{1}')
        
        # Add a stop to the regexp, to avoid getting EC100s and the like
        current_ECx = paste0(current_ECx, '$')
        
        ECOTOX$endpoint = ifelse(grepl(ECOTOX$endpoint, pattern = current_ECx),
                                          'NOEC',
                                          ECOTOX$endpoint)
        
        
        head(sort(table(ECOTOX$endpoint), decreasing = T), n = 200)
        
        
      } else if(length(ECx_to_NOEC) > 1){
        
        # Multivalue, loop over values
        for(i in 1:length(ECx_to_NOEC)){
          
          current_ECx = ECx_to_NOEC[i]
          
          # Make a regexp out of the ECx (including LCx, ICx)
          current_ECx = str_replace(current_ECx, pattern = '^[EIL]{1}[DC]{1}', replacement = '[EIL]{1}[DC]{1}')
          
          # Add a stop to the regexp, to avoid getting EC100s and the like
          current_ECx = paste0(current_ECx, '$')
          
          ECOTOX$endpoint = ifelse(grepl(ECOTOX$endpoint, pattern = current_ECx),
                                   'NOEC',
                                   ECOTOX$endpoint)
          
          
        }
      }
    }
    
    
    
    
    
    
    
    
    # Initial filtering using filter function (to reduce computational load of further steps)
    print('Running ECOTOX_filter_function for the first time')
    ECOTOX_filtered_initial = ECOTOX_filter_function(ECOTOX_database = ECOTOX, 
                                                     filters = filters)
    
    # reporting loss of data due to filtering etc, this is the second report
    print('Number of data points in ECOTOX after initial filters:')
    print(nrow(ECOTOX_filtered_initial))
    print('Number of unique CAS in ECOTOX after initial filters:')
    print(length(unique(ECOTOX_filtered_initial$cas_number)))
    
    # Clean up ECOTOX with function
    print('Running ECOTOX_cleanup_function')
    ECOTOX_clean = ECOTOX_cleanup_function(ECOTOX_database = ECOTOX_filtered_initial, 
                                           molweights = molweight_lookup,
                                           ECx_to_NOEC = ECx_to_NOEC,
                                           settings = c('fix species', 'fix molweights'))
    
    
    # reporting loss of data due to filtering etc, this is the third report
    print('Number of data points in ECOTOX after cleanup:')
    print(nrow(ECOTOX_clean))
    print('Number of unique CAS in ECOTOX after cleanup:')
    print(length(unique(ECOTOX_clean$cas_number)))
    
    # Extract identifiers from ECOTOX
    ECOTOX_unique_cas = unique(ECOTOX_clean$cas_number)
    ECOTOX_unique_cas_fixed = as.cas(ECOTOX_unique_cas)
    ECOTOX_identifiers = data.frame(original_CAS = ECOTOX_unique_cas_fixed)
    
    
    # Set load of lookupfile or not
    if(rerun){
      
      smiles_function_settings = c('no load', 'query cir')
      
    } else {
      
      smiles_function_settings = c('query cir')
      
    }
    
    # Add SMILES from cir/pubchem using add_SMILES function
    ECOTOX_identifiers = QSAR_add_smiles_function(ECOTOX_identifiers, 
                                                  local_lookupfile = ecotox_cir_lookupfile,
                                                  settings = smiles_function_settings)
    
    
    # Add smiles from the identifiers frame to the original ECOTOX frame
    ECOTOX_clean = merge(ECOTOX_clean, ECOTOX_identifiers, by.x = 'cas_number', by.y = 'original_CAS', all.x = T)
    
    # Remove entries without SMILES
    ECOTOX_clean = ECOTOX_clean[!ECOTOX_clean$original_SMILES == '' & !is.na(ECOTOX_clean$original_SMILES),]
    
    # Filter ECOTOX with function again, this time removing some additional
    print('Running ECOTOX_filter_function for the second time')
    ECOTOX_filtered = ECOTOX_filter_function(ECOTOX_clean, filters = list(salts = '\\.'))
    
    # reporting loss of data due to filtering etc, this is the fourth report
    print('Number of data points in ECOTOX after second filter:')
    print(nrow(ECOTOX_filtered))
    print('Number of unique CAS in ECOTOX after second filter:')
    print(length(unique(ECOTOX_filtered$cas_number)))
    
    # Add pesticide class as NA to match with EFSA (Not needed anymore)
    # ECOTOX_filtered$pesticide_class = NA
    
    # Add Database = 'ECOTOX'
    ECOTOX_filtered$Database = 'ECOTOX'
    
    # Save ECOTOX_filteredfiltered, for work on other systems
    save(ECOTOX_filtered, file = paste0(intermediate_directory, '/ECOTOX_filtered_', version, '.Rda'))
    save(ECOTOX_identifiers, file = paste0(intermediate_directory, '/ECOTOX_identifiers_', version, '.Rda'))
    
  } else {
    
    load(file = paste0(intermediate_directory, '/ECOTOX_filtered_', version, '.Rda'))
    load(file = paste0(intermediate_directory, '/ECOTOX_identifiers_', version, '.Rda'))
    
  }
  
  return(list('ECOTOX_data' = ECOTOX_filtered,
              'ECOTOX_identifiers' = ECOTOX_identifiers))
  
}