################################################################################
#    A function for importing and cleaning up EFSA for QSAR_main               #
################################################################################
#
# Original author: Patrik Svedberg
#
# Based on the script 
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

# Made minor changes to facilitate separation of data collection from analysis


################################################################################
# 1. Initiating function and parsing arguments
################################################################################

function(EFSA_filepath,
         version = NULL,
         rerun = FALSE,
         filters = list('operator' = '=',
                        'reliability' = '1|2',
                        'active substance' = 'active substance',
                        'salts' = '\\.',
                        'medium' = 'freshwater'),
         efsa_cir_dumpfile = paste0(intermediate_directory, '/EFSA_CIR_dump.Rda'),
         settings = NULL){
  
  # Get functions
  EFSA_cleanup_function = dget('Functions/EFSA_cleanup_function.R')
  EFSA_filter_function = dget('Functions/EFSA_filter_function.R')
  QSAR_add_smiles_function = dget('Functions/QSAR_add_smiles_function.R')
  
  # Check that the input makes sense
  if(is.null(EFSA_filepath)){
    print('EFSA_filepath not provided')
    print('Rerun function with correct arguments')
    return(FALSE)
  }
  
  print('Parsing settings')
  if(is.null(settings)){
    print('No additional settings provided')
  }

  # Handling of cached/dump-files
  load_cache_file = F
  if(rerun){
    
    print('Rerunning function, this might take some time due to CIR_query')
    
  } else if(!is.null(version)) {
    
    load_cache_file = T
    print('Loading cache_file')
    
  } else {
    
    print('No rerun set, but no cache file version provided. Exiting function.')
    return(FALSE)
    
  }
  
  
  ################################################################################
  # 2. function body
  ################################################################################
  
  # If we are not loading a cached version
  if(!load_cache_file){
    
    print('Loading EFSA file')
    
    # The actual reading of the file
    EFSA = as.data.frame(read_xlsx(EFSA_filepath))
    
    # Check colnames in EFSA
    
    EFSA_expected_colnames = c("documentReferencePK", "substanceKey", "substanceDesc", "CASNO", 
                               "ECNO", "studyName", "endpointKind", "qualifier", "name", "deviation", 
                               "type", "year", "title", "reportNo", "company", "companyStudyNo", 
                               "resutType", "reliability", "rationalReliability", "GLPCompliance", 
                               "testMaterial", "testMaterialIndicator", "applicationMethod", 
                               "analyticalMonitoring", "vehicle", "testType", "medium", "testDuration", 
                               "limitTest", "exposureDurationValue", "exposureDurationUnit", 
                               "organism", "expDurationValue", "expDurationUnit", "endpointType", 
                               "minQualifier", "minValue", "maxQualifier", "maxValue", "endpointUnit", 
                               "basisConc", "effectConcType", "basisEffect", "remarks")
    
    if(any(!colnames(EFSA) == EFSA_expected_colnames)){
      
      print('EFSA frame has other columns than expected, might cause errors')
      
      # Print missing colnames
      print('Missing colnames:')
      print(EFSA_expected_colnames[!EFSA_expected_colnames %in% colnames(EFSA)])
      
      # Print not expcted colnames
      print('Unrecognized colnames:')
      print(colnames(EFSA)[!colnames(EFSA) %in% EFSA_expected_colnames])
      
      
    }
    
    # reporting loss of data due to filtering etc, this is the first report
    print('Number of data points in EFSA raw:')
    print(nrow(EFSA))
    print('Number of unique CAS in EFSA raw:')
    print(length(unique(EFSA$CASNO)))
    
    # Run EFSA cleanup function (Replaces NA CAS with CAS manually collected from PubChem 2022)
    # With the setting "fix species" to translate all old species names to new and add species_group
    print('Running EFSA_cleanup_function')
    EFSA_clean = EFSA_cleanup_function(EFSA,
                                       settings = c('fix species'))
    
    ## We drop entries that has no CAS after cleanup, since it will introduce issues downstream
    EFSA_clean = EFSA_clean[!is.na(EFSA_clean$CASNO),]
    
    # reporting loss of data due to filtering etc, this is the second report
    print('Number of data points in EFSA after cleanup:')
    print(nrow(EFSA_clean))
    print('Number of unique CAS in EFSA after cleanup:')
    print(length(unique(EFSA_clean$CASNO)))
    
    
    # Get all unique CAS
    EFSA_unique_cas = data.frame('original_CAS' = unique(EFSA_clean$CASNO))
    
    # Set load of dumpfile or not
    if(rerun){
      
      smiles_function_settings = c('no load', 'query cir')
      
    } else {
      
      smiles_function_settings = c('query cir')
      
    }
    
    # Run function to get SMILES for EFSA
    # and make new dataframe with identifiers
    EFSA_identifiers = QSAR_add_smiles_function(identifiers_database = EFSA_unique_cas, 
                                                local_dumpfile = efsa_cir_dumpfile,
                                                settings = smiles_function_settings)
    
    # Add smiles from the identifiers frame to the original EFSA frame
    EFSA_clean = merge(EFSA_clean, EFSA_identifiers, by.x = 'CASNO', by.y = 'original_CAS', all.x = T)
    EFSA_clean = EFSA_clean %>% 
      rename(original_CAS = CASNO)
    
    EFSA_clean = EFSA_clean %>% 
      relocate(original_SMILES, .after = original_CAS)
    
    
    # Add database identifier
    EFSA_clean$Database = 'EFSA'
    
    
    # filter EFSA based on function, with an initial general filter
    print('Running EFSA_filter_function')
    EFSA_filtered = EFSA_filter_function(EFSA_clean, 
                                         filters = filters)
    
    # Remove columns with only NAs
    EFSA_filtered = Filter(function(x)!all(is.na(x)), EFSA_filtered)
    
    # reporting loss of data due to filtering etc, this is the third report
    print('Number of data points in EFSA after filter:')
    print(nrow(EFSA_filtered))
    print('Number of unique CAS in EFSA after filter:')
    print(length(unique(EFSA_filtered$CASNO)))
    
    
    # Save EFSA_filtered, for work on other systems
    print('Saving output')
    save(EFSA_filtered, file = paste0(intermediate_directory, '/EFSA_filtered_', version, '.Rda'))
    save(EFSA_identifiers, file = paste0(intermediate_directory, '/EFSA_identifiers_', version, '.Rda'))
    
  } else {
    
    # If we are loading cached version
    load(file = paste0(intermediate_directory, '/EFSA_filtered_', version, '.Rda'))
    load(file = paste0(intermediate_directory, '/EFSA_identifiers_', version, '.Rda'))
    
  }
  
  print('Finishing EFSA import')
  return(list('EFSA_data' = EFSA_filtered,
              'EFSA_identifiers' = EFSA_identifiers))
  
}