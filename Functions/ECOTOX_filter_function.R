################################################################################
#    A function for filtering ECOTOX based on predefined filters                 #
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
# 2. Function body


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

function(ECOTOX_database,
         filters = NULL,
         settings = NULL){
  
  # Check that the input makes sense
  if(is.null(ECOTOX_database)|!is.data.frame(ECOTOX_database)){
    print('ECOTOX_database not provided or not a dataframe')
    print('Rerun function with correct arguments')
    return(FALSE)
  }
  
  # Check if additional settings provided
  if(is.null(settings)){
    print('No additional settings provided, only fixing CAS-numbers')
    
    # Reset settings to empty string (to allow tolower to work)
    settings = ''
    
  } 
  
  # Check if filters provided
  if(is.null(filters)){
    print('No filters set, returning frame')
    return(ECOTOX_database)
    
  } 
  
  #check if filters is a list, and if it is named
  if(!is.list(filters)){
    print('Filters not a list')
    return(FALSE)
  } else if(is.null(names(filters))){
    print('Filters not named list')
    return(FALSE)
  }
  
  ################################################################################
  # 2. function body
  ################################################################################
  
  # Operator filter (if used)
  if('operator' %in% names(filters)& 'conc1_mean_op' %in% colnames(ECOTOX_database)){
    # filter out limit values
    ECOTOX_database = ECOTOX_database[ECOTOX_database$conc1_mean_op %in% filters[['operator']],]
    
  }
  
  # formulation filter (if used)
  if('active substance' %in% names(filters) & 'test_formulation' %in% colnames(ECOTOX_database)){
    # filter out formulations
    
    ## filter out formulations (ECOTOX handles this differently than EFSA and so this part is a bit wonky)
    
    # Get terms for non-formulations (using regexp to match descriptions or definitions in ecotox appendix)
    ecotox_chemical_formulations = read.csv(file = 'C:/Users/admin/Downloads/ECOTOX-Term-Appendix-C.csv')
    ecotox_chemical_formulations_no_formulation = 
      ecotox_chemical_formulations[
        !grepl(ecotox_chemical_formulations$DEFINITION, 
               pattern = '([Ff]ormulation)|([Cc]ommercial)|([Ss]uspoemulsion)|([Tt]ablet)') &
          !grepl(ecotox_chemical_formulations$DESCRIPTION, 
                 pattern = '([Ff]ormulation)|([Cc]ommercial)|([Ss]uspoemulsion)|([Tt]ablet)'),]
    
    # Get all entries which match non-formulations
    ECOTOX_database = ECOTOX_database[ECOTOX_database$test_formulation %in% ecotox_chemical_formulations_no_formulation$TERM,]
    
  }
  
  # salts filter (if used)
  if('salts' %in% names(filters) & 'original_SMILES' %in% colnames(ECOTOX_database)){
    # filter out salts
    ECOTOX_database = ECOTOX_database[!grepl(ECOTOX_database$original_SMILES, pattern = filters[['salts']]),]
  }
  
  # medium filter (if used)
  if('medium' %in% names(filters)& 'media_type' %in% colnames(ECOTOX_database)){
    # Filter for freshwater
    ECOTOX_database = ECOTOX_database[grepl(ECOTOX_database$media_type, pattern = filters[['medium']]),]
  }
  
  # species filter (if used)
  if('species' %in% names(filters)& 'latin_name' %in% colnames(ECOTOX_database)){
    
    ECOTOX_database = ECOTOX_database[ECOTOX_database$latin_name %in% filters[['species']] ,]
  }
  
  # duration filter (if used)
  if('duration' %in% names(filters)& 'obs_duration_mean' %in% colnames(ECOTOX_database)){
    # Durations come in different forms, check if this is the case
    if(grepl(filters[['duration']], pattern = '>')){
      ECOTOX_database = ECOTOX_database[ECOTOX_database$obs_duration_mean > str_extract(filters[['duration']], pattern = '[0-9]') ,]
    } else if(grepl(filters[['duration']], pattern = '<')){
      ECOTOX_database = ECOTOX_database[ECOTOX_database$obs_duration_mean < str_extract(filters[['duration']], pattern = '[0-9]') ,]
    } else if(grepl(filters[['duration']], pattern = '-')){
      min_dur = str_extract(filters[['duration']], pattern = '^[0-9]')
      max_dur = str_extract(filters[['duration']], pattern = '[0-9]$')
      ECOTOX_database = ECOTOX_database[ECOTOX_database$obs_duration_mean >= min_dur  & ECOTOX_database$obs_duration_mean <= max_dur,]
    } else {
      ECOTOX_database = ECOTOX_database[ECOTOX_database$obs_duration_mean %in% filters[['duration']] ,]
    }
  }
  
  # endpoint filter (if used)
  if('endpoint' %in% names(filters) & 'endpoint' %in% colnames(ECOTOX_database)){
    
    ECOTOX_database = ECOTOX_database[ECOTOX_database$endpoint %in% filters[['endpoint']] ,]
  }
  
  # effect filter (if used)
  if('effect' %in% names(filters) & 'effect' %in% colnames(ECOTOX_database)){
    
    ECOTOX_database = ECOTOX_database[grepl(ECOTOX_database$effect, pattern = filters[['effect']]) ,]
  }
  
  return(ECOTOX_database) 
  
  
}