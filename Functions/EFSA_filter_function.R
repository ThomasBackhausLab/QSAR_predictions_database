################################################################################
#    A function for filtering EFSA based on predefined filters                 #
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

################################################################################
# 1. Initiating function and parsing arguments
################################################################################

function(EFSA_database,
         filters = NULL,
         settings = NULL){
  
  # Check that the input makes sense
  if(is.null(EFSA_database)|!is.data.frame(EFSA_database)){
    print('EFSA_database not provided or not a dataframe')
    print('Rerun function with correct arguments')
    return(FALSE)
  }
  
  # Check if additional settings provided
  if(is.null(settings)){
    print('No additional settings provided')
    
    # Reset settings to empty string (to allow tolower to work)
    settings = ''
    
  } 
  
  # Check if filters provided
  if(is.null(filters)){
    print('No filters set, returning frame')
    return(EFSA_database)
    
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
  if('operator' %in% names(filters) & 'expValueOp' %in% colnames(EFSA_database)){
    # filter out limit values
    EFSA_database = EFSA_database[EFSA_database$expValueOp %in% filters[['operator']],]
    
  }
  
  # Reliability filter (if used)
  if('reliability' %in% names(filters)& 'reliability' %in% colnames(EFSA_database)){
    # filter for reliability scores in filter_vals
    EFSA_database = EFSA_database[grepl(EFSA_database$reliability, pattern = filters[['reliability']]),]
  }
  
  # formulation filter (if used)
  if('active substance' %in% names(filters)& 'substanceDesc' %in% colnames(EFSA_database)){
    # filter out formulations
    EFSA_database = EFSA_database[grepl(EFSA_database$substanceDesc, pattern = filters[['active substance']], ignore.case = T),]
  }
  
  # salts filter (if used)
  if('salts' %in% names(filters)& 'original_SMILES' %in% colnames(EFSA_database)){
    # filter out salts
    EFSA_database = EFSA_database[!grepl(EFSA_database$original_SMILES, pattern = filters[['salts']]),]
  }
  
  # medium filter (if used)
  if('medium' %in% names(filters)& 'medium' %in% colnames(EFSA_database)){
    # Filter for freshwater
    EFSA_database = EFSA_database[grepl(EFSA_database$medium, pattern = filters[['medium']]),]
  }
  
  # species filter (if used)
  if('species' %in% names(filters)& 'organism' %in% colnames(EFSA_database)){
  
    EFSA_database = EFSA_database[EFSA_database$organism %in% filters[['species']] ,]
  }
  
  # duration filter (if used)
  if('duration' %in% names(filters)& 'Duration_hour' %in% colnames(EFSA_database)){
    # Durations come in different forms, check if this is the case
    if(grepl(filters[['duration']], pattern = '>')){
      EFSA_database = EFSA_database[EFSA_database$Duration_hour > str_extract(filters[['duration']], pattern = '[0-9]') ,]
    } else if(grepl(filters[['duration']], pattern = '<')){
      EFSA_database = EFSA_database[EFSA_database$Duration_hour < str_extract(filters[['duration']], pattern = '[0-9]') ,]
    } else if(grepl(filters[['duration']], pattern = '-')){
      min_dur = str_extract(filters[['duration']], pattern = '^[0-9]')
      max_dur = str_extract(filters[['duration']], pattern = '[0-9]$')
      EFSA_database = EFSA_database[EFSA_database$obs_duration_mean >= min_dur  & EFSA_database$obs_duration_mean <= max_dur,]
    } else {
      EFSA_database = EFSA_database[EFSA_database$Duration_hour %in% filters[['duration']] ,]
    }
  }
  
  # endpoint filter (if used)
  if('endpoint' %in% names(filters)& 'endpointType' %in% colnames(EFSA_database)){
    
    EFSA_database = EFSA_database[EFSA_database$endpointType %in% filters[['endpoint']] ,]
  }
  
  # effect filter (if used)
  if('effect' %in% names(filters)& 'basisEffect' %in% colnames(EFSA_database)){
    
    EFSA_database = EFSA_database[grepl(EFSA_database$basisEffect, pattern = filters[['effect']]) ,]
  }
  
  return(EFSA_database) 
  
  
}