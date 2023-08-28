################################################################################
#    A function for adding InChIKeys cir_query and                             #
#           a local lookup file                                                  #
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


function(identifiers_database,
         local_lookupfile = NULL,
         settings = NULL){
  
  
  # Check if "original_SMILES" is a column in input
  input_smiles = F
  if('original_SMILES' %in% colnames(identifiers_database)){
    print('original_SMILES detected in input, using as a base')
    input_smiles = T
  }
  
  # Check if additional settings provided
  if(is.null(settings)){
    print('No additional settings provided')
    
    # Reset settings to empty string (to allow tolower to work)
    settings = ''
    
  } 
  
  
  ### FIX lookupFILE STUFFFFS
  
  
  # How to handle local lookupfile loading
  if('ignore lookupfile' %in% tolower(settings) | tolower(settings) == 'ignore lookupfile'){
    
    print('Ignoring local lookupfile, no load')
    load_lookupfile = F
    
  } else if('no load' %in% tolower(settings) | tolower(settings) == 'no load'){
    
    print('Not loading local lookupfile')
    load_lookupfile = F
    
  } else if(is.null(local_lookupfile)){
    
    print('No local lookupfile provided')
    load_lookupfile = F
    
  } else if(!is.null(local_lookupfile) & !file.exists(local_lookupfile)){
    
    print('Local lookupfile doesnt exist')
    load_lookupfile = F
    
  } else if(!is.null(local_lookupfile) & file.exists(local_lookupfile)){
    
    print('Local lookupfile exists, loading')
    load_lookupfile = T
    
  } else {
    
    print('Local lookupfile error, check code/coder and retry')
    return(FALSE)
  }
  
  # How to handle local lookupfile saving
  if('ignore lookupfile' %in% tolower(settings) | tolower(settings) == 'ignore lookupfile'){
    
    save_lookupfile = F
    
  } else if('no save' %in% tolower(settings) | tolower(settings) == 'no save'){
  
    print('Not saving local lookupfile')
    save_lookupfile = F
  
  } else if(is.null(local_lookupfile)){
    
    save_lookupfile = F
      
  } else if(!is.null(local_lookupfile) & !file.exists(local_lookupfile)){
    
    print(paste0('Local lookupfile will be saved as: ', local_lookupfile))
    
    save_lookupfile = T
    
  
    
  } else if(!is.null(local_lookupfile) & file.exists(local_lookupfile)){
    
    print(paste0('Local lookupfile will be saved as: ', local_lookupfile))
    save_lookupfile = T
    
  } else {
    
    print('Local lookupfile error, check code/coder and retry')
    return(FALSE)
    
  } 
  
  
  
  
  ################################################################################
  # 2. function body
  ################################################################################
  
  
  if(!load_lookupfile){
    
    # Add InChIKey if there is none, using cir_query
    
    identifiers_database$InChIKey = NA
    
    for(i in 1:nrow(identifiers_database)){
      
      if(!is.na(identifiers_database[i, "InChIKey"])){
        next
      }
      
      current_smiles = identifiers_database[i, "original_SMILES"]
      
      temp_inchi = cir_query(current_smiles, representation = 'stdinchikey')
      
      # If we didnt get inchikeys from SMILES we do another run of cir_query with CAS
      if(is.na(temp_inchi$stdinchikey[1])){
        
        current_cas = identifiers_database[i, "original_CAS"]
        
        temp_inchi = cir_query(current_cas, representation = 'stdinchikey')
        
      }
      
      
      identifiers_database[i, 'InChIKey'] = temp_inchi$stdinchikey[1]
      
    }
  
    # # Fix some errors that disambigous smiles and/or cas may cause (getting Inchikeys from pubchem), not applicable with new version
    # 
    # identifiers_database[identifiers_database$original_CAS == '7440-48-4', "InChIKey"]  = "InChIKey=GUTLYIVDDKVIGB-UHFFFAOYSA-N"
    # identifiers_database[identifiers_database$original_CAS == '64-17-5', "InChIKey"]    = "InChIKey=LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
    # identifiers_database[identifiers_database$original_CAS == '74-89-5', "InChIKey"]    = "InChIKey=BAVYZALUXZFZLV-UHFFFAOYSA-N"
    # identifiers_database[identifiers_database$original_CAS == '302-01-2', "InChIKey"]   = "InChIKey=OAKJQQAXSVQMHS-UHFFFAOYSA-N"
    # identifiers_database[identifiers_database$original_CAS == '67-56-1', "InChIKey"]    =  "InChIKey=OKKJLVBELUTLKV-UHFFFAOYSA-N"
    # identifiers_database[identifiers_database$original_CAS == '75-08-1', "InChIKey"]    =  "InChIKey=DNJIEGIFACGWOD-UHFFFAOYSA-N"
    # identifiers_database[identifiers_database$original_CAS == '107-15-3', "InChIKey"]   =  "InChIKey=PIICEJLVQHRZGT-UHFFFAOYSA-N"
    # identifiers_database[identifiers_database$original_CAS == '60-34-4', "InChIKey"]    =  "InChIKey=HDZGCSFEDULWCS-UHFFFAOYSA-N"
    # identifiers_database[identifiers_database$original_CAS == '74-87-3', "InChIKey"]    =  "InChIKey=NEHMKBQYUWJMIP-UHFFFAOYSA-N"
    # 
    
    # Add molecule number (for re-matching)
    identifiers_database$molecule_number <- seq(1,nrow(identifiers_database), 1)
    
    
  } else {
    
    # Apparently my workaround didnt work, so here is one by ricardo at: https://stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file
    loadRData <- function(fileName){
      #loads an RData file, and returns it
      load(fileName)
      get(ls()[ls() != "fileName"])
    }
    
    lookup_file = loadRData(local_lookupfile)
    
    identifiers_database = merge(identifiers_database, lookup_file[,c('original_CAS', 'InChIKey')], by = 'original_CAS', all.x = T)
    
    identifiers_database = identifiers_database[rowSums(is.na(identifiers_database)) != ncol(identifiers_database), ]
    
    for(i in 1:nrow(identifiers_database)){
      
      if(!is.na(identifiers_database[i, "InChIKey"])){
        next
      }
      
      current_smiles = identifiers_database[i, "original_SMILES"]
      
      print('InChIKey not in lookupfile, querying CIR')
      
      # Clear tempicnhi to avoid duplicates being saved
      temp_inchi = NA
      cir_error = F
      temp_inchi <- tryCatch(
        {
          withTimeout({cir_query(current_smiles, representation = 'stdinchikey')}, 
                      timeout = 10)
        },
        TimeoutException = function(ex) cat("Timeout. Skipping.\n"),
        error=function(cond) {
          message(cond)
          # Choose a return value in case of error
          cir_error = T
          return(NA)
        },
        warning=function(cond) {
          message(cond)
          # Choose a return value in case of warning
          cir_error = T
          return(NA)
        }      
        )   
      
      
      # If we didnt get inchikeys from SMILES we do another run of cir_query with CAS
      if(is.na(temp_inchi$stdinchikey[1])){
        
        current_cas = identifiers_database[i, "original_CAS"]
        
        temp_inchi = cir_query(current_cas, representation = 'stdinchikey')
        
      }
      
      if(!cir_error){
        identifiers_database[i, 'InChIKey'] = temp_inchi$stdinchikey[1]
      } else {
        identifiers_database[i, 'InChIKey'] = NA
      }
      
      
    }
    
  }
  
  if(save_lookupfile){
    
    save(identifiers_database, file = local_lookupfile)
    
  }
  
  return(identifiers_database)
  
}
  