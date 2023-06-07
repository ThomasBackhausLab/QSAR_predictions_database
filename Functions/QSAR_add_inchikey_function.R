################################################################################
#    A function for adding InChIKeys cir_query and                             #
#           a local dump file                                                  #
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

################################################################################
# 1. Initiating function and parsing arguments
################################################################################


function(identifiers_database,
         local_dumpfile = NULL,
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
  
  
  ### FIX DUMPFILE STUFFFFS
  
  
  # How to handle local dumpfile loading
  if('ignore dumpfile' %in% tolower(settings) | tolower(settings) == 'ignore dumpfile'){
    
    print('Ignoring local dumpfile, no load')
    load_dumpfile = F
    
  } else if('no load' %in% tolower(settings) | tolower(settings) == 'no load'){
    
    print('Not loading local dumpfile')
    load_dumpfile = F
    
  } else if(is.null(local_dumpfile)){
    
    print('No local dumpfile provided')
    load_dumpfile = F
    
  } else if(!is.null(local_dumpfile) & !file.exists(local_dumpfile)){
    
    print('Local dumpfile doesnt exist')
    load_dumpfile = F
    
  } else if(!is.null(local_dumpfile) & file.exists(local_dumpfile)){
    
    print('Local dumpfile exists, loading')
    load_dumpfile = T
    
  } else {
    
    print('Local dumpfile error, check code/coder and retry')
    return(FALSE)
  }
  
  # How to handle local dumpfile saving
  if('ignore dumpfile' %in% tolower(settings) | tolower(settings) == 'ignore dumpfile'){
    
    save_dumpfile = F
    
  } else if('no save' %in% tolower(settings) | tolower(settings) == 'no save'){
  
    print('Not saving local dumpfile')
    save_dumpfile = F
  
  } else if(is.null(local_dumpfile)){
    
    save_dumpfile = F
      
  } else if(!is.null(local_dumpfile) & !file.exists(local_dumpfile)){
    
    print(paste0('Local dumpfile will be saved as: ', local_dumpfile))
    
    save_dumpfile = T
    
  
    
  } else if(!is.null(local_dumpfile) & file.exists(local_dumpfile)){
    
    print(paste0('Local dumpfile will be saved as: ', local_dumpfile))
    save_dumpfile = T
    
  } else {
    
    print('Local dumpfile error, check code/coder and retry')
    return(FALSE)
    
  } 
  
  
  
  
  ################################################################################
  # 2. function body
  ################################################################################
  
  
  if(!load_dumpfile){
    
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
    
    dump_file = loadRData(local_dumpfile)
    
    identifiers_database = merge(identifiers_database, dump_file[,c('original_CAS', 'InChIKey')], by = 'original_CAS', all.x = T)
    
    identifiers_database = identifiers_database[rowSums(is.na(identifiers_database)) != ncol(identifiers_database), ]
    
    for(i in 1:nrow(identifiers_database)){
      
      if(!is.na(identifiers_database[i, "InChIKey"])){
        next
      }
      
      current_smiles = identifiers_database[i, "original_SMILES"]
      
      print('InChIKey not in dumpfile, querying CIR')
      
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
  
  if(save_dumpfile){
    
    save(identifiers_database, file = local_dumpfile)
    
  }
  
  return(identifiers_database)
  
}
  