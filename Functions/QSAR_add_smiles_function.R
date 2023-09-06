################################################################################
#    A function for adding SMILES from CAS using cir_query and                 #
#           a local lookup file                                                  #
################################################################################
#
# Original author: Patrik Svedberg
#
# Contact email: patrik.svedberg@bioenv.gu.se
# (if the above doesnt work (i.e. years down the line) try p.a.svedberg@gmail.com)
#
# Based on the script :
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

# Minor changes to feedback and documentation


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
  } else {
    
    print('original_SMILES not detected in input, using CAS as a base')
    
  }
  
  # Check if additional settings provided
  if(is.null(settings)){
    print('No additional settings provided')
    
    # Reset settings to empty string (to allow tolower to work)
    settings = ''
    
  } 
  
  # Check if no manual fixes in settings
  manual_fix = T
  if('no manual fix' %in% tolower(settings) | tolower(settings) == 'no manual fix' | 'no manual fixes' %in% tolower(settings) | tolower(settings) == 'no manual fixes'){
    print('Settings: No manual fix')
    
    manual_fix = F
  }
  
  query_cir = F
  if('query cir' %in% tolower(settings) | tolower(settings) == 'query cir'){
    print('Settings: Query CIR')
    
    query_cir = T
    
    library(webchem)
  }
  
  # How to handle local lookupfile loading
  if('ignore lookupfile' %in% tolower(settings) | tolower(settings) == 'ignore lookupfile'){
    
    print('Settings: Ignoring local lookupfile, no load')
    load_lookupfile = F
    
  } else if('no load' %in% tolower(settings) | tolower(settings) == 'no load'){
    
    print('Settings: Not loading local lookupfile')
    load_lookupfile = F
    
  } else if(!is.null(local_lookupfile) & !file.exists(local_lookupfile)){
    
    print('Local lookupfile doesnt exist')
    load_lookupfile = F
    
  } else if(is.null(local_lookupfile)){
    
    print('Settings: No local lookupfile provided')
    load_lookupfile = F
    
  } else if(!is.null(local_lookupfile) & file.exists(local_lookupfile)){
    
    print('Settings: Local lookupfile exists, loading')
    load_lookupfile = T
    
  } else {
    
    print('Local lookupfile error, check code/coder and retry')
    return(FALSE)
  }
  
  # How to handle local lookupfile saving
  if('ignore lookupfile' %in% tolower(settings) | tolower(settings) == 'ignore lookupfile'){
    
    save_lookupfile = F
    
  } else if('no save' %in% tolower(settings) | tolower(settings) == 'no save'){
    
    print('Settings: Not saving local lookupfile')
    save_lookupfile = F
    
  } else if(!is.null(local_lookupfile) & !file.exists(local_lookupfile)){
    
    print(paste0('Settings: Local lookupfile will be saved as: ', local_lookupfile))
    
    save_lookupfile = T
    
  } else if(is.null(local_lookupfile)){
    
    save_lookupfile = F
    
  } else if(!is.null(local_lookupfile) & file.exists(local_lookupfile)){
    
    print(paste0('Settings: Local lookupfile will be saved as: ', local_lookupfile))
    save_lookupfile = T
    
  } else {
    
    print('Local lookupfile error, check code/coder and retry')
    return(FALSE)
    
  } 

  
  
  
  ################################################################################
  # 2. function body
  ################################################################################

  ################################################################################
  # 2.1 processing lookupfile, input smiles etc
  ################################################################################
  
  
  # If we load a lookup_file, do that first
  if(load_lookupfile){
    
    print('Loading local lookupfile')
    lookup_file = get(load(local_lookupfile))
  
  } else {
    
    lookup_file = NULL
    
  }
  
  # Make an outfile frame
  out_file = data.frame('original_CAS' = unique(identifiers_database$original_CAS))
  
  # if SMILES provided in input frame
  if(input_smiles){
    
    out_file = merge(out_file, identifiers_database[,c('original_CAS', 'original_SMILES')], by = 'original_CAS', all.x = T)  
    
  }
  
  # Merge if there is a lookupfile, else add NAs
  if(!is.null(lookup_file)){
  
    out_file = merge(out_file, lookup_file[,c('original_CAS', 'original_SMILES')], by = 'original_CAS', all.x = T)  
    
  } 
  
  # If both lookupfile and input smiles were used, fix the merger issues
  if("original_SMILES.x" %in% colnames(out_file) & "original_SMILES.y" %in% colnames(out_file)){
    
    out_file$original_SMILES = out_file$original_SMILES.y
    
    # Use x (input SMILES) if there are both input and local_lookupfile SMILES
    out_file$original_SMILES[!is.na(out_file$original_SMILES.x)] = out_file$original_SMILES.x[!is.na(out_file$original_SMILES.x)]
    
    # Remove the .x and .y columns
    out_file$original_SMILES.x = NULL
    out_file$original_SMILES.y = NULL
    
  }
  
  # If there still is no original_SMILES in out_file, add NAs
  if(!"original_SMILES" %in% colnames(out_file)){
    
    out_file$original_SMILES = NA
    
  }
  
  
  ################################################################################
  # 2.2 Query cir
  ################################################################################
  
  # Get all CAS missing SMILES
  temp_missing_SMILES = data.frame('original_CAS' = out_file[is.na(out_file$original_SMILES), 'original_CAS'])
  
  # query CIR if some are missing
  if(nrow(temp_missing_SMILES)>0 & query_cir){
    
    print('SMILES missing, querying CIR')
    
    # Message number of CAS missing SMILES
    message(paste0('Number of CAS missing SMILES: ', nrow(temp_missing_SMILES), ', starting query...'))
    
    # Get smiles for CAS missing SMILES
    temp_missing_SMILES$original_SMILES = as.data.frame(cir_query(temp_missing_SMILES$original_CAS, representation = 'smiles', match = 'first'))[,2]
    
    # Get all still missing
    temp_still_missing_SMILES = temp_missing_SMILES[is.na(temp_missing_SMILES$original_SMILES), 'original_CAS']
    
    # Message number of still missing
    message(paste0('Number of CAS still missing SMILES after cir_query: ', length(temp_still_missing_SMILES)))
    
    # Rebind the query results with original frame
    out_file = merge(out_file[,c('original_CAS', 'original_SMILES')], temp_missing_SMILES[,c('original_CAS', 'original_SMILES')], by = 'original_CAS', all = T)
    
    
    # coalesce columns
    if("original_SMILES.x" %in% colnames(out_file)){
      out_file = out_file %>% mutate(original_SMILES = coalesce(original_SMILES.x, original_SMILES.y)) %>% dplyr::select(original_CAS, original_SMILES)
    }
    
  } else {
    
    print('No SMILES missing, skipping CIR query')
    
  }
  
  # Remove any NA rows that may have been created
  out_file = out_file[rowSums(is.na(out_file)) != ncol(out_file), ]
  
  # And remove dupicate rows
  out_file = distinct(out_file)
  
  
  ################################################################################
  # 2.3 Manually add some SMILES from PubChem
  ################################################################################
  
  if(manual_fix){
    
    # Add SMILES from Pubchem
    # collected in 2022
    out_file[out_file$original_CAS == '203313-25-1', "original_SMILES"] = 'CCOC(=O)OC1=C(C(=O)NC12CCC(CC2)OC)C3=C(C=CC(=C3)C)C'
    out_file[out_file$original_CAS == '130561-48-7', "original_SMILES"] = 'COCCOC1=CC=CC2=C1C(=O)C(=NN2C3=CC=C(C=C3)Cl)C(=O)O'
    out_file[out_file$original_CAS == '106325-08-0', "original_SMILES"] = 'C1=CC=C(C(=C1)C2C(O2)(CN3C=NC=N3)C4=CC=C(C=C4)F)Cl'
    out_file[out_file$original_CAS == '500008-45-7', "original_SMILES"] = 'CC1=CC(=CC(=C1NC(=O)C2=CC(=NN2C3=C(C=CC=N3)Cl)Br)C(=O)NC)Cl'
    out_file[out_file$original_CAS == '145701-23-1', "original_SMILES"] = 'COC1=NC=C(C2=NC(=NN21)S(=O)(=O)NC3=C(C=CC=C3F)F)F'
    out_file[out_file$original_CAS == '149961-52-4', "original_SMILES"] = 'CC1=CC(=C(C=C1)C)OCC2=CC=CC=C2C(=NOC)C(=O)NC'
    out_file[out_file$original_CAS == '98243-83-5', "original_SMILES"] = 'CC1=C(C(=CC=C1)C)N(C(C)C(=O)OC)C(=O)CC2=CC=CC=C2'
    out_file[out_file$original_CAS == '865318-97-4', "original_SMILES"] = 'CCCCCCCCC1=C(N2C(=NC=N2)N=C1CC)N'
    out_file[out_file$original_CAS == '317815-83-1', "original_SMILES"] = 'CC1=C(C(=CS1)C(=O)OC)S(=O)(=O)NC(=O)N2C(=O)N(C(=N2)OC)C'
    out_file[out_file$original_CAS == '881685-58-1', "original_SMILES"] = 'CC(C)C1C2CCC1C3=C2C=CC=C3NC(=O)C4=CN(N=C4C(F)F)C'
    out_file[out_file$original_CAS == '374726-62-2', "original_SMILES"] = 'COC1=C(C=CC(=C1)CCNC(=O)C(C2=CC=C(C=C2)Cl)OCC#C)OCC#C'
    out_file[out_file$original_CAS == '658066-35-4', "original_SMILES"] = 'C1=CC=C(C(=C1)C(=O)NCCC2=C(C=C(C=N2)C(F)(F)F)Cl)C(F)(F)F'
    out_file[out_file$original_CAS == '210631-68-8', "original_SMILES"] = 'CC1=C(C=CC(=C1C2=NOCC2)S(=O)(=O)C)C(=O)C3=CNN(C3=O)C'
    out_file[out_file$original_CAS == '39300-45-3', "original_SMILES"] = 'CCCCCCC(C)C1=C(C(=CC(=C1)[N+](=O)[O-])[N+](=O)[O-])OC(=O)C=CC'
    out_file[out_file$original_CAS == '139968-49-3', "original_SMILES"] = 'C1=CC(=CC(=C1)C(F)(F)F)C(=NNC(=O)NC2=CC=C(C=C2)OC(F)(F)F)CC3=CC=C(C=C3)C#N'
    out_file[out_file$original_CAS == '95977-29-0', "original_SMILES"] = 'CC(C(=O)O)OC1=CC=C(C=C1)OC2=C(C=C(C=N2)C(F)(F)F)Cl'
    out_file[out_file$original_CAS == '189278-12-4', "original_SMILES"] = 'CCCN1C(=O)C2=C(C=CC(=C2)I)N=C1OCCC'
    
    out_file[out_file$original_CAS == '935545-74-7', "original_SMILES"] = 'CCC1CCCC(C(C(=O)C2=CC3C(C2CC(=O)O1)CCC4C3CC(C4)OC5C(C(C(C(O5)C)OC)OCC)OC)C)OC6CCC(C(O6)C)N(C)C'
    out_file[out_file$original_CAS == '400882-07-7', "original_SMILES"] = 'CC(C)(C)C1=CC=C(C=C1)C(C#N)(C(=O)C2=CC=CC=C2C(F)(F)F)C(=O)OCCOC'
    out_file[out_file$original_CAS == '494793-67-8', "original_SMILES"] = 'CC1=NN(C(=C1C(=O)NC2=CC=CC=C2C(C)CC(C)C)F)C'
    out_file[out_file$original_CAS == '335104-84-2', "original_SMILES"] = 'CS(=O)(=O)C1=C(C(=C(C=C1)C(=O)C2C(=O)CCCC2=O)Cl)COCC(F)(F)F'
    out_file[out_file$original_CAS == '874967-67-6', "original_SMILES"] = 'CN1C=C(C(=N1)C(F)F)C(=O)NC2=CC=CC=C2C3CC3C4CC4'
    out_file[out_file$original_CAS == '581809-46-3', "original_SMILES"] = 'CN1C=C(C(=N1)C(F)F)C(=O)NC2=C(C=C(C=C2)F)C3=CC(=C(C=C3)Cl)Cl'
    out_file[out_file$original_CAS == '348635-87-0', "original_SMILES"] = 'CC1=C(C2=C(N1S(=O)(=O)C3=NN(C=N3)S(=O)(=O)N(C)C)C=C(C=C2)F)Br'
    out_file[out_file$original_CAS == '183675-82-3', "original_SMILES"] = 'CC(C)CC(C)C1=C(C=CS1)NC(=O)C2=CN(N=C2C(F)(F)F)C'
    out_file[out_file$original_CAS == '8007-46-3', "original_SMILES"] = 'CC1=CCC(CC1)C(C)(C)O.CC1=CC=C(C=C1)C(C)C.CC1=CC(=C(C=C1)C(C)C)O.CC(C)C12CCC(C1C2)(C)O.CC(=CCCC(C)(C=C)O)C'
    out_file[out_file$original_CAS == '283159-90-0', "original_SMILES"] = 'CC(C)C(C(=O)NC(CC(=O)OC)C1=CC=C(C=C1)Cl)NC(=O)OC(C)C'
    out_file[out_file$original_CAS == '422556-08-9', "original_SMILES"] = 'COC1=CC(=NC2=NC(=NN12)NS(=O)(=O)C3=C(C=CN=C3OC)C(F)(F)F)OC'
    out_file[out_file$original_CAS == '213464-77-8', "original_SMILES"] = 'CN(C)C(=O)C1=CC=CC=C1NS(=O)(=O)NC(=O)NC2=NC(=CC(=N2)OC)OC'
    out_file[out_file$original_CAS == '473798-59-3', "original_SMILES"] = 'CC1=CC=CC=C1C2=C(N(N(C2=O)C(C)C)C(=O)SCC=C)N'
    out_file[out_file$original_CAS == '32588-36-6', "original_SMILES"] = 'C1=CC=C2C(=C1)C=C(N2)CC(=O)O'
    
    out_file[out_file$original_CAS == '153233-91-1', "original_SMILES"] = 'CCOC1=C(C=CC(=C1)C(C)(C)C)C2COC(=N2)C3=C(C=CC=C3F)F'
  
    # collected in
    out_file[out_file$original_CAS == '898-84-0' ,'original_SMILES'] = 'CC12CC(C3C(C1CCC2=O)CCC4=CC(=O)C=CC34C)O'
    out_file[out_file$original_CAS == '74610-55-2' ,'original_SMILES'] = 'CCC1C(C=C(C=CC(=O)C(CC(C(C(C(CC(=O)O1)O)C)OC2C(C(C(C(O2)C)OC3CC(C(C(O3)C)O)(C)O)N(C)C)O)CC=O)C)C)COC4C(C(C(C(O4)C)O)OC)OC.C(C(C(=O)O)O)(C(=O)O)O'
    out_file[out_file$original_CAS == '30125-65-6' ,'original_SMILES'] = 'CC(C)(C)NC1=NC(=NC(=N1)N)SC'
    out_file[out_file$original_CAS == '2610-61-9' ,'original_SMILES'] = 'CC(C)OC(=O)NC1=CC(=CC=C1)O'
    out_file[out_file$original_CAS == '812-70-4' ,'original_SMILES'] = 'C(CC(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)C(=O)O'
    out_file[out_file$original_CAS == '53826-12-3' ,'original_SMILES'] = 'C(C(=O)O)C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F'
    out_file[out_file$original_CAS == '27854-31-5' ,'original_SMILES'] = 'C(C(=O)O)C(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F'
    out_file[out_file$original_CAS == '914637-49-3' ,'original_SMILES'] = 'C(CC(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)C(=O)O'
    out_file[out_file$original_CAS == '26952-21-6' ,'original_SMILES'] = 'CC(C)CCCCCO'
    out_file[out_file$original_CAS == '114247-09-5' ,'original_SMILES'] = 'CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F.Cl'
    out_file[out_file$original_CAS == '25155-30-0' ,'original_SMILES'] = 'CCCCCCCCCC(CC)C1=CC=C(C=C1)S(=O)(=O)[O-].[Na+]'
    out_file[out_file$original_CAS == '70887-84-2' ,'original_SMILES'] = 'C(=C(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)F)C(=O)O'
    out_file[out_file$original_CAS == '70887-88-6' ,'original_SMILES'] = 'C(=C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)F)C(=O)O'
    out_file[out_file$original_CAS == '166812-75-5' ,'original_SMILES'] = 'CCCCC=CCCCCCCCCCCC(=O)C(F)(F)F'
    out_file[out_file$original_CAS == '7722-64-7' ,'original_SMILES'] = '[O-][Mn](=O)(=O)=O.[K+]'
    out_file[out_file$original_CAS == '53826-13-4' ,'original_SMILES'] = 'C(C(=O)O)C(C(C(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F'
    out_file[out_file$original_CAS == '70887-94-4' ,'original_SMILES'] = 'C(=C(C(C(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)F)C(=O)O'
    out_file[out_file$original_CAS == '100473-08-3' ,'original_SMILES'] = 'CC(C)(C)C1=CC(=C(C=C1)OP(=O)(OC2=CC=CC=C2)OC3=CC=CC=C3)C(C)(C)C'
    out_file[out_file$original_CAS == '225789-38-8' ,'original_SMILES'] = 'CCP(=O)(CC)[O-].CCP(=O)(CC)[O-].CCP(=O)(CC)[O-].[Al+3]'
    out_file[out_file$original_CAS == '1203791-41-6' ,'original_SMILES'] = 'C1CC2N3CCN(C3=C(C1O2)[N+](=O)[O-])CC4=CN=C(C=C4)Cl'
    out_file[out_file$original_CAS == '145307-19-3' ,'original_SMILES'] = 'CCOC(=O)CC(C(=O)OCC)SP(=O)(OC)OC'
    out_file[out_file$original_CAS == '308068-56-6' ,'original_SMILES'] = '[C]'
    out_file[out_file$original_CAS == '3663-24-9' ,'original_SMILES'] = 'CCCCC1=CC2=NNN=C2C=C1'
    out_file[out_file$original_CAS == '8001-54-5' ,'original_SMILES'] = 'CCCCCCCCCCCC[N+](C)(C)CC1=CC=CC=C1.[Cl-]'
    out_file[out_file$original_CAS == '45298-90-6' ,'original_SMILES'] = 'C(C(C(C(C(F)(F)S(=O)(=O)[O-])(F)F)(F)F)(F)F)(C(C(C(F)(F)F)(F)F)(F)F)(F)F'
    out_file[out_file$original_CAS == '116498-46-5' ,'original_SMILES'] = 'CCCC1COC(O1)(CN2C=NC=N2)C3=C(C=C(C=C3)Cl)Cl'
    out_file[out_file$original_CAS == '25085-02-3' ,'original_SMILES'] = 'C=CC(=O)N.C=CC(=O)O.[Na]'
    out_file[out_file$original_CAS == '37337-13-6' ,'original_SMILES'] = 'O=[Cr](=O)=O.O=[Cu].O=[As](=O)O[As](=O)=O'
    out_file[out_file$original_CAS == '141318-03-8' ,'original_SMILES'] = 'CCOC(=O)CC(C(=O)OCC)SP(=S)(OC)OC'
    out_file[out_file$original_CAS == '8003-34-7' ,'original_SMILES'] = 'CC1=C(C(=O)CC1OC(=O)C2C(C2(C)C)C=C(C)C)CC=CC=C.CC1=C(C(=O)CC1OC(=O)C2C(C2(C)C)C=C(C)C(=O)OC)CC=CC=C'
    out_file[out_file$original_CAS == '1034343-98-0' ,'original_SMILES'] = '[C]'
    out_file[out_file$original_CAS == '1161016-80-3' ,'original_SMILES'] = 'CC1COC(O1)(CN2C=NC=N2)C3=C(C=C(C=C3)OC4=CC=C(C=C4)Cl)Cl'
    out_file[out_file$original_CAS == '70630-17-0' ,'original_SMILES'] = 'CC1=C(C(=CC=C1)C)N(C(C)C(=O)OC)C(=O)COC'
    out_file[out_file$original_CAS == '15067-52-4' ,'original_SMILES'] = 'CCCCOC(C)COC(=O)C(C)OC1=CC(=C(C=C1Cl)Cl)Cl'
    out_file[out_file$original_CAS == '145307-18-2' ,'original_SMILES'] = 'CCOC(=O)CC(C(=O)OCC)SP(=O)(OC)OC'
    out_file[out_file$original_CAS == '1344-28-1' ,'original_SMILES'] = '[O-2].[O-2].[O-2].[Al+3].[Al+3]'
    out_file[out_file$original_CAS == '736994-63-1' ,'original_SMILES'] = 'CC1=CC(=CC(=C1NC(=O)C2=CC(=NN2C3=C(C=CC=N3)Cl)Br)C(=O)NC)C#N'
    out_file[out_file$original_CAS == '41318-75-6' ,'original_SMILES'] = 'C1=CC(=CC=C1OC2=C(C=C(C=C2)Br)Br)Br'
    out_file[out_file$original_CAS == '298-07-7' ,'original_SMILES'] = 'CCCCC(CC)COP(=O)(O)OCC(CC)CCCC'
    out_file[out_file$original_CAS == '9012-54-8' ,'original_SMILES'] = 'C(C1C(C(C(C(O1)OC2C(OC(C(C2O)O)OC3C(OC(C(C3O)O)O)CO)CO)O)O)O)O'
    out_file[out_file$original_CAS == '1333-86-4' ,'original_SMILES'] = '[C]'
    out_file[out_file$original_CAS == '1161016-82-5' ,'original_SMILES'] = 'CC1COC(O1)(CN2C=NC=N2)C3=C(C=C(C=C3)OC4=CC=C(C=C4)Cl)Cl'
    out_file[out_file$original_CAS == '116498-43-2' ,'original_SMILES'] = 'CCCC1COC(O1)(CN2C=NC=N2)C3=C(C=C(C=C3)Cl)Cl'
    out_file[out_file$original_CAS == '114247-06-2' ,'original_SMILES'] = 'CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F.Cl'
    out_file[out_file$original_CAS == '15096-52-3' ,'original_SMILES'] = 'F[Al-3](F)(F)(F)(F)F.[Na+].[Na+].[Na+]'
    out_file[out_file$original_CAS == '65589-70-0' ,'original_SMILES'] = 'C[N+]1=C2C=C(C=CC2=CC3=C1C=C(C=C3)N)N.C1=CC(=CC2=NC3=C(C=CC(=C3)N)C=C21)N.[Cl-]'
    out_file[out_file$original_CAS == '15773-35-0' ,'original_SMILES'] = 'C1(=C(C(=C(C(=C1Cl)Cl)Cl)Cl)Cl)[O-].C1(=C(C(=C(C(=C1Cl)Cl)Cl)Cl)Cl)[O-].[Cu+2]'
    out_file[out_file$original_CAS == '174501-64-5' ,'original_SMILES'] = 'CCCCN1C=C[N+](=C1)C.F[P-](F)(F)(F)(F)F'
    out_file[out_file$original_CAS == '1161016-84-7' ,'original_SMILES'] = 'CC1COC(O1)(CN2C=NC=N2)C3=C(C=C(C=C3)OC4=CC=C(C=C4)Cl)Cl'
    out_file[out_file$original_CAS == '500008-45-7' ,'original_SMILES'] = 'CC1=CC(=CC(=C1NC(=O)C2=CC(=NN2C3=C(C=CC=N3)Cl)Br)C(=O)NC)Cl'
    out_file[out_file$original_CAS == '1330-69-4' ,'original_SMILES'] = 'CCCCCCCCCCCC.C1=CC=CC=C1.[O-]S(=O)[O-]'
    out_file[out_file$original_CAS == '174501-65-6' ,'original_SMILES'] = '[B-](F)(F)(F)F.CCCCN1C=C[N+](=C1)C'
    out_file[out_file$original_CAS == '71662-11-8' ,'original_SMILES'] = 'C1COCCN1CC2CN(C(=O)O2)N=CC3=CC=C(O3)[N+](=O)[O-].C(C(C(=O)O)O)(C(=O)O)O'
    out_file[out_file$original_CAS == '85100-77-2' ,'original_SMILES'] = 'CCCCN1C=C[N+](=C1)C.[Br-]'
    out_file[out_file$original_CAS == '11132-78-8' ,'original_SMILES'] = '[Cl-].[Cl-].[Mn+2]'
    out_file[out_file$original_CAS == '8000-48-4' ,'original_SMILES'] = 'CC1(C2CCC(O1)(CC2)C)C'
    out_file[out_file$original_CAS == '1161016-86-9' ,'original_SMILES'] = 'CC1COC(O1)(CN2C=NC=N2)C3=C(C=C(C=C3)OC4=CC=C(C=C4)Cl)Cl'
    out_file[out_file$original_CAS == '960003-92-3' ,'original_SMILES'] = 'CCC(C)SP(=O)(N1CCSC1=O)OCC'
    out_file[out_file$original_CAS == '20762-60-1' ,'original_SMILES'] = '[N-]=[N+]=[N-].[K+]'
    out_file[out_file$original_CAS == '68333-79-9' ,'original_SMILES'] = '[NH4+].[NH4+].[NH4+].[O-]P(=O)([O-])[O-]'
    out_file[out_file$original_CAS == '37841-25-1' ,'original_SMILES'] = 'C1=C(C=C(C(=C1[N+](=O)[O-])C#N)[N+](=O)[O-])[N+](=O)[O-]'
    out_file[out_file$original_CAS == '31512-74-0' ,'original_SMILES'] = 'CCC[N+](C)(C)CC[N+](C)(C)CCOC.[Cl-].[Cl-]'
    out_file[out_file$original_CAS == '960003-90-1' ,'original_SMILES'] = 'CCC(C)SP(=O)(N1CCSC1=O)OCC'
    out_file[out_file$original_CAS == '960003-91-2' ,'original_SMILES'] = 'CCC(C)SP(=O)(N1CCSC1=O)OCC'
    out_file[out_file$original_CAS == '960003-93-4' ,'original_SMILES'] = 'CCC(C)SP(=O)(N1CCSC1=O)OCC' # Stereo isomers
    out_file[out_file$original_CAS == '5138-93-2' ,'original_SMILES'] = 'C1=CC(=C(C=C1Cl)S(=O)(=O)[O-])Cl.[Na+]'
    out_file[out_file$original_CAS == '18402-10-3' ,'original_SMILES'] = 'CCCCCCCCCCO[Si](C)(C)C'
    out_file[out_file$original_CAS == '439680-76-9' ,'original_SMILES'] = 'CC1=C(C=CC=C1C2=CC=CC=C2)COC(=O)C3C(C3(C)C)C=C(C(F)(F)F)Cl'
    out_file[out_file$original_CAS == '87648-90-6' ,'original_SMILES'] = 'CC1=C(C=CC=C1C2=CC=CC=C2)COC(=O)C3C(C3(C)C)C=C(C(F)(F)F)Cl'
    out_file[out_file$original_CAS == '27342-88-7' ,'original_SMILES'] = 'CCCCCCCCCCCCO'
    out_file[out_file$original_CAS == '1318-02-1' ,'original_SMILES'] = 'O=[Al]O[Al]=O.O=[Si]=O'
    out_file[out_file$original_CAS == '70133-86-7' ,'original_SMILES'] = 'CN1CCN(CC1)CCN2C3=CC=CC=C3OCOC4=C2C=C(C=C4)Cl.Cl.Cl'
    out_file[out_file$original_CAS == '96989-24-1' ,'original_SMILES'] = 'C1=CC=C2C3C(N=CN=[N+]3C=CC2=C1)(C4=CC=C(C=C4)Cl)O.[Br-]'
    out_file[out_file$original_CAS == '70133-82-3' ,'original_SMILES'] = 'CN(C)CCCN1C2=CC=CC=C2OCOC3=C1C=C(C=C3)Cl.C(=CC(=O)O)C(=O)O'
    out_file[out_file$original_CAS == '11118-72-2' ,'original_SMILES'] = 'CC1C(OC=C2C1=C(C(=C(C2=O)C(=O)O)O)C)C'
    out_file[out_file$original_CAS == '256412-89-2' ,'original_SMILES'] = 'CC(C(=O)N(C)C1=CC=CC=C1F)OC2=CC=C(C=C2)OC3=NC4=C(O3)C=C(C=C4)Cl'
    out_file[out_file$original_CAS == '37228-47-0' ,'original_SMILES'] = 'C(CO)O.OP(O)O'
    out_file[out_file$original_CAS == '17439-94-0' ,'original_SMILES'] = 'C1CC2C(C(C1O2)C(=O)O)C(=O)O.N.N'
    out_file[out_file$original_CAS == '149961-52-4' ,'original_SMILES'] = 'CC1=CC(=C(C=C1)C)OCC2=CC=CC=C2C(=NOC)C(=O)NC'
    out_file[out_file$original_CAS == '144171-61-9' ,'original_SMILES'] = 'COC(=O)C12CC3=C(C1=NN(CO2)C(=O)N(C4=CC=C(C=C4)OC(F)(F)F)C(=O)OC)C=CC(=C3)Cl'
    out_file[out_file$original_CAS == '73606-19-6' ,'original_SMILES'] = 'C(C(C(C(F)(F)Cl)(F)F)(F)F)(C(C(OC(C(F)(F)S(=O)(=O)[O-])(F)F)(F)F)(F)F)(F)F.[K+]'
    out_file[out_file$original_CAS == '71526-07-3' ,'original_SMILES'] = 'C1CCC2(CC1)N(CCO2)C(=O)C(Cl)Cl'
    out_file[out_file$original_CAS == '205650-65-3' ,'original_SMILES'] = 'C1=C(C=C(C(=C1Cl)N2C(=C(C(=N2)C#N)C(F)(F)F)N)Cl)C(F)(F)F'
    out_file[out_file$original_CAS == '2274-99-9' ,'original_SMILES'] = 'COP(=O)(OC)SC(C1=CC=CC=C1)SP(=O)(OC)OC'
    out_file[out_file$original_CAS == '74051-80-2' ,'original_SMILES'] = 'CCCC(=NOCC)C1=C(CC(CC1=O)CC(C)SCC)O'
    out_file[out_file$original_CAS == '633-96-5' ,'original_SMILES'] = 'C1=CC=C2C(=C1)C=CC(=C2N=NC3=CC=C(C=C3)S(=O)(=O)[O-])O.[Na+]'
    out_file[out_file$original_CAS == '302578-97-8' ,'original_SMILES'] = 'C1=C(C=C(C(=C1Cl)N2C(=C(C(=N2)C#N)S(=O)C(F)(F)F)N)Cl)C(F)(F)F'
    out_file[out_file$original_CAS == '1300-21-6' ,'original_SMILES'] = 'CC(Cl)Cl'
    out_file[out_file$original_CAS == '8007-46-3' ,'original_SMILES'] = 'CC1=CCC(CC1)C(C)(C)O.CC1=CC=C(C=C1)C(C)C.CC1=CC(=C(C=C1)C(C)C)O.CC(C)C12CCC(C1C2)(C)O.CC(=CCCC(C)(C=C)O)C'
    out_file[out_file$original_CAS == '66-76-2' ,'original_SMILES'] = 'C1=CC=C2C(=C1)C(=C(C(=O)O2)CC3=C(C4=CC=CC=C4OC3=O)O)O'
    out_file[out_file$original_CAS == '71706-07-5' ,'original_SMILES'] = 'C1C(CN(CC1([N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-])([N+](=O)[O-])[N+](=O)[O-]'
    out_file[out_file$original_CAS == '302578-96-7' ,'original_SMILES'] = 'C1=C(C=C(C(=C1Cl)N2C(=C(C(=N2)C#N)S(=O)C(F)(F)F)N)Cl)C(F)(F)F'
    out_file[out_file$original_CAS == '1315501-18-8' ,'original_SMILES'] = 'CC1(C(C1C(=O)OC(C#N)C2=CC(=CC=C2)OC3=CC=CC=C3)C=C(Cl)Cl)C'
    out_file[out_file$original_CAS == '8003-19-8' ,'original_SMILES'] = 'CC(CCl)Cl.C(C=CCl)Cl'
    out_file[out_file$original_CAS == '13814-96-5' ,'original_SMILES'] = '[B-](F)(F)(F)F.[B-](F)(F)(F)F.[Pb+2]'
    out_file[out_file$original_CAS == '1686-59-5' ,'original_SMILES'] = 'CC1(CCC2C(=C1)CCC3C2(CCCC3(C)CO)C)C=C'
    out_file[out_file$original_CAS == '23060-14-2' ,'original_SMILES'] = 'CCOC(=O)CC(C(=O)OCC)S'
    out_file[out_file$original_CAS == '1686-64-2' ,'original_SMILES'] = 'CC1(CCC2C(=CCC3C2(CCCC3(C)CO)C)C1)C=C'
    out_file[out_file$original_CAS == '8046-53-5' ,'original_SMILES'] = 'C1=CC=C(C=C1)S(=O)(=O)[O-].[Na+]'
    out_file[out_file$original_CAS == '120067-83-6' ,'original_SMILES'] = 'C1=C(C=C(C(=C1Cl)N2C(=C(C(=N2)C#N)SC(F)(F)F)N)Cl)C(F)(F)F'
    out_file[out_file$original_CAS == '87392-12-9' ,'original_SMILES'] = 'CCC1=CC=CC(=C1N(C(C)COC)C(=O)CCl)C'
    out_file[out_file$original_CAS == '1112-38-5' ,'original_SMILES'] = 'COP(=S)(O)OC'
    out_file[out_file$original_CAS == '11138-66-2' ,'original_SMILES'] = 'C1=CC(=C(C=C1N)N)OCCO.Cl.Cl'
    out_file[out_file$original_CAS == '502-39-6' ,'original_SMILES'] = 'C[Hg]N=C(N)NC#N'
    out_file[out_file$original_CAS == '1207727-04-5' ,'original_SMILES'] = 'CN(C1=CC=CC(=C1F)C(=O)NC2=C(C=C(C=C2Br)C(C(F)(F)F)(C(F)(F)F)F)C(F)(F)F)C(=O)C3=CC=CC=C3'
    out_file[out_file$original_CAS == '62037-80-3' ,'original_SMILES'] = 'C(=O)(C(C(F)(F)F)(OC(C(C(F)(F)F)(F)F)(F)F)F)[O-].[NH4+]'
    out_file[out_file$original_CAS == "25474-41-3" ,'original_SMILES'] = 'CCC(C)C1=CC(=CC=C1)OC(=O)N(C)SC2=CC=CC=C2'
    out_file[out_file$original_CAS == "39300-45-3" ,'original_SMILES'] = 'CCCCCCC(C)C1=C(C(=CC(=C1)[N+](=O)[O-])[N+](=O)[O-])OC(=O)C=CC'
    out_file[out_file$original_CAS == '45187-15-3' ,'original_SMILES'] = 'C(C(C(F)(F)S(=O)(=O)[O-])(F)F)(C(F)(F)F)(F)F'
    out_file[out_file$original_CAS == '45285-51-6' ,'original_SMILES'] = 'C(=O)(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)[O-]'
    out_file[out_file$original_CAS == '70901-12-1' ,'original_SMILES'] = 'C(C(=O)O)NCP(=O)(O)[O-].[K+]'
    out_file[out_file$original_CAS == '77501-87-2' ,'original_SMILES'] = 'CC(C(=O)O)OC(=O)C1=C(C=CC(=C1)OC2=C(C=C(C=C2)C(F)(F)F)Cl)[N+](=O)[O-]'
    out_file[out_file$original_CAS == '8067-49-0' ,'original_SMILES'] = 'CCCN(CCC)C1=C(C=C(C=C1[N+](=O)[O-])C(F)(F)F)[N+](=O)[O-].CN(C)C(=O)C(C1=CC=CC=C1)C2=CC=CC=C2'
    out_file[out_file$original_CAS == "117-10-2" ,'original_SMILES'] = 'C1=CC2=C(C(=C1)O)C(=O)C3=C(C2=O)C=CC=C3O'
    out_file[out_file$original_CAS == "131-53-3" ,'original_SMILES'] = 'COC1=CC(=C(C=C1)C(=O)C2=CC=CC=C2O)O'
    out_file[out_file$original_CAS == "131-56-6" ,'original_SMILES'] = 'C1=CC=C(C=C1)C(=O)C2=C(C=C(C=C2)O)O'
    out_file[out_file$original_CAS == "14866-68-3" ,'original_SMILES'] = '[O-]Cl(=O)=O'
    out_file[out_file$original_CAS == "189-55-9" ,'original_SMILES'] = 'C1=CC=C2C3=C4C(=CC2=C1)C=CC5=CC6=CC=CC=C6C(=C54)C=C3'
    out_file[out_file$original_CAS == "191-24-2" ,'original_SMILES'] = 'C1=CC2=C3C(=C1)C4=CC=CC5=C4C6=C(C=C5)C=CC(=C36)C=C2'
    out_file[out_file$original_CAS == "192-97-2" ,'original_SMILES'] = 'C1=CC=C2C(=C1)C3=CC=CC4=C3C5=C(C=CC=C25)C=C4'
    out_file[out_file$original_CAS == "2498-66-0" ,'original_SMILES'] = 'C1=CC=C2C(=C1)C=CC3=C2C(=O)C4=CC=CC=C4C3=O'
    out_file[out_file$original_CAS == "30777-19-6" ,'original_SMILES'] = 'C1C2=CC=CC=C2C3=CC4=CC=CC=C4C=C31'
    out_file[out_file$original_CAS == "62-74-8" ,'original_SMILES'] = 'C(C(=O)[O-])F.[Na+]'
    out_file[out_file$original_CAS == "635-12-1" ,'original_SMILES'] = 'C1=CC=C2C=C3C(=O)C=CC(=O)C3=CC2=C1'
    out_file[out_file$original_CAS == "72-48-0" ,'original_SMILES'] = 'C1=CC=C2C(=C1)C(=O)C3=C(C2=O)C(=C(C=C3)O)O'
    out_file[out_file$original_CAS == "81-54-9" ,'original_SMILES'] = 'C1=CC=C2C(=C1)C(=O)C3=C(C2=O)C(=C(C=C3O)O)O'
    out_file[out_file$original_CAS == "81-61-8" ,'original_SMILES'] = 'C1=CC(=C(C2=C1C(=O)C3=C(C=CC(=C3C2=O)O)O)O)O'
    out_file[out_file$original_CAS == "81-64-1" ,'original_SMILES'] = 'C1=CC=C2C(=C1)C(=O)C3=C(C=CC(=C3C2=O)O)O'
    out_file[out_file$original_CAS == "92-24-0" ,'original_SMILES'] = 'C1=CC=C2C=C3C=C4C=CC=CC4=CC3=CC2=C1'
    
    # collected in sept 2023
    out_file[out_file$original_CAS == "104390-56-9" ,'original_SMILES'] = 'CC1=C(C(=CC=C1)C(=O)O)N(C(C)C(=O)O)C(=O)COC'
    out_file[out_file$original_CAS == "107-66-4" ,'original_SMILES'] = 'CCCCOP(=O)(O)OCCCC'
    out_file[out_file$original_CAS == "1173021-76-5" ,'original_SMILES'] = 'CCC1=CC=CC(=C1NC(=O)CS(=O)(=O)[O-])C.[Na+]'
    out_file[out_file$original_CAS == "1185255-09-7" ,'original_SMILES'] = 'COC=C(C1=CC=CC=C1OC2=NC=NC(=C2)OC3=CC=CC=C3C#N)C(=O)O'
    out_file[out_file$original_CAS == "1217465-10-5" ,'original_SMILES'] = 'CCC1=CC=CC(=C1N(C(C)C(=O)O)C(=O)C(=O)O)C'
    out_file[out_file$original_CAS == "1231244-60-2" ,'original_SMILES'] = 'CC1=C(C(=CC=C1)C)N(CN2C=CC=N2)C(=O)C(=O)O'
    out_file[out_file$original_CAS == "1246215-97-3" ,'original_SMILES'] = 'CC1=C(C(=CC=C1)C)N(CN2C=CC=N2)C(=O)CS(=O)CC(=O)O'
    out_file[out_file$original_CAS == "137361-04-7" ,'original_SMILES'] = 'C1CCC(CC1)OC(=O)CC(C(=O)OC2CCCCC2)S(=O)(=O)O'
    out_file[out_file$original_CAS == "152019-73-3" ,'original_SMILES'] = 'CCC1=CC=CC(=C1N(C(C)COC)C(=O)C(=O)O)C'
    out_file[out_file$original_CAS == "171118-09-5" ,'original_SMILES'] = 'CCC1=CC=CC(=C1N(C(C)COC)C(=O)CS(=O)(=O)O)C'
    out_file[out_file$original_CAS == "172960-62-2" ,'original_SMILES'] = 'CC1=C(C(=CC=C1)C)N(CN2C=CC=N2)C(=O)CS(=O)(=O)O'
    out_file[out_file$original_CAS == "205939-58-8" ,'original_SMILES'] = 'CC1=CSC(=C1N(C(C)COC)C(=O)CS(=O)(=O)O)C'
    out_file[out_file$original_CAS == "216667-08-2" ,'original_SMILES'] = 'CCCCCCCCCCCCCC(=O)NCCC[N+](C)(C)CCCS(=O)(=O)[O-]'
    out_file[out_file$original_CAS == "252913-85-2" ,'original_SMILES'] = 'CC(=NOCC1=CC=CC=C1C(=NOC)C(=O)O)C2=CC(=CC=C2)C(F)(F)F'
    out_file[out_file$original_CAS == "26725-51-9" ,'original_SMILES'] = 'C1=CC2=NNN=C2C(=C1)O'
    out_file[out_file$original_CAS == "28291-75-0" ,'original_SMILES'] = 'C1CCC(CC1)NC2=NC3=CC=CC=C3S2'
    out_file[out_file$original_CAS == "31468-12-9" ,'original_SMILES'] = 'CN(C)C(=O)NC1CCCCC1'
    out_file[out_file$original_CAS == "47324-98-1" ,'original_SMILES'] = 'CC[N+](CC)(CC1=CC=CC=C1)CC(=O)NC2=C(C=CC=C2C)C'

    # get a vector with CAS missing smiles
    temp_missing_SMILES = out_file[is.na(out_file$original_SMILES), 'original_CAS']
    
    # Check that all have been handled
    nrow(out_file[is.na(out_file$original_SMILES),])
    
    # Remake empty smiles NA for later (so we get the correct format on the outfile)
    out_file[ !is.na(out_file$original_SMILES) & out_file$original_SMILES == '','original_SMILES'] = NA
    
    # Get all still missing
    temp_still_missing_SMILES = out_file[is.na(out_file$original_SMILES), 'original_CAS']
    
    # Message number of still missing
    message(paste0('Number of CAS still missing SMILES after manual fix: ', length(temp_still_missing_SMILES)))
    
  }
  
  
  # Fix the SMILES column to be character
  class(out_file$original_SMILES) = 'character'
  
  
  
  ################################################################################
  # 2.4 Save (if set) and output
  ################################################################################
  
  # Remove duplicate rows
  out_file = distinct(out_file)
  
  # Make sure all output is present in the input
  out_file = out_file[out_file$original_CAS %in% identifiers_database$original_CAS,]
  
  # If we save a lookupfile, do that now
  if(save_lookupfile){
    
    # Add new entries to lookup_file
    lookup_file = distinct(rbind(lookup_file[,c("original_CAS", "original_SMILES")], out_file))
    
    # Remove all NA smiles from lookupfile
    lookup_file = lookup_file[!is.na(lookup_file$original_SMILES),]
    
    # save lookup_file
    save(lookup_file, file = local_lookupfile)
    
  } 
  
  print('Finishing add SMILES function')
  return(out_file)  
    
  
}