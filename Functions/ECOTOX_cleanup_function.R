################################################################################
#    A function for fixing errors in ECOTOX                                      #
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

function(ECOTOX_database,
         molweights = NULL,
         settings = NULL){
  
  # Check that the input makes sense
  if(is.null(ECOTOX_database)|!is.data.frame(ECOTOX_database)){
    print('ECOTOX_database not provided or not a dataframe')
    print('Rerun function with correct arguments')
    return(FALSE)
  }
  
  # Handling of molweights reference
  mol_type = NA
  if(is.null(molweights)){
    print('No molweights provided')
    print('Do you want to continue without reference (will take time) - answer "yes" or "no"')
    mol_type = 'none'
    temp_answer = readline()
    if(tolower(temp_answer) == 'no'){
      return(FALSE)
    } else {
      print('Continuing without reference - will search for molweights')
    }
  } else if(is.character(molweights)){
    
    print('Molweights provided as reference file (location should be working directory)')
    mol_type = 'file'
    
  } else if(is.data.frame(molweights)){
    
    print('Molweights provided as dataframe')
    mol_type = 'dataframe'
    
  }
  
  print('Parsing settings')
  if(is.null(settings)){
    print('No additional settings provided, only fixing CAS-numbers')
  }
  
  fix_species = F
  if(ifelse(is.vector(settings), 'fix species' %in% tolower(settings), tolower(settings) == 'fix species')){
    fix_species = T
    print('Fixing species names (old names -> new)')
  }
  
  fix_molweights = F
  if(ifelse(is.vector(settings), 'fix molweights' %in% tolower(settings), tolower(settings) == 'fix molweights')){
    fix_molweights = T
    print('Adding molweights, and calculating molar units into mg per l')
  }
  
  molweight_search = F
  if(ifelse(is.vector(settings), 'search molweights' %in% tolower(settings), tolower(settings) == 'search molweights')){
    molweight_search = T
    print('Searching for molweights')
  }
  
  ################################################################################
  # 2. function body
  ################################################################################
  
  print('Performing function tasks')
  
  # Remove NA CAS 
  ECOTOX_database = ECOTOX_database[!is.na(ECOTOX_database$cas_number),]
  
  # Remove '*' in conc1_mean field
  ECOTOX_database$conc1_mean_original = str_replace(ECOTOX_database$conc1_mean, pattern = '\\*', replacement = '')
  
  ## Transform con1_mean into mg/L
  # Some units have ben commented out either due to incompatible conversions or lack of molecular weight data
  ECOTOX_database$mgperL = NA
  
  ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'g/L', as.numeric(ECOTOX_database$conc1_mean) * 1000, ECOTOX_database$mgperL)
  ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'mg/L', as.numeric(ECOTOX_database$conc1_mean), ECOTOX_database$mgperL)
  ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'ug/L', as.numeric(ECOTOX_database$conc1_mean) / 1000, ECOTOX_database$mgperL)
  ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'ng/L', as.numeric(ECOTOX_database$conc1_mean) / 1000^2, ECOTOX_database$mgperL)
  ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'ug/mL', as.numeric(ECOTOX_database$conc1_mean), ECOTOX_database$mgperL)
  
  ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'g/l', as.numeric(ECOTOX_database$conc1_mean) * 1000, ECOTOX_database$mgperL)
  ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'mg/l', as.numeric(ECOTOX_database$conc1_mean), ECOTOX_database$mgperL)
  ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'ug/l', as.numeric(ECOTOX_database$conc1_mean) / 1000, ECOTOX_database$mgperL)
  ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'ng/l', as.numeric(ECOTOX_database$conc1_mean) / 1000^2, ECOTOX_database$mgperL)
  ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'ug/ml', as.numeric(ECOTOX_database$conc1_mean), ECOTOX_database$mgperL)
  
  ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'AI g/L', as.numeric(ECOTOX_database$conc1_mean) * 1000, ECOTOX_database$mgperL)
  ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'AI mg/L', as.numeric(ECOTOX_database$conc1_mean), ECOTOX_database$mgperL)
  ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'AI ug/L', as.numeric(ECOTOX_database$conc1_mean) / 1000, ECOTOX_database$mgperL)
  ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'AI ng/L', as.numeric(ECOTOX_database$conc1_mean) / 1000^2, ECOTOX_database$mgperL)
  ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'AI ug/mL', as.numeric(ECOTOX_database$conc1_mean), ECOTOX_database$mgperL)
  
  ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'AI g/l', as.numeric(ECOTOX_database$conc1_mean) * 1000, ECOTOX_database$mgperL)
  ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'AI mg/l', as.numeric(ECOTOX_database$conc1_mean), ECOTOX_database$mgperL)
  ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'AI ug/l', as.numeric(ECOTOX_database$conc1_mean) / 1000, ECOTOX_database$mgperL)
  ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'AI ng/l', as.numeric(ECOTOX_database$conc1_mean) / 1000^2, ECOTOX_database$mgperL)
  ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'AI ug/ml', as.numeric(ECOTOX_database$conc1_mean), ECOTOX_database$mgperL)
  
  
  # mole-based units conversion requires molweights
  if(fix_molweights){
    
    # if conc1_unit contain any molar unit, collect molweights from webchem
    ECOTOX_unique_cas = unique(ECOTOX_database$cas_number)
    
    # Add molweights (skipped if we have a molweight dump file, ran if we force a rerun)
    if(mol_type != 'file'){
      
      # If we have manually loaded the molweight dump file we use that as a basis to avoid overusing cir_query
      if(mol_type == 'dataframe'){
        
        ECOTOX_database$molweight = NULL
        
        ECOTOX_database = merge(ECOTOX_database, molweights, by = 'cas_number', all.x = T)
        
      } else {
        
        ECOTOX_database$molweight = NA
        
      }
      
      if(molweight_search){
      
        print('Sreaching for molweights with cir_query')
        
        for(i in 1:length(ECOTOX_unique_cas)){
          
          current_cas = ECOTOX_unique_cas[i]
          
          current_cas_unit_subset = ECOTOX_database[ECOTOX_database$cas_number == current_cas & !is.na(ECOTOX_database$cas_number), 'conc1_unit']
          
          
          if(any(grepl(current_cas_unit_subset, pattern = 'M', ignore.case = F)) & any(is.na(ECOTOX_database[ECOTOX_database$cas_number == current_cas & !is.na(ECOTOX_database$cas_number), 'molweight']))){
            
            current_molweight = cir_query(identifier = current_cas, representation = 'mw', match = 'first', resolver = 'cas_number')
            
            ECOTOX_database[ECOTOX_database$cas_number == current_cas, 'molweight'] = current_molweight[1,2]
            
          }
          
          if(i%%100 == 0){
            
            print(paste0(i, ' unique cas processed out of ', length(ECOTOX_unique_cas)))
          }
          
        }
        
        molweight_dump = distinct(ECOTOX_database[!is.na(ECOTOX_database$molweight), c("cas_number", "molweight")])
        
        print('Save molweights reference? - "no" or filename')
        temp_input = readlines()
        
        if(tolower(temp_input) != 'no'){
          
          print(paste0('Saving molweight reference as ', temp_input))
          save(molweight_dump, file = temp_input)
          
        } else {
          
          print('Not saving molweight reference')
          
        }  
        
      } else {
        
        print('Not searching for additional molweights')
        
      }
      
      
      
      
    } else {
      
      load(molweights)
      
      ECOTOX_database$molweight = NULL
      
      ECOTOX_database = merge(ECOTOX_database, molweight_dump, by = 'cas_number', all.x = T)
      
      if(molweight_search){
        
        print('Sreaching for molweights with cir_query')
        
        for(i in 1:length(ECOTOX_unique_cas)){
          
          current_cas = ECOTOX_unique_cas[i]
          
          current_cas_unit_subset = ECOTOX_database[ECOTOX_database$cas_number == current_cas & !is.na(ECOTOX_database$cas_number), 'conc1_unit']
          
          
          if(any(grepl(current_cas_unit_subset, pattern = 'M', ignore.case = F)) & any(is.na(ECOTOX_database[ECOTOX_database$cas_number == current_cas & !is.na(ECOTOX_database$cas_number), 'molweight']))){
            
            current_molweight = cir_query(identifier = current_cas, representation = 'mw', match = 'first', resolver = 'cas_number')
            
            ECOTOX_database[ECOTOX_database$cas_number == current_cas, 'molweight'] = current_molweight[1,2]
            
          }
          
          if(i%%100 == 0){
            
            print(paste0(i, ' unique cas processed out of ', length(ECOTOX_unique_cas)))
          }
          
        }
        
        molweight_dump = distinct(ECOTOX_database[!is.na(ECOTOX_database$molweight), c("cas_number", "molweight")])
        
        save(molweight_dump, file = molweights)
        
      } else {
        
        print('Not searching for additional molweights')
        
      }
      
      
    }
    
    # Molar based units
    ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'M', as.numeric(ECOTOX_database$conc1_mean) * as.numeric(ECOTOX_database$molweight) * 1000, ECOTOX_database$mgperL)
    ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'mol/L', as.numeric(ECOTOX_database$conc1_mean) * as.numeric(ECOTOX_database$molweight) * 1000, ECOTOX_database$mgperL)
    ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'mol/l', as.numeric(ECOTOX_database$conc1_mean) * as.numeric(ECOTOX_database$molweight) * 1000, ECOTOX_database$mgperL)
    ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'mM', as.numeric(ECOTOX_database$conc1_mean) * as.numeric(ECOTOX_database$molweight), ECOTOX_database$mgperL)
    ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'mmol/L', as.numeric(ECOTOX_database$conc1_mean) * as.numeric(ECOTOX_database$molweight), ECOTOX_database$mgperL)
    ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'mmol/l', as.numeric(ECOTOX_database$conc1_mean) * as.numeric(ECOTOX_database$molweight), ECOTOX_database$mgperL)
    ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'mmol/mL', as.numeric(ECOTOX_database$conc1_mean) * as.numeric(ECOTOX_database$molweight) / 1000, ECOTOX_database$mgperL)
    ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'mmol/ml', as.numeric(ECOTOX_database$conc1_mean) * as.numeric(ECOTOX_database$molweight) / 1000, ECOTOX_database$mgperL)
    ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'uM', as.numeric(ECOTOX_database$conc1_mean) * as.numeric(ECOTOX_database$molweight) / 1000, ECOTOX_database$mgperL)
    ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'umol/L', as.numeric(ECOTOX_database$conc1_mean) * as.numeric(ECOTOX_database$molweight) / 1000, ECOTOX_database$mgperL)
    ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'umol/l', as.numeric(ECOTOX_database$conc1_mean) * as.numeric(ECOTOX_database$molweight) / 1000, ECOTOX_database$mgperL)
    ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'umol/mL', as.numeric(ECOTOX_database$conc1_mean) * as.numeric(ECOTOX_database$molweight), ECOTOX_database$mgperL)
    ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'umol/ml', as.numeric(ECOTOX_database$conc1_mean) * as.numeric(ECOTOX_database$molweight), ECOTOX_database$mgperL)
    ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'nM', as.numeric(ECOTOX_database$conc1_mean) * as.numeric(ECOTOX_database$molweight) / 1000^2, ECOTOX_database$mgperL)
    
    ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'AI M', as.numeric(ECOTOX_database$conc1_mean) * as.numeric(ECOTOX_database$molweight) * 1000, ECOTOX_database$mgperL)
    ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'AI mol/L', as.numeric(ECOTOX_database$conc1_mean) * as.numeric(ECOTOX_database$molweight) * 1000, ECOTOX_database$mgperL)
    ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'AI mM', as.numeric(ECOTOX_database$conc1_mean) * as.numeric(ECOTOX_database$molweight), ECOTOX_database$mgperL)
    ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'AI mmol/L', as.numeric(ECOTOX_database$conc1_mean) * as.numeric(ECOTOX_database$molweight), ECOTOX_database$mgperL)
    ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'AI mmol/mL', as.numeric(ECOTOX_database$conc1_mean) * as.numeric(ECOTOX_database$molweight) / 1000, ECOTOX_database$mgperL)
    ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'AI uM', as.numeric(ECOTOX_database$conc1_mean) * as.numeric(ECOTOX_database$molweight) / 1000, ECOTOX_database$mgperL)
    ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'AI umol/L', as.numeric(ECOTOX_database$conc1_mean) * as.numeric(ECOTOX_database$molweight) / 1000, ECOTOX_database$mgperL)
    ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'AI umol/mL', as.numeric(ECOTOX_database$conc1_mean) * as.numeric(ECOTOX_database$molweight), ECOTOX_database$mgperL)
    ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'AI nM', as.numeric(ECOTOX_database$conc1_mean) * as.numeric(ECOTOX_database$molweight) / 1000^2, ECOTOX_database$mgperL)
    
  }
  
  
  
  ## Percentage based units
  #ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == '%', as.numeric(ECOTOX_database$conc1_mean) * 10000, ECOTOX_database$mgperL)
  #ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'ppm', as.numeric(ECOTOX_database$conc1_mean), ECOTOX_database$mgperL)
  #ECOTOX_database$mgperL <- ifelse(ECOTOX_database$conc1_unit == 'ppb', as.numeric(ECOTOX_database$conc1_mean) / 1000, ECOTOX_database$mgperL)
  
  # Remove rows with NA mgperL
  ECOTOX_database = ECOTOX_database[!is.na(ECOTOX_database$mgperL),]
  
  
  
  
  if(fix_species){
    
    # Vectors with old names
    old_raphidocelis =  c('Pseudokirchneriella subcapitata', 'Selenastrum capricornutum', 'Kirchneriella subcapitata')
    old_desmodesmus =   c('Scenedesmus gutwinskii', 'Scenedesmus spicatus', 'Scenedesmus subspicatus')
    old_danio_rerio =   c('Cyprinus rerio', 'Brachydanio rerio', 'Barilius rerio', 'Nuria rerio', 'Cyprinus chapalio', 'Perilampus striatus', 'Danio lineatus', 'Brachydanio frankei', 'Danio frankei')
    old_oryzias =       c('Poecilia latipes', 'Aplocheilus latipes', 'Oryzias latipes latipes')
    old_poecilia =      c('Acanthophacelus reticulata', 'Poecilia (Acanthophacelus) reticulata', 'Poecilia latipinna reticulata', 'Lebistes poecilioides', 'Girardinus guppii')
    old_mykiss =        c('Salmo mykiss', 'Oncorhynchus nerka mykiss', 'Parasalmo mykiss', 'Salmo purpuratus', 'Salmo penshinensis', 'Salmo gairdnerii', 'Salmo rivularis', 'Salmo kamloops whitehousei', 'Salmo irideus argentatus', 'Salmo nelsoni')
    
    # Replace old names with new
    ECOTOX_database$latin_name = ifelse(ECOTOX_database$latin_name %in% old_raphidocelis,
                                      'Raphidocelis subcapitata',
                                      ifelse(ECOTOX_database$latin_name %in% old_desmodesmus,
                                             'Desmodesmus subspicatus',
                                             ifelse(ECOTOX_database$latin_name %in% old_danio_rerio,
                                                    'Danio rerio',
                                                    ifelse(ECOTOX_database$latin_name %in% old_oryzias,
                                                           'Oryzias latipes',
                                                           ifelse(ECOTOX_database$latin_name %in% old_poecilia,
                                                                  'Poecilia reticulata',
                                                                  ifelse(ECOTOX_database$latin_name %in% old_mykiss,
                                                                         'Oncorhynchus mykiss',
                                                                         ECOTOX_database$latin_name
                                                                  )
                                                           )
                                                    )
                                             )
                                      )
    )
    
    
  }
  
  # calculate durations into hours (if days etc)
  ECOTOX_database$Duration_hour = NA
  
  ECOTOX_database$obs_duration_mean = as.numeric(ECOTOX_database$obs_duration_mean)
  
  ECOTOX_database$Duration_hour = ifelse(ECOTOX_database$obs_duration_unit == 'h',
                                         ECOTOX_database$obs_duration_mean,
                                         ifelse(ECOTOX_database$obs_duration_unit == 'd',
                                                ECOTOX_database$obs_duration_mean * 24,
                                                ifelse(ECOTOX_database$obs_duration_unit == 'wk',
                                                       ECOTOX_database$obs_duration_mean * 24 * 7,
                                                       ifelse(ECOTOX_database$obs_duration_unit == 'mo',
                                                              ECOTOX_database$obs_duration_mean * 24 * 7 * 30,
                                                              ECOTOX_database$Duration_hour))))
  
  
  
  

  return(ECOTOX_database)
  
}
