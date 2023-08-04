################################################################################
#    A function for fixing errors in EFSA                                      #
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

function(EFSA_database,
         settings = NULL){
  
  # Check that the input makes sense
  if(is.null(EFSA_database)|!is.data.frame(EFSA_database)){
    print('EFSA_database not provided or not a dataframe')
    print('Rerun function with correct arguments')
    return(FALSE)
  }
  
  if(is.null(settings)){
    print('No additional settings provided')
  }
  
  fix_species = F
  if(ifelse(is.vector(settings), 'fix species' %in% tolower(settings), tolower(settings) == 'fix species')){
    fix_species = T
    print('Fixing species setting (old names -> new), adding species group and removing anything in parenthesis')
  }
  
  ################################################################################
  # 2. function body
  ################################################################################
  
  print('Fixing EFSA CAS numbers')
  
  # Replace CAS with CAS collected from PubChem in 2022
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Ss]pinetoram'),'CASNO']            = '935545-74-7'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Cc]yflumetofen'),'CASNO']          = '400882-07-7'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Aa]scorbic ?[Aa]cid'),'CASNO']     = '50-81-7'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Gg]amma-?[Cc]yhalothrin'),'CASNO'] = '76703-62-3'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Pp]enflufen'),'CASNO']             = '494793-67-8'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Bb]itertanol'),'CASNO']            = '55179-31-2'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Pp]yridaben'),'CASNO']             = '96489-71-3'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Tt]riflumuron'),'CASNO']           = '64628-44-0'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Tt]embotrione'),'CASNO']           = '335104-84-2'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Pp]otassium ?[Tt]hiocyanate'),'CASNO']           = '333-20-0'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Tt]hiobencarb'),'CASNO']           = '28249-77-6'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Cc]yproconazole'),'CASNO']         = '94361-06-5'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Pp]inoxaden'),'CASNO']             = '243973-20-8'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Ss]edaxane'),'CASNO']              = '874967-67-6'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Bb]ixafen'),'CASNO']               = '581809-46-3'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Oo]range ?[Oo]il'),'CASNO']        = '60066-88-8'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Aa]minopyralid'),'CASNO']          = '150114-71-9'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Oo]ryzalin'),'CASNO']              = '19044-88-3'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Aa]misulbrom'),'CASNO']            = '348635-87-0'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Hh]alosulfuron ?[Mm]ethyl'),'CASNO']            = '100784-20-1'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Ii]pconazole'),'CASNO']            = '125225-28-7'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Pp]enthiopyrad'),'CASNO']          = '183675-82-3'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Tt]hyme ?[Oo]il'),'CASNO']         = '8007-46-3'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Cc]hromafenozide'),'CASNO']        = '143807-66-3'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Vv]aliphenal'),'CASNO']            = '283159-90-0'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Tt]riadimenol'),'CASNO']           = '55219-65-3'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Cc]hlorpyrifos'),'CASNO']          = '2921-88-2'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Ff]lumioxazine'),'CASNO']          = '103361-09-7'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Pp]yroxsulam'),'CASNO']            = '422556-08-9'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Oo]rthosulfamuron'),'CASNO']       = '213464-77-8'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Dd]isodium ?[Pp]hosphonate'),'CASNO']       = '13708-85-5'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Ff]enpyrazamine'),'CASNO']         = '473798-59-3'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Cc]yazofamid'),'CASNO']            = '120116-88-3'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '1,4-? ?[Dd]imethylnaphthalene'),'CASNO']            = '571-58-4'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Dd]iuron'),'CASNO']                = '330-54-1'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Ee]thametsulfuron ?[Mm]ethyl'),'CASNO']            = '97780-06-8'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Pp]hosphane'),'CASNO']             = '7803-51-2'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Ff]luopyram'),'CASNO']             = '658066-35-4'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Pp]otassium ?[Ii]odide'),'CASNO']  = '7681-11-0'
  EFSA_database[grepl(EFSA_database$substanceDesc, pattern = '[Ii]ndolacetic ?[Aa]cid'),'CASNO']  = '32588-36-6'
  
  #### No CAS for
  known_missing_cas = 
    c("Silver thiosulphate",
      "Potassium iodide + thiocyanate",
      "Adoxophyes orana granulovirus Strain BV-0001",
      "Bacillus firmus",
      "Potassium phosphite",
      "Kresoxim methyl BAS 494 02 F",
      "Kresoxim methyl BAS 494 04 F",
      "Quassia",
      "Carfentrazone-ethyl + IPUWG",
      "Tagetes Oil",
      "Tagetes oil",
      "Paecilomyces fumosoroseus strain FE9901",
      "Paecilomyces fumosoroseus strain Fe9901",
      "helicoverpa armigera nucleopolyhedrovirus",
      "Carfentrazone-ethyl MCPPp",
      "Candida oleophila",
      "Spinetonam")
  
  missing_cas_substDesc = unique(str_trim(str_extract(string = unique(EFSA_database[is.na(EFSA_database$CASNO), "substanceDesc"]), pattern = '^[A-Za-z0-9 \\-\\+]*(?= ?[\\(/])')))
  
  if(any(!missing_cas_substDesc %in% known_missing_cas)){
    
    message(paste0('Missing CAS handling for ', missing_cas_substDesc[!missing_cas_substDesc %in% known_missing_cas]))
    
    return(FALSE)
    
  }
  
  # Remove NA CAS (it should be the compounds above)
  EFSA_database = EFSA_database[!is.na(EFSA_database$CASNO),]
  
  
  # EFSA has the wrong CAS for Etoxazole 153233-91-1 (they had a 2 too much and called it 1523233-91-1)
  EFSA_database[EFSA_database$CASNO == '1523233-91-1', 'CASNO'] = '153233-91-1'
  
  
  
  if(fix_species){
    
    print('Fixing EFSA species names and group')
    
    #Fix EFSA_database$organism to remove all parenthesises (they contain new names for same species, which we handle in a different way)
    EFSA_database$organism = str_remove_all(EFSA_database$organism, pattern = ' ?\\(.*\\)$')
    
    # Add species group
    EFSA_database$species_group = ifelse(EFSA_database$organism %in% oecd_fish_species,
                                'Fish',
                                ifelse(EFSA_database$organism %in% oecd_algae_species,
                                       'Algae',
                                       ifelse(EFSA_database$organism %in% daphnia_species,
                                              'Crustaceans',
                                              NA)
                                )
    )
    
    # Vectors with old names
    old_raphidocelis =  c('Pseudokirchneriella subcapitata', 'Selenastrum capricornutum', 'Kirchneriella subcapitata')
    old_desmodesmus =   c('Scenedesmus gutwinskii', 'Scenedesmus spicatus', 'Scenedesmus subspicatus')
    old_danio_rerio =   c('Cyprinus rerio', 'Brachydanio rerio', 'Barilius rerio', 'Nuria rerio', 'Cyprinus chapalio', 'Perilampus striatus', 'Danio lineatus', 'Brachydanio frankei', 'Danio frankei')
    old_oryzias =       c('Poecilia latipes', 'Aplocheilus latipes', 'Oryzias latipes latipes')
    old_poecilia =      c('Acanthophacelus reticulata', 'Poecilia (Acanthophacelus) reticulata', 'Poecilia latipinna reticulata', 'Lebistes poecilioides', 'Girardinus guppii')
    old_mykiss =        c('Salmo mykiss', 'Oncorhynchus nerka mykiss', 'Parasalmo mykiss', 'Salmo purpuratus', 'Salmo penshinensis', 'Salmo gairdnerii', 'Salmo rivularis', 'Salmo kamloops whitehousei', 'Salmo irideus argentatus', 'Salmo nelsoni')
    
    # Replace old names with new
    EFSA_database$organism = ifelse(EFSA_database$organism %in% old_raphidocelis,
                                    'Raphidocelis subcapitata',
                                    ifelse(EFSA_database$organism %in% old_desmodesmus,
                                           'Desmodesmus subspicatus',
                                           ifelse(EFSA_database$organism %in% old_danio_rerio,
                                                  'Danio rerio',
                                                  ifelse(EFSA_database$organism %in% old_oryzias,
                                                         'Oryzias latipes',
                                                         ifelse(EFSA_database$organism %in% old_poecilia,
                                                                'Poecilia reticulata',
                                                                ifelse(EFSA_database$organism %in% old_mykiss,
                                                                       'Oncorhynchus mykiss',
                                                                       EFSA_database$organism
                                                                )
                                                         )
                                                  )
                                           )
                                    )
    )
    
    
  }
  
  
  print('Fixing EFSA durations')
  
  # Make all durations numeric
  EFSA_database$expDurationValue = as.numeric(EFSA_database$expDurationValue)
  EFSA_database$exposureDurationValue = as.numeric(EFSA_database$exposureDurationValue)
  
  
  ## This is the old version, in new version we just use exposureDuration in all cases
  # # Get specific duration per experiment in temp column
  # EFSA$temp_duration = ifelse(!is.na(EFSA$expDurationValue),
  #                             EFSA$expDurationValue,
  #                             EFSA$exposureDurationValue)
  # 
  # # Get specific duration unit per experiment in temp column
  # EFSA$temp_durationUnit = ifelse(!is.na(EFSA$expDurationValue),
  #                                 EFSA$expDurationUnit,
  #                                 EFSA$exposureDurationUnit)
  
  # Get duration
  EFSA_database$temp_duration = as.numeric(EFSA_database$exposureDurationValue)
  EFSA_database$temp_durationUnit = EFSA_database$exposureDurationUnit
  
  # Convert EFSA duration to hours
  EFSA_database$Duration_hour = ifelse(EFSA_database$temp_durationUnit == 'h',
                                    EFSA_database$temp_duration,
                                    ifelse(EFSA_database$temp_durationUnit == 'd',
                                           EFSA_database$temp_duration * 24,
                                           ifelse(EFSA_database$temp_durationUnit == 'wk',
                                                  EFSA_database$temp_duration * 24 * 7,
                                                  ifelse(EFSA_database$temp_durationUnit == 'mo',
                                                         EFSA_database$temp_duration * 24 * 30,
                                                         ifelse(EFSA_database$temp_durationUnit == 'min',
                                                                EFSA_database$temp_duration / 60,
                                                                EFSA_database$temp_duration
                                                         )
                                                  )
                                           )
                                    )
  )
  
  # Remove temporary columns
  EFSA_database = EFSA_database[,!grepl(colnames(EFSA_database), pattern = 'temp')]
  
  
  print('Fixing EFSA experimental values')
  
  # Make min and max values numeric
  EFSA_database$minValue = as.numeric(EFSA_database$minValue)
  EFSA_database$maxValue = as.numeric(EFSA_database$maxValue)
  
  # Calculate single conc from range fields
  EFSA_database$expValue = ifelse(!is.na(EFSA_database$minValue) & is.na(EFSA_database$maxValue),
                               EFSA_database$minValue,
                               ifelse(!is.na(EFSA_database$minValue) & !is.na(EFSA_database$maxValue),
                                      rowSums(EFSA_database[,c("minValue", "maxValue")])/2,
                                      ifelse(is.na(EFSA_database$minValue) & !is.na(EFSA_database$maxValue),
                                             EFSA_database$maxValue,
                                             NA)
                               )
  )
  
  # Get qualifier/operation of conc
  EFSA_database$expValueOp = ifelse(!is.na(EFSA_database$minValue) & is.na(EFSA_database$maxValue) & !is.na(EFSA_database$minQualifier),
                                 EFSA_database$minQualifier,
                                 ifelse(!is.na(EFSA_database$minValue) & is.na(EFSA_database$maxValue) & is.na(EFSA_database$minQualifier),
                                        '=',
                                        ifelse(is.na(EFSA_database$minValue) & !is.na(EFSA_database$maxValue) & !is.na(EFSA_database$maxQualifier),
                                               EFSA_database$maxQualifier,
                                               ifelse(!is.na(EFSA_database$minValue) & !is.na(EFSA_database$maxValue),
                                                      'range',
                                                      'UNHANDLED'
                                               )
                                        )
                                 )
  )
  
  
  
  
  
  
  
  print('Finishing EFSA cleanup')
  
  return(EFSA_database)
  
}
  