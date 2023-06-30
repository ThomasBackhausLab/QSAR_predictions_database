
################################################################################
#    A function for reducing output from VEGA, ECOSAR and T.E.S.T              #
#   to species specific subsets with reduced model outputs (depending on QSAR) #
################################################################################
#
# Original author: Patrik Svedberg
#
# Based on the output of the function QSAR_processing_function.R
# 
#
################################################################################
#                  Table of contents                                           #
# 
# 0. Planned changes and WIP
# 1. Initiating function and parsing arguments
# 2. Subsetting based on filters
# 3. Wide in_format processing
#   3.1. VEGA processing
#   3.2. ECOSAR processing
#   3.3. Fixing some TEST naming issues
#   3.4. Return data.frame
# 4. Long in_format processing
#   4.1. VEGA processing
#   4.2. ECOSAR processing
#   4.3. Return long data.frame

################################################################################
# 0. Planned changes and WIP
################################################################################
#
# Add species handling
#
# Add method alternatives for ECOSAR
#
# Make the function able to handle different TEST outputs


################################################################################
#         0.1. Versions and changelog                                          #
################################################################################

## Version 1.0

# While some changes has been done prior, versions were not controlled.
# Version 1.0 only handles the output of all QSARs, both wide and long format output

## Version 2.0 changes

# Added option to only process specific QSARs based upon which QSAR data is present in the indata
# But only for long format

## Version 3.0 changes

# Added model description differences for calculations - i.e. geo_mean_no_baseline for ecosar
# Added option to do second calculation methods for vega predictions
# changed how filter works to improve customization (and some minor following changes)



################################################################################
# 1. Initiating function and parsing arguments
################################################################################


function(QSAR_output, 
         species, 
         in_format = "long",
         out_format = "long",
         endpoint = 'acute',
         vega_filter = c('good', 'moderate'),
         vega_filter2 = NULL,
         vega_method = c('geo_mean'), 
         vega_method2 = NULL,
         ecosar_method = c('low'),
         ecosar_method2 = NULL,
         data_filter = 'all'){
  #' A helping function to run and compile QSAR models
  #'
  #' @Description This function reduces the output of the QSAR_processing_function
  #' to a single prediction per model for specified species 
  #' (currently only working for fish and daphnia)
  #' 
  #' @param QSAR_output vector of strings. A vector with paths to the 
  #' executables of the QSAR models
  #' 
  #' @param species string specifying the target species
  #' 
  #' @param in_format string, either "long" or "wide" depending on how the information
  #' is arranged in the input frame
  #' 
  #' @param out_format string, either "long" or "wide" depending on how the information
  #' should be used. Long format allows keeping qualty of predictions (vega) and
  #' multiple prediction classes (ecosar), and is the correct format for ggplot. 
  #' Wide format is more readable for the human eye and is good when the 
  #' predictions should be used as is.
  #' 
  #' @param endpoint string.  'acute' or 'chronic' or 'all', , specifying which endpoints
  #' should be included in the filter, default is 'acute'.
  #' 
  #' @param vega_method vector of strings specifying which reliability of the output of VEGA to keep.
  #' filter options 'experimental', 'low', 'moderate', 'good'
  #' calculation options 'mean', 'geo_mean' and 'lowest'.
  #' Default is c('good', 'moderate', 'geo_mean') 
  #' 
  #' @param vega_method2 vector of strings. The same options as vega_method, 
  #' and used if you want mutliple handlings of VEGA calculations in the same output frame
  #' 
  #' @param ecosar_method vector of strings. Specifying how to handle the output
  #' of ECOSAR. 
  #' Calculation options 'low', 'mean', 'geo_mean', 'drop baseline'
  #' Default is 'low'
  #' 
  #' @param ecosar_method2 vector of strings. NOT IMPLEMENTED RIGHT NOW.
  #' Same options as ecosar_method, used if you want multiple handlings of ecosar calculations.
  #' 
  #' @param data_filter string. How to handle filtered-out data (quality filter),
  #' or raw data. Options: 'all', 'quality', 'calculated'
  #' 
  #' COMPLETE THIS HELP TEXT!
  
  ## Print time of function start
  start_time = Sys.time()
  print(paste0('Starting script at ', start_time))
  
  print('Processing input arguments')
  
  ## Set function version
  version = 2
  
  # Ensure lower case in formats arguments
  in_format = tolower(in_format)
  out_format = tolower(out_format)
  
  if (in_format != out_format){
    print('####################################################################')
    print('I would recommend using the same informat as out format')
    print('Format swapping in the processing function may lead to loss of information')
    print('####################################################################')
  }
  
  if(in_format == 'long'){
    # Get which QSARs are present in the data
    QSAR_vector = as.vector(unique(QSAR_output$QSAR_tool))
    
    # Set boolean depending on presence
    if('vega' %in% tolower(QSAR_vector)){
      process_vega = T
    } else {
      process_vega = F
    }
    
    if('ecosar' %in% tolower(QSAR_vector)){
      process_ecosar = T
    } else {
      process_ecosar = F
    }
    
    if('test' %in% tolower(QSAR_vector) | 't.e.s.t.' %in% tolower(QSAR_vector)){
      process_test = T
    } else {
      process_test = F
    }
  }
  
  # Parse data_filter
  if(!tolower(data_filter) %in% c('all', 'quality', 'calculated')){
    
    message('data_filter not processable. Use default ("all") or specify: "all", "quality" or "calculated"')
    
  }
  
  # Double run vega or ecosar
  vega_dual_run = F
  if(!is.null(vega_filter2)){
    vega_dual_run = T
  }
  ecosar_dual_run = F
  if(!is.null(ecosar_method2)){
    ecosar_dual_run = T
  }
  
  
  # Set column order (they can be out of order based on identifiers fed to the QSAR_processing script)
  
  column_order = c("original_CAS", 
                   "ecosar_CAS", 
                   "test_CAS", 
                   "original_SMILES",
                   "ecosar_SMILES",
                   "vega_SMILES",
                   "test_SMILES",
                   "InChIKey",
                   "model",
                   "value",
                   "reliability",
                   "duration",
                   "model_organism",
                   "model_endpoint",
                   "QSAR_tool",
                   "Chemical",
                   "ECOSAR Class",
                   "Max Log Kow",
                   "error_warning",
                   "Molecular Weight (g/mol)",
                   "Log Kow",
                   "Water Solubility (mg/L)",
                   "Melting Point (C)",
                   "ecosar_CAS_fixed",
                   "original_CAS_short",
                   "model_duration (h)")
  
  
  
  QSAR_output[,column_order[!column_order %in% colnames(QSAR_output)]] = NA
  
  QSAR_output = QSAR_output[column_order]
  
  ################################################################################
  # 2. Subsetting based on filters
  ################################################################################
  
  print('Subsetting input frame based on filter settings')
  
  # Set endpoint filters
  current_endpoint_filter = '(FAILSAFE)'
  if(grepl(x = endpoint, pattern = '[Aa]cute')){
    
    current_endpoint_filter = paste0(current_endpoint_filter, '|(EC50)')
    
  } else if(grepl(x = endpoint, pattern = '[Cc]hronic')){
    
    current_endpoint_filter = paste0(current_endpoint_filter, '|(NOEC)')
    
  } else if(grepl(x = endpoint, pattern = '([Aa]ll)|([Bb]oth)')){
    
    current_endpoint_filter = paste0(current_endpoint_filter, '|(EC50)|(NOEC)')
  
  }
  
  current_species = species
  
  # Special salt water/marine considerations
  salt_water = F
  if(grepl(species, pattern = '(SW)|([Mm]arine)')){
    salt_water = T
  }
  
  current_duration_filter = '(FAILSAFE)'
  
  # Set durations for acute tests
  if(grepl(current_endpoint_filter, pattern = 'EC50')){
    
    if(grepl(x = current_species, pattern = '[Dd]aphn')){
      
      current_duration_filter = paste0(current_duration_filter, '|(48h)')
      
    } else if(grepl(x = current_species, pattern = '[Ff]ish')){
      
      current_duration_filter = paste0(current_duration_filter, '|(96h)')
      
    } else if(grepl(x = current_species, pattern = '[Aa]lgae')){
      
      current_duration_filter = paste0(current_duration_filter, '|(72h)|(96h)')
      
    } else if(grepl(x = current_species, pattern = '[Mm]ysid')){
      
      current_duration_filter = paste0(current_duration_filter, '|(96h)')
      
    }
    
  }
  
  # Set durations for chronic tests (ECOSAR does not have durations for chonic tests)
  if(grepl(current_endpoint_filter, pattern = 'NOEC')){
    
    if(grepl(x = current_species, pattern = '[Dd]aphn')){
      
      current_duration_filter = paste0(current_duration_filter, '|(21d)')
      
    } else if(grepl(x = current_species, pattern = '[Ff]ish')){
      
      current_duration_filter = paste0(current_duration_filter, '|(ELS)|(NA)')
      
    } else if(grepl(x = current_species, pattern = '[Aa]lgae')){
      
      current_duration_filter = paste0(current_duration_filter, '|(72h)')
      
    }
    
  }
  
  if(in_format == 'wide'){
    
    # This isnt updated, so we drop it...
    
    print('Wide format input not supported in current version, exiting function')
    
    return(FALSE)
    
    
    ## This is not run
    
    current_subset = QSAR_output[,(grepl(colnames(QSAR_output),pattern = paste0(tolower(substr(current_species, 1, 4)), '|', str_to_sentence(substr(current_species, 1, 4)))) & 
                                     grepl(colnames(QSAR_output), pattern = current_duration_filter) &
                                     grepl(colnames(QSAR_output), pattern = paste0(current_endpoint_filter,'|(ECOSAR)|(TEST)') )) | 
                                   colnames(QSAR_output) %in% c("original_CAS", 
                                                                "original_SMILES",
                                                                "InChIKey",
                                                                "ecosar_CAS",
                                                                "test_CAS", 
                                                                "ecosar_SMILES",
                                                                "vega_SMILES",
                                                                "test_SMILES" )]
    
    
  } else if (in_format == 'long'){
    
    if(grepl(current_endpoint_filter, pattern = 'NOEC')){
      
      current_subset = QSAR_output[grepl(tolower(QSAR_output$model_organism), pattern = tolower(substr(current_species, 1, 4))) &
                                     if(salt_water) grepl(tolower(QSAR_output$model_organism), pattern = tolower('SW')) else !grepl(tolower(QSAR_output$model_organism), pattern = tolower('SW')) &
                                     grepl(QSAR_output$model_endpoint, pattern = current_endpoint_filter) &
                                     (grepl(QSAR_output$duration, pattern = current_duration_filter)|
                                        is.na(QSAR_output$duration)),]
      
    } else {
    
      current_subset = QSAR_output[grepl(tolower(QSAR_output$model_organism), pattern = tolower(substr(current_species, 1, 4))) &
                                     if(salt_water) grepl(tolower(QSAR_output$model_organism), pattern = tolower('SW')) else !grepl(tolower(QSAR_output$model_organism), pattern = tolower('SW')) &
                                     grepl(QSAR_output$model_endpoint, pattern = current_endpoint_filter) &
                                     grepl(QSAR_output$duration, pattern = current_duration_filter),]
      
      
    }
    
    
    
    
  }
  
  # Simplify current duration and current endpoint
  current_duration = 
    str_replace_all(
      str_replace_all(
        str_replace_all(
          current_duration_filter, 
          pattern = '\\|',
          replacement = ' '),
        pattern = '[\\)\\(]', 
        replacement = ''),
      pattern = 'FAILSAFE',
      replacement = '')
  
  current_endpoint = 
    str_trim(
      str_replace_all(
        str_replace_all(
          str_replace_all(
            current_endpoint_filter, 
            pattern = '\\|',
            replacement = ' '),
          pattern = '[\\)\\(]', 
          replacement = ''),
        pattern = 'FAILSAFE',
        replacement = '')
    )
  
  # Make sure we have data left, otherwise return NULL
  if(nrow(current_subset) == 0){
    message('Every data point filtered out! Tweak filters or check input file. Exiting function.')
    return(NULL)
  }
  
  ################################################################################
  # 3. Wide in_format processing
  ################################################################################
  # 
  # if(in_format == 'wide'){
  #   
  #   print('Wide format processing')
  #   
  #   ################################################################################
  #   # 3.1. VEGA processing
  #   ################################################################################
  #   
  #   print('Processing vega first method')
  #   
  #   subset_filter = '(FAILSAFE)'
  #   
  #   # Invert the filters (translate from "what to keep" to "what to remove")
  #   if('all' %in% vega_method){
  #     inverted_filter = c()
  #   } else {
  #     filter_inverter_vector = c('experimental', 'low', 'moderate', 'good')
  #     inverted_filter = filter_inverter_vector[!filter_inverter_vector %in% tolower(vega_method)]
  #   }
  #   
  #   # Produce a filter string regexp
  #   if('experimental' %in% inverted_filter){
  #     subset_filter = paste(subset_filter, '(EXPERIMENTAL)', sep = '|') 
  #   }
  #   if('low' %in% inverted_filter){
  #     subset_filter = paste(subset_filter, '(low)', sep = '|') 
  #   }
  #   if('moderate' %in% inverted_filter){
  #     subset_filter = paste(subset_filter, '(moderate)', sep = '|') 
  #   }
  #   if('good' %in% inverted_filter){
  #     subset_filter = paste(subset_filter, '(good)', sep = '|')
  #   }
  #   
  #   # Get logical with TRUE for experimental or low quality
  #   temp_replacement_matrix = mapply(current_subset, FUN = function(x){grepl(x = x, pattern = subset_filter) & !is.na(x)})
  #   
  #   # Remove colnames otherwise it's weird
  #   colnames(temp_replacement_matrix) = NULL
  #   
  #   # remove first column and add last (shift the logicals one step)
  #   temp_replacement_matrix = temp_replacement_matrix[,-1]
  #   temp_replacement_matrix = cbind(temp_replacement_matrix, matrix(data = FALSE, nrow = nrow(temp_replacement_matrix)))
  #   
  #   # Replace true positions with NA
  #   current_subset[temp_replacement_matrix] = NA
  #   
  #   # Save all columns with reliability and drop them (keeping only predicted values)
  #   current_vega_reliabilities = current_subset[,grepl(colnames(current_subset),pattern = 'reliability') | colnames(current_subset) %in% c("original_CAS", 
  #                                                                                                                                          "original_SMILES",
  #                                                                                                                                          "InChIKey",
  #                                                                                                                                          "ecosar_CAS",
  #                                                                                                                                          "test_CAS", 
  #                                                                                                                                          "ecosar_SMILES",
  #                                                                                                                                          "vega_SMILES",
  #                                                                                                                                          "test_SMILES" )]
  #   current_subset = current_subset[,!grepl(colnames(current_subset),pattern = 'reliability')]
  #   
  #   # Make all predictions numeric
  #   current_subset[,9:ncol(current_subset)] = sapply(current_subset[,9:ncol(current_subset)], as.numeric)
  #   
  #   
  #   # define number of vega predictions used for mean
  #   current_subset$VEGA_number_of_quality_predictions = rowSums(!is.na(current_subset[,grepl(x = colnames(current_subset), pattern = 'VEGA')]))
  #   
  #   current_subset = current_subset %>% relocate(VEGA_number_of_quality_predictions, .after = InChIKey)
  #   
  #   colnames(current_subset)[colnames(current_subset) == 'VEGA_number_of_quality_predictions'] = paste0('VEGA_number_of_quality_predictions_', current_species, '_', current_duration)
  #   
  #   #If geo_mean chosen as method
  #   if('geo_mean' %in% vega_method){
  #     current_subset$`VEGA_mean_prediction [mg/l]` = apply(
  #       current_subset[,grepl(x = colnames(current_subset), pattern = 'VEGA_[^Nn]')], 
  #       FUN = function(x){
  #         geom_mean = exp(
  #           sum(
  #             log(
  #               x[x > 0 & !is.na(x)]
  #             ), 
  #             na.rm = T
  #           ) / length(x[x > 0 & !is.na(x)])
  #         )
  #         if(all(is.na(x) | x == 0)){
  #           return(NA)
  #         } else {
  #           return(geom_mean)
  #         }
  #         
  #       }, 
  #       MARGIN = 1)
  #     
  #     current_subset = current_subset %>% relocate(`VEGA_mean_prediction [mg/l]`, .after = InChIKey)
  #     
  #     colnames(current_subset)[colnames(current_subset) == 'VEGA_mean_prediction [mg/l]'] = paste0('VEGA_mean_prediction_', current_species, '_', current_duration, ' [mg/l]')
  #     
  #   } else if('mean' %in% vega_method){
  #     
  #     current_subset$`VEGA_mean_prediction [mg/l]` = apply(
  #       current_subset[,grepl(x = colnames(current_subset), pattern = 'VEGA_[^Nn]')], 
  #       FUN = function(x){
  #         mean_calc = mean(x[x > 0 & !is.na(x)])
  #         if(all(is.na(x) | is.infinite(x) | x == 0)){
  #           return(NA)
  #         } else {
  #           return(mean_calc)
  #         }
  #       }, 
  #       MARGIN = 1)
  #     
  #     current_subset = current_subset %>% relocate(`VEGA_mean_prediction [mg/l]`, .after = InChIKey)
  #     
  #     colnames(current_subset)[colnames(current_subset) == 'VEGA_mean_prediction [mg/l]'] = paste0('VEGA_mean_prediction_', current_species, '_', current_duration, ' [mg/l]')
  #     
  #   } else if('lowest' %in% vega_method){
  #     
  #     current_subset$`VEGA_lowest_prediction [mg/l]` = apply(
  #       current_subset[,grepl(x = colnames(current_subset), pattern = 'VEGA_[^Nn]')], 
  #       FUN = function(x){
  #         min_pred = suppressWarnings(min(x[x > 0 & !is.na(x)]))
  #         if(all(is.na(x) | x == 0)){
  #           return(NA)
  #         } else {
  #           return(min_pred)
  #         }
  #       }, 
  #       MARGIN = 1)
  #     
  #     current_subset = current_subset %>% relocate(`VEGA_lowest_prediction [mg/l]`, .after = InChIKey)
  #     
  #     colnames(current_subset)[colnames(current_subset) == 'VEGA_lowest_prediction [mg/l]'] = paste0('VEGA_lowest_prediction_', current_species, '_', current_duration, ' [mg/l]')
  #     
  #   }
  #   
  #   if(data_filter == 'calculated'){
  #     
  #     current_subset = current_subset[,!grepl(x = colnames(current_subset), pattern = 'VEGA_[^nml]')]
  #     
  #   }
  #   
  #   ################################################################################
  #   # 3.2. ECOSAR processing
  #   ################################################################################
  #   
  #   print('processing ecosar first setting')
  #   
  #   # define number of ECOSAR predictions
  #   current_subset$ECOSAR_number_of_prediction_classes = rowSums(!is.na(current_subset[,grepl(x = colnames(current_subset), pattern = 'ECOSAR')]))
  #   
  #   current_subset = current_subset %>% relocate(ECOSAR_number_of_prediction_classes, .after = InChIKey)
  #   
  #   colnames(current_subset)[colnames(current_subset) == 'ECOSAR_number_of_prediction_classes'] = paste0('ECOSAR_number_of_prediction_classes_', current_species, '_', current_duration)
  #   
  #   
  #   # ECOSAR states in their user manual that "Traditionally, in the absence of adequate measured data, the most conservative effect level is
  #   # used when predictions are identified from multiple classes", so pick out the most conservative (lowest) value if there are multiple prediction classes
  #   
  #   
  #   if('low' %in% ecosar_method){
  #     
  #     current_subset$`ECOSAR_low [mg/l]` = apply(
  #       current_subset, 
  #       MARGIN = 1, 
  #       FUN = function(x){
  #         
  #         temp = x[grepl(x = colnames(current_subset), pattern = 'ECOSAR_[^n]')]
  #         
  #         temp2 = suppressWarnings(min(as.numeric(temp), na.rm = T))
  #         
  #         if(!is.infinite(temp2)){
  #           return(temp2)
  #         } else {
  #           return(NA)
  #         }
  #         
  #       }
  #     )
  #     
  #     current_subset = current_subset %>% relocate(`ECOSAR_low [mg/l]`, .after = InChIKey)
  #     
  #     colnames(current_subset)[colnames(current_subset) == 'ECOSAR_low [mg/l]'] = paste0('ECOSAR_low_', current_species, '_', current_duration, ' [mg/l]')
  #     
  #   } else if('mean' %in% ecosar_method){
  #     
  #     print('mean not currently supported for ECOSAR, rerun with "low" ')
  #     
  #   }
  #   
  #   
  #   if(data_filter == 'calculated'){
  #     
  #     current_subset = current_subset[,!grepl(x = colnames(current_subset), pattern = 'ECOSAR_[^nml]')]
  #     
  #   }
  #   
  #   if(ecosar_dual_run){
  #     
  #     print('processing ecosar second setting')
  #     
  #     print('Short format needs some work, you doofus!')
  #     
  #     
  #   }
  #   
  #   ################################################################################
  #   # 3.3 Fixing some TEST naming issues
  #   ################################################################################
  #   
  #   colnames(current_subset)[grepl(colnames(current_subset), pattern = 'TEST')] = str_replace(colnames(current_subset)[grepl(colnames(current_subset), pattern = 'TEST')], pattern = 'Fish', replacement = 'fish')
  #   
  #   
  #   ################################################################################
  #   # 3.4 Return data.frame
  #   ################################################################################
  #   
  #   print('Finalizing output')
  #   
  #   if(out_format == 'wide'){
  #     
  #     return(current_subset)
  #     
  #   } else if(out_format == 'long') {
  #     
  #     current_value_columns = 
  #       colnames(current_subset)[!colnames(current_subset) %in% c("original_CAS", 
  #                                                                 "original_SMILES",
  #                                                                 "InChIKey",
  #                                                                 "ecosar_CAS",
  #                                                                 "test_CAS", 
  #                                                                 "ecosar_SMILES",
  #                                                                 "vega_SMILES",
  #                                                                 "test_SMILES" ) &
  #                                  !grepl(colnames(current_subset), pattern = 'number_of')]
  #     
  #     
  #     
  #     current_subset_long = gather(current_subset,
  #                                  key = model_name,
  #                                  value = value,
  #                                  all_of(current_value_columns),
  #                                  factor_key = T,
  #                                  na.rm = F)
  #     
  #     # Add information about vega predictions
  #     
  #     current_vega_reliability_columns = 
  #       colnames(current_vega_reliabilities)[!colnames(current_vega_reliabilities) %in% c("original_CAS", 
  #                                                                                         "original_SMILES",
  #                                                                                         "InChIKey",
  #                                                                                         "ecosar_CAS",
  #                                                                                         "test_CAS", 
  #                                                                                         "ecosar_SMILES",
  #                                                                                         "vega_SMILES",
  #                                                                                         "test_SMILES" ) &
  #                                            grepl(colnames(current_vega_reliabilities), pattern = 'reliability')]
  #     
  #     
  #     current_vega_reliabilities_long = gather(current_vega_reliabilities,
  #                                            key = model_name,
  #                                            value = reliability,
  #                                            all_of(current_vega_reliability_columns),
  #                                            factor_key = T,
  #                                            na.rm = F)
  #     
  #     # Prepare merging of reliability values by renaming model_name in both sets
  #     current_vega_reliabilities_long$model_name = str_replace(current_vega_reliabilities_long$model_name, 
  #                                                              pattern = '(?<=assessment).+$',
  #                                                              replacement = '')
  #     current_subset_long$model_name = str_replace(current_subset_long$model_name, 
  #                                                  pattern = '(?<=assessment).+$',
  #                                                  replacement = '')
  #     
  #     
  #     current_subset_long = merge(current_subset_long, current_vega_reliabilities_long, all.x = T, by = c("original_CAS", 
  #                                                                                          "original_SMILES",
  #                                                                                          "InChIKey",
  #                                                                                          "ecosar_CAS",
  #                                                                                          "test_CAS", 
  #                                                                                          "ecosar_SMILES",
  #                                                                                          "vega_SMILES",
  #                                                                                          "test_SMILES",
  #                                                                                          'model_name'))
  #     
  #     # Remove the number_of columns
  #     current_subset_long = current_subset_long[,!grepl(colnames(current_subset_long), pattern = 'number_of')]
  #     
  #     # Remove all rows with NA in value
  #     current_subset_long = current_subset_long[!is.na(current_subset_long$value),]
  #     
  #     # Add information columns on endpoints, durations, species etc
  #     
  #     
  #     
  #     
  #   }
  #   
  #   
  # }
  
  ################################################################################
  # 4. Long in_format processing
  ################################################################################
  
  if(in_format == 'long'){
   
    print('Long format processing')
    
    # Make the value column numerical
    current_subset$value = as.numeric(current_subset$value)
    
    
    
    ################################################################################
    # 4.1. VEGA processing
    ################################################################################
    
    if(process_vega){
        
      print('Parsing vega first filter (all methods)')
      
      subset_filter = '(FAILSAFE)'
      
      # Invert the filters (translate from "what to keep" to "what to remove")
      if('all' %in% vega_filter){
        inverted_filter = c()
      } else {
        filter_inverter_vector = c('experimental', 'low', 'moderate', 'good')
        inverted_filter = filter_inverter_vector[!filter_inverter_vector %in% tolower(vega_filter)]
      }
      
      # Produce a filter string regexp
      if('experimental' %in% inverted_filter){
        subset_filter = paste(subset_filter, '(EXPERIMENTAL)', sep = '|') 
      }
      if('low' %in% inverted_filter){
        subset_filter = paste(subset_filter, '(low)', sep = '|') 
      }
      if('moderate' %in% inverted_filter){
        subset_filter = paste(subset_filter, '(moderate)', sep = '|') 
      }
      if('good' %in% inverted_filter){
        subset_filter = paste(subset_filter, '(good)', sep = '|')
      }
      
      # Loop over CAS numbers represented in vega to get mean and number of predictions (make sure there are VEGA predictions left)
      for(i in 1:ifelse(length(unique(current_subset[current_subset$QSAR_tool == 'VEGA' & 
                                              !is.na(current_subset$value), 'original_CAS'])) > 0,
                        length(unique(current_subset[current_subset$QSAR_tool == 'VEGA' & 
                                                     !is.na(current_subset$value), 'original_CAS'])),
                        1)){
        
        if(length(unique(current_subset[current_subset$QSAR_tool == 'VEGA' & 
                                        !is.na(current_subset$value), 'original_CAS'])) == 0){
          
          print('All vega predictions filtered out with current settings')
          break
        }
        
        # Get current CAS
        current_cas = unique(current_subset[current_subset$QSAR_tool == 'VEGA' &
                                              !is.na(current_subset$value), 'original_CAS'])[i]
        
        
        # Get current number of predictions for this CAS
        current_vega_number = nrow(current_subset[current_subset$QSAR_tool == 'VEGA' &
                                                    !is.na(current_subset$value) &
                                                    current_subset$original_CAS == current_cas &
                                                    !grepl(current_subset$reliability, pattern = subset_filter),])
        
        # Skip if number is 0
        if(current_vega_number == 0){
          next
        }
        
        ## Calculate mean/low etc based on calculation method
        
        # First we get "current_values"
        current_values = current_subset[current_subset$QSAR_tool == 'VEGA' &
                                          !is.na(current_subset$value) &
                                          current_subset$original_CAS == current_cas &
                                          !grepl(current_subset$reliability, pattern = subset_filter) &
                                          !is.na(current_subset$reliability) &
                                          current_subset$reliability != 'calculated', 'value']
        
        current_values = as.numeric(current_values)
        
        #If geo_mean chosen as method
        if('geo_mean' %in% vega_method){
          
          geo_mean = ifelse(
            length(current_values) == 1,
            current_values,
            ifelse(
              all(is.na(current_values) | current_values == 0),
              NA,
              exp(
                sum(
                  log(
                    current_values[current_values > 0 & !is.na(current_values)]
                    ), 
                  na.rm = T
                  ) / length(current_values[current_values > 0 & !is.na(current_values)])
                )
              )
            )
          
          # Set model name based on method and filters (to be compatible with running calculations twice)
          model_name = paste0('VEGA_',
                              species,
                              '_',
                              endpoint,
                              '_geo_mean_',
                              ifelse(grepl(subset_filter, pattern = 'low'),
                                     'no_low_',
                                     ''),
                              ifelse(grepl(subset_filter, pattern = 'moderate'),
                                     'no_moderate_',
                                     ''),
                              ifelse(grepl(subset_filter, pattern = 'good'),
                                     'no_good_',
                                     ''),
                              ifelse(grepl(subset_filter, pattern = 'EXPERIMENTAL'),
                                     'no_exp_',
                                     ''),
                              'prediction[mg/l]')
          
          
          # Make sure no spaces from input arguments follows
          model_name = str_replace_all(model_name, pattern = ' ', replacement = '_')
          
          
          new_row = c("original_CAS" = current_cas, 
                      "ecosar_CAS" = NA, 
                      "test_CAS" = NA, 
                      "original_SMILES" = current_subset[current_subset$original_CAS == current_cas, "original_SMILES"][1],
                      "ecosar_SMILES" = NA,
                      "vega_SMILES" = current_subset[current_subset$original_CAS == current_cas, "vega_SMILES"][1],
                      "test_SMILES" = NA,
                      "InChIKey" = current_subset[current_subset$original_CAS == current_cas, "InChIKey"][1],
                      "model" = model_name,
                      "value" = geo_mean,
                      "reliability" = 'calculated',
                      "duration" = current_duration,
                      "model_organism" = current_species,
                      "model_endpoint" = current_endpoint,
                      "QSAR_tool" = 'VEGA',
                      "Chemical" = NA,
                      "ECOSAR Class" = NA,
                      "Max Log Kow" = NA,
                      "error_warning" = NA,
                      "Molecular Weight (g/mol)" = NA,
                      "Log Kow" = NA,
                      "Water Solubility (mg/L)" = NA,
                      "Melting Point (C)" = NA,
                      "ecosar_CAS_fixed" = NA,
                      "original_CAS_short" = NA,
                      "model_duration (h)" = NA)
          
          
          if(new_row['value'] != 0 & !is.na(new_row['value'])){
            current_subset = rbind(current_subset, new_row)
          }
          
          
        } 
        
        if('mean' %in% vega_method){
          
          vega_mean = ifelse(
            length(current_values) == 1,
            current_values,
            ifelse(all(is.na(current_values) | is.infinite(current_values) | current_values == 0),
                   NA,
                   mean(current_values[current_values > 0 & !is.na(current_values)]
                        )
                   )
            )
          model_name = paste0('VEGA_',
                              species,
                              '_',
                              endpoint,
                              '_mean_',
                              ifelse(grepl(subset_filter, pattern = 'low'),
                                     'no_low_',
                                     ''),
                              ifelse(grepl(subset_filter, pattern = 'moderate'),
                                     'no_moderate_',
                                     ''),
                              ifelse(grepl(subset_filter, pattern = 'good'),
                                     'no_good_',
                                     ''),
                              ifelse(grepl(subset_filter, pattern = 'EXPERIMENTAL'),
                                     'no_exp_',
                                     ''),
                              'prediction[mg/l]')
          
          # Make sure no spaces from input arguments follows
          model_name = str_replace_all(model_name, pattern = ' ', replacement = '_')
          
          new_row = c("original_CAS" = current_cas, 
                      "ecosar_CAS" = NA, 
                      "test_CAS" = NA, 
                      "original_SMILES" = current_subset[current_subset$original_CAS == current_cas, "original_SMILES"][1],
                      "ecosar_SMILES" = NA,
                      "vega_SMILES" = current_subset[current_subset$original_CAS == current_cas, "vega_SMILES"][1],
                      "test_SMILES" = NA,
                      "InChIKey" = current_subset[current_subset$original_CAS == current_cas, "InChIKey"][1],
                      "model" = model_name,
                      "value" = vega_mean,
                      "reliability" = 'calculated',
                      "duration" = current_duration,
                      "model_organism" = current_species,
                      "model_endpoint" = current_endpoint,
                      "QSAR_tool" = 'VEGA',
                      "Chemical" = NA,
                      "ECOSAR Class" = NA,
                      "Max Log Kow" = NA,
                      "error_warning" = NA,
                      "Molecular Weight (g/mol)" = NA,
                      "Log Kow" = NA,
                      "Water Solubility (mg/L)" = NA,
                      "Melting Point (C)" = NA,
                      "ecosar_CAS_fixed" = NA,
                      "original_CAS_short" = NA,
                      "model_duration (h)" = NA)
          
          
          if(new_row['value'] != 0 & !is.na(new_row['value'])){
            current_subset = rbind(current_subset, new_row)
          }
          
          
        } 
          
          
        if('low' %in% vega_method){
          
          vega_low = ifelse(
            length(current_values) == 1,
            current_values,
            ifelse(all(is.na(current_values) | is.infinite(current_values) | current_values == 0),
                   NA,
                   suppressWarnings(min(current_values[current_values > 0 & !is.na(current_values)]
                                        )
                                    )
                   )
            )
          model_name = paste0('VEGA_',
                              species,
                              '_',
                              endpoint,
                              '_lowest_',
                              ifelse(grepl(subset_filter, pattern = 'low'),
                                     'no_low_',
                                     ''),
                              ifelse(grepl(subset_filter, pattern = 'moderate'),
                                     'no_moderate_',
                                     ''),
                              ifelse(grepl(subset_filter, pattern = 'good'),
                                     'no_good_',
                                     ''),
                              ifelse(grepl(subset_filter, pattern = 'EXPERIMENTAL'),
                                     'no_exp_',
                                     ''),
                              'prediction[mg/l]')
          
          # Make sure no spaces from input arguments follows
          model_name = str_replace_all(model_name, pattern = ' ', replacement = '_')
          
          new_row = c("original_CAS" = current_cas, 
                      "ecosar_CAS" = NA, 
                      "test_CAS" = NA, 
                      "original_SMILES" = current_subset[current_subset$original_CAS == current_cas, "original_SMILES"][1],
                      "ecosar_SMILES" = NA,
                      "vega_SMILES" = current_subset[current_subset$original_CAS == current_cas, "vega_SMILES"][1],
                      "test_SMILES" = NA,
                      "InChIKey" = current_subset[current_subset$original_CAS == current_cas, "InChIKey"][1],
                      "model" = model_name,
                      "value" = vega_low,
                      "reliability" = 'calculated',
                      "duration" = current_duration,
                      "model_organism" = current_species,
                      "model_endpoint" = current_endpoint,
                      "QSAR_tool" = 'VEGA',
                      "Chemical" = NA,
                      "ECOSAR Class" = NA,
                      "Max Log Kow" = NA,
                      "error_warning" = NA,
                      "Molecular Weight (g/mol)" = NA,
                      "Log Kow" = NA,
                      "Water Solubility (mg/L)" = NA,
                      "Melting Point (C)" = NA,
                      "ecosar_CAS_fixed" = NA,
                      "original_CAS_short" = NA,
                      "model_duration (h)" = NA)
          
          
          if(new_row['value'] != 0 & !is.na(new_row['value'])){
            current_subset = rbind(current_subset, new_row)
          }
          
        }
        
      }
      
      # Remove vega entries with NA values
      current_subset = current_subset[(current_subset$QSAR_tool == 'VEGA' &
                                        !is.na(current_subset$value)) |
                                        current_subset$QSAR_tool != 'VEGA',]
      
      
      if(vega_dual_run){
        
        print('Parsing vega second method')
        
        subset_filter2 = '(FAILSAFE)'
        
        # Invert the filters (translate from "what to keep" to "what to remove")
        if('all' %in% vega_filter2){
          inverted_filter = c()
        } else {
          filter_inverter_vector = c('experimental', 'low', 'moderate', 'good')
          inverted_filter = filter_inverter_vector[!filter_inverter_vector %in% tolower(vega_filter2)]
        }
        
        # Produce a filter string regexp
        if('experimental' %in% inverted_filter){
          subset_filter2 = paste(subset_filter2, '(EXPERIMENTAL)', sep = '|') 
        }
        if('low' %in% inverted_filter){
          subset_filter2 = paste(subset_filter2, '(low)', sep = '|') 
        }
        if('moderate' %in% inverted_filter){
          subset_filter2 = paste(subset_filter2, '(moderate)', sep = '|') 
        }
        if('good' %in% inverted_filter){
          subset_filter2 = paste(subset_filter2, '(good)', sep = '|')
        }
        
        
        # Loop over CAS numbers represented in vega to get mean and number of predictions (make sure there are VEGA predictions left)
        for(i in 1:ifelse(length(unique(current_subset[current_subset$QSAR_tool == 'VEGA' & 
                                                       !is.na(current_subset$value), 'original_CAS'])) > 0,
                          length(unique(current_subset[current_subset$QSAR_tool == 'VEGA' & 
                                                       !is.na(current_subset$value), 'original_CAS'])),
                          1)){
          
          if(length(unique(current_subset[current_subset$QSAR_tool == 'VEGA' & 
                                          !is.na(current_subset$value), 'original_CAS'])) == 0){
            
            print('All vega predictions filtered out with current settings')
            break
          }
          
          # Get current CAS
          current_cas = unique(current_subset[current_subset$QSAR_tool == 'VEGA' &
                                                !is.na(current_subset$value), 'original_CAS'])[i]
          
          
          # Get current number of predictions for this CAS
          current_vega_number = nrow(current_subset[current_subset$QSAR_tool == 'VEGA' &
                                                      !is.na(current_subset$value) &
                                                      current_subset$original_CAS == current_cas &
                                                      !grepl(current_subset$reliability, pattern = subset_filter2),])
          
          # Skip if number is 0
          if(current_vega_number == 0){
            next
          }
          
          ## Calculate mean/low etc based on calculation method
          
          # First we get "current_values"
          current_values = current_subset[current_subset$QSAR_tool == 'VEGA' &
                                            !is.na(current_subset$value) &
                                            current_subset$original_CAS == current_cas &
                                            !grepl(current_subset$reliability, pattern = subset_filter2) &
                                            !is.na(current_subset$reliability) &
                                            current_subset$reliability != 'calculated', 'value']
          
          current_values = as.numeric(current_values)
          
          #If geo_mean chosen as method
          if('geo_mean' %in% vega_method2){
            
            geo_mean = ifelse(
              length(current_values) == 1,
              current_values,
              ifelse(
                all(is.na(current_values) | current_values == 0),
                NA,
                exp(
                  sum(
                    log(
                      current_values[current_values > 0 & !is.na(current_values)]
                    ), 
                    na.rm = T
                  ) / length(current_values[current_values > 0 & !is.na(current_values)])
                )
              )
            )
            model_name = paste0('VEGA_',
                                species,
                                '_',
                                endpoint,
                                '_geo_mean_',
                                ifelse(grepl(subset_filter2, pattern = 'low'),
                                       'no_low_',
                                       ''),
                                ifelse(grepl(subset_filter2, pattern = 'moderate'),
                                       'no_moderate_',
                                       ''),
                                ifelse(grepl(subset_filter2, pattern = 'good'),
                                       'no_good_',
                                       ''),
                                ifelse(grepl(subset_filter2, pattern = 'EXPERIMENTAL'),
                                       'no_exp_',
                                       ''),
                                'prediction[mg/l]')
            
            # Make sure no spaces from input arguments follows
            model_name = str_replace_all(model_name, pattern = ' ', replacement = '_')
            
            new_row = c("original_CAS" = current_cas, 
                        "ecosar_CAS" = NA, 
                        "test_CAS" = NA, 
                        "original_SMILES" = current_subset[current_subset$original_CAS == current_cas, "original_SMILES"][1],
                        "ecosar_SMILES" = NA,
                        "vega_SMILES" = current_subset[current_subset$original_CAS == current_cas, "vega_SMILES"][1],
                        "test_SMILES" = NA,
                        "InChIKey" = current_subset[current_subset$original_CAS == current_cas, "InChIKey"][1],
                        "model" = model_name,
                        "value" = geo_mean,
                        "reliability" = 'calculated',
                        "duration" = current_duration,
                        "model_organism" = current_species,
                        "model_endpoint" = current_endpoint,
                        "QSAR_tool" = 'VEGA',
                        "Chemical" = NA,
                        "ECOSAR Class" = NA,
                        "Max Log Kow" = NA,
                        "error_warning" = NA,
                        "Molecular Weight (g/mol)" = NA,
                        "Log Kow" = NA,
                        "Water Solubility (mg/L)" = NA,
                        "Melting Point (C)" = NA,
                        "ecosar_CAS_fixed" = NA,
                        "original_CAS_short" = NA,
                        "model_duration (h)" = NA)
            
            
            # If new_row has a value, save it. Otherwise we move on
            if(new_row['value'] != 0 & !is.na(new_row['value'])){
              current_subset = rbind(current_subset, new_row)
            }
            
            
          } 
          
          if('mean' %in% vega_method2){
            
            vega_mean = ifelse(
              length(current_values) == 1,
              current_values,
              ifelse(all(is.na(current_values) | is.infinite(current_values) | current_values == 0),
                     NA,
                     mean(current_values[current_values > 0 & !is.na(current_values)]
                     )
              )
            )
            
            model_name = paste0('VEGA_',
                                species,
                                '_',
                                endpoint,
                                '_mean_',
                                ifelse(grepl(subset_filter2, pattern = 'low'),
                                       'no_low_',
                                       ''),
                                ifelse(grepl(subset_filter2, pattern = 'moderate'),
                                       'no_moderate_',
                                       ''),
                                ifelse(grepl(subset_filter2, pattern = 'good'),
                                       'no_good_',
                                       ''),
                                ifelse(grepl(subset_filter2, pattern = 'EXPERIMENTAL'),
                                       'no_exp_',
                                       ''),
                                'prediction[mg/l]')
            
            # Make sure no spaces from input arguments follows
            model_name = str_replace_all(model_name, pattern = ' ', replacement = '_')
            
            new_row = c("original_CAS" = current_cas, 
                        "ecosar_CAS" = NA, 
                        "test_CAS" = NA, 
                        "original_SMILES" = current_subset[current_subset$original_CAS == current_cas, "original_SMILES"][1],
                        "ecosar_SMILES" = NA,
                        "vega_SMILES" = current_subset[current_subset$original_CAS == current_cas, "vega_SMILES"][1],
                        "test_SMILES" = NA,
                        "InChIKey" = current_subset[current_subset$original_CAS == current_cas, "InChIKey"][1],
                        "model" = model_name,
                        "value" = vega_mean,
                        "reliability" = 'calculated',
                        "duration" = current_duration,
                        "model_organism" = current_species,
                        "model_endpoint" = current_endpoint,
                        "QSAR_tool" = 'VEGA',
                        "Chemical" = NA,
                        "ECOSAR Class" = NA,
                        "Max Log Kow" = NA,
                        "error_warning" = NA,
                        "Molecular Weight (g/mol)" = NA,
                        "Log Kow" = NA,
                        "Water Solubility (mg/L)" = NA,
                        "Melting Point (C)" = NA,
                        "ecosar_CAS_fixed" = NA,
                        "original_CAS_short" = NA,
                        "model_duration (h)" = NA)
            
            
            if(new_row['value'] != 0 & !is.na(new_row['value'])){
              current_subset = rbind(current_subset, new_row)
            }
            
            
          } 
          
          if('low' %in% vega_method2){
            
            vega_low = ifelse(
              length(current_values) == 1,
              current_values,
              ifelse(all(is.na(current_values) | is.infinite(current_values) | current_values == 0),
                     NA,
                     suppressWarnings(min(current_values[current_values > 0 & !is.na(current_values)]
                     )
                     )
              )
            )
            model_name = paste0('VEGA_',
                                species,
                                '_',
                                endpoint,
                                '_lowest_',
                                ifelse(grepl(subset_filter2, pattern = 'low'),
                                       'no_low_',
                                       ''),
                                ifelse(grepl(subset_filter2, pattern = 'moderate'),
                                       'no_moderate_',
                                       ''),
                                ifelse(grepl(subset_filter2, pattern = 'good'),
                                       'no_good_',
                                       ''),
                                ifelse(grepl(subset_filter2, pattern = 'EXPERIMENTAL'),
                                       'no_exp_',
                                       ''),
                                'prediction[mg/l]')
            
            # Make sure no spaces from input arguments follows
            model_name = str_replace_all(model_name, pattern = ' ', replacement = '_')
            
            new_row = c("original_CAS" = current_cas, 
                        "ecosar_CAS" = NA, 
                        "test_CAS" = NA, 
                        "original_SMILES" = current_subset[current_subset$original_CAS == current_cas, "original_SMILES"][1],
                        "ecosar_SMILES" = NA,
                        "vega_SMILES" = current_subset[current_subset$original_CAS == current_cas, "vega_SMILES"][1],
                        "test_SMILES" = NA,
                        "InChIKey" = current_subset[current_subset$original_CAS == current_cas, "InChIKey"][1],
                        "model" = model_name,
                        "value" = vega_low,
                        "reliability" = 'calculated',
                        "duration" = current_duration,
                        "model_organism" = current_species,
                        "model_endpoint" = current_endpoint,
                        "QSAR_tool" = 'VEGA',
                        "Chemical" = NA,
                        "ECOSAR Class" = NA,
                        "Max Log Kow" = NA,
                        "error_warning" = NA,
                        "Molecular Weight (g/mol)" = NA,
                        "Log Kow" = NA,
                        "Water Solubility (mg/L)" = NA,
                        "Melting Point (C)" = NA,
                        "ecosar_CAS_fixed" = NA,
                        "original_CAS_short" = NA,
                        "model_duration (h)" = NA)
            
            
            if(new_row['value'] != 0 & !is.na(new_row['value'])){
              current_subset = rbind(current_subset, new_row)
            }
            
          }
          
        }
        
        # Remove vega entries with NA values
        current_subset = current_subset[(current_subset$QSAR_tool == 'VEGA' &
                                           !is.na(current_subset$value)) |
                                          current_subset$QSAR_tool != 'VEGA',]
        
        
        
      }
      
      # If we dont keep all data, remove predictions using the quality filter
      if(data_filter == 'calculated'){
        
        current_subset = current_subset[(current_subset$QSAR_tool == 'VEGA' &
                                           current_subset$reliability == 'calculated' &
                                           !is.na(current_subset$reliability)) |
                                          current_subset$QSAR_tool != 'VEGA',]
        
      }
    }
    
    ################################################################################
    #   4.2. ECOSAR processing    
    ################################################################################
    
    if(process_ecosar){
      
      print('Parsing ecosar first method')
      
      
      # Loop over CAS numbers present in ECOSAR entries
      for(i in 1:length(unique(current_subset[current_subset$QSAR_tool == 'ECOSAR' &
                                              !is.na(current_subset$value), 'original_CAS']))){
      
        
        current_cas = unique(current_subset[current_subset$QSAR_tool == 'ECOSAR' &
                                              !is.na(current_subset$value), 'original_CAS'])[i]
        
        
        # ECOSAR states in their user manual that "Traditionally, in the absence of adequate measured data, the most conservative effect level is
        # used when predictions are identified from multiple classes", so pick out the most conservative (lowest) value if there are multiple prediction classes
        
        # First we get "current_values"
        current_values = current_subset[current_subset$QSAR_tool == 'ECOSAR' &
                                          !is.na(current_subset$value) &
                                          current_subset$original_CAS == current_cas, 'value']
        
        current_values = as.numeric(current_values)
        
        
        # The drop baseline script (Dropping neutral organics if there are other values)
        baseline = ''
        
        if('drop baseline' %in% ecosar_method){
          
          if(sum(grepl(current_subset[current_subset$QSAR_tool == 'ECOSAR' &
                                      !is.na(current_subset$value) &
                                      current_subset$original_CAS == current_cas, 'model'], pattern = 'Neutral_Organics')) < length(current_values)){
            
            # Set model name to notify non-baseline predictions
            baseline = 'no_baseline_'
            
            current_values = current_subset[current_subset$QSAR_tool == 'ECOSAR' &
                                              !is.na(current_subset$value) &
                                              current_subset$original_CAS == current_cas &
                                              !grepl(current_subset$model, pattern = 'Neutral_Organics'), 'value']
            
            current_values = as.numeric(current_values)
            
          } else {
            
            # Set model name to notify only baseline predictions
            baseline = 'baseline_only_'
            
          }
        }
        
        
        if('low' %in% ecosar_method){
          
          current_model_name = paste0('ECOSAR_',
                                      species,
                                      '_',
                                      endpoint,
                                      '_lowest_', 
                                      baseline, 
                                      'prediction [mg/l]')
          
          # Make sure no spaces from input arguments follows
          current_model_name = str_replace_all(current_model_name, pattern = ' ', replacement = '_')
          
          ecosar_low = ifelse(
            length(current_values) == 1,
            current_values,
            ifelse(all(is.na(current_values) | is.infinite(current_values) | current_values == 0),
                   NA,
                   suppressWarnings(min(current_values[current_values > 0 & !is.na(current_values)]
                   )
                   )
            )
          )
          
          
          new_row = c("original_CAS" = current_cas, 
                      "ecosar_CAS" = current_subset[current_subset$original_CAS == current_cas, "ecosar_CAS"][1],
                      "test_CAS" = NA, 
                      "original_SMILES" = current_subset[current_subset$original_CAS == current_cas, "original_SMILES"][1],
                      "ecosar_SMILES" = NA,
                      "vega_SMILES" = NA,
                      "test_SMILES" = NA,
                      "InChIKey" = current_subset[current_subset$original_CAS == current_cas, "InChIKey"][1],
                      "model" = current_model_name,
                      "value" = ecosar_low,
                      "reliability" = 'calculated',
                      "duration" = current_duration,
                      "model_organism" = current_species,
                      "model_endpoint" = current_endpoint,
                      "QSAR_tool" = 'ECOSAR',
                      "Chemical" = current_subset[current_subset$original_CAS == current_cas, "Chemical"][1],
                      "ECOSAR Class" = NA,
                      "Max Log Kow" = NA,
                      "error_warning" = NA,
                      "Molecular Weight (g/mol)" = current_subset[current_subset$original_CAS == current_cas, "Molecular Weight (g/mol)"][1],
                      "Log Kow" = NA,
                      "Water Solubility (mg/L)" = current_subset[current_subset$original_CAS == current_cas, "Water Solubility (mg/L)"][1],
                      "Melting Point (C)" = current_subset[current_subset$original_CAS == current_cas, "Melting Point (C)"][1],
                      "ecosar_CAS_fixed" = NA,
                      "original_CAS_short" = NA,
                      "model_duration (h)" = NA)
          
          
          if(new_row['value'] != 0 & !is.na(new_row['value'])){
            current_subset = rbind(current_subset, new_row)
          }
          
        } else if('mean' %in% ecosar_method){
          
          print('mean not currently supported for ECOSAR, rerun with "low" or "geo_mean" ')
          break
          
        } else if('geo_mean' %in% ecosar_method ){
          
          current_model_name = paste0('ECOSAR_',
                                      species,
                                      '_',
                                      endpoint,
                                      '_geo_mean_', 
                                      baseline, 
                                      'prediction [mg/l]')
          
          
          # Make sure no spaces from input arguments follows
          current_model_name = str_replace_all(current_model_name, pattern = ' ', replacement = '_')
          
          
          geo_mean = ifelse(
            length(current_values) == 1,
            current_values,
            ifelse(
              all(is.na(current_values) | current_values == 0),
              NA,
              exp(
                sum(
                  log(
                    current_values[current_values > 0 & !is.na(current_values)]
                  ), 
                  na.rm = T
                ) / length(current_values[current_values > 0 & !is.na(current_values)])
              )
            )
          )
          
          
          new_row = c("original_CAS" = current_cas, 
                      "ecosar_CAS" = current_subset[current_subset$original_CAS == current_cas, "ecosar_CAS"][1],
                      "test_CAS" = NA, 
                      "original_SMILES" = current_subset[current_subset$original_CAS == current_cas, "original_SMILES"][1],
                      "ecosar_SMILES" = NA,
                      "vega_SMILES" = NA,
                      "test_SMILES" = NA,
                      "InChIKey" = current_subset[current_subset$original_CAS == current_cas, "InChIKey"][1],
                      "model" = current_model_name,
                      "value" = geo_mean,
                      "reliability" = 'calculated',
                      "duration" = current_duration,
                      "model_organism" = current_species,
                      "model_endpoint" = current_endpoint,
                      "QSAR_tool" = 'ECOSAR',
                      "Chemical" = current_subset[current_subset$original_CAS == current_cas, "Chemical"][1],
                      "ECOSAR Class" = NA,
                      "Max Log Kow" = NA,
                      "error_warning" = NA,
                      "Molecular Weight (g/mol)" = current_subset[current_subset$original_CAS == current_cas, "Molecular Weight (g/mol)"][1],
                      "Log Kow" = NA,
                      "Water Solubility (mg/L)" = current_subset[current_subset$original_CAS == current_cas, "Water Solubility (mg/L)"][1],
                      "Melting Point (C)" = current_subset[current_subset$original_CAS == current_cas, "Melting Point (C)"][1],
                      "ecosar_CAS_fixed" = NA,
                      "original_CAS_short" = NA,
                      "model_duration (h)" = NA)
          
          if(new_row['value'] != 0 & !is.na(new_row['value'])){
            current_subset = rbind(current_subset, new_row)
          }
          
          
        }
        
          
      }
      
      
      if(ecosar_dual_run & ecosar_method2[1] != ecosar_method[1]){
        
        print('Parsing ecosar second method')
        
        # Loop over CAS numbers present in ECOSAR entries
        for(i in 1:length(unique(current_subset[current_subset$QSAR_tool == 'ECOSAR' &
                                                !is.na(current_subset$value), 'original_CAS']))){
          
          
          current_cas = unique(current_subset[current_subset$QSAR_tool == 'ECOSAR' &
                                                !is.na(current_subset$value), 'original_CAS'])[i]
          
          
          # ECOSAR states in their user manual that "Traditionally, in the absence of adequate measured data, the most conservative effect level is
          # used when predictions are identified from multiple classes", so pick out the most conservative (lowest) value if there are multiple prediction classes
          
          # First we get "current_values"
          current_values = current_subset[current_subset$QSAR_tool == 'ECOSAR' &
                                            is.na(current_subset$reliability) &
                                            !is.na(current_subset$value) &
                                            current_subset$original_CAS == current_cas, 'value']
          
          current_values = as.numeric(current_values)
          
          
          # The drop baseline script (Dropping neutral organics if there are other values)
          baseline = ''
          
          if('drop baseline' %in% ecosar_method2){
            
            if(sum(grepl(current_subset[current_subset$QSAR_tool == 'ECOSAR' &
                                        is.na(current_subset$reliability) &
                                        !is.na(current_subset$value) &
                                        current_subset$original_CAS == current_cas, 'model'], pattern = 'Neutral_Organics')) < length(current_values)){
              
              # # Set model name to notify non-baseline predictions
              # current_model_name = paste0('ECOSAR_',
              #                             species,
              #                             '_',
              #                             endpoint,
              #                             '_mean_no_baseline_prediction [mg/l]')
              # 
              #
              # # Make sure no spaces from input arguments follows
              # current_model_name = str_replace_all(current_model_name, pattern = ' ', replacement = '_')
              
              baseline = 'no_baseline_'
              
              current_values = current_subset[current_subset$QSAR_tool == 'ECOSAR' &
                                                is.na(current_subset$reliability) &
                                                !is.na(current_subset$value) &
                                                current_subset$original_CAS == current_cas &
                                                !grepl(current_subset$model, pattern = 'Neutral_Organics'), 'value']
              
              current_values = as.numeric(current_values)
              
            } else {
              
              # Set model name to notify only baseline predictions
              baseline = 'baseline_only_'
              
              
            }
          }
          
          
          if('low' %in% ecosar_method2){
            
            current_model_name = paste0('ECOSAR_',
                                        species,
                                        '_',
                                        endpoint,
                                        '_lowest_', baseline, 'prediction [mg/l]')
            
            
            # Make sure no spaces from input arguments follows
            current_model_name = str_replace_all(current_model_name, pattern = ' ', replacement = '_')
            
            
            ecosar_low = ifelse(
              length(current_values) == 1,
              current_values,
              ifelse(all(is.na(current_values) | is.infinite(current_values) | current_values == 0),
                     NA,
                     suppressWarnings(min(current_values[current_values > 0 & !is.na(current_values)]
                     )
                     )
              )
            )
            
            
            new_row = c("original_CAS" = current_cas, 
                        "ecosar_CAS" = current_subset[current_subset$original_CAS == current_cas, "ecosar_CAS"][1],
                        "test_CAS" = NA, 
                        "original_SMILES" = current_subset[current_subset$original_CAS == current_cas, "original_SMILES"][1],
                        "ecosar_SMILES" = NA,
                        "vega_SMILES" = NA,
                        "test_SMILES" = NA,
                        "InChIKey" = current_subset[current_subset$original_CAS == current_cas, "InChIKey"][1],
                        "model" = current_model_name,
                        "value" = ecosar_low,
                        "reliability" = 'calculated',
                        "duration" = current_duration,
                        "model_organism" = current_species,
                        "model_endpoint" = current_endpoint,
                        "QSAR_tool" = 'ECOSAR',
                        "Chemical" = current_subset[current_subset$original_CAS == current_cas, "Chemical"][1],
                        "ECOSAR Class" = NA,
                        "Max Log Kow" = NA,
                        "error_warning" = NA,
                        "Molecular Weight (g/mol)" = current_subset[current_subset$original_CAS == current_cas, "Molecular Weight (g/mol)"][1],
                        "Log Kow" = NA,
                        "Water Solubility (mg/L)" = current_subset[current_subset$original_CAS == current_cas, "Water Solubility (mg/L)"][1],
                        "Melting Point (C)" = current_subset[current_subset$original_CAS == current_cas, "Melting Point (C)"][1],
                        "ecosar_CAS_fixed" = NA,
                        "original_CAS_short" = NA,
                        "model_duration (h)" = NA)
            
            
            if(new_row['value'] != 0 & !is.na(new_row['value'])){
              current_subset = rbind(current_subset, new_row)
            }
            
          } else if('mean' %in% ecosar_method2){
            
            print('mean not currently supported for ECOSAR, rerun with "low" ')
            break
            
          } else if('geo_mean' %in% ecosar_method2 ){
            
            current_model_name = paste0('ECOSAR_',
                                        species,
                                        '_',
                                        endpoint,
                                        '_geo_mean_', baseline, 'prediction [mg/l]')
            
            
            # Make sure no spaces from input arguments follows
            current_model_name = str_replace_all(current_model_name, pattern = ' ', replacement = '_')
            
            
            geo_mean = ifelse(
              length(current_values) == 1,
              current_values,
              ifelse(
                all(is.na(current_values) | current_values == 0),
                NA,
                exp(
                  sum(
                    log(
                      current_values[current_values > 0 & !is.na(current_values)]
                    ), 
                    na.rm = T
                  ) / length(current_values[current_values > 0 & !is.na(current_values)])
                )
              )
            )
            
            
            new_row = c("original_CAS" = current_cas, 
                        "ecosar_CAS" = current_subset[current_subset$original_CAS == current_cas, "ecosar_CAS"][1],
                        "test_CAS" = NA, 
                        "original_SMILES" = current_subset[current_subset$original_CAS == current_cas, "original_SMILES"][1],
                        "ecosar_SMILES" = NA,
                        "vega_SMILES" = NA,
                        "test_SMILES" = NA,
                        "InChIKey" = current_subset[current_subset$original_CAS == current_cas, "InChIKey"][1],
                        "model" = current_model_name,
                        "value" = geo_mean,
                        "reliability" = 'calculated',
                        "duration" = current_duration,
                        "model_organism" = current_species,
                        "model_endpoint" = current_endpoint,
                        "QSAR_tool" = 'ECOSAR',
                        "Chemical" = current_subset[current_subset$original_CAS == current_cas, "Chemical"][1],
                        "ECOSAR Class" = NA,
                        "Max Log Kow" = NA,
                        "error_warning" = NA,
                        "Molecular Weight (g/mol)" = current_subset[current_subset$original_CAS == current_cas, "Molecular Weight (g/mol)"][1],
                        "Log Kow" = NA,
                        "Water Solubility (mg/L)" = current_subset[current_subset$original_CAS == current_cas, "Water Solubility (mg/L)"][1],
                        "Melting Point (C)" = current_subset[current_subset$original_CAS == current_cas, "Melting Point (C)"][1],
                        "ecosar_CAS_fixed" = NA,
                        "original_CAS_short" = NA,
                        "model_duration (h)" = NA)
            
            if(new_row['value'] != 0 & !is.na(new_row['value'])){
              current_subset = rbind(current_subset, new_row)
            }
            
            
          }
          
          
        }
        
      }
      
      
      
      
      if(data_filter == 'calculated'){
        
        current_subset = current_subset[(current_subset$QSAR_tool == 'ECOSAR' &
                                           current_subset$reliability == 'calculated' &
                                          !is.na(current_subset$reliability)) |
                                          current_subset$QSAR_tool != 'ECOSAR',]
                                          
        
      }
    
    }
    
    
    ################################################################################
    # 4.3. Return long data.frame
    ################################################################################
    
    print('Finalizing output')
    
    end_time = Sys.time()
    
    print(paste0('Ending script at: ', end_time, ' after ', hms::as_hms(difftime(end_time, start_time))))
    
    return(current_subset)
     
  }
  
}
