################################################################################
#    A function for calculating distances (and other calculated metrics)       #
#       in the QSAR merger script                                              #
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

################################################################################
# 1. Initiating function and parsing arguments
################################################################################

function(experimental_averages_df,
         QSAR_df,
         identifiers,
         settings = NULL){
  
  # Check that the input makes sense
  if(is.null(experimental_averages_df)|!is.data.frame(experimental_averages_df)){
    print('Experimental dataframe not provided or wrong format')
    print('Rerun function with correct arguments')
    return(FALSE)
  }
  
  if(is.null(QSAR_df)|!is.data.frame(QSAR_df)){
    print('QSAR dataframe not provided or wrong format')
    print('Rerun function with correct arguments')
    return(FALSE)
  }
  
  if(is.null(settings)){
    print('No additional settings provided, running base script')
  }
  
  ################################################################################
  # 2. function body
  ################################################################################
  
  
  ## Make a dataframe for endpoint calculated metrics
  
  df_calculated = experimental_averages_df[,c('original_CAS', "experimental_mean", "EFSA_mean", "ECOTOX_mean")]
  
  df_calculated = merge(df_calculated, QSAR_df[QSAR_df$reliability == 'calculated' &
                                                 QSAR_df$QSAR_tool == 'VEGA' &
                                                 grepl(QSAR_df$model, pattern = 'geo_mean_no_low', ignore.case = T), 
                                               c('original_CAS', 'value')],
                        all = T, by = 'original_CAS')
  
  df_calculated = df_calculated %>%
    rename(vega_geo_mean_no_low = value)
  
  df_calculated = merge(df_calculated, QSAR_df[QSAR_df$reliability == 'calculated' &
                                                 QSAR_df$QSAR_tool == 'VEGA' &
                                                 grepl(QSAR_df$model, pattern = 'lowest_no_low', ignore.case = T), 
                                               c('original_CAS', 'value')],
                        all = T, by = 'original_CAS')
  
  df_calculated = df_calculated %>%
    rename(vega_lowest_no_low = value)
  
  
  df_calculated = merge(df_calculated, QSAR_df[QSAR_df$reliability == 'calculated' &
                                                 QSAR_df$QSAR_tool == 'VEGA' &
                                                 grepl(QSAR_df$model, pattern = 'geo_mean_no_moderate_no_good', ignore.case = T), 
                                               c('original_CAS', 'value')],
                        all = T, by = 'original_CAS')
  
  df_calculated = df_calculated %>%
    rename(vega_geo_mean_only_low = value)
  
  
  df_calculated = merge(df_calculated, QSAR_df[QSAR_df$reliability == 'calculated' &
                                                 QSAR_df$QSAR_tool == 'VEGA' &
                                                 grepl(QSAR_df$model, pattern = 'lowest_no_moderate_no_good', ignore.case = T), 
                                               c('original_CAS', 'value')],
                        all = T, by = 'original_CAS')
  
  df_calculated = df_calculated %>%
    rename(vega_lowest_only_low = value)
  
  df_calculated = merge(df_calculated, QSAR_df[QSAR_df$reliability == 'calculated' &
                                                 !is.na(QSAR_df$reliability) &
                                                 grepl(QSAR_df$model, pattern = 'mean', ignore.case = T) &
                                                 QSAR_df$QSAR_tool == 'ECOSAR', 
                                               c('original_CAS', 'value', 'model')],
                        all = T, by = 'original_CAS')
  
  
  df_calculated = df_calculated %>%
    rename(ecosar_mean = value)
  
  df_calculated = merge(df_calculated, QSAR_df[QSAR_df$reliability == 'calculated' &
                                                 !is.na(QSAR_df$reliability) &
                                                 grepl(QSAR_df$model, pattern = 'low', ignore.case = T) &
                                                 QSAR_df$QSAR_tool == 'ECOSAR', 
                                               c('original_CAS', 'value')],
                        all = T, by = 'original_CAS')
  
  
  df_calculated = df_calculated %>%
    rename(ecosar_lowest = value)
  
  # Also get which ecosar predictions are baseline only from the model-field
  df_calculated$ecosar_model = ifelse(grepl(df_calculated$model, pattern = 'no_baseline'), 'no_baseline', 'baseline')
  
  # Drop model column (since we are adding another one)
  df_calculated$model = NULL
  
  df_calculated = merge(df_calculated, QSAR_df[QSAR_df$QSAR_tool == 'T.E.S.T.' &
                                                 QSAR_df$reliability != 'EXPERIMENTAL', c('original_CAS', 'value')],
                                   all = T, by = 'original_CAS')
  
  df_calculated = df_calculated %>%
    rename(test_consensus = value)
  
  # Remove NAs
  df_calculated = df_calculated[!is.na(df_calculated$original_CAS),]
  
  # Make numeric
  df_calculated$experimental_mean = as.numeric(df_calculated$experimental_mean)
  df_calculated$vega_geo_mean_no_low = as.numeric(df_calculated$vega_geo_mean_no_low)
  df_calculated$vega_lowest_no_low = as.numeric(df_calculated$vega_lowest_no_low)
  df_calculated$vega_geo_mean_only_low = as.numeric(df_calculated$vega_geo_mean_only_low)
  df_calculated$vega_lowest_only_low = as.numeric(df_calculated$vega_lowest_only_low)
  df_calculated$ecosar_mean = as.numeric(df_calculated$ecosar_mean)
  df_calculated$ecosar_low = as.numeric(df_calculated$ecosar_lowest)
  df_calculated$test_consensus = as.numeric(df_calculated$test_consensus)
  
  
  # Add log columns and log_distance and distance_abs
  df_calculated$experimental_log = log10(df_calculated$experimental_mean)
  df_calculated$vega_geo_mean_no_low_log = log10(df_calculated$vega_geo_mean_no_low)
  df_calculated$vega_lowest_no_low_log = log10(df_calculated$vega_lowest_no_low)
  df_calculated$vega_geo_mean_only_low_log = log10(df_calculated$vega_geo_mean_only_low)
  df_calculated$vega_lowest_only_low_log = log10(df_calculated$vega_lowest_only_low)
  df_calculated$ecosar_mean_log = log10(df_calculated$ecosar_mean)
  df_calculated$ecosar_low_log = log10(df_calculated$ecosar_low)
  df_calculated$test_log = log10(df_calculated$test_consensus)
  
  
  df_calculated$vega_geo_mean_no_low_log_distance = df_calculated$vega_geo_mean_no_low_log - df_calculated$experimental_log
  df_calculated$vega_lowest_no_low_log_distance = df_calculated$vega_lowest_no_low_log - df_calculated$experimental_log
  df_calculated$vega_geo_mean_only_low_log_distance = df_calculated$vega_geo_mean_only_low_log - df_calculated$experimental_log
  df_calculated$vega_lowest_only_low_log_distance = df_calculated$vega_lowest_only_low_log - df_calculated$experimental_log
  df_calculated$ecosar_mean_log_distance = df_calculated$ecosar_mean_log - df_calculated$experimental_log
  df_calculated$ecosar_low_log_distance = df_calculated$ecosar_low_log - df_calculated$experimental_log
  df_calculated$test_log_distance = df_calculated$test_log - df_calculated$experimental_log
  
  
  df_calculated$vega_geo_mean_no_low_log_distance_abs = abs(df_calculated$vega_geo_mean_no_low_log_distance)
  df_calculated$vega_lowest_no_low_log_distance_abs = abs(df_calculated$vega_lowest_no_low_log_distance)
  df_calculated$vega_geo_mean_only_low_log_distance_abs = abs(df_calculated$vega_geo_mean_only_low_log_distance)
  df_calculated$vega_lowest_only_low_log_distance_abs = abs(df_calculated$vega_lowest_only_low_log_distance)
  df_calculated$ecosar_mean_log_distance_abs = abs(df_calculated$ecosar_mean_log_distance)
  df_calculated$ecosar_low_log_distance_abs = abs(df_calculated$ecosar_low_log_distance)
  df_calculated$test_log_distance_abs = abs(df_calculated$test_log_distance)
  
  ## Also do distance to separate indatasets
  # EFSA
  df_calculated$vega_geo_mean_no_low_efsa_log_distance = df_calculated$vega_geo_mean_no_low_log - log10(df_calculated$EFSA_mean)
  df_calculated$vega_lowest_no_low_efsa_log_distance = df_calculated$vega_lowest_no_low_log - log10(df_calculated$EFSA_mean)
  df_calculated$vega_geo_mean_only_low_efsa_log_distance = df_calculated$vega_geo_mean_only_low_log - log10(df_calculated$EFSA_mean)
  df_calculated$vega_lowest_only_low_efsa_log_distance = df_calculated$vega_lowest_only_low_log - log10(df_calculated$EFSA_mean)
  df_calculated$ecosar_mean_efsa_log_distance = df_calculated$ecosar_mean_log - log10(df_calculated$EFSA_mean)
  df_calculated$ecosar_low_efsa_log_distance = df_calculated$ecosar_low_log - log10(df_calculated$EFSA_mean)
  df_calculated$test_efsa_log_distance = df_calculated$test_log - log10(df_calculated$EFSA_mean)
  
  df_calculated$vega_geo_mean_no_low_efsa_log_distance_abs = abs(df_calculated$vega_geo_mean_no_low_log - log10(df_calculated$EFSA_mean))
  df_calculated$vega_lowest_no_low_efsa_log_distance_abs = abs(df_calculated$vega_lowest_no_low_log - log10(df_calculated$EFSA_mean))
  df_calculated$vega_geo_mean_only_low_efsa_log_distance_abs = abs(df_calculated$vega_geo_mean_only_low_log - log10(df_calculated$EFSA_mean))
  df_calculated$vega_lowest_only_low_efsa_log_distance_abs = abs(df_calculated$vega_lowest_only_low_log - log10(df_calculated$EFSA_mean))
  df_calculated$ecosar_mean_efsa_log_distance_abs = abs(df_calculated$ecosar_mean_efsa_log_distance)
  df_calculated$ecosar_low_efsa_log_distance_abs = abs(df_calculated$ecosar_low_efsa_log_distance)
  df_calculated$test_efsa_log_distance_abs = abs(df_calculated$test_efsa_log_distance)
  
  # ECOTOX
  df_calculated$vega_geo_mean_no_low_ecotox_log_distance = df_calculated$vega_geo_mean_no_low_log - log10(df_calculated$ECOTOX_mean)
  df_calculated$vega_lowest_no_low_ecotox_log_distance = df_calculated$vega_lowest_no_low_log - log10(df_calculated$ECOTOX_mean)
  df_calculated$vega_geo_mean_only_low_ecotox_log_distance = df_calculated$vega_geo_mean_only_low_log - log10(df_calculated$ECOTOX_mean)
  df_calculated$vega_lowest_only_low_ecotox_log_distance = df_calculated$vega_lowest_only_low_log - log10(df_calculated$ECOTOX_mean)
  df_calculated$ecosar_mean_ecotox_log_distance = df_calculated$ecosar_mean_log - log10(df_calculated$ECOTOX_mean)
  df_calculated$ecosar_low_ecotox_log_distance = df_calculated$ecosar_low_log - log10(df_calculated$ECOTOX_mean)
  df_calculated$test_ecotox_log_distance = df_calculated$test_log - log10(df_calculated$ECOTOX_mean)
  
  df_calculated$vega_geo_mean_no_low_ecotox_log_distance_abs = abs(df_calculated$vega_geo_mean_no_low_log - log10(df_calculated$ECOTOX_mean))
  df_calculated$vega_lowest_no_low_ecotox_log_distance_abs = abs(df_calculated$vega_lowest_no_low_log - log10(df_calculated$ECOTOX_mean))
  df_calculated$vega_geo_mean_only_low_ecotox_log_distance_abs = abs(df_calculated$vega_geo_mean_only_low_log - log10(df_calculated$ECOTOX_mean))
  df_calculated$vega_lowest_only_low_ecotox_log_distance_abs = abs(df_calculated$vega_lowest_only_low_log - log10(df_calculated$ECOTOX_mean))
  df_calculated$ecosar_mean_ecotox_log_distance_abs = abs(df_calculated$ecosar_mean_ecotox_log_distance)
  df_calculated$ecosar_low_ecotox_log_distance_abs = abs(df_calculated$ecosar_low_ecotox_log_distance)
  df_calculated$test_ecotox_log_distance_abs = abs(df_calculated$test_ecotox_log_distance)
  
  
  
  
  
  
  
  
  
  # Inaccurate above
  
  # Add inaccurate true/false for one order, two orders and three orders above (all data)
  
  df_calculated$vega_geo_mean_no_low_inaccurate_one_order_above = ifelse(df_calculated$vega_geo_mean_no_low_log_distance > 1, 1, 0)
  df_calculated$vega_lowest_no_low_inaccurate_one_order_above = ifelse(df_calculated$vega_lowest_no_low_log_distance > 1, 1, 0)
  df_calculated$vega_geo_mean_only_low_inaccurate_one_order_above = ifelse(df_calculated$vega_geo_mean_only_low_log_distance > 1, 1, 0)
  df_calculated$vega_lowest_only_low_inaccurate_one_order_above = ifelse(df_calculated$vega_lowest_only_low_log_distance > 1, 1, 0)
  df_calculated$ecosar_mean_inaccurate_one_order_above = ifelse(df_calculated$ecosar_mean_log_distance > 1, 1, 0)
  df_calculated$ecosar_low_inaccurate_one_order_above = ifelse(df_calculated$ecosar_low_log_distance > 1, 1, 0)
  df_calculated$test_inaccurate_one_order_above = ifelse(df_calculated$test_log_distance > 1, 1, 0)
  
  df_calculated$vega_geo_mean_no_low_inaccurate_two_order_above = ifelse(df_calculated$vega_geo_mean_no_low_log_distance > 2, 1, 0)
  df_calculated$vega_lowest_no_low_inaccurate_two_order_above = ifelse(df_calculated$vega_lowest_no_low_log_distance > 2, 1, 0)
  df_calculated$vega_geo_mean_only_low_inaccurate_two_order_above = ifelse(df_calculated$vega_geo_mean_only_low_log_distance > 2, 1, 0)
  df_calculated$vega_lowest_only_low_inaccurate_two_order_above = ifelse(df_calculated$vega_lowest_only_low_log_distance > 2, 1, 0)
  df_calculated$ecosar_mean_inaccurate_two_order_above = ifelse(df_calculated$ecosar_mean_log_distance > 2, 1, 0)
  df_calculated$ecosar_low_inaccurate_two_order_above = ifelse(df_calculated$ecosar_low_log_distance > 2, 1, 0)
  df_calculated$test_inaccurate_two_order_above = ifelse(df_calculated$test_log_distance > 2, 1, 0)
  
  df_calculated$vega_geo_mean_no_low_inaccurate_three_order_above = ifelse(df_calculated$vega_geo_mean_no_low_log_distance > 3, 1, 0)
  df_calculated$vega_lowest_no_low_inaccurate_three_order_above = ifelse(df_calculated$vega_lowest_no_low_log_distance > 3, 1, 0)
  df_calculated$vega_geo_mean_only_low_inaccurate_three_order_above = ifelse(df_calculated$vega_geo_mean_only_low_log_distance > 3, 1, 0)
  df_calculated$vega_lowest_only_low_inaccurate_three_order_above = ifelse(df_calculated$vega_lowest_only_low_log_distance > 3, 1, 0)
  df_calculated$ecosar_mean_inaccurate_three_order_above = ifelse(df_calculated$ecosar_mean_log_distance > 3, 1, 0)
  df_calculated$ecosar_low_inaccurate_three_order_above = ifelse(df_calculated$ecosar_low_log_distance > 3, 1, 0)
  df_calculated$test_inaccurate_three_order_above = ifelse(df_calculated$test_log_distance > 3, 1, 0)
  
  
  # Add inaccurate true/false for one order, two orders and three orders (EFSA)
  df_calculated$vega_geo_mean_no_low_efsa_inaccurate_one_order_above = ifelse(df_calculated$vega_geo_mean_no_low_efsa_log_distance > 1, 1, 0)
  df_calculated$vega_lowest_no_low_efsa_inaccurate_one_order_above = ifelse(df_calculated$vega_lowest_no_low_efsa_log_distance > 1, 1, 0)
  df_calculated$vega_geo_mean_only_low_efsa_inaccurate_one_order_above = ifelse(df_calculated$vega_geo_mean_only_low_efsa_log_distance > 1, 1, 0)
  df_calculated$vega_lowest_only_low_efsa_inaccurate_one_order_above = ifelse(df_calculated$vega_lowest_only_low_efsa_log_distance > 1, 1, 0)
  df_calculated$ecosar_mean_efsa_inaccurate_one_order_above = ifelse(df_calculated$ecosar_mean_efsa_log_distance > 1, 1, 0)
  df_calculated$ecosar_low_efsa_inaccurate_one_order_above = ifelse(df_calculated$ecosar_low_efsa_log_distance > 1, 1, 0)
  df_calculated$test_efsa_inaccurate_one_order_above = ifelse(df_calculated$test_efsa_log_distance > 1, 1, 0)
  
  df_calculated$vega_geo_mean_no_low_efsa_inaccurate_two_order_above = ifelse(df_calculated$vega_geo_mean_no_low_efsa_log_distance > 2, 1, 0)
  df_calculated$vega_lowest_no_low_efsa_inaccurate_two_order_above = ifelse(df_calculated$vega_lowest_no_low_efsa_log_distance > 2, 1, 0)
  df_calculated$vega_geo_mean_only_low_efsa_inaccurate_two_order_above = ifelse(df_calculated$vega_geo_mean_only_low_efsa_log_distance > 2, 1, 0)
  df_calculated$vega_lowest_only_low_efsa_inaccurate_two_order_above = ifelse(df_calculated$vega_lowest_only_low_efsa_log_distance > 2, 1, 0)
  df_calculated$ecosar_mean_efsa_inaccurate_two_order_above = ifelse(df_calculated$ecosar_mean_efsa_log_distance > 2, 1, 0)
  df_calculated$ecosar_low_efsa_inaccurate_two_order_above = ifelse(df_calculated$ecosar_low_efsa_log_distance > 2, 1, 0)
  df_calculated$test_efsa_inaccurate_two_order_above = ifelse(df_calculated$test_efsa_log_distance > 2, 1, 0)
  
  df_calculated$vega_geo_mean_no_low_efsa_inaccurate_three_order_above = ifelse(df_calculated$vega_geo_mean_no_low_efsa_log_distance > 3, 1, 0)
  df_calculated$vega_lowest_no_low_efsa_inaccurate_three_order_above = ifelse(df_calculated$vega_lowest_no_low_efsa_log_distance > 3, 1, 0)
  df_calculated$vega_geo_mean_only_low_efsa_inaccurate_three_order_above = ifelse(df_calculated$vega_geo_mean_only_low_efsa_log_distance > 3, 1, 0)
  df_calculated$vega_lowest_only_low_efsa_inaccurate_three_order_above = ifelse(df_calculated$vega_lowest_only_low_efsa_log_distance > 3, 1, 0)
  df_calculated$ecosar_mean_efsa_inaccurate_three_order_above = ifelse(df_calculated$ecosar_mean_efsa_log_distance > 3, 1, 0)
  df_calculated$ecosar_low_efsa_inaccurate_three_order_above = ifelse(df_calculated$ecosar_low_efsa_log_distance > 3, 1, 0)
  df_calculated$test_efsa_inaccurate_three_order_above = ifelse(df_calculated$test_efsa_log_distance > 3, 1, 0)
  
  
  # Add inaccurate true/false for one order, two orders and three orders (ECOTOX)
  df_calculated$vega_geo_mean_no_low_ecotox_inaccurate_one_order_above = ifelse(df_calculated$vega_geo_mean_no_low_ecotox_log_distance > 1, 1, 0)
  df_calculated$vega_lowest_no_low_ecotox_inaccurate_one_order_above = ifelse(df_calculated$vega_lowest_no_low_ecotox_log_distance > 1, 1, 0)
  df_calculated$vega_geo_mean_only_low_ecotox_inaccurate_one_order_above = ifelse(df_calculated$vega_geo_mean_only_low_ecotox_log_distance > 1, 1, 0)
  df_calculated$vega_lowest_only_low_ecotox_inaccurate_one_order_above = ifelse(df_calculated$vega_lowest_only_low_ecotox_log_distance > 1, 1, 0)
  df_calculated$ecosar_mean_ecotox_inaccurate_one_order_above = ifelse(df_calculated$ecosar_mean_ecotox_log_distance > 1, 1, 0)
  df_calculated$ecosar_low_ecotox_inaccurate_one_order_above = ifelse(df_calculated$ecosar_low_ecotox_log_distance > 1, 1, 0)
  df_calculated$test_ecotox_inaccurate_one_order_above = ifelse(df_calculated$test_ecotox_log_distance > 1, 1, 0)
  
  df_calculated$vega_geo_mean_no_low_ecotox_inaccurate_two_order_above = ifelse(df_calculated$vega_geo_mean_no_low_ecotox_log_distance > 2, 1, 0)
  df_calculated$vega_lowest_no_low_ecotox_inaccurate_two_order_above = ifelse(df_calculated$vega_lowest_no_low_ecotox_log_distance > 2, 1, 0)
  df_calculated$vega_geo_mean_only_low_ecotox_inaccurate_two_order_above = ifelse(df_calculated$vega_geo_mean_only_low_ecotox_log_distance > 2, 1, 0)
  df_calculated$vega_lowest_only_low_ecotox_inaccurate_two_order_above = ifelse(df_calculated$vega_lowest_only_low_ecotox_log_distance > 2, 1, 0)
  df_calculated$ecosar_mean_ecotox_inaccurate_two_order_above = ifelse(df_calculated$ecosar_mean_ecotox_log_distance > 2, 1, 0)
  df_calculated$ecosar_low_ecotox_inaccurate_two_order_above = ifelse(df_calculated$ecosar_low_ecotox_log_distance > 2, 1, 0)
  df_calculated$test_ecotox_inaccurate_two_order_above = ifelse(df_calculated$test_ecotox_log_distance > 2, 1, 0)
  
  df_calculated$vega_geo_mean_no_low_ecotox_inaccurate_three_order_above = ifelse(df_calculated$vega_geo_mean_no_low_ecotox_log_distance > 3, 1, 0)
  df_calculated$vega_lowest_no_low_ecotox_inaccurate_three_order_above = ifelse(df_calculated$vega_lowest_no_low_ecotox_log_distance > 3, 1, 0)
  df_calculated$vega_geo_mean_only_low_ecotox_inaccurate_three_order_above = ifelse(df_calculated$vega_geo_mean_only_low_ecotox_log_distance > 3, 1, 0)
  df_calculated$vega_lowest_only_low_ecotox_inaccurate_three_order_above = ifelse(df_calculated$vega_lowest_only_low_ecotox_log_distance > 3, 1, 0)
  df_calculated$ecosar_mean_ecotox_inaccurate_three_order_above = ifelse(df_calculated$ecosar_mean_ecotox_log_distance > 3, 1, 0)
  df_calculated$ecosar_low_ecotox_inaccurate_three_order_above = ifelse(df_calculated$ecosar_low_ecotox_log_distance > 3, 1, 0)
  df_calculated$test_ecotox_inaccurate_three_order_above = ifelse(df_calculated$test_ecotox_log_distance > 3, 1, 0)
  
  
  # Inaccurate below
  
  # Add inaccurate true/false for one order, two orders and three orders above (all data)
  
  df_calculated$vega_geo_mean_no_low_inaccurate_one_order_below = ifelse(df_calculated$vega_geo_mean_no_low_log_distance < -1, 1, 0)
  df_calculated$vega_lowest_no_low_inaccurate_one_order_below = ifelse(df_calculated$vega_lowest_no_low_log_distance < -1, 1, 0)
  df_calculated$vega_geo_mean_only_low_inaccurate_one_order_below = ifelse(df_calculated$vega_geo_mean_only_low_log_distance < -1, 1, 0)
  df_calculated$vega_lowest_only_low_inaccurate_one_order_below = ifelse(df_calculated$vega_lowest_only_low_log_distance < -1, 1, 0)
  df_calculated$ecosar_mean_inaccurate_one_order_below = ifelse(df_calculated$ecosar_mean_log_distance < -1, 1, 0)
  df_calculated$ecosar_low_inaccurate_one_order_below = ifelse(df_calculated$ecosar_low_log_distance < -1, 1, 0)
  df_calculated$test_inaccurate_one_order_below = ifelse(df_calculated$test_log_distance < -1, 1, 0)
  
  df_calculated$vega_geo_mean_no_low_inaccurate_two_order_below = ifelse(df_calculated$vega_geo_mean_no_low_log_distance < -2, 1, 0)
  df_calculated$vega_lowest_no_low_inaccurate_two_order_below = ifelse(df_calculated$vega_lowest_no_low_log_distance < -2, 1, 0)
  df_calculated$vega_geo_mean_only_low_inaccurate_two_order_below = ifelse(df_calculated$vega_geo_mean_only_low_log_distance < -2, 1, 0)
  df_calculated$vega_lowest_only_low_inaccurate_two_order_below = ifelse(df_calculated$vega_lowest_only_low_log_distance < -2, 1, 0)
  df_calculated$ecosar_mean_inaccurate_two_order_below = ifelse(df_calculated$ecosar_mean_log_distance < -2, 1, 0)
  df_calculated$ecosar_low_inaccurate_two_order_below = ifelse(df_calculated$ecosar_low_log_distance < -2, 1, 0)
  df_calculated$test_inaccurate_two_order_below = ifelse(df_calculated$test_log_distance < -2, 1, 0)
  
  df_calculated$vega_geo_mean_no_low_inaccurate_three_order_below = ifelse(df_calculated$vega_geo_mean_no_low_log_distance < -3, 1, 0)
  df_calculated$vega_lowest_no_low_inaccurate_three_order_below = ifelse(df_calculated$vega_lowest_no_low_log_distance < -3, 1, 0)
  df_calculated$vega_geo_mean_only_low_inaccurate_three_order_below = ifelse(df_calculated$vega_geo_mean_only_low_log_distance < -3, 1, 0)
  df_calculated$vega_lowest_only_low_inaccurate_three_order_below = ifelse(df_calculated$vega_lowest_only_low_log_distance < -3, 1, 0)
  df_calculated$ecosar_mean_inaccurate_three_order_below = ifelse(df_calculated$ecosar_mean_log_distance < -3, 1, 0)
  df_calculated$ecosar_low_inaccurate_three_order_below = ifelse(df_calculated$ecosar_low_log_distance < -3, 1, 0)
  df_calculated$test_inaccurate_three_order_below = ifelse(df_calculated$test_log_distance < -3, 1, 0)
  
  
  # Add inaccurate true/false for one order, two orders and three orders (EFSA)
  df_calculated$vega_geo_mean_no_low_efsa_inaccurate_one_order_below = ifelse(df_calculated$vega_geo_mean_no_low_efsa_log_distance < -1, 1, 0)
  df_calculated$vega_lowest_no_low_efsa_inaccurate_one_order_below = ifelse(df_calculated$vega_lowest_no_low_efsa_log_distance < -1, 1, 0)
  df_calculated$vega_geo_mean_only_low_efsa_inaccurate_one_order_below = ifelse(df_calculated$vega_geo_mean_only_low_efsa_log_distance < -1, 1, 0)
  df_calculated$vega_lowest_only_low_efsa_inaccurate_one_order_below = ifelse(df_calculated$vega_lowest_only_low_efsa_log_distance < -1, 1, 0)
  df_calculated$ecosar_mean_efsa_inaccurate_one_order_below = ifelse(df_calculated$ecosar_mean_efsa_log_distance < -1, 1, 0)
  df_calculated$ecosar_low_efsa_inaccurate_one_order_below = ifelse(df_calculated$ecosar_low_efsa_log_distance < -1, 1, 0)
  df_calculated$test_efsa_inaccurate_one_order_below = ifelse(df_calculated$test_efsa_log_distance < -1, 1, 0)
  
  df_calculated$vega_geo_mean_no_low_efsa_inaccurate_two_order_below = ifelse(df_calculated$vega_geo_mean_no_low_efsa_log_distance < -2, 1, 0)
  df_calculated$vega_lowest_no_low_efsa_inaccurate_two_order_below = ifelse(df_calculated$vega_lowest_no_low_efsa_log_distance < -2, 1, 0)
  df_calculated$vega_geo_mean_only_low_efsa_inaccurate_two_order_below = ifelse(df_calculated$vega_geo_mean_only_low_efsa_log_distance < -2, 1, 0)
  df_calculated$vega_lowest_only_low_efsa_inaccurate_two_order_below = ifelse(df_calculated$vega_lowest_only_low_efsa_log_distance < -2, 1, 0)
  df_calculated$ecosar_mean_efsa_inaccurate_two_order_below = ifelse(df_calculated$ecosar_mean_efsa_log_distance < -2, 1, 0)
  df_calculated$ecosar_low_efsa_inaccurate_two_order_below = ifelse(df_calculated$ecosar_low_efsa_log_distance < -2, 1, 0)
  df_calculated$test_efsa_inaccurate_two_order_below = ifelse(df_calculated$test_efsa_log_distance < -2, 1, 0)
  
  df_calculated$vega_geo_mean_no_low_efsa_inaccurate_three_order_below = ifelse(df_calculated$vega_geo_mean_no_low_efsa_log_distance < -3, 1, 0)
  df_calculated$vega_lowest_no_low_efsa_inaccurate_three_order_below = ifelse(df_calculated$vega_lowest_no_low_efsa_log_distance < -3, 1, 0)
  df_calculated$vega_geo_mean_only_low_efsa_inaccurate_three_order_below = ifelse(df_calculated$vega_geo_mean_only_low_efsa_log_distance < -3, 1, 0)
  df_calculated$vega_lowest_only_low_efsa_inaccurate_three_order_below = ifelse(df_calculated$vega_lowest_only_low_efsa_log_distance < -3, 1, 0)
  df_calculated$ecosar_mean_efsa_inaccurate_three_order_below = ifelse(df_calculated$ecosar_mean_efsa_log_distance < -3, 1, 0)
  df_calculated$ecosar_low_efsa_inaccurate_three_order_below = ifelse(df_calculated$ecosar_low_efsa_log_distance < -3, 1, 0)
  df_calculated$test_efsa_inaccurate_three_order_below = ifelse(df_calculated$test_efsa_log_distance < -3, 1, 0)
  
  
  # Add inaccurate true/false for one order, two orders and three orders (ECOTOX)
  df_calculated$vega_geo_mean_no_low_ecotox_inaccurate_one_order_below = ifelse(df_calculated$vega_geo_mean_no_low_ecotox_log_distance < -1, 1, 0)
  df_calculated$vega_lowest_no_low_ecotox_inaccurate_one_order_below = ifelse(df_calculated$vega_lowest_no_low_ecotox_log_distance < -1, 1, 0)
  df_calculated$vega_geo_mean_only_low_ecotox_inaccurate_one_order_below = ifelse(df_calculated$vega_geo_mean_only_low_ecotox_log_distance < -1, 1, 0)
  df_calculated$vega_lowest_only_low_ecotox_inaccurate_one_order_below = ifelse(df_calculated$vega_lowest_only_low_ecotox_log_distance < -1, 1, 0)
  df_calculated$ecosar_mean_ecotox_inaccurate_one_order_below = ifelse(df_calculated$ecosar_mean_ecotox_log_distance < -1, 1, 0)
  df_calculated$ecosar_low_ecotox_inaccurate_one_order_below = ifelse(df_calculated$ecosar_low_ecotox_log_distance < -1, 1, 0)
  df_calculated$test_ecotox_inaccurate_one_order_below = ifelse(df_calculated$test_ecotox_log_distance < -1, 1, 0)
  
  df_calculated$vega_geo_mean_no_low_ecotox_inaccurate_two_order_below = ifelse(df_calculated$vega_geo_mean_no_low_ecotox_log_distance < -2, 1, 0)
  df_calculated$vega_lowest_no_low_ecotox_inaccurate_two_order_below = ifelse(df_calculated$vega_lowest_no_low_ecotox_log_distance < -2, 1, 0)
  df_calculated$vega_geo_mean_only_low_ecotox_inaccurate_two_order_below = ifelse(df_calculated$vega_geo_mean_only_low_ecotox_log_distance < -2, 1, 0)
  df_calculated$vega_lowest_only_low_ecotox_inaccurate_two_order_below = ifelse(df_calculated$vega_lowest_only_low_ecotox_log_distance < -2, 1, 0)
  df_calculated$ecosar_mean_ecotox_inaccurate_two_order_below = ifelse(df_calculated$ecosar_mean_ecotox_log_distance < -2, 1, 0)
  df_calculated$ecosar_low_ecotox_inaccurate_two_order_below = ifelse(df_calculated$ecosar_low_ecotox_log_distance < -2, 1, 0)
  df_calculated$test_ecotox_inaccurate_two_order_below = ifelse(df_calculated$test_ecotox_log_distance < -2, 1, 0)
  
  df_calculated$vega_geo_mean_no_low_ecotox_inaccurate_three_order_below = ifelse(df_calculated$vega_geo_mean_no_low_ecotox_log_distance < -3, 1, 0)
  df_calculated$vega_lowest_no_low_ecotox_inaccurate_three_order_below = ifelse(df_calculated$vega_lowest_no_low_ecotox_log_distance < -3, 1, 0)
  df_calculated$vega_geo_mean_only_low_ecotox_inaccurate_three_order_below = ifelse(df_calculated$vega_geo_mean_only_low_ecotox_log_distance < -3, 1, 0)
  df_calculated$vega_lowest_only_low_ecotox_inaccurate_three_order_below = ifelse(df_calculated$vega_lowest_only_low_ecotox_log_distance < -3, 1, 0)
  df_calculated$ecosar_mean_ecotox_inaccurate_three_order_below = ifelse(df_calculated$ecosar_mean_ecotox_log_distance < -3, 1, 0)
  df_calculated$ecosar_low_ecotox_inaccurate_three_order_below = ifelse(df_calculated$ecosar_low_ecotox_log_distance < -3, 1, 0)
  df_calculated$test_ecotox_inaccurate_three_order_below = ifelse(df_calculated$test_ecotox_log_distance < -3, 1, 0)
  
  
  # Add logkow and pka
  df_calculated = merge(df_calculated, identifiers[,c("original_CAS", "logkow", "pka")], by = 'original_CAS', all.x = T)
  
  
  # Make sure it is numeric
  df_calculated$logkow = as.numeric(df_calculated$logkow)
  
  return(df_calculated)
}

