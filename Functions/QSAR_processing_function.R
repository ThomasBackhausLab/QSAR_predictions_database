################################################################################
#    A function for processing QSAR data from VEGA, ECOSAR and T.E.S.T         #
#                                                                              #
################################################################################
#
# Original author: Patrik Svedberg
#
# Contact email: patrik.svedberg@bioenv.gu.se
# (if the above doesnt work (i.e. years down the line) try p.a.svedberg@gmail.com)
#
# Based on the output of QSAR models:
# 
# Vega 1.1.5 LC50/EC50/NOEC models
# ECOSAR 2.0 (2.2)
# T.E.S.T 5.1.1.0
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

# While some changes has been done prior, versions were not controlled.
# Version 1.0 always runs all QSARs, both wide and long format output

## Version 2.0 changes

# Added option to only run specific QSARs, specified by QSAR_chosen vector
# default set to run all three


################################################################################
# 1. Initiating function and parsing arguments
################################################################################


function(identifiers,
         QSAR_chosen = c('vega', 'ecosar', 'test'),
         QSAR_paths,
         working_directory, 
         run_QSARs = T,
         format = 'wide',
         ecosar_fix = T,
         ecosar_marine = F){
  #' A helping function to run and compile QSAR models
  #'
  #' @Description This function opens and helps to run three QSAR
  #' models; VEGA, ECOSAR and T.E.S.T. It then compiles the results
  #' into a single frame, with mean values for daphnia 48h EC50
  #' and fish 96h EC50. 
  #' 
  #' @param identifiers data.frame. A frame with CAS, SMILES and InChIKeys
  #' for the substances to be processed with the QSAR models
  #' 
  #' @param QSAR_chosen A vector of strings. Each string presents one of the QSAR 
  #' models the script should run. Default is running all of them.
  #' 
  #' @param QSAR_paths vector of strings. A vector with paths to the 
  #' executables of the QSAR models
  #' 
  #' @param working_directory string. 
  #' 
  #' @param run_QSARs boolean. If TRUE, the QSARS will be run. Recommended to 
  #' set to FALSE if rerunning on pre-existing files
  #' 
  #' @param format string, "wide" or "long". Wide format better readability for
  #' the human eye. Long format is better for ggplot and computer readability
  #'
  #' 
  #'
  #' COMPLETE THIS HELP TEXT!
  
  
  ################################################################################
  # 2. Packages and housekeeping
  ################################################################################
  
  # The package handling within functions is from an older version of the script
  # It has been kept for debugging purposes and if any error occurs
  
  #### Prepare loading packages if this is the first time ####
  #library(BiocManager)
  #BiocManager::install(c('ChemmineOB', 'ChemmineR'))
  #
  #### Load packages ####
  # packages <- c('readxl', 'dplyr', 'data.table', 'readr', 'stringr', 'webchem', "ChemmineR", "ChemmineOB", 'tidyr', 'data.table', 'rcdk')
  # 
  # tryCatch(
  #   expr = {lapply(packages, library, character.only = TRUE)
  #     },
  #   error = function(e){
  #     message('Error loading packages, make sure all the following packages are loaded:')
  #     print(packages)
  #   },
  #   warning = function(w){
  #     message('Warnings issued while loading packages. Continuing function.')
  #   }
  # )
  # 
  ## Set function version
  version = 2
  
  #######################
  
  ## Check provided arguments
  
  # Identifiers
  if(missing(identifiers)){ 
    message('Identifiers not provided. Exiting function.')
    return(FALSE)
  }
  
  # Qsars chosen
  if('vega' %in% tolower(QSAR_chosen)){
    run_vega = T
  } else {
    run_vega = F
  }
  
  if('ecosar' %in% tolower(QSAR_chosen)){
    run_ecosar = T
  } else {
    run_ecosar = F
  }
  
  if('test' %in% tolower(QSAR_chosen) | 't.e.s.t.' %in% tolower(QSAR_chosen)){
    run_test = T
  } else {
    run_test = F
  }
  
  if(!any(run_vega, run_ecosar, run_test)){
    message('No QSARs chosen, either run default (all) or specify QSARs. Exiting function')
    return(FALSE)
  }
  
  # Qsar paths
  if(missing(QSAR_paths)){
    message('QSAR_paths not provided. Using defaults.')
    QSAR_paths = c('C:/Program Files/Vega 1.1.5/vega-1.1.5-b48/VEGA.jar',
                   'C:/Program Files/Ecosar 2.2/ecosarapplication/bin/ecosarapplication64.exe',
                   'C:/Program Files/TEST 5.1.1.0/TEST.exe')
  }
  
  if(missing(working_directory)){
    message('Identifiers not provided. Exiting function.')
    return(FALSE)
  }
  
  message(paste0('Arguments passed successfully. Running the follwoing QSARs: ', ifelse(run_vega, 'vega ',''), ifelse(run_ecosar, 'ecosar ',''), ifelse(run_test, 't.e.s.t. ','')))
  
  
  ## Create and set working directories:
  
  # identifier directory
  identifierwd <- paste0(working_directory, '/identifiers')
  if(!file.exists(identifierwd)){
    
    dir.create(identifierwd)
    
  }
  
  # QSAR-specific directories
  if(run_vega){
    vegawd <- paste0(working_directory, '/VEGAoutput')
    if(!file.exists(vegawd)){
      
      dir.create(vegawd)
      
    }
  }
  
  if(run_ecosar){
    ecosarwd <- paste0(working_directory, '/ECOSARoutput')
    if(!file.exists(ecosarwd)){
      
      dir.create(ecosarwd)
      
    }
  }
  
  if(run_test){
    testwd1 <- paste0(working_directory, '/TESToutput')
    if(!file.exists(testwd1)){
      
      dir.create(testwd1)
      
    }
    testwd2 <- paste0(testwd1, '/MyToxicity')
    if(!file.exists(testwd2)){
      
      dir.create(testwd2)
      
    }
  }
  
  
  # General outwd
  outwd <- paste0(working_directory, '/mergedoutput')
  if(!file.exists(outwd)){
    
    dir.create(outwd)
    
  }
  
  # Auxiliary wd for external databases (experimental, monitoring etc)
  auxwd <- paste0(working_directory, '/auxiliaryfiles')
  if(!file.exists(auxwd)){
    
    dir.create(auxwd)
    
  }
  
  
  ## Set paths to your QSAR models (for opening them remotely)
  
  if(run_vega){
    if(any(grepl(QSAR_paths, pattern = 'VEGA.jar'))){
      vegapath <- QSAR_paths[grepl(QSAR_paths, pattern = 'VEGA.jar')]
    } else {
      message('VEGA.jar not supplied. Exiting function.')
      return(FALSE)
    }
  }
  
  if(run_ecosar){
    if(any(grepl(QSAR_paths, pattern = 'ecosarapplication64.exe'))){
      ecosarpath <- QSAR_paths[grepl(QSAR_paths, pattern = 'ecosarapplication64.exe')]
    } else {
      message('ecosarapplication64.exe not supplied. Exiting function.')
      return(FALSE)
    }
  }
  
  if(run_test){
    if(any(grepl(QSAR_paths, pattern = 'TEST.exe'))){
      testpath <- QSAR_paths[grepl(QSAR_paths, pattern = 'TEST.exe')]
    } else {
      message('TEST.exe not supplied. Exiting function.')
      return(FALSE)
    }
  }
  
  # Set filenames of the QSAR output files (these will have to be supplied manually into the QSAR applications, and then identified by the script)
  vegafile = 'report_summary.txt'
  ecosarfile = 'ecosar_output.xlsx'
  testfiles = c('Batch_Daphnia_magna_LC50_(48_hr)_Consensus.csv',
                'Batch_Fathead_minnow_LC50_(96_hr)_Consensus.csv')
  


  ################################################################################
  #                               3. Running QSARs                               #
  ################################################################################
  if(run_QSARs){
    
    ##  Run VEGA
    if(run_vega){
      print('Starting the Vega application')
      
      shell.exec(vegapath)
      
      temp = readline(prompt = 'Do you need instructions on VEGA handling? (yes/no):')
      
      if(temp == 'yes'){
        
        cat(
          'Vega will open in a separate window. Click "import file" on the right of that window.
          Navigate to your SMILES-file as exported above (should be in', paste0(working_directory, '/identifiers/'), ') and click open
          On the left of the Vega-window click "select" to select which models should be applied
          Open the "Ecotox"-tab (at the top, second from the left) and select the following models:

          "Fish Acute (LC50) Toxicity model (KNN/Read-Across) (version 1.0.0)"
          "Fish Acute (LC50) Toxicity model (NIC) (version 1.0.0)"
          "Fish Acute (LC50) Toxicity model (IRFMN) (version 1.0.0)"
          "Fish Acute (LC50) Toxicity model (IRFMN/Combase) (version 1.0.0)"
          "Fish Chronic (NOEC) Toxicity model (IRFMN) (version 1.0.0)"
          "Fathead Minnow LC50 96h (EPA) (version 1.0.7)"
          "Fathead Minnow LC50 model (KNN/IRFMN) (version 1.1.0)"
          "Guppy LC50 model (KNN/IRFMN) (version 1.1.0)"
          "Daphnia Magna LC50 48h (EPA) (version 1.0.7)"
          "Daphnia Magna LC50 48h (DEMETRA) (version 1.0.4)"
          "Daphnia Magna Acute (EC50) Toxicity model (IRFMN) (version 1.0.0)"
          "Daphnia Magna Acute (EC50) Toxicity model (IRFMN/Combase) (version 1.0.0)"
          "Daphnia Magna Chronic (NOEC) Toxicity model (IRFMN) (version 1.0.0)"
          "Algae Acute (EC50) Toxicity model (IRFMN) (version 1.0.0)"
          "Algae Acute (EC50) Toxicity model (ProtoQSAR/Combase) (version 1.0.0)"
          "Algae Chronic (NOEC) Toxicity model (IRFMN) (version 1.0.0)"
          "Sludge (EC50) Toxicity model (ProtoQSAR/Combase) (version 1.0.0)"

          Click the "EXPORT"-button on the left and click to checkmark "summary (plain text file)" on the right of that tab.
          Enter, or use the button on the right of the textbox to navigate to, your vega-output folder
          Click predict and wait for results')
        
        readline("Have you read and understood the instructions? (press ENTER to continue)")
        
      } else if(temp == 'no'){
        
        print('No instructions needed. Continuing function')
        
      } else {
        
        print('Answer not recognized, interpreting as "no". Continuing function')
        
      }
      
      print(paste0('Set the report folder to the following: ', vegawd))
      
      question1 <- readline("Has VEGA finished its process? (press ENTER to continue)")
    }
    
    
    ## Run ECOSAR
    if(run_ecosar){
      print('Starting the ECOSAR application')
      
      shell.exec(ecosarpath)
      
      temp = readline(prompt = 'Do you need instructions on ECOSAR handling? (yes/no):')
      
      if(temp == 'yes'){
        
        cat(
        'Ecosar will open in a separate window. The first page will show an EULA-type agreement.
        Read the short agreement and click accept (if you do accept the agreement, otherwise this script may not be applied properly). 
        You will end up in the tab "Organic Molecule".
        Click the "Batch" button on the right side of the ECOSAR window. Click "Load" on the right side.
        A file browser will open. Navigate to your identifiers-file (should be in', paste0(working_directory, '/identifiers/'), 
        'ECOSAR accepts both SMILES and CAS, but here we use SMILES.
        Note that if you are running this script for more than 1000 compoinds, you will have SMILES-files with names like SMILES_1to1000 and so on until SMILES_X001toEnd.
        ECOSAR is very slow with large datasets, and speed is improved by running smaller batches - you will need to run this for each of the files.
        Select the file you want to use, and click "open".
        ECOSAR may take some time to identify the substances entered. Usually in the range of 1-10 minutes. 
        You will now see the CAS, name and SMILES of all identified substances (Unidentified will lack name and/or SMILES or CAS depending on which identifier you supplied)
        
        The following step is important:
          click "Report" and navigate to your ECOSAR-output folder, name your output file and click save.
          When naming, use "ecosar_output" if you have less than 1000 substances, else put ecosar_output1 for the first file, and 2 for the second etc.
        !!!DO NOT click "submit" as this will run all substances individually in the ecosar window and present results visually instead of an output file.!!!
        
        Wait for results. This can take hours even with splitting up the compounds.')
        
        readline("Have you read and understood the instructions? (press ENTER to continue)")
        
      } else if(temp == 'no'){
        
        print('No instructions needed. Continuing function')
        
      } else {
        
        print('Answer not recognized, interpreting as "no". Continuing function')
        
      }
      
      print(paste0('Set the report folder to the following: ', ecosarwd))
      print(paste0('Input the following filename: ', ecosarfile, ', or if multiple split SMILES-lists add a 1, 2 or 3 etc before the file-extension depending on which SMILES-file you used'))
      
      question2 <- readline("Has ECOSAR finished its process? (press ENTER to continue)")
      
    }
    
    ## Run T.E.S.T.
    if(run_test){
      print('Starting the T.E.S.T. application')
      
      shell.exec(testpath)
      
      temp = readline(prompt = 'Do you need instructions on T.E.S.T. handling? (yes/no):')
      
      if(tolower(temp) == 'yes'){
        
        cat(
          'T.E.S.T. will open in a separate window.
          From the first page, find the turquoise "Switch to Batch Mode"-button on the bottom right and click it.
          In Batch Mode, there will be an empty text-box in the upper left part of the window. 
          Manually copy and paste the identifiers (use SMILES here) from your identifiers-file (you should find your file in ', 
          paste0(working_directory, '/identifiers/'),') 
          Make sure that the drop down menu to the top left says "Smiles" and not automatic.
          Click "Search" to load in the identifiers into the table on the right. 
          In this table you can see if any compounds are unidentified or exhibits any other error.
          On the left you can chose endpoint and method. 
          For this script we will run fathead minnow LC50 and Daphnia magna LC50, but these will each be run separately. 
          For now, pick the Daphnia magna endpoint. 
          For method, use the default "Concensus"-method.
          Further down is a field for output-folder; click "Browse..." and navigate to your test-ouptut-folder.
          Click the green "Calculate!"-button on the bottom right to start calculations for your selected model (a results-window opens)
          In the results window wait for the predictions to finish, then click "save to text (.csv)", which will open a file browser (which should already be in your chosen output folder)
          Click "save". Here the file may be opened automatically by your default CSV-reader (e.g. Excell), just close it.
          Close the results window and switch "Endpoint" (if you follow this guide, it will be Fathead minnow LC50 now) and
          click the green "Calculate!"-button once again. And once more save the results by clicking "save to text (.csv)".
          
          The results are now in and ready to be imported!')
          
        readline("Have you read and understood the instructions? (press ENTER to continue)")
        
      } else if(tolower(temp) == 'no'){
        
        print('No instructions needed. Continuing function')
        
      } else {
        
        print('Answer not recognized, interpreting as "no". Continuing function')
        
      }
      
      question3 <- readline("Has T.E.S.T. finished its process? (press ENTER to continue)")
    }
    
  }
  
  
  ################################################################################
  #               4. Importing and standardizing VEGA model data                 #
  ################################################################################
  if(run_vega){
    
    print('Processing VEGA output')
    
    ## Vega import from summary text file
    
    
    vega_dataframe = as.data.frame(fread(paste0(vegawd, '/', vegafile), na.strings = c('N/A', '[ERRROR]')))
    
    
    # Remove columns with unitless info
    vega_dataframe = vega_dataframe[,colnames(vega_dataframe)[colnames(vega_dataframe) %in% c('Id', "SMILES") | grepl(colnames(vega_dataframe), pattern = 'assessment')]]
    
    
    # Split prediction quality from numeric value (and remove unit, it is all mg/l)
    for(i in seq(3, 2*ncol(vega_dataframe)-3, by=2)){
      
      # Get current column
      current_colname = colnames(vega_dataframe)[i]
      
      vega_dataframe = separate(
        vega_dataframe,
        current_colname,
        c(current_colname, paste0(current_colname, ' - reliability')),
        sep = " mg/L \\("
      )
      
      # Remove trailing value) or reliability) from reliability column
      vega_dataframe[,c(paste0(current_colname, ' - reliability'))] = str_replace(vega_dataframe[,c(paste0(current_colname, ' - reliability'))], pattern = ' value\\)$| reliability\\)$', replacement = '')
      
    }
    
    # Rename SMILES to vega_SMILES
    vega_dataframe = vega_dataframe %>% rename(vega_SMILES = SMILES)
    
    # if no original_CAS supplied
    if(all(is.na(identifiers$original_CAS))){
      
      # Repack Id to molecule number
      vega_dataframe$molecule_number = as.numeric(str_extract(vega_dataframe$Id, pattern = ' [0-9]*$'))
      
      # Add original smiles and Inchikeys
      vega_dataframe = merge(vega_dataframe, identifiers[,c("original_CAS", "original_SMILES", "InChIKey", "molecule_number")], by = 'molecule_number', all.x = T)
      
      # Drop molecule number and ID to fit with the rest of the script
      vega_dataframe$molecule_number = NULL
      vega_dataframe$Id = NULL
      
      # And move original_CAS to the front (even if it is empty)
      vega_dataframe = vega_dataframe %>% 
        relocate(original_CAS)
      
      
    } else {
      
      # Add original smiles and Inchikeys
      vega_dataframe = merge(vega_dataframe, identifiers[,c("original_CAS", "original_SMILES", "InChIKey")], by.x = 'Id', by.y = 'original_CAS', all.x = T)
      
      # Rename Id to original_CAS
      vega_dataframe = vega_dataframe %>% 
        rename(original_CAS = Id)
      
    }
    
    # Relocate identifiers to beginning
    vega_dataframe = vega_dataframe %>%
      relocate(original_SMILES, .before = vega_SMILES) %>%
      relocate(InChIKey, .after = original_CAS)
    
    # Add duration in all columns
    vega_durations = data.frame('models' = c("Fish Acute (LC50) Toxicity model (KNN/Read-Across) (version 1.0.0)",
                                             "Fish Acute (LC50) Toxicity model (NIC) (version 1.0.0)",
                                             "Fish Acute (LC50) Toxicity model (IRFMN) (version 1.0.0)",
                                             "Fish Acute (LC50) Toxicity model (IRFMN/Combase) (version 1.0.0)",
                                             "Fish Chronic (NOEC) Toxicity model (IRFMN) (version 1.0.0)",
                                             "Fathead Minnow LC50 96h (EPA) (version 1.0.7)",
                                             "Fathead Minnow LC50 model (KNN/IRFMN) (version 1.1.0)",
                                             "Guppy LC50 model (KNN/IRFMN) (version 1.1.0)",
                                             "Daphnia Magna LC50 48h (EPA) (version 1.0.7)",
                                             "Daphnia Magna LC50 48h (DEMETRA) (version 1.0.4)",
                                             "Daphnia Magna Acute (EC50) Toxicity model (IRFMN) (version 1.0.0)",
                                             "Daphnia Magna Acute (EC50) Toxicity model (IRFMN/Combase) (version 1.0.0)",
                                             "Daphnia Magna Chronic (NOEC) Toxicity model (IRFMN) (version 1.0.0)",
                                             "Algae Acute (EC50) Toxicity model (IRFMN) (version 1.0.0)",
                                             "Algae Acute (EC50) Toxicity model (ProtoQSAR/Combase) (version 1.0.0)",
                                             "Algae Chronic (NOEC) Toxicity model (IRFMN) (version 1.0.0)",
                                             "Sludge (EC50) Toxicity model (ProtoQSAR/Combase) (version 1.0.0)"),
                                'durations' = c('96h', 
                                                '96h',
                                                '96h',
                                                '96h',
                                                'ELS',
                                                '96h',
                                                '96h',
                                                '96h', 
                                                '48h',
                                                '48h',
                                                '48h',
                                                '48h',
                                                '21d',
                                                '72h',
                                                '72h',
                                                '72h',
                                                '3h'
                                )
    )
    
    # Add durations to model names lacking it where we have duration information
    vega_durations$model_withduration = apply(vega_durations, FUN = function(x){
      if(is.na(x[2]) | grepl(x[1], pattern = x[2])){
        return(x[1])
      } else {
        temp = str_replace(x[1], pattern = '(?<=\\) )(?=\\(version)', replacement = paste0(x[2], ' '))
        return(temp)  
      }
    },
    MARGIN = 1)
    
    # Drop version numbers (to match with vega summary file output)
    vega_durations$models = sapply(vega_durations$models, FUN = function(x){str_remove_all(x, pattern = ' \\(version [0-9\\.]{5}\\)')})
    vega_durations$model_withduration = sapply(vega_durations$model_withduration, FUN = function(x){str_remove_all(x, pattern = ' \\(version [0-9\\.]{5}\\)')})
    
    # Then we replace model names in the columns of vega_dataframe with the new versions
    for(i in 1:nrow(vega_durations)){
      colnames(vega_dataframe) = sapply(colnames(vega_dataframe), FUN = function(x){gsub(pattern = vega_durations[i,1], replacement = vega_durations[i,3], x = x, fixed = T)})
    }
    
    
    # Add [mg/l] unit to all prediction columns
    colnames(vega_dataframe) = ifelse(grepl(colnames(vega_dataframe), pattern = 'reliability'),
                                      colnames(vega_dataframe),
                                      str_replace(colnames(vega_dataframe), pattern = '(?<=t)(?=$)', replacement = ' [mg/l]'))
    
    # Add platform name to all non-identifier columns
    colnames(vega_dataframe)[5:length(colnames(vega_dataframe))] = str_replace(colnames(vega_dataframe)[5:length(colnames(vega_dataframe))], pattern = '(?=^)', replacement = 'VEGA_')
    
    if(format == 'wide'){
      
      vega_dataframe_wide = vega_dataframe
      
    } else if(format == 'long'){
      
      # Repack vega into long format
      identifier_columns = c('original_CAS',
                             'InChIKey',
                             'original_SMILES',
                             'vega_SMILES')
      
      # We split on two steps for values and reliabilities
      value_columns = c("VEGA_Fish Acute (LC50) Toxicity model (KNN/Read-Across) 96h - assessment [mg/l]",
                        "VEGA_Fish Acute (LC50) Toxicity model (NIC) 96h - assessment [mg/l]",
                        "VEGA_Fish Acute (LC50) Toxicity model (IRFMN) 96h - assessment [mg/l]",
                        "VEGA_Fish Acute (LC50) Toxicity model (IRFMN/Combase) 96h - assessment [mg/l]",
                        "VEGA_Fish Chronic (NOEC) Toxicity model (IRFMN) ELS - assessment [mg/l]",
                        "VEGA_Fathead Minnow LC50 96h (EPA) - assessment [mg/l]",
                        "VEGA_Fathead Minnow LC50 model (KNN/IRFMN) 96h - assessment [mg/l]",
                        "VEGA_Guppy LC50 model (KNN/IRFMN) 96h - assessment [mg/l]",
                        "VEGA_Daphnia Magna LC50 48h (EPA) - assessment [mg/l]",
                        "VEGA_Daphnia Magna LC50 48h (DEMETRA) - assessment [mg/l]",
                        "VEGA_Daphnia Magna Acute (EC50) Toxicity model (IRFMN) 48h - assessment [mg/l]",
                        "VEGA_Daphnia Magna Acute (EC50) Toxicity model (IRFMN/Combase) 48h - assessment [mg/l]",
                        "VEGA_Daphnia Magna Chronic (NOEC) Toxicity model (IRFMN) 21d - assessment [mg/l]",
                        "VEGA_Algae Acute (EC50) Toxicity model (IRFMN) 72h - assessment [mg/l]",
                        "VEGA_Algae Acute (EC50) Toxicity model (ProtoQSAR/Combase) 72h - assessment [mg/l]",
                        "VEGA_Algae Chronic (NOEC) Toxicity model (IRFMN) 72h - assessment [mg/l]",
                        "VEGA_Sludge (EC50) Toxicity model (ProtoQSAR/Combase) 3h - assessment [mg/l]")
      
      reliability_columns = c("VEGA_Fish Acute (LC50) Toxicity model (KNN/Read-Across) 96h - assessment - reliability",
                              "VEGA_Fish Acute (LC50) Toxicity model (NIC) 96h - assessment - reliability",
                              "VEGA_Fish Acute (LC50) Toxicity model (IRFMN) 96h - assessment - reliability",
                              "VEGA_Fish Acute (LC50) Toxicity model (IRFMN/Combase) 96h - assessment - reliability",
                              "VEGA_Fish Chronic (NOEC) Toxicity model (IRFMN) ELS - assessment - reliability",
                              "VEGA_Fathead Minnow LC50 96h (EPA) - assessment - reliability",
                              "VEGA_Fathead Minnow LC50 model (KNN/IRFMN) 96h - assessment - reliability",
                              "VEGA_Guppy LC50 model (KNN/IRFMN) 96h - assessment - reliability",
                              "VEGA_Daphnia Magna LC50 48h (EPA) - assessment - reliability",
                              "VEGA_Daphnia Magna LC50 48h (DEMETRA) - assessment - reliability",
                              "VEGA_Daphnia Magna Acute (EC50) Toxicity model (IRFMN) 48h - assessment - reliability",
                              "VEGA_Daphnia Magna Acute (EC50) Toxicity model (IRFMN/Combase) 48h - assessment - reliability",
                              "VEGA_Daphnia Magna Chronic (NOEC) Toxicity model (IRFMN) 21d - assessment - reliability",
                              "VEGA_Algae Acute (EC50) Toxicity model (IRFMN) 72h - assessment - reliability",
                              "VEGA_Algae Acute (EC50) Toxicity model (ProtoQSAR/Combase) 72h - assessment - reliability",
                              "VEGA_Algae Chronic (NOEC) Toxicity model (IRFMN) 72h - assessment - reliability",
                              "VEGA_Sludge (EC50) Toxicity model (ProtoQSAR/Combase) 3h - assessment - reliability")
      
      vega_values = vega_dataframe[,colnames(vega_dataframe) %in% append(identifier_columns, value_columns)]
      
      vega_reliabilities = vega_dataframe[,colnames(vega_dataframe) %in% append(identifier_columns, reliability_columns)]
      
      vega_values_long = gather(vega_values, 
                                key = model_name,
                                value = value,
                                all_of(value_columns),
                                factor_key = T, 
                                na.rm = F)
      
      vega_values_long[,"model_name"] = str_replace(vega_values_long[,"model_name"], pattern = ' - assessment \\[mg/l]$', replacement = '')
      
      vega_reliabilities_long = gather(vega_reliabilities, 
                                key = model_name,
                                value = reliability,
                                all_of(reliability_columns),
                                factor_key = T, 
                                na.rm = F)
      
      
      vega_dataframe_long = cbind(vega_values_long, vega_reliabilities_long[,"reliability"])
      
      vega_dataframe_long = vega_dataframe_long %>%
        rename(reliability = 'vega_reliabilities_long[, \"reliability\"]') %>%
        rename(model = model_name)
      
      # Fix duration, model organism and model endpoint
      vega_dataframe_long$duration = apply(vega_dataframe_long,MARGIN = 1, FUN = function(x){str_extract(x[5], pattern = '([0-9]{2}h)|(ELS)|([0-9]{2}d)')})
      
      vega_dataframe_long$model_organism = apply(vega_dataframe_long,MARGIN = 1, FUN = function(x){str_extract(x[5], pattern = '([Aa]lgae)|([Dd]aphnia [Mm]agna)|([Ff]athead [Mm]innow)|([Ff]ish)|([Gg]uppy)|([Ss]ludge)')})
      
      # Change model organism to "fish" if guppy or fathead minnow
      vega_dataframe_long$model_organism = ifelse(grepl(vega_dataframe_long$model_organism, pattern = '([Gg]uppy)|([Ff]athead [Mm]innow)'),
                                                  'fish',
                                                  vega_dataframe_long$model_organism)
      
      vega_dataframe_long$model_endpoint = apply(vega_dataframe_long,MARGIN = 1, FUN = function(x){str_extract(x[5], pattern = '(EC50)|(LC50)|(NOEC)')})
      
      # Set QSAR_tool to VEGA
      vega_dataframe_long$QSAR_tool = 'VEGA'
      
      # Cleanup in aisle 5 aka remove NA predictions
      vega_dataframe_long = vega_dataframe_long[!is.na(vega_dataframe_long$value),]
      
    }
    
  }  
  
  ################################################################################
  #               5. Importing and standardizing ECOSAR model data               #
  ################################################################################
  
  if(run_ecosar){
    
    print('Processing ECOSAR output')
    
    print('Before reading file')
    
    ## Import ECOSAR data (if single file <1001 compounds, else try to dir the files and order them)
    if(file.exists(paste0(ecosarwd, '/', ecosarfile))){
      ecosar_dataframe <- suppressWarnings(as.data.frame(readxl::read_xlsx(paste0(ecosarwd, '/', ecosarfile), col_types = c("numeric", "text", "text", "text", "text", "text", "text", "text", "numeric", "numeric", "text", "numeric", "numeric", "numeric", "numeric"))))
      
      # Ecosar 2.2 adds flag info at the bottom, remove it and save it in a separate frame (we dont use it but good to have)
      ecosar_flaginfo = ecosar_dataframe[is.na(ecosar_dataframe$Number),]
      ecosar_dataframe = ecosar_dataframe[!is.na(ecosar_dataframe$Number),]
    } else if(file.exists(str_replace(paste0(ecosarwd, '/', ecosarfile), pattern = '\\.', replacement = '1.'))){
      # If ecosar file doesnt exist, see if there are multiple files
      file_list = dir(paste0(ecosarwd, '/'))[grepl(x = dir(paste0(ecosarwd, '/')), pattern = '[0-9]\\.xlsx')]
      ecosar_dataframe = data.frame()
      
      # loop over files and bind them into a single frame
      for(file in file_list){
        
        current_filenumber = as.numeric(str_extract(file, pattern = '[0-9]{1,2}'))
        
        temp_frame = suppressWarnings(as.data.frame(readxl::read_xlsx(paste0(ecosarwd, '/', file), col_types = c("numeric", "text", "text", "text", "text", "text", "text", "text", "numeric", "numeric", "text", "numeric", "numeric", "numeric", "numeric"))))
        
        # Remove flag info for each file
        ecosar_flaginfo = temp_frame[is.na(temp_frame$Number),]
        temp_frame = temp_frame[!is.na(temp_frame$Number),]
        
        # fix molecule number to fit number of files
        temp_frame$Number = temp_frame$Number + (current_filenumber - 1) * 1000
        
        ecosar_dataframe = rbind(ecosar_dataframe, temp_frame)
        
      }
      
    }
    
    # In case of duplicate rows (happens sometimes if the output file has been tampered with) remove duplicates
    ecosar_dataframe = distinct(ecosar_dataframe)
    
    # Add 1 to the ecosar dataframe number field (ecosar starts at 0, our numbering starts at 1)
    ecosar_dataframe$Number = ecosar_dataframe$Number + 1
    
    # Check that we have all the molecules processed in ecosar
    ecosar_unprocessed = identifiers[!identifiers$molecule_number %in% ecosar_dataframe$Number,]
    
    rerun_ecosar = FALSE
    if(nrow(ecosar_unprocessed) > 0){
      
      rerun_ecosar = TRUE
      
      # Add new numbering to ecosar_unprocessed (order number for import)
      ecosar_unprocessed$molecule_number_rerun = seq(1,nrow(ecosar_unprocessed), 1)
      
    }
    
    if(nrow(ecosar_unprocessed) > 0 & !file.exists(paste0(ecosarwd, '/ecosar_rerun_output.xlsx'))){
      
      rerun_ecosar = TRUE
      
      message('Ecosar unable to process SMILES for all molecules!')
      message('rerun ecosar with CAS for molecules: ')
      print(ecosar_unprocessed[,c("molecule_number", "original_CAS")])
      
      write.table(ecosar_unprocessed[,"original_CAS"],
                  file = paste0(identifierwd, '/ECOSARrerunCAS.txt'),
                  col.names = F,
                  row.names = F,
                  quote = F)
      
      message(paste0('Using the following identifier file ', paste0(identifierwd, '/ECOSARrerunCAS.txt')))
      message(paste0('and the following output file ', paste0(ecosarwd, '/ecosar_rerun_output.xlsx')))
      
      # Starting ECOSAR in case you closed it
      shell.exec(ecosarpath)
      
      question2 <- readline("Has ECOSAR finished its rerun process? (press ENTER to continue)")
      
    }
    
    if(rerun_ecosar & file.exists(paste0(ecosarwd, '/ecosar_rerun_output.xlsx'))){
      
      rerun_results = suppressWarnings(as.data.frame(readxl::read_xlsx(paste0(ecosarwd, '/ecosar_rerun_output.xlsx'), col_types = c("numeric", "text", "text", "text", "text", "text", "text", "text", "numeric", "numeric", "text", "numeric", "numeric", "numeric", "numeric"))))
      
      # Ecosar 2.2 adds flag info at the bottom, remove it
      rerun_results = rerun_results[!is.na(rerun_results$Number),]
      
      # In case of duplicate rows (happens sometimes if the output file has been tampered with) remove duplicates
      rerun_results = distinct(rerun_results)
      
      # Add 1 to the ecosar dataframe number field (ecosar starts at 0, our numbering starts at 1)
      rerun_results$Number = rerun_results$Number + 1
      
      
      # Add the old molecule number instead of the new
      rerun_results = merge(rerun_results, ecosar_unprocessed[,c("molecule_number", "molecule_number_rerun")], by.x = 'Number', by.y = 'molecule_number_rerun')
      
      # Replace Number with molecule number
      rerun_results$Number = rerun_results$molecule_number
      
      # And remove molecule number
      rerun_results = rerun_results[,!colnames(rerun_results) %in% c('molecule_number')]
      
      # Add these rows to the other ECOSAR data
      ecosar_dataframe = rbind(ecosar_dataframe, rerun_results)
      
      # Double check if we still have missing molecules
      ecosar_unprocessed_2 = identifiers[!identifiers$molecule_number %in% ecosar_dataframe$Number,]
      
      message('Rerun of ECOSAR missing entries done')
      if(nrow(ecosar_unprocessed_2) >0){
        
        message('Still missing the following chemicals: ')
        print(ecosar_unprocessed_2[,c("molecule_number", "original_CAS")])
        message('They are probably not available in ECOSAR')
        
      }
      
    } else if(rerun_ecosar){
      
      message('No reprocessing of ECOSAR done, please do and rerun this part of the script')
      
    }
    
    ## Make a CSV version of the data (for other applications)
    # translate name from .xlsx to .csv
    ecosar_csvfile <- stringr::str_replace(string = ecosarfile, pattern = '.xlsx$', '.csv')
    # export as csv
    write.table(ecosar_dataframe, file = paste0(ecosarwd, '/', ecosar_csvfile), sep = '\t', row.names = F)
    
    
    # Rename SMILES to ecosar_SMILES to keep track of source of SMILES
    ecosar_dataframe = ecosar_dataframe %>%
      rename(ecosar_SMILES = SMILES) %>%
      rename(ecosar_CAS = CAS) %>%
      rename(duration = Duration)
    
    # Make a CAS_fixed version of ecosar_CAS 
    ecosar_dataframe$ecosar_CAS_fixed = as.cas(ecosar_dataframe$ecosar_CAS)
    
    
    # Merge ecosar predictions with original identifiers
    ecosar_dataframe = merge(ecosar_dataframe, identifiers, by.x = 'Number', by.y = 'molecule_number', all.x = T)
    
    
    
    # Add ECOSAR in the QSAR_tool field of the frame
    ecosar_dataframe$QSAR_tool <- 'ECOSAR'
    
    # Translate durations into hours and make new column
    ecosar_dataframe$'model_duration (h)' <- ifelse(grepl(ecosar_dataframe$duration, pattern = 'd'),
                                                    as.numeric(str_extract(ecosar_dataframe$duration, 
                                                                           pattern = '[0-9]{1,3}(?=d)')) * 24,
                                                    as.numeric(str_extract(ecosar_dataframe$duration, 
                                                                           pattern = '[0-9]{1,3}(?=h)')))
    
    
    
    # Fix some more column names in the ECOSAR output
    
    ecosar_dataframe <- 
      ecosar_dataframe %>%
      rename(model_endpoint = `End Point`) %>%
      rename(model_organism = `Organism`) %>%
      rename(`Predicted value [mg/L]` = `Concentration (mg/L)`)
    
    
    if(ecosar_marine){
      
      # Remove >1000 g/mol substances (AD)
      ecosar_outside_AD = ecosar_dataframe[ecosar_dataframe$`Molecular Weight (g/mol)` >= 1000 |
                                             grepl(x = ecosar_dataframe$Alert, pattern = 'LogKowCutOff'),]
      
      # Remove these entries from the resulting data
      ecosar_dataframe <- ecosar_dataframe[ecosar_dataframe$`Molecular Weight (g/mol)` < 1000 &
                                             !grepl(x = ecosar_dataframe$Alert, pattern = 'LogKowCutOff'),]
      
      
      
    } else {
      
      ## Extract and all seawater and substances not within the AD (based on molweight)
      
      ecosar_seawater <- ecosar_dataframe[grepl(ecosar_dataframe$model_organism, pattern = '\\(SW\\)'),]
      
      ecosar_outside_AD = ecosar_dataframe[ecosar_dataframe$`Molecular Weight (g/mol)` >= 1000 |
                                             grepl(x = ecosar_dataframe$Alert, pattern = 'LogKowCutOff'),]
      
      # Remove these entries from the resulting data
      ecosar_dataframe <- ecosar_dataframe[!grepl(ecosar_dataframe$model_organism, pattern = '\\(SW\\)') &
                                             ecosar_dataframe$`Molecular Weight (g/mol)` < 1000 &
                                             !grepl(x = ecosar_dataframe$Alert, pattern = 'LogKowCutOff'),]
      
      
    }
    
    
    
    # Deduplicate the dataframe (If you dont clear the memory before running there might be duplicate rows. Remove them)
    ecosar_dataframe = distinct(ecosar_dataframe)
    
    if(format == 'long'){
      
      ecosar_dataframe_long = ecosar_dataframe
      
      # Fix some columns to improve merging
      ecosar_dataframe_long = ecosar_dataframe_long %>%
        rename(error_warning = Alert) %>%
        rename(value = `Predicted value [mg/L]`)
      
      # Add "model" column
      ecosar_dataframe_long$model = apply(ecosar_dataframe_long, MARGIN = 1, FUN = function(x){paste('ECOSAR', x[5], x[6], x[8], x[7], sep = '_')})
      
      ecosar_dataframe_long$model = str_replace_all(ecosar_dataframe_long$model, pattern = ' ', replacement = '_')
      
      # Remove number column
      ecosar_dataframe_long = ecosar_dataframe_long[,!colnames(ecosar_dataframe_long) %in% c('Number')]
      
      # We sometimes get multiple predictions from the same SMILES (ECOSAR thinks it might be different substances)
      # We fix this if ecosar_fix is in settings.
      if(ecosar_fix){
        
        ECOSAR_CAS = unique(ecosar_dataframe_long$original_CAS)
        
        # Loop over CAS to fix them
        for(i in 1:length(ECOSAR_CAS)){
          
          # Extract CAS
          current_CAS = ECOSAR_CAS[i]
          
          # Subset non-calculated values for aid CAS
          current_subset = ecosar_dataframe_long[ecosar_dataframe_long$original_CAS == current_CAS,]
          
          # Distinct at input SMILES, model and value fields
          temp_frame = 
            distinct_at(
              current_subset,
              .vars = c('original_SMILES', 
                        'model', 
                        'value'),
              .keep_all = T)
          
          # Check if we have duplicates
          if(any(duplicated(temp_frame$model))){
            
            #print(paste0('CAS with an issue: ', current_CAS))
            
            # If duplicates fixed by CAS identification (most of the time)
            if(current_CAS %in% as.cas(unique(temp_frame$ecosar_CAS))){
              
              #print(paste0('CAS matched, options: ', paste(unique(as.cas(unique(temp_frame$ecosar_CAS))), collapse = ', ')))
              
              # Remove old rows
              ecosar_dataframe_long = ecosar_dataframe_long[!ecosar_dataframe_long$original_CAS == current_CAS,]
              
              # Add new rows
              ecosar_dataframe_long = rbind(ecosar_dataframe_long, temp_frame[temp_frame$original_CAS == as.cas(temp_frame$ecosar_CAS),])
              
            } else{
              
              print(paste0('CAS mismatch, current_SMILES: ', unique(unique(temp_frame$original_SMILES))))
              print(paste0('Matching metals, SMILES options: ', paste(unique(unique(temp_frame$ecosar_SMILES)), collapse = ', ')))
              
              # Special cases where depreciated CAS are used
              
              # Find which metal (if any) is causing the issue and match it with ecosar SMILES
              current_metal = str_extract(unique(temp_frame$original_SMILES), pattern = '\\[[A-Za-z]{1,2}\\]')
              
              if(!is.na(current_metal)){
                
                print(paste0('Chosen SMILES: ', unique(temp_frame[grepl(temp_frame$ecosar_SMILES, pattern = current_metal, fixed = T), 'ecosar_SMILES'])))
                
                # Use this metal to get the right rows from temp_frame
                
                # Remove old rows
                ecosar_dataframe_long = ecosar_dataframe_long[!(ecosar_dataframe_long$original_CAS == current_CAS),]
                
                # Add new rows
                ecosar_dataframe_long = rbind(ecosar_dataframe_long, temp_frame[grepl(temp_frame$ecosar_SMILES, pattern = current_metal, fixed = T),])
                
                
              } else if(unique(temp_frame$original_SMILES) %in% unique(temp_frame$ecosar_SMILES)){
                
                # Then we check if we have matching SMILES instead
                print(paste0('Chosen SMILES: ', unique(temp_frame[temp_frame$ecosar_SMILES == unique(temp_frame$original_SMILES), 'ecosar_SMILES'])))
                
                # Remove old rows
                ecosar_dataframe_long = ecosar_dataframe_long[!(ecosar_dataframe_long$original_CAS == current_CAS),]
                
                # Add new rows
                ecosar_dataframe_long = rbind(ecosar_dataframe_long, temp_frame[temp_frame$ecosar_SMILES == unique(temp_frame$original_SMILES),])
                
              } else if(current_CAS == '72-20-8'){
                
                # Then we check if we have matching SMILES instead
                print(paste0('Manually fixing HEXADRIN/ENDRIN mismatch'))
                
                # Remove old rows
                ecosar_dataframe_long = ecosar_dataframe_long[!(ecosar_dataframe_long$original_CAS == current_CAS),]
                
                # Add new rows
                ecosar_dataframe_long = rbind(ecosar_dataframe_long, temp_frame[temp_frame$ecosar_CAS == '128109',])
                
              } else {
                print('Something else')
                break
                
              }

            }
            
          }
          
        }
        
      }
      
      
      
    }else if(format == 'wide'){
      
      ##### Repack ecosar into one line per substance (wide format) acute predictions
      
      ecosar_dataframe_wide = data.frame(original_CAS = unique(ecosar_dataframe$original_CAS[!is.na(ecosar_dataframe$original_CAS)]))
      
      # Add empty columns to ecosar_dataframe_wide based on organisms
      
      columns_to_add = c('ecosar_CAS', 'original_SMILES', 'ecosar_SMILES', 'InChIKey', 
                         'Daphnid_EC50_48h_1 [mg/l]', 'Daphnid_EC50_1_prediction_class', 'Daphnid_EC50_48h_2 [mg/l]', 'Daphnid_EC50_2_prediction_class', 'Daphnid_EC50_48h_3 [mg/l]', 'Daphnid_EC50_3_prediction_class', 'Daphnid_EC50_48h_4 [mg/l]', 'Daphnid_EC50_4_prediction_class', 'Daphnid_EC50_48h_5 [mg/l]', 'Daphnid_EC50_5_prediction_class', 'Daphnid_EC50_48h_6 [mg/l]', 'Daphnid_EC50_6_prediction_class', 
                         'Fish_EC50_96h_1 [mg/l]', 'Fish_EC50_1_prediction_class', 'Fish_EC50_96h_2 [mg/l]', 'Fish_EC50_2_prediction_class', 'Fish_EC50_96h_3 [mg/l]', 'Fish_EC50_3_prediction_class', 'Fish_EC50_96h_4 [mg/l]', 'Fish_EC50_4_prediction_class', 'Fish_EC50_96h_5 [mg/l]', 'Fish_EC50_5_prediction_class', 'Fish_EC50_96h_6 [mg/l]', 'Fish_EC50_6_prediction_class', 
                         'Fish14d_1 [mg/l]', 'Fish14d_1_prediction_class', 'Fish14d_2 [mg/l]', 'Fish14d_2_prediction_class',
                         'Green_Algae_EC50_96h_1 [mg/l]', 'Green_Algae_EC50_1_prediction_class', 'Green_Algae_EC50_96h_2 [mg/l]', 'Green_Algae_EC50_2_prediction_class', 'Green_Algae_EC50_96h_3 [mg/l]', 'Green_Algae_EC50_3_prediction_class', 'Green_Algae_EC50_96h_4 [mg/l]', 'Green_Algae_EC50_4_prediction_class', 'Green_Algae_EC50_96h_5 [mg/l]', 'Green_Algae_EC50_5_prediction_class', 'Green_Algae_EC50_96h_6 [mg/l]', 'Green_Algae_EC50_6_prediction_class',
                         'Mysid_EC50_96h_1 [mg/l]', 'Mysid_EC50_1_prediction_class', 'Mysid_EC50_96h_2 [mg/l]', 'Mysid_EC50_2_prediction_class', 'Mysid_EC50_96h_3 [mg/l]', 'Mysid_EC50_3_prediction_class', 'Mysid_EC50_96h_4 [mg/l]', 'Mysid_EC50_4_prediction_class', 'Mysid_EC50_96h_5 [mg/l]', 'Mysid_EC50_5_prediction_class', 'Mysid_EC50_96h_6 [mg/l]', 'Mysid_EC50_6_prediction_class',
                         'Earthworm_EC50_14d_1 [mg/l]', 'Earthworm_EC50_1_prediction_class', 'Earthworm_EC50_14d_2 [mg/l]', 'Earthworm_EC50_2_prediction_class', 'Earthworm_EC50_14d_3 [mg/l]', 'Earthworm_EC50_3_prediction_class', 'Earthworm_EC50_14d_4 [mg/l]', 'Earthworm_EC50_4_prediction_class', 'Earthworm_EC50_14d_5 [mg/l]', 'Earthworm_EC50_5_prediction_class', 'Earthworm_EC50_14d_6 [mg/l]', 'Earthworm_EC50_6_prediction_class',
                         'Lemna_gibba_EC50_7d_1 [mg/l]', 'Lemna_gibba_EC50_1_prediction_class', 'Lemna_gibba_EC50_7d_2 [mg/l]', 'Lemna_gibba_EC50_2_prediction_class', 'Lemna_gibba_EC50_7d_3 [mg/l]', 'Lemna_gibba_EC50_3_prediction_class', 'Lemna_gibba_EC50_7d_4 [mg/l]', 'Lemna_gibba_EC50_4_prediction_class', 'Lemna_gibba_EC50_7d_5 [mg/l]', 'Lemna_gibba_EC50_5_prediction_class', 'Lemna_gibba_EC50_7d_6 [mg/l]', 'Lemna_gibba_EC50_6_prediction_class')
      
      ecosar_dataframe_wide[columns_to_add] = NA
      
      # Make a vector of species used in ECOSAR (this can be updated along with the clumns above if ECOSAR adds more species)
      ecosar_organimsms = c('Daphnid', 'Fish', 'Green_Algae', 'Mysid', 'Earthworm', 'Lemna_gibba')
      
      # Loop over unique cas and add data for that cas to the wide frame (we will possibly add unmatchable entries later - the entries with NA original_cas)
      for(i in 1:nrow(ecosar_dataframe_wide)){
        
        # Get current CAS
        current_cas = ecosar_dataframe_wide$original_CAS[i]
        
        # Subset ECOSAR frame with CAS (ignore NAs since we will handle those later)
        current_subset = ecosar_dataframe[ecosar_dataframe$original_CAS == current_cas & !is.na(ecosar_dataframe$original_CAS) & !is.na(ecosar_dataframe$duration),]
        
        # First we check if we have multiple ecosar_CAS (ecosar has identified two different substances with the SMILES we predicted for)
        if(length(unique(current_subset$ecosar_CAS)) > 1){
          
          # If we have multiple ecosar_CAS we check if they have the same prediction (the expected scenario) (We are checking if the number of prediction values are a multiple of number of non-duplicated prediction values)
          if(nrow(current_subset) %% sum(!duplicated(current_subset$`Predicted value [mg/L]`)) == 0 ){
            
            # We discard the predictions for the CAS that does not match the original CAS
            temp_subset = current_subset[current_subset$original_CAS == current_subset$ecosar_CAS_fixed & !is.na(current_subset$ecosar_CAS_fixed),]
            
            # Make sure we didnt delete all points because of no matching CAS, if so we just pick the first CAS
            if(nrow(temp_subset) == 0){
              
              temp_subset = current_subset[current_subset$ecosar_CAS_fixed == unique(current_subset$ecosar_CAS_fixed)[1] & !is.na(current_subset$ecosar_CAS_fixed),]
              
            }
            
            current_subset = temp_subset
            
          }
          
          # Check that the issue is solved, otherwise report error
          if(length(unique(current_subset$ecosar_CAS)) > 1){
            print('Multiple different predictions for the same original CAS')
            print(paste0('Please check the original CAS ', current_cas, ' for mismatching!'))
          }
          
        }
        
        # Add all the identifiers
        ecosar_dataframe_wide[ecosar_dataframe_wide$original_CAS == current_cas ,c('ecosar_CAS', 'original_SMILES', 'ecosar_SMILES', 'InChIKey')] = current_subset[1, c("ecosar_CAS_fixed", "original_SMILES", "ecosar_SMILES", "InChIKey")]
        
        
        
        
        # Then we loop over organisms
        for (j in 1:length(ecosar_organimsms)){
          
          current_organism = ecosar_organimsms[j]
          
          
          # Check if the organism is not in the subset, if so skip
          if(!str_replace(current_organism, pattern = '_', replacement = ' ') %in% current_subset$model_organism){
            next
          }
          
          # Subset the subset with organism
          current_organism_subset = current_subset[current_subset$model_organism == str_replace(current_organism, pattern = '_', replacement = ' '),]
          
          
          # Fish do have a 14d model that does not fit the pattern, and needs special care in some cases
          special_case = FALSE
          if(current_organism == 'Fish' & '14d' %in% current_organism_subset$duration){
            special_case = TRUE
          }
          
          # Check if special case (with fish 14d - there are only two predictionclasses for this)
          if(!special_case){
            # If the subset does not contain 6 prediction classes, fill up with NAs
            if(nrow(current_organism_subset) < 6){
              current_organism_subset[seq(nrow(current_organism_subset)+1, 6, 1),] = NA
            } else if(nrow(current_organism_subset) > 6){
              
              print(paste0('More than 6 prediction classes for organism - ', current_organism, ' -  CAS number - ', current_cas, ' - Please return to code to add more columns for storing data'))
              
            }
            
            # Then we save the predictions and prediction classes into the proper columns (the paste0 adding underscore is to avoid filling in fish14d data if only fish 96h)
            ecosar_dataframe_wide[ecosar_dataframe_wide$original_CAS == current_cas ,grepl(colnames(ecosar_dataframe_wide), pattern = paste0(current_organism, '_EC50_'))] = c(rbind(current_organism_subset$`Predicted value [mg/L]`, current_organism_subset$`ECOSAR Class`))
            
          } else {
            
            
            # first we handle fish96h
            current_fish96h = current_organism_subset[current_organism_subset$duration == '96h',]
            
            if(nrow(current_fish96h) < 6){
              current_fish96h[seq(nrow(current_fish96h)+1, 6, 1),] = NA
            } else if(nrow(current_fish96h) > 6){
              
              print(paste0('More than 6 prediction classes for organism - ', current_organism, ' -  CAS number - ', current_cas, ' - Please return to code to add more columns for storing data'))
              
            }
            
            # Then we save the predictions and prediction classes into the proper columns (the paste0 adding underscore is to avoid filling in fish14d data if only fish 96h)
            ecosar_dataframe_wide[ecosar_dataframe_wide$original_CAS == current_cas ,grepl(colnames(ecosar_dataframe_wide), pattern = paste0(current_organism, '_EC50_'))] = c(rbind(current_fish96h$`Predicted value [mg/L]`, current_fish96h$`ECOSAR Class`))
            
            
            
            # then handle fish14d
            current_fish14d = current_organism_subset[current_organism_subset$duration == '14d',]
            
            if(nrow(current_fish14d) < 2){
              current_fish14d[seq(nrow(current_fish14d)+1, 2, 1),] = NA
            } else if(nrow(current_fish14d) > 2){
              
              print(paste0('More than 2 prediction classes (which should be impossible) for organism - fish_14d - for CAS number - ', current_cas, ' - Please return to code to add more columns for storing data'))
              
            }
            
            # Then we save the predictions and prediction classes into the proper columns 
            ecosar_dataframe_wide[ecosar_dataframe_wide$original_CAS == current_cas ,grepl(colnames(ecosar_dataframe_wide), pattern = 'Fish14d')] = c(rbind(current_fish14d$`Predicted value [mg/L]`, current_fish14d$`ECOSAR Class`))
            
          }
        }
        
      }
      
      
      
      
      ##### Repack ecosar into one line per substance (wide format) chronic predictions
      
      # Add empty columns to ecosar_dataframe_wide based on organisms
      
      columns_to_add = c('Daphnid_ChV_1 [mg/l]', 'Daphnid_ChV_1_prediction_class', 'Daphnid_ChV_2 [mg/l]', 'Daphnid_ChV_2_prediction_class', 'Daphnid_ChV_3 [mg/l]', 'Daphnid_ChV_3_prediction_class', 'Daphnid_ChV_4 [mg/l]', 'Daphnid_ChV_4_prediction_class', 'Daphnid_ChV_5 [mg/l]', 'Daphnid_ChV_5_prediction_class', 'Daphnid_ChV_6 [mg/l]', 'Daphnid_ChV_6_prediction_class', 
                         'Fish_ChV_1 [mg/l]', 'Fish_ChV_1_prediction_class', 'Fish_ChV_2 [mg/l]', 'Fish_ChV_2_prediction_class', 'Fish_ChV_3 [mg/l]', 'Fish_ChV_3_prediction_class', 'Fish_ChV_4 [mg/l]', 'Fish_ChV_4_prediction_class', 'Fish_ChV_5 [mg/l]', 'Fish_ChV_5_prediction_class', 'Fish_ChV_6 [mg/l]', 'Fish_ChV_6_prediction_class', 
                         'Green_Algae_ChV_1 [mg/l]', 'Green_Algae_ChV_1_prediction_class', 'Green_Algae_ChV_2 [mg/l]', 'Green_Algae_ChV_2_prediction_class', 'Green_Algae_ChV_3 [mg/l]', 'Green_Algae_ChV_3_prediction_class', 'Green_Algae_ChV_4 [mg/l]', 'Green_Algae_ChV_4_prediction_class', 'Green_Algae_ChV_5 [mg/l]', 'Green_Algae_ChV_5_prediction_class', 'Green_Algae_ChV_6 [mg/l]', 'Green_Algae_ChV_6_prediction_class',
                         'Mysid_ChV_1 [mg/l]', 'Mysid_ChV_1_prediction_class', 'Mysid_ChV_2 [mg/l]', 'Mysid_ChV_2_prediction_class', 'Mysid_ChV_3 [mg/l]', 'Mysid_ChV_3_prediction_class', 'Mysid_ChV_4 [mg/l]', 'Mysid_ChV_4_prediction_class', 'Mysid_ChV_5 [mg/l]', 'Mysid_ChV_5_prediction_class', 'Mysid_ChV_6 [mg/l]', 'Mysid_ChV_6_prediction_class',
                         'Earthworm_ChV_1 [mg/l]', 'Earthworm_ChV_1_prediction_class', 'Earthworm_ChV_2 [mg/l]', 'Earthworm_ChV_2_prediction_class', 'Earthworm_ChV_3 [mg/l]', 'Earthworm_ChV_3_prediction_class', 'Earthworm_ChV_4 [mg/l]', 'Earthworm_ChV_4_prediction_class', 'Earthworm_ChV_5 [mg/l]', 'Earthworm_ChV_5_prediction_class', 'Earthworm_ChV_6 [mg/l]', 'Earthworm_ChV_6_prediction_class',
                         'Lemna_gibba_ChV_1 [mg/l]', 'Lemna_gibba_ChV_1_prediction_class', 'Lemna_gibba_ChV_2 [mg/l]', 'Lemna_gibba_ChV_2_prediction_class', 'Lemna_gibba_ChV_3 [mg/l]', 'Lemna_gibba_ChV_3_prediction_class', 'Lemna_gibba_ChV_4 [mg/l]', 'Lemna_gibba_ChV_4_prediction_class', 'Lemna_gibba_ChV_5 [mg/l]', 'Lemna_gibba_ChV_5_prediction_class', 'Lemna_gibba_ChV_6 [mg/l]', 'Lemna_gibba_ChV_6_prediction_class')
      
      ecosar_dataframe_wide[columns_to_add] = NA
      
      # Make a vector of species used in ECOSAR (this can be updated along with the clumns above if ECOSAR adds more species)
      ecosar_organimsms = c('Daphnid', 'Fish', 'Green_Algae', 'Mysid', 'Earthworm', 'Lemna_gibba')
      
      # Loop over unique cas and add data for that cas to the wide frame (we will possibly add unmatchable entries later - the entries with NA original_cas)
      for(i in 1:nrow(ecosar_dataframe_wide)){
        
        # Get current CAS
        current_cas = ecosar_dataframe_wide$original_CAS[i]
        
        # Subset ECOSAR frame with CAS (ignore NAs since we will handle those later)
        current_subset = ecosar_dataframe[ecosar_dataframe$original_CAS == current_cas & !is.na(ecosar_dataframe$original_CAS) & is.na(ecosar_dataframe$duration),]
        
        # First we check if we have multiple ecosar_CAS (ecosar has identified two different substances with the SMILES we predicted for)
        if(length(unique(current_subset$ecosar_CAS)) > 1){
          
          # If we have multiple ecosar_CAS we check if they have the same prediction (the expected scenario) (We are checking if the number of prediction values are a multiple of number of non-duplicated prediction values)
          if(nrow(current_subset) %% sum(!duplicated(current_subset$`Predicted value [mg/L]`)) == 0 ){
            
            # We discard the predictions for the CAS that does not match the original CAS
            temp_subset = current_subset[current_subset$original_CAS == current_subset$ecosar_CAS_fixed & !is.na(current_subset$ecosar_CAS_fixed),]
            
            # Make sure we didnt delete all points because of no matching CAS, if so we just pick the first CAS
            if(nrow(temp_subset) == 0){
              
              temp_subset = current_subset[current_subset$ecosar_CAS_fixed == unique(current_subset$ecosar_CAS_fixed)[1] & !is.na(current_subset$ecosar_CAS_fixed),]
              
            }
            
            current_subset = temp_subset
            
          }
          
          # Check that the issue is solved, otherwise report error
          if(length(unique(current_subset$ecosar_CAS)) > 1){
            print('Multiple different predictions for the same original CAS')
            print(paste0('Please check the original CAS ', current_cas, ' for mismatching!'))
          }
          
        }
        
        # Add all the identifiers
        ecosar_dataframe_wide[ecosar_dataframe_wide$original_CAS == current_cas ,c('ecosar_CAS', 'original_SMILES', 'ecosar_SMILES', 'InChIKey')] = current_subset[1, c("ecosar_CAS_fixed", "original_SMILES", "ecosar_SMILES", "InChIKey")]
        
        
        
        
        # Then we loop over organisms
        for (j in 1:length(ecosar_organimsms)){
          
          current_organism = ecosar_organimsms[j]
          
          
          # Check if the organism is not in the subset, if so skip
          if(!str_replace(current_organism, pattern = '_', replacement = ' ') %in% current_subset$model_organism){
            next
          }
          
          # Subset the subset with organism
          current_organism_subset = current_subset[current_subset$model_organism == str_replace(current_organism, pattern = '_', replacement = ' '),]
          
          
          # Fish do have a 14d model that does not fit the pattern, and needs special care in some cases
          special_case = FALSE
          if(current_organism == 'Fish' & '14d' %in% current_organism_subset$duration){
            special_case = TRUE
          }
          
          # Check if special case (with fish 14d - there are only two predictionclasses for this)
          if(!special_case){
            # If the subset does not contain 6 prediction classes, fill up with NAs
            if(nrow(current_organism_subset) < 6){
              current_organism_subset[seq(nrow(current_organism_subset)+1, 6, 1),] = NA
            } else if(nrow(current_organism_subset) > 6){
              
              print(paste0('More than 6 prediction classes for organism - ', current_organism, ' -  CAS number - ', current_cas, ' - Please return to code to add more columns for storing data'))
              
            }
            
            # Then we save the predictions and prediction classes into the proper columns (the paste0 adding underscore is to avoid filling in fish14d data if only fish 96h)
            ecosar_dataframe_wide[ecosar_dataframe_wide$original_CAS == current_cas ,grepl(colnames(ecosar_dataframe_wide), pattern = paste0(current_organism, '_ChV_'))] = c(rbind(current_organism_subset$`Predicted value [mg/L]`, current_organism_subset$`ECOSAR Class`))
          } else {
            
            
            # first we handle fish96h
            current_fish96h = current_organism_subset[current_organism_subset$duration == '96h',]
            
            if(nrow(current_fish96h) < 6){
              current_fish96h[seq(nrow(current_fish96h)+1, 6, 1),] = NA
            } else if(nrow(current_fish96h) > 6){
              
              print(paste0('More than 6 prediction classes for organism - ', current_organism, ' -  CAS number - ', current_cas, ' - Please return to code to add more columns for storing data'))
              
            }
            
            # Then we save the predictions and prediction classes into the proper columns (the paste0 adding underscore is to avoid filling in fish14d data if only fish 96h)
            ecosar_dataframe_wide[ecosar_dataframe_wide$original_CAS == current_cas ,grepl(colnames(ecosar_dataframe_wide), pattern = paste0(current_organism, '_ChV_'))] = c(rbind(current_fish96h$`Predicted value [mg/L]`, current_fish96h$`ECOSAR Class`))
            
            
            
            # then handle fish14d
            current_fish14d = current_organism_subset[current_organism_subset$duration == '14d',]
            
            if(nrow(current_fish14d) < 2){
              current_fish14d[seq(nrow(current_fish14d)+1, 2, 1),] = NA
            } else if(nrow(current_fish14d) > 2){
              
              print(paste0('More than 2 prediction classes (which should be impossible) for organism - fish_14d - for CAS number - ', current_cas, ' - Please return to code to add more columns for storing data'))
              
            }
            
            # Then we save the predictions and prediction classes into the proper columns 
            ecosar_dataframe_wide[ecosar_dataframe_wide$original_CAS == current_cas ,grepl(colnames(ecosar_dataframe_wide), pattern = 'Fish14d')] = c(rbind(current_fish14d$`Predicted value [mg/L]`, current_fish14d$`ECOSAR Class`))
            
          }
        }
        
      }
      
      # Then we filter out all columns with all NA values
      test = Filter(function(x)!all(is.na(x)), ecosar_dataframe_wide)
      
      # Add platform name to all non-identifier columns
      colnames(ecosar_dataframe_wide)[6:length(colnames(ecosar_dataframe_wide))] = str_replace(colnames(ecosar_dataframe_wide)[6:length(colnames(ecosar_dataframe_wide))], pattern = '(?=^)', replacement = 'ECOSAR_')
      
      
      
    }
  }
  
  ################################################################################
  #             6. Importing and standardizing T.E.S.T. model data               #
  ################################################################################
  
  if(run_test){
    
    print('Processing T.E.S.T. output')
    
    # Read all csv-files from the TEST output folder
    test_file_list <- list.files(path = testwd2, all.files = F, pattern= '.*\\.csv$', ignore.case = T)
    
    test_dataframe <- NULL
    
    for(i in 1:length(test_file_list)){
      
      filename <- test_file_list[i]
      
      current_organism <- str_extract(filename, pattern = '(?<=Batch_)[A-Za-z]{0,100}_[A-Za-z]{0,100}(?=_)')
      
      current_endpoint <- str_extract(filename, pattern = '(?<=Batch_[A-Za-z]{0,100}_[A-Za-z]{0,100}_)[A-Z0-9]{0,6}(?=_)')
      
      current_duration <- str_extract(filename, pattern = '(?<=Batch_[A-Za-z]{0,100}_[A-Za-z]{0,100}_[A-Z0-9]{0,6}_\\()[0-9]{2}_hr(?=\\)_)')
      
      current_frame <- as.data.frame(read.csv(file = paste0(testwd2,'/', filename), header = T, stringsAsFactors = F, dec = ",", sep = ",", quote = '"', na.strings = 'N/A'))
      
      
      
      
      current_frame$model_organism <- str_replace_all(current_organism, pattern = '_', replacement = ' ')
      
      current_frame$model_endpoint <- current_endpoint
      
      current_frame$'model_duration (h)' <- str_split(current_duration, pattern = '_')[[1]][1]
      
      ## Merge to ecosar frame
      if(!exists('test_dataframe')){
        
        test_dataframe <- current_frame
        
      } else {
        
        # Save and deduplicate (Index may be ignored if the original entry contained)
        test_dataframe <- as.data.frame(distinct(rbind(test_dataframe, current_frame), across(-'Index'),  .keep_all = T))
        
      }
      
    }
    
    
    
    
    # Add T.E.S.T. in the QSAR_tool field of the frame
    test_dataframe$QSAR_tool <- 'T.E.S.T.'
    
    # Add duration column
    test_dataframe$duration = ifelse(!is.na(test_dataframe$`model_duration (h)`),
                                     paste0(test_dataframe$`model_duration (h)`, 'h'),
                                     test_dataframe$`model_duration (h)`)
    
    # Rename columns to match standard
    test_dataframe <- 
      test_dataframe %>%
      rename(test_CAS = ID) %>%
      rename(`Predicted value [mg/L]` = `Pred_Value._mg.L`) %>%
      rename(test_SMILES = SmilesRan) %>%
      rename(original_SMILES = Query)
    
    # Add inchi and original CAS from orignial data
    test_dataframe <- merge(test_dataframe, identifiers[,c("original_SMILES", "original_CAS", "InChIKey")], by = 'original_SMILES', all.x = T)
    
    # Move columns around to make more sense
    test_dataframe <- 
      test_dataframe %>%
      relocate(original_CAS) %>%
      relocate(test_CAS, .after = original_CAS) %>%
      relocate(original_SMILES, .after =test_CAS) %>%
      relocate(test_SMILES, .after =original_SMILES) %>%
      relocate(InChIKey, .after =test_SMILES)
    
    
    if(format == 'long'){
      
      test_dataframe_long = test_dataframe
      
      # Remove some columns, and fix some
      test_dataframe_long = test_dataframe_long[,!grepl(colnames(test_dataframe_long), pattern = 'Log10')]
      
      test_dataframe_long = test_dataframe_long %>%
        rename(prediction = 'Predicted value [mg/L]') %>%
        rename(error_warning = 'Error') %>%
        rename(EXPERIMENTAL = 'Exp_Value._mg.L')
  
      # Fix model column
      test_dataframe_long$model_organism = ifelse(test_dataframe_long$model_organism == 'Fathead minnow', 'Fish', test_dataframe_long$model_organism)
      test_dataframe_long$model_organism = str_replace_all(test_dataframe_long$model_organism, pattern = ' ', replacement = '_')
      test_dataframe_long$model = apply(test_dataframe_long,MARGIN = 1, FUN = function(x){paste('TEST', x[10],x[11],x[12],sep = '_')})
      test_dataframe_long$model = str_replace_all(test_dataframe_long$model, pattern = '(?=$)', replacement = 'h')
      
      # Gather the experimentl values with the predicted
      test_dataframe_long = gather(test_dataframe_long, key = reliability, value = value, c('EXPERIMENTAL', 'prediction'))
      
      # Remove rows with NA or empty values
      test_dataframe_long = test_dataframe_long[!is.na(test_dataframe_long$value) & test_dataframe_long$value != '',]
      
      # Remove index column
      test_dataframe_long = test_dataframe_long[, !colnames(test_dataframe_long) %in% c('Index')]
      
      
    } else if(format == 'wide'){
    
      ### Make test_dataframe on wide format
      
      test_dataframe_wide = data.frame(original_CAS = unique(test_dataframe$original_CAS[!is.na(test_dataframe$original_CAS)]))
      
      # Add columns for identifiers and predictions
      
      columns_to_add = c('test_CAS', 'original_SMILES', 'test_SMILES', 'InChIKey', 
                         'Daphnia_magna_48h [mg/l]',
                         'Fish_96h [mg/l]')
      
      test_dataframe_wide[columns_to_add] = NA
      
      for(i in 1:nrow(test_dataframe_wide)){
        
        current_cas = test_dataframe_wide$original_CAS[i]
        
        current_subset = test_dataframe[test_dataframe$original_CAS == current_cas & !is.na(test_dataframe$original_CAS),]
        
        test_dataframe_wide[test_dataframe_wide$original_CAS == current_cas, c("test_CAS", "original_SMILES", "test_SMILES", "InChIKey")] = current_subset[1, c("test_CAS", "original_SMILES", "test_SMILES", "InChIKey")]
        
        test_dataframe_wide[test_dataframe_wide$original_CAS == current_cas, c("Daphnia_magna_48h [mg/l]", "Fish_96h [mg/l]")] = 
          c(current_subset[current_subset$model_organism == 'Daphnia magna',"Predicted value [mg/L]"], current_subset[current_subset$model_organism == 'Fathead minnow',"Predicted value [mg/L]"])
        
      }
      
      
      # Add paltform name to all non-identifier columns
      colnames(test_dataframe_wide)[6:length(colnames(test_dataframe_wide))] = str_replace(colnames(test_dataframe_wide)[6:length(colnames(test_dataframe_wide))], pattern = '(?=^)', replacement = 'TEST_')
      
      
    }
    
  }
  
  ################################################################################
  #                     7. Merging model data                                    #
  ################################################################################
  
  print('Finalizing function output')
  
  if(format == 'wide'){
    
    # If all three QSARs were run
    if(all(run_vega, run_ecosar, run_test)){
      # Merge Vega with ecosar, skipping the identifiers that overlap. Doing the same with test.
      
      merged_output_wide = merge(vega_dataframe_wide, ecosar_dataframe_wide[,colnames(ecosar_dataframe_wide)[!colnames(ecosar_dataframe_wide) %in% c('original_SMILES', "InChIKey")]], by = 'original_CAS', all = T)
      merged_output_wide = merge(merged_output_wide, test_dataframe_wide[colnames(test_dataframe_wide)[!colnames(test_dataframe_wide) %in% c('original_SMILES', 'InChIKey')]], all = T)
      
    }
    
    # If two QSARs were run
    if(sum(c(run_vega, run_ecosar, run_test)) == 2){
      
      if(run_vega){
        df1 = vega_dataframe_wide
      } else {
        df1 = test_dataframe_wide
      }
      
      if(run_ecosar){
        df2 = ecosar_dataframe_wide 
      } else {
        df2 = test_dataframe_wide
      }
      
      merged_output_wide = merge(df1, df2[,colnames(df2)[!colnames(df2) %in% c('original_SMILES', "InChIKey")]], by = 'original_CAS', all = T)
      
    }
    
    # If only one QSAR was run
    if(sum(c(run_vega, run_ecosar, run_test)) == 1){
      
      if(run_vega){
        df1 = vega_dataframe_wide
      } else if(run_ecosar){
        df1 = ecosar_dataframe_wide
      } else {
        df1 = test_dataframe_wide
      }
        
      merged_output_wide = df1
      
      }
    
    # Relocate columns to get all identifiers at the beginning
    if(!run_vega){
      merged_output_wide$vega_SMILES = NA
    }
    if(!run_ecosar){
      merged_output_wide$ecosar_CAS = NA
      merged_output_wide$ecosar_SMILES = NA
    }
    if(!run_test){
      merged_output_wide$test_CAS = NA
      merged_output_wide$test_SMILES = NA
    }
    
    merged_output_wide = 
      merged_output_wide %>%
      relocate(original_CAS) %>%
      relocate(ecosar_CAS, .after = original_CAS) %>%
      relocate(test_CAS, .after = ecosar_CAS) %>%
      relocate(original_SMILES, .after = test_CAS) %>%
      relocate(vega_SMILES, .after = original_SMILES) %>%
      relocate(ecosar_SMILES, .after = vega_SMILES) %>%
      relocate(test_SMILES, .after = ecosar_SMILES) %>%
      relocate(InChIKey, .after = test_SMILES)
    
    
    # Return the full output with added mean calculations
    return(merged_output_wide)
    
    
  } else if (format == 'long'){
    
    # If all QSARs are run
    if(all(run_vega, run_ecosar, run_test)){
      
      # Add empty columns so we can rbind dataframes, 3 frames version
      vega_dataframe_long[,colnames(ecosar_dataframe_long)[!colnames(ecosar_dataframe_long) %in% colnames(vega_dataframe_long)]] <- NA
      vega_dataframe_long[,colnames(test_dataframe_long)[!colnames(test_dataframe_long) %in% colnames(vega_dataframe_long)]] <- NA
      ecosar_dataframe_long[,colnames(vega_dataframe_long)[!colnames(vega_dataframe_long) %in% colnames(ecosar_dataframe_long)]] <- NA
      test_dataframe_long[,colnames(vega_dataframe_long)[!colnames(vega_dataframe_long) %in% colnames(test_dataframe_long)]] <- NA
      
      merged_output_long = rbind(vega_dataframe_long, ecosar_dataframe_long)
      merged_output_long = rbind(merged_output_long, test_dataframe_long)
      
    }
    
    # If two QSARs were run
    if(sum(c(run_vega, run_ecosar, run_test)) == 2){
      
      if(run_vega){
        df1 = vega_dataframe_long
      } else {
        df1 = test_dataframe_long
      }
      
      if(run_ecosar){
        df2 = ecosar_dataframe_long 
      } else {
        df2 = test_dataframe_long
      }
      
      df1[,colnames(df2)[!colnames(df2) %in% colnames(df1)]] <- NA
      df2[,colnames(df1)[!colnames(df1) %in% colnames(df2)]] <- NA
      
      merged_output_long = rbind(df1, df2)
      
    }
    
    # If only one QSAR was run
    if(sum(c(run_vega, run_ecosar, run_test)) == 1){
      
      if(run_vega){
        df1 = vega_dataframe_long
      } else if(run_ecosar){
        df1 = ecosar_dataframe_long
      } else {
        df1 = test_dataframe_long
      }
      
      merged_output_long = df1
      
    }
    
    
    if(!run_vega){
      merged_output_long$vega_SMILES = NA
    }
    if(!run_ecosar){
      merged_output_long$ecosar_CAS = NA
      merged_output_long$ecosar_SMILES = NA
    }
    if(!run_test){
      merged_output_long$test_CAS = NA
      merged_output_long$test_SMILES = NA
    }
    
    # Move some columns around to improve readability
    merged_output_long = merged_output_long %>%
      relocate(ecosar_CAS, .after = original_CAS) %>%
      relocate(test_CAS, .after = ecosar_CAS) %>%
      relocate(original_SMILES, .after = test_CAS) %>%
      relocate(ecosar_SMILES, .after = original_SMILES) %>%
      relocate(vega_SMILES, .after = ecosar_SMILES) %>%
      relocate(test_SMILES, .after = vega_SMILES) 
    
    # We change LC50 to EC50
    merged_output_long$model_endpoint = ifelse(merged_output_long$model_endpoint == 'LC50', 'EC50', merged_output_long$model_endpoint)
    
    # Fix ECOSAR ChV -> NOEC - according to ECHA REACH Guidance R.10, the ChV can be divided by √2 to derive a NOEC
    merged_output_long$value = ifelse(merged_output_long$model_endpoint == 'ChV', as.numeric(merged_output_long$value) / sqrt(2), merged_output_long$value)
    merged_output_long$model_endpoint = ifelse(merged_output_long$model_endpoint == 'ChV', 'NOEC', merged_output_long$model_endpoint)
    
    return(merged_output_long)
    
  }
  
  
  
  

}