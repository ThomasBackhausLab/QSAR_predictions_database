################################################################################
#    A function for building ECOTOX from ASCII-files                                                  #
################################################################################
#
# Original author: Thomas Backhaus, Francis Spilsbury
#
# Adapted by: Patrik Svedberg
#
# Contact email: patrik.svedberg@bioenv.gu.se
# (if the above doesnt work (i.e. years down the line) try p.a.svedberg@gmail.com)
#
# Based on the script: AquireConverter_v.8.R
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

# Adapted from AquireConverter_v.8.R into function

################################################################################
# 1. Initiating function and parsing arguments
################################################################################

function(ECOTOX_ascii_path,
         settings = NULL){
  
  
  # packages, imports, functions
  
  library(readxl)
  library(writexl)
  library(tidyverse)
  library(stringr)
  library(data.table)
  
  # Set the current workdirectories
  inwd <- ECOTOX_ascii_path
  
  
  # This script workis with wd, so we save the old wd, set a new one and then return to the original at the end
  original_wd = getwd()
  
  setwd(inwd)
  
  on.exit(setwd(original_wd))
  
  save_lookup = F
  if(is.null(settings)){
    
    print('No additional settings provided')
    
  } else if(!is.null(settings) & (settings == 'save' | 'save' %in% settings)){
    
    save_lookup = T
    
    outwd <- "C:/Git/QSAR_predictions_database/Output"
  }
  
  ################################################################################
  # 2. function body
  ################################################################################
  
  
  # timestamp of the Ecotox database
  database_timestamp <- str_extract(ECOTOX_ascii_path, pattern = "[0-9]{2}_15_[0-9]{4}")
  
  
  # Read tables from the ascii
  SpeciesLookup =  read.table("validation/species.txt", 
                              comment.char = "", 
                              quote = "", 
                              sep="|", 
                              stringsAsFactors=FALSE ,
                              header = TRUE,
                              strip.white = TRUE, 
                              check.names = FALSE,
                              fill = TRUE)
  
  ResultsLookup =  read.table("results.txt", 
                              comment.char = "", 
                              quote = "",
                              sep="|", 
                              stringsAsFactors=FALSE ,
                              header = TRUE,
                              strip.white = TRUE, 
                              check.names = FALSE,
                              fill = TRUE)
  
  TestsLookup =  read.table("tests.txt", 
                            comment.char = "", 
                            quote = "",
                            sep="|", 
                            stringsAsFactors=FALSE ,
                            header = TRUE,
                            strip.white = TRUE, 
                            check.names = FALSE,
                            fill = TRUE,
                            colClasses = c("test_cas"="character")) # forces to import "test_cas" as character (not int)
  
  ReferencesLookup =  read.table("validation/references.txt", 
                                 comment.char = "", 
                                 quote = "", 
                                 sep="|", 
                                 stringsAsFactors=FALSE ,
                                 header = TRUE,
                                 strip.white = TRUE, 
                                 check.names = FALSE,
                                 fill = TRUE,
                                 na.strings = c("","NA"))
  
  EffectLookup =  read.table("validation/effect_codes.txt", 
                             comment.char = "", 
                             quote = "",
                             sep="|", 
                             stringsAsFactors=FALSE ,
                             header = TRUE,
                             strip.white = TRUE, 
                             check.names = FALSE,
                             fill = TRUE)
  
  MeasurementLookup =  read.table("validation/measurement_codes.txt", 
                                  comment.char = "", 
                                  quote = "", 
                                  sep="|", 
                                  stringsAsFactors=FALSE ,
                                  header = TRUE,
                                  strip.white = TRUE, 
                                  check.names = FALSE,
                                  fill = TRUE)
  
  # Timeunits werden an 4 Stellen ben?tigt
  TimeUnitLookup_1 = read.table("validation/duration_unit_codes.txt", 
                                comment.char = "", 
                                quote = "",
                                sep="|", 
                                stringsAsFactors=FALSE ,
                                header = TRUE,
                                strip.white = TRUE, 
                                check.names = FALSE,
                                fill = TRUE)
  
  TimeUnitLookup_2 = read.table("validation/duration_unit_codes.txt", 
                                comment.char = "", 
                                quote = "",
                                sep="|", 
                                stringsAsFactors=FALSE ,
                                header = TRUE,
                                strip.white = TRUE, 
                                check.names = FALSE,
                                fill = TRUE)
  
  TimeUnitLookup_3 = read.table("validation/duration_unit_codes.txt", 
                                comment.char = "", 
                                quote = "",
                                sep="|", 
                                stringsAsFactors=FALSE ,
                                header = TRUE,
                                strip.white = TRUE, 
                                check.names = FALSE,
                                fill = TRUE)
  
  TimeUnitLookup_4 = read.table("validation/duration_unit_codes.txt", 
                                comment.char = "", 
                                quote = "",
                                sep="|", 
                                stringsAsFactors=FALSE ,
                                header = TRUE,
                                strip.white = TRUE, 
                                check.names = FALSE,
                                fill = TRUE)    
  
  
  CASLookup = read.table("validation/chemicals.txt", 
                         comment.char = "", 
                         quote = "",
                         sep="|", 
                         stringsAsFactors=FALSE ,
                         header = TRUE,
                         strip.white = TRUE, 
                         check.names = FALSE,
                         fill = TRUE,
                         colClasses = c("character","character","character")) # forces to import everything as character
  
  LifeStageLookup = read.table("validation/lifestage_codes.txt", 
                               comment.char = "", 
                               quote = "", 
                               sep="|", 
                               stringsAsFactors=FALSE ,
                               header = TRUE,
                               strip.white = TRUE, 
                               check.names = FALSE,
                               fill = TRUE)
  
  EndpointLookup = read.table("validation/endpoint_codes.txt", 
                              comment.char = "", 
                              quote = "", 
                              sep="|", 
                              stringsAsFactors=FALSE ,
                              header = TRUE,
                              strip.white = TRUE, 
                              check.names = FALSE,
                              fill = TRUE)  
  
  ########
  # Convert CAS number to Form with Hyphens
  ########
  
  # in CASLookup
  CASLookup$dummy<-str_sub(CASLookup$cas_number,nchar(CASLookup$cas_number),nchar(CASLookup$cas_number))
  CASLookup$dummy2<-str_sub(CASLookup$cas_number,nchar(CASLookup$cas_number)-2,nchar(CASLookup$cas_number)-1)
  CASLookup$dummy3<-str_sub(CASLookup$cas_number,1,nchar(CASLookup$cas_number)-3)
  
  CASLookup$original_cas_number<-CASLookup$cas_number
  CASLookup$cas_number<-paste0(CASLookup$dummy3,"-",CASLookup$dummy2,"-",CASLookup$dummy)
  
  CASLookup$dummy<-NULL    
  CASLookup$dummy2<-NULL
  CASLookup$dummy3<-NULL
  CASLookup$dummy4<-NULL
  
  # in TestsLookup
  TestsLookup$dummy<-str_sub(TestsLookup$test_cas,nchar(TestsLookup$test_cas),nchar(TestsLookup$test_cas))
  TestsLookup$dummy2<-str_sub(TestsLookup$test_cas,nchar(TestsLookup$test_cas)-2,nchar(TestsLookup$test_cas)-1)
  TestsLookup$dummy3<-str_sub(TestsLookup$test_cas,1,nchar(TestsLookup$test_cas)-3)
  
  TestsLookup$original_test_cas<-TestsLookup$test_cas
  TestsLookup$test_cas<-paste0(TestsLookup$dummy3,"-",TestsLookup$dummy2,"-",TestsLookup$dummy)
  
  TestsLookup$dummy<-NULL    
  TestsLookup$dummy2<-NULL
  TestsLookup$dummy3<-NULL
  TestsLookup$dummy4<-NULL
  
  ########
  # Rename Variables that occur under the same name in various importfiles
  # although they mean different things
  ########
  setnames(ResultsLookup, "additional_comments", "results_additional_comments")
  setnames(TestsLookup, "additional_comments", "tests_additional_comments")
  setnames(ResultsLookup, "created_date", "results_created_date")
  setnames(TestsLookup, "created_date", "tests_created_date")
  setnames(ResultsLookup, "modified_date", "results_modified_date")
  setnames(TestsLookup, "modified_date", "tests_modified_date")
  
  setnames(SpeciesLookup, "ecotox_group", "species_group")
  setnames(CASLookup, "ecotox_group", "chemical_group")
  
  setnames(EffectLookup, "description", "effect_description")
  setnames(EndpointLookup, "description", "endpoint_description")
  setnames(MeasurementLookup, "description", "measurement_description")
  
  setnames(TimeUnitLookup_1, "description", "duration_unit_description_organism")
  setnames(TimeUnitLookup_2, "description", "duration_unit_description_study")
  setnames(TimeUnitLookup_3, "description", "duration_unit_description_exposure")
  setnames(TimeUnitLookup_4, "description", "duration_unit_description_obs")
  
  
  ########
  # Rename Variables that have different names but are in fact keys needed for
  # merging
  ########
  setnames(TestsLookup, "test_cas", "cas_number")
  setnames(EffectLookup, "code", "effect_code")
  setnames(EndpointLookup, "code", "endpoint_code")
  setnames(MeasurementLookup, "code", "measurement_code")
  setnames(TimeUnitLookup_1, "code", "organism_age_unit")
  setnames(TimeUnitLookup_2, "code", "study_duration_unit")
  setnames(TimeUnitLookup_3, "code", "exposure_duration_unit")
  setnames(TimeUnitLookup_4, "code", "obs_duration_unit")
  
  Backup1<-ResultsLookup
  # ResultsLookup<-Backup1
  
  #######
  # Generate keys from existing fields that are polluted 
  # with / ~ *
  # keep the original effect, measurement and endpoint code
  # will be used later.
  #######
  ResultsLookup$effect_code<-gsub("/","",ResultsLookup$effect)
  ResultsLookup$effect_code<-gsub("~","",ResultsLookup$effect_code)
  ResultsLookup$effect_code<-gsub("\\*","",ResultsLookup$effect_code)
  
  ResultsLookup$measurement_code<-gsub("/","",ResultsLookup$measurement)
  ResultsLookup$measurement_code<-gsub("~","",ResultsLookup$measurement_code)
  ResultsLookup$measurement_code<-gsub("\\*","",ResultsLookup$measurement_code)
  
  ResultsLookup$endpoint_code<-gsub("/","",ResultsLookup$endpoint)
  ResultsLookup$endpoint_code<-gsub("~","",ResultsLookup$endpoint_code)
  ResultsLookup$endpoint_code<-gsub("\\*","",ResultsLookup$endpoint_code)  
  
  
  ########
  # Mergers
  ########
  
  # merge with results table
  ECOTOX <- merge(ResultsLookup,TestsLookup, by="test_id")
  
  # merge with References
  ECOTOX <- merge(ECOTOX, ReferencesLookup, by="reference_number")
  
  # merge with Chemical Names
  ECOTOX <- merge(ECOTOX, CASLookup, by= "cas_number")
  
  # merge with Species Names
  ECOTOX <- merge(ECOTOX, SpeciesLookup, by="species_number")
  
  # merge with effect codes
  ECOTOX <- merge(ECOTOX, EffectLookup, by="effect_code")
  
  Backup2<-ECOTOX  
  
  # merge with endpoint codes
  ECOTOX <- merge(ECOTOX,EndpointLookup,by="endpoint_code")
  
  # by measurement
  ECOTOX<-merge(ECOTOX,MeasurementLookup, by="measurement_code", all.x=TRUE)
  
  # by organism_age_unit
  ECOTOX<-merge(ECOTOX,TimeUnitLookup_1, by="organism_age_unit", all.x=TRUE)
  
  # by study_duration_unit
  ECOTOX<-merge(ECOTOX,TimeUnitLookup_2, by="study_duration_unit", all.x=TRUE)
  
  # by exposure_duration_unit
  ECOTOX<-merge(ECOTOX,TimeUnitLookup_3, by="exposure_duration_unit", all.x=TRUE)
  
  # by obs_duration_unit
  ECOTOX<-merge(ECOTOX,TimeUnitLookup_4, by="obs_duration_unit", all.x=TRUE)
  
  if(save_lookup){
    
    #save dataframe
    save(ECOTOX,file=paste0(outwd,"/ECOTOX", " ", database_timestamp," v",vers,".Rda"))
    
    
  }
  
  return(ECOTOX)
  
}
