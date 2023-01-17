# Function version of C3 index condition discovery/scoring
# James Stanley, james.stanley@otago.ac.nz
# Adapted from M3 code
# v1.0 -- 17/01/2023: first version on GitHub

c3index_scoring <- function(in_df, in_df_precancer,
                            cancersite,
                            clinical_code_col="CLIN_CD", 
                            id_cols=c("MASTER_ENC", "snz_uid", "snz_moh_uid"),
                            return_condition_cols = TRUE,
                            return_score = TRUE){
  
  # in_df is a long dataframe with the ICD-10 codes for all events in the lookback period
    # this needs to have a per-person idenitifier(s) on it -- see below
    # default for datalab is snz_uid and snz_moh_uid
    # do NOT use per-event identifiers if you want to capture all
      # past admission diagnoses for an individual
  
  # in_df_precancer is the same dataset restricted to prior to the cancer-date admission
    # this is to rule out coding conditions that might be sequelae of the cancer
    # including e.g. mental health conditions arising from the cancer diagnosis
  
  # cancersite is a character string (single-element) that identifies the cancer 
  # site being used in analysis. This is used to select "other malignancies" to
  # detect in the analysis coding stage.
  
  # clinical_code_col is a character string (single-element) 
    # giving the name of the column that contains the ICD codes themselves
  
  # id_cols is a vector of character strings (can be length 1 or contain multiple identifiers)
    # giving the name of the columns that contain the unique identifiers
      # default for datalab is snz_uid and snz_moh_uid, from MoH often MASTER_ENC
    # DO NOT use per-event identifiers if you want to capture all past admission diagnoses for an individual
      # which is the way C3 is intended to work
  
  # Start with a message: if you ask for no results to be returned,
  # then don't bother running the function!
  if(!return_condition_cols & !return_score){
    cat("Function not run: you should ask for condition columns and/or score\n")
    return(NULL)
  }
  
  # Convert cancersite specification to uppercase (only used internally)
  cancersite <- str_to_upper(cancersite)
  
  if(!(cancersite %in% c("BLADDER", "BREAST", "COLON", "KIDNEY", "OVARIAN" ,
                         "UTERINE", "RECTAL", "STOMACH", "LIVER", 
                         "LUNG", "PROSTATE", "HEADNECK"))){
    cat("Function not run: cancersite argument not in list of supported cancer sites\n")
    cat("Contact james.stanley@otago.ac.nz for assistance\n")
    return(NULL)
  }

  # These next elements are used to "rule out" certain conditions that might
  # actually be symptoms of the cancer 
  # (e.g. anemia with both lower and upper GI cancers)
  
  # Set all potential condition overscales to 1 by default (which keeps that element)
  anemi_overscale <- coagu_overscale <- intes_overscale <- 
    liver_overscale <- renal_overscale <- uppgi_overscale <- urina_overscale <- 1
    
  if (cancersite %in% c("COLON", "RECTAL")){
    anemi_overscale <- intes_overscale <- 0
  }
  if (cancersite %in% c("LIVER", "STOMACH")){
    anemi_overscale <- coagu_overscale <- uppgi_overscale <- liver_overscale <- 0
  }
  if (cancersite %in% c("KIDNEY", "BLADDER")){  
    renal_overscale <- urina_overscale <- 0
  }

  # Now work out primary malignancy ICD codes
  PrimaryMalignancy <- case_when(
                         cancersite=="BLADDER" ~ "C67",
                         cancersite=="BREAST"  ~ "C50",
                         cancersite=="COLON"   ~ "C18\\|C19",
                         cancersite=="KIDNEY"  ~ "C64\\|C65",
                         cancersite=="OVARIAN" ~ "C56",
                         cancersite=="UTERINE" ~ "C53\\|C54\\|C55",
                         cancersite=="RECTAL"  ~ "C20",
                         cancersite=="STOMACH" ~ "C16",
                         cancersite=="LIVER"   ~ "C22",
                         
                         # Following are "bolt on" as not
                         # in the original C3 index
                         
                         cancersite=="LUNG" ~ "C34", # note that lung generally NOT taken as another malignancy...
                         cancersite=="PROSTATE" ~ "C61",
                         # Note slightly more complex specification here 
                         # to select two non-adjacent stretches for deletion
                         cancersite=="HEADNECK" ~ paste0("C0\\|C10\\|C11\\|C12\\|C13\\|C14\\|",
                                                         "|C30\\|C31\\|C32")
                         )
  
  # Master list of all commonly counted Other Malignancies
  OtherMalignancy <- paste0("C0|",
                            "C10|C11|C12|C13|C14|C15|C16|C17|C18|C19|",
                            "C20|C21|C23|C24|C25|C26|",
                            "C30|C31|C32|C33|C37|C38|C39|",
                            "C43|C45|C46|C47|C48|C49|",
                            "C50|C51|C52|C53|C54|C55|C56|C57|C58|",
                            "C60|C61|C62|C63|C64|C65|C66|C67|C68|C69|",
                            "C70|C72|C73|C74|C75|",
                            "C81|C82|C83|C84|C85|C88|",
                            "C90|C91|C92|C93|C94|C95"
                            )
  
  # And now delete primary malignancy list from OtherMalignancy
  # (so not counted as a second pre-existing cancer!)
  OtherMalignancy <- str_replace_all(OtherMalignancy, 
                                     paste0(PrimaryMalignancy, "\\|"), "")

  # Prelims section: Code to set up ICD 10 code lists --------------------------------------
  
  # This section defines the ICD codes for each C3 condition
  # (as regular expressions, separated by | delimiter)
  
  # Then codes exclusions/complications
  # using a similar process.
  
  # Prelims 1: ICD code lists for C3 conditions --------------------------------------------
  
  # CAN_003 Peripheral vascular
  f_CAN_003 <- "I70|I71|I720|I731|I738|I739|I771|K551|K552|K558|K559"

  # CAN_004 Venous insufficiency
  f_CAN_004 <- "I830|I832|I872"

  # CAN_005 Cerebrovascular disease
  f_CAN_005 <- "G45|G46|I60|I61|I62|I63|I64|I65|I66|I67|I69"

  # CAN_006 Dementia
  f_CAN_006 <- "F00|F01|F020|F021|F022|F023|F03|F051|G30|G310|G311"

  # CAN_008 Chronic pulmonary # Long line!
  f_CAN_008 <- paste0("E84|",
                      "J40|J41|J42|J43|J44|J45|J46|J47|",
                      "J60|J61|J62|J63|J64|J65|J66|J67|J684|",
                      "J701|J703|",
                      "J84|",
                      "J961|J980|J982|J983|J984")
  
  # CAN_009 Connective tissue
  f_CAN_009 <- paste0("L93|M05|M06|M08|M120|M123|",
                      "M30|M31|M32|M33|M34|M350|M351|M352|M353|M354|M355|M356|M358|M359")
  
  # CAN_010 Upper GI
  f_CAN_010 <- "K220|K221|K224|K225|K228|K229|K25|K26|K27|K28|K311|K312|K314|K316"

  # CAN_011 Diabetes uncomplicated
  f_CAN_011 <- "E100|E101|E109|E110|E111|E119|E120|E121|E129|E130|E131|E139|E140|E141|E149"
  
  # CAN_012 Diabetes complicated
  f_CAN_012 <- paste0("E102|E103|E104|E105|E106|E107|E108|",
                      "E112|E113|E114|E115|E116|E117|E118|",
                      "E122|E123|E124|E125|E126|E127|E128|",
                      "E132|E133|E134|E135|E136|E137|E138|",
                      "E142|E143|E144|E145|E146|E147|E148")

  # CAN_013 Paralysis  
  f_CAN_013 <- "G041|G114|G800|G801|G802|G81|G82|G830|G831|G832|G833|G834|G839"
  
  # CAN_014 Chronic renal
  f_CAN_014 <- paste0("I120|I131|",
                      "N032|N033|N034|N035|N036|N037|N038|N039|",
                      "N042|N043|N044|N045|N046|N047|N048|N049|",
                      "N052|N053|N054|N055|N056|N057|N058|N059|",
                      "N11|N18|N19|N250|N258|N259|",
                      "Z49|Z940|Z992")

  # CAN_017 Liver disease: moderate or severe
  f_CAN_017 <- paste0("I85|I864|I982|",
                      "K70|K711|K713|K714|K715|K717|K721|K729|K73|K74|",
                      "K760|K762|K763|K764|K765|K766|K767|K768|K769|",
                      "Z944")

  # CAN_019 Angina
  f_CAN_019 <- "I20"
  
  # CAN_023 Cardiac valve
  f_CAN_023 <- paste0("I05|I06|I07|I08|I091|I098|I34|I35|I36|I37|I38|",
                      "Q230|Q231|Q232|Q233|Q238|Q239|",
                      "T820|",
                      "Z952|Z953|Z954")

  # CAN_024 Bowel disease inflammatory
  f_CAN_024 <- "K50|K51|K522|K528|K529"
  
  # CAN_025 Other neurological disorders exc epilepsy
  f_CAN_025 <- paste0("G10|G110|G111|G112|G113|G118|G119|G12|G13|",
                      "G20|G21|G23|G255|",
                      "G312|G318|G319|G35|G36|G37|",
                      "G90|G934|",
                      "R470")
  
  # CAN_026 Epilepsy
  f_CAN_026 <- "G400|G401|G402|G403|G404|G406|G407|G408|G409|G41"
  
  # CAN_027 Muscular peripheral nerve disorder
  f_CAN_027 <- "G60|G61|G620|G621|G622|G628|G629|G64|G70|G71|G720|G721|G722|G723|G724|G728|G729|G731"

  # CAN_028 Major psychiatric disorder
  f_CAN_028 <- "F20|F22|F25|F28|F29|F302|F31|F321|F322|F323|F328|F329|F33|F39"

  # CAN_030 Coagulopathy and other blood disorder
  f_CAN_030 <- paste0("D55|D56|D57|D58|D590|D591|D592|D593|D594|D598|D599|",
                      "D60|D61|D64|D66|D67|D680|D681|D682|D688|D689|",
                      "D691|D692|D693|D694|D696|D698|D699|",
                      "D70|D71|D72|D74|D750|D752|D758|D759")
  
  # CAN_032 Obesity
  f_CAN_032 <- "E66"
  
  # CAN_033 Alcohol abuse
  f_CAN_033 <- paste0("F101|F102|F103|F104|F105|F106|F107|F108|F109|",
                      "K292|",
                      "Z502|Z714")
  
  # CAN_036 Endocrine disorder
  f_CAN_036 <- paste0("E01|E02|E03|E05|E062|E063|E065|E07|",
                      "E163|E164|E168|E169|",
                      "E20|E210|E212|E213|E214|E215|E22|E230|E232|E233|E236|E237|",
                      "E240|E241|E243|E244|E248|E249|E25|E26|E27|",
                      "E31|E32|E345|E348|E349")
  
  # CAN_037 Urinary tract problem chronic
  f_CAN_037 <- "N301|N302|N31|N32|N35|N36"

  # CAN_039 Osteoporosis and bone disorders
  f_CAN_039 <- paste0("M80|M810|M811|M815|M818|M819|",
                      "M831|M832|M833|M834|M835|M838|M839|",
                      "M85|M863|M864|M865|M866|M88")

  # CAN_041 Metabolic disorder
  f_CAN_041 <- "E70|E71|E72|E74|E75|E76|E77|E78|E791|E798|E799|E80|E83|E85|E88"

  # CAN_043 Hepatitis Chronic viral
  f_CAN_043 <- "B18|B942|Z225"
  
  # CAN_044 Sleep disorder
  f_CAN_044 <- "F51|G470|G471|G472|G473"
  
  # CAN_045 Inner ear disorder
  f_CAN_045 <- "H80|H81|H83|H90|H910|H911|H913|H918|H919|H930|H931|H932|H933"

  # CAN_048 Eye problem long term
  f_CAN_048 <- paste0("H16|H181|H184|H185|H186|",
                      "H201|H212|",
                      "H301|H311|H312|H313|H314|H330|H332|H333|H334|H335|H34|H35|",
                      "H43|H46|H47|H49|",
                      "H50|H51|H530|H531|H532|H533|H534|H536|H538|H539|H54|",
                      "Q12|Q13|Q14|Q15")
  
  # CAN_049 Cardiac disease other
  f_CAN_049 <- "I248|I249|I250|I251|I253|I254|I256|I258|I259|I310|I311|I421|I422|I424"

  # CAN_050 Intestinal disorder
  f_CAN_050 <- "K57|K592|K593|K90"
  
  # CAN_051 Joint spinal disorder
  f_CAN_051 <- paste0("M07|M13|M150|M151|M152|M154|M158|M159|",
                      "M400|M402|M403|M404|M405|M41|M42|M43|M45|M460|M461|M462|",
                      "M47|M480|M481|M482|M485|M488|M489|",
                      "G950|G951")
  
  
  # Pre-admission event diagnoses  --------------------------------------------
  
  # CAP_001 Myocardial infarction
  f_CAP_001 <- "I21|I22|I23|I241|I252"
  
  # CAP_002 Congestive heart failure
  f_CAP_002 <- "I099|I110|I130|I132|I255|I420|I425|I426|I427|I428|I429|I43|I50"

  # CAP_015 Other malignancy
  f_CAP_015 <- OtherMalignancy

  # CAP_020 Hypertension uncomplicated
  f_CAP_020 <- "I10|I119|I129|I139" 
  
  # CAP_021 Cardiac arrhythmia
  f_CAP_021 <- "I441|I442|I443|I456|I459|I47|I48|I49|T821|Z450|Z950"
  
  # CAP_022 Pulmonary circulation disorder
  f_CAP_022 <- "I26|I27|I280|I281|I288|I289"

  # CAP_029 Anxiety and Behavioral disorders
  f_CAP_029 <- paste0("F40|F41|F42|F44|F45|F48|",
                      "F50|F55|F59|F60|F61|F63|F64|F65|F66|F68|F69")

  # CAP_031 Anemia deficiency
  f_CAP_031 <- "D50|D51|D52|D53"

  # CAP_047 Malnutrition nutritional
  f_CAP_047 <- paste0("E40|E41|E42|E43|E44|E45|E46|",
                      "E50|E51|E52|E53|E54|E55|E56|E58|E59|",
                      "E60|E61|E63|E64")

  # Note on codes not used in C3  --------------------------------------------
  
  # I don't have notes on the provenance of their removal, but these are listed
  # here so that the full enumeration is clearer...
  
  # 7='Mental etc disorders from brain damage'
  # 16='Primary malignancy as in title'
  # 18='AIDS'
  # 34='Drug abuse'
  # 35='Pancreatitis (alcohol and not alcohol induced)'
  # 38='TB'
  # 40='Immune system disorders'
  # 42='Mental retardation'
  # 46='Chronic infection NOS'
  
  # Prelims 2: Make vector of ICD codes for C3 conditions ------------------------------

  f_CAP_codes <- c(f_CAP_001,
                   f_CAP_002,
                   f_CAP_015,
                   f_CAP_020,
                   f_CAP_021,
                   f_CAP_022,
                   f_CAP_029,
                   f_CAP_031,
                   f_CAP_047)
  
  f_CAN_codes <- c(# f_CAN_001, # in pre-period only
                   # f_CAN_002, # in pre-period only
                   f_CAN_003,
                   f_CAN_004,
                   f_CAN_005,
                   f_CAN_006,
                   # f_CAN_007, # not used in C3
                   f_CAN_008,
                   f_CAN_009,
                   f_CAN_010,
                   f_CAN_011,
                   f_CAN_012,
                   f_CAN_013,
                   f_CAN_014,
                   # f_CAN_015, # in pre-period only
                   # f_CAN_016, # not used in C3
                   f_CAN_017,
                   # f_CAN_018, # not used in C3
                   f_CAN_019,
                   # f_CAN_020, # in pre-period only
                   # f_CAN_021, # in pre-period only
                   # f_CAN_022, # in pre-period only
                   f_CAN_023,
                   f_CAN_024,
                   f_CAN_025,
                   f_CAN_026,
                   f_CAN_027,
                   f_CAN_028,
                   # f_CAN_029, # in pre-period only
                   f_CAN_030,
                   # f_CAN_031, # in pre-period only
                   f_CAN_032,
                   f_CAN_033,
                   # f_CAN_034, # not used in C3
                   # f_CAN_035, # not used in C3
                   f_CAN_036,
                   f_CAN_037,
                   # f_CAN_038, # not used in C3
                   f_CAN_039,
                   # f_CAN_040, # not used in C3
                   f_CAN_041,
                   # f_CAN_042, # not used in C3
                   f_CAN_043,
                   f_CAN_044,
                   f_CAN_045,
                   # f_CAN_046, # not used in C3
                   # f_CAN_047, # in pre-period only
                   f_CAN_048,
                   f_CAN_049,
                   f_CAN_050,
                   f_CAN_051)
  
  # Prelims 3: Make vector of names for C3 conditions ------------------------------
  
  # Labels in same order as the conditions above
  # End with a semi-colon for easy replacement at next step
  # otherwise using str_replace a partial match for J42 as "MY DISEASE" 
  # on a stored data code of J421 will result in "MY DISEASE1"
  f_CAN_label <- c("CAN_003 Peripheral vascular disease;",
                   "CAN_004 Venous insufficiency;",
                   "CAN_005 Cerebrovascular disease;",
                   "CAN_006 Dementia;",
                   "CAN_008 COPD and asthma;",
                   "CAN_009 Connective tissue disorders;",
                   "CAN_010 Upper GI disorders;",
                   "CAN_011 Diabetes no complications;",
                   "CAN_012 Diabetes with complications;",
                   "CAN_013 Paralysis;",
                   "CAN_014 Renal disease;",
                   "CAN_017 Liver disease moderate or severe;",
                   "CAN_019 Angina;",
                   "CAN_023 Cardiac valve disorders;",
                   "CAN_024 Inflammatory bowel disease;",
                   "CAN_025 Neurological excl epilepsy;",
                   "CAN_026 Epilepsy;",
                   "CAN_027 Peripheral nerve or muscular disorder;",
                   "CAN_028 Major psychiatric disorders;",
                   "CAN_030 Coagulopathies and other blood disorders;",
                   "CAN_032 Obesity;",
                   "CAN_033 Alcohol abuse;",                   
                   "CAN_036 Endocrine disorders;",
                   "CAN_037 Urinary tract disorder;",
                   "CAN_039 Osteoporosis and bone disorders;",
                   "CAN_041 Metabolic disorder;",
                   "CAN_043 Chronic viral hepatitis;",
                   "CAN_044 Sleep disorder;",
                   "CAN_045 Inner ear disorders;",
                   "CAN_048 Eye problems;",
                   "CAN_049 Other cardiac conditions;",
                   "CAN_060 Intestinal disorders;",
                   "CAN_061 Joint and spinal disorders;"
                   )
  
  f_CAP_label <- c("CAP_001 Myocardial infarction;",
                   "CAP_002 Congestive heart failure;",
                   "CAP_015 Other malignancy;",
                   "CAP_020 Hypertension;",
                   "CAP_021 Cardiac arrhythmia;",
                   "CAP_022 Pulmonary circulation disorder;",
                   "CAP_029 Anxiety and behavioral disorders;",
                   "CAP_031 Anemia;",
                   "CAP_047 Nutritional disorders;")
  
  # Prelims 4: Make a named vector of conditions and codes ----------------------------
  
  # Now take the LABELS as a vector
  f_CAN_coding_vector <- f_CAN_label
  # and name these labels with the matching text expressions
  # (counterintuitive... but this is how it works!)
  names(f_CAN_coding_vector) <- f_CAN_codes
  
  # Now take the LABELS as a vector
  f_CAP_coding_vector <- f_CAP_label
  # and name these labels with the matching text expressions
  # (counterintuitive... but this is how it works!)
  names(f_CAP_coding_vector) <- f_CAP_codes
  
  # C3 Step 1: Code for main condition set ----------------------------------

  # Discover and code the ICD codes across to the C3 index conditions
  # This looks at pre-admission diagnostic codes (prior to cancer admission)
  
  coded_step1_C3 <- in_df %>% 
    # Two steps here: str_replace_all uses the f_CAN_coding_vector to search and replace C3 codes
    # with their corresponding C3 variable name
    # Wrapped inside str_split_fixed, which was an inelegant way to strip trailing numbers.
    mutate(C3_code = str_split_fixed(str_replace_all(!!sym(clinical_code_col), f_CAN_coding_vector), ";", 2)[, 1]) %>%
    # If not an C3 code: then recode text as z_Exclude
    mutate(C3_code = ifelse(str_sub(C3_code, 1, 4)=="CAN_", C3_code, "z_Exclude")) %>%
    # Filter out exclude rows (can turn off line for testing coding works)
    filter(C3_code != "z_Exclude") %>%
    # Set C3_indi value to 1 for spread function to make a binary variable for each condition
    mutate(C3_indi = 1) %>%
    # Arrange by study ID
    arrange(across(all_of(id_cols))) %>%
    # Just keep study IDs and C3_code variable
    select(all_of(id_cols), C3_code, C3_indi) %>%
    # get rid of duplicate rows (i.e. reduce to one row per condition per patient)
    distinct() %>%
    # UPDATE Nov 2022: to make sure all codes are incorporated as columns
    bind_rows(data.frame(C3_code = str_replace_all(f_CAN_label,  ";", ""), remove_me = 1)) %>% 
    # now pivot_wider (convert LONG table to WIDE)
    # and put a zero in place for not-recorded conditions
    pivot_wider(names_from = C3_code, values_from=C3_indi, values_fill=0, names_sort = TRUE) %>% 
    # remove dummy rows from above
    filter(is.na(remove_me)) %>% select(-remove_me) %>% 
    # Set Diabetes uncomplicated to zero if 012 diabetes complicated recorded
    mutate(`CAN_011 Diabetes no complications` = ifelse(`CAN_012 Diabetes with complications`==1, 0, `CAN_011 Diabetes no complications`))
  
  # C3 Step 2: Code for conditions prior to cancer admission  ----------------------------------
  coded_step2_preadmission <- in_df_precancer %>%
    # See annotations as per step 1
    mutate(C3_code = str_split_fixed(str_replace_all(!!sym(clinical_code_col), f_CAP_coding_vector), ";", 2)[, 1]) %>%
    mutate(C3_code = ifelse(str_sub(C3_code, 1, 4)=="CAP_", C3_code, "z_Exclude")) %>%
    filter(C3_code != "z_Exclude") %>%
    mutate(C3_indi = 1) %>%
    arrange(across(all_of(id_cols))) %>%
    select(all_of(id_cols), C3_code, C3_indi) %>%
    distinct() %>%
    bind_rows(data.frame(C3_code = str_replace_all(f_CAP_label,  ";", ""), remove_me = 1)) %>% 
    pivot_wider(names_from = C3_code, values_from=C3_indi, values_fill=0, names_sort = TRUE) %>% 
    filter(is.na(remove_me)) %>% select(-remove_me)

  # C3 Step 3: bring together all diagnoses for both reviewed periods --------------
  coded_step3_joinalldiagnoses <-  
    full_join(coded_step1_C3, coded_step2_preadmission) %>% 
    mutate(across(starts_with("CAN"), ~ifelse(is.na(.), 0, .)),
           across(starts_with("CAP"), ~ifelse(is.na(.), 0, .))) 
  
  # Reorder columns based on numbering system here (not currently in use)
  # coded_step3_joinalldiagnoses <- coded_step3_joinalldiagnoses %>% 
  #   select(any_of(id_cols), order(colnames(coded_step3_joinalldiagnoses)))

  # C3 Step 4: rename some columns, set a maximum length of 32 chars --------------  
  coded_step4_rename <- coded_step3_joinalldiagnoses %>% 
    # Change names from CAN_001 to C3_
    rename_with(~paste0("C3_", str_sub(.x, 9, 99)), -any_of(id_cols)) %>% 
    # Replace annoying characters
    rename_with(~str_replace_all(.x, " |:", "_"), -any_of(id_cols)) %>% 
    # Also trim names to 32 characters (or 31 characters if ends in _)
    rename_with(~ifelse(str_sub(.x, 32, 32)=="_",
                        str_sub(.x, 1, 31),
                        str_sub(.x, 1, 32)), -any_of(id_cols))

  # Alphabetical order reordering of columns
  coded_step4_rename <- coded_step4_rename %>%
    select(any_of(id_cols), order(colnames(coded_step4_rename)))
  
  # C3 Step 5: Make final scored dataset, and prepare for return dataframe --------------
  
  final_C3_scored <- coded_step4_rename %>% 
    mutate(C3score_allsites = 
           C3_Alcohol_abuse                    * 1.08     +
           C3_Anemia                           * 0.59  * anemi_overscale +
           C3_Angina                           * 0.51     +
           C3_Anxiety_and_behavioral_disord    * 0.57     +
           C3_Cardiac_arrhythmia               * 0.77     +
           C3_Cardiac_valve_disorders          * 1.10     +
           C3_Cerebrovascular_disease          * 1.09     +
           C3_Chronic_viral_hepatitis          * 0.39     +
           C3_Coagulopathies_and_other_bloo    * 0.75  * coagu_overscale +
           C3_Congestive_heart_failure         * 1.26     +
           C3_Connective_tissue_disorders      * 0.51     +
           C3_COPD_and_asthma                  * 1.09     +
           C3_Dementia                         * 1.35     +
           C3_Diabetes_no_complications        *-0.03     +
           C3_Diabetes_with_complications      * 0.88     +
           C3_Endocrine_disorders              * 0.77     +
           C3_Epilepsy                         * 1.04     +
           C3_Eye_problems                     * 0.63     +
           C3_Hypertension                     * 0.72     +     
           C3_Inflammatory_bowel_disease       * 0.52     +
           C3_Inner_ear_disorders              * 0.54     +
           C3_Intestinal_disorders             * 0.11  * intes_overscale +
           C3_Joint_and_spinal_disorders       * 0.69     +
           C3_Liver_disease_moderate_or_sev    * 0.92  * liver_overscale +
           C3_Major_psychiatric_disorders      * 0.79     +
           C3_Metabolic_disorder               * 0.61     +
           C3_Myocardial_infarction            * 0.93     +   
           C3_Neurological_excl_epilepsy       * 1.06     +
           C3_Nutritional_disorders            * 1.16     +
           C3_Obesity                          * 0.83     +
           C3_Osteoporosis_and_bone_disorde    * 0.49     +
           C3_Other_cardiac_conditions         * 0.62     +
           C3_Other_malignancy                 * 0.17     +
           C3_Paralysis                        * 1.03     +
           C3_Peripheral_nerve_or_muscular     * 1.20     +
           C3_Peripheral_vascular_disease      * 0.98     +
           C3_Pulmonary_circulation_disorde    * 0.95     +
           C3_Renal_disease                    * 1.38  * renal_overscale +    
           C3_Sleep_disorder                   * 1.41     +
           C3_Upper_GI_disorders               * 0.11  * uppgi_overscale +
           C3_Urinary_tract_disorder           * 0.12  * urina_overscale +
           C3_Venous_insufficiency             * 0.70
    ) %>% 

  # Added in categorised version for C3 descriptive statistics
  mutate(C3cat_allsites = case_when(C3score_allsites <= 0 ~ 0,
                                    C3score_allsites <= 1 ~ 1,
                                    C3score_allsites <= 2 ~ 2,
                                    C3score_allsites >  2 ~ 3,
                                    is.na(C3score_allsites) ~ as.numeric(NA))
         )

  # Finally: remove unwanted columns from C3 process as per arguments
  # If didn't want the individual-column-conditions, then drop those
  if(!return_condition_cols){
    final_C3_scored <- select(final_C3_scored, -starts_with("C3_"))
  }
  # If didn't want the C3Score, then drop it
  if(!return_score){
    final_C3_scored <- select(final_C3_scored, -C3score_allsites, -C3cat_allsites)
  }
    
  # Return this final object
  final_C3_scored
   
}