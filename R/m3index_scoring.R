# Function version of M3 index condition discovery/scoring
# James Stanley, james.stanley@otago.ac.nz
# Adapted from standalone code written in 2018
# v0.1 -- initial non-function version, was tested outside of Datalab against the original SAS macro
# v0.7 -- first version of function Nov 23 2022, built in IDI Datalab environment
# v0.8 -- refined and tested outside datalab on National Collections data

# v1.0 -- 17/02/2023. Tested and ready for release on website.

m3index_scoring <- function(in_df, 
                            clinical_code_col="CLIN_CD", 
                            id_cols=c("MASTER_ENC", "snz_uid", "snz_moh_uid"),
                            return_condition_cols = TRUE,
                            return_score = TRUE){
  
  # in_df is a long dataframe with the ICD-10 codes for all events in the lookback period
    # this needs to have a per-person idenitifier(s) on it -- see below
    # default for datalab is snz_uid and snz_moh_uid
    # do NOT use per-event identifiers if you want to capture all
      # past admission diagnoses for an individual
  
  # clinical_code_col is a character string (single-element) 
    # giving the name of the column that contains the ICD codes themselves
  
  # id_cols is a vector of character strings (can be length 1 or contain multiple identifiers)
    # giving the name of the columns that contain the unique identifiers
      # default for datalab is snz_uid and snz_moh_uid
    # DO NOT use per-event identifiers if you want to capture all past admission diagnoses for an individual
      # which is the way M3 is intended to work
  
  # Start with a message: if you ask for no results to be returned,
  # then don't bother running the function!
  if(!return_condition_cols & !return_score){
    cat("Function not run: you should ask for condition columns and/or score\n")
    return(NULL)
  }
  
  # Prelims section: Code to set up ICD 10 code lists --------------------------------------
  
  # This section defines the ICD codes for each M3 condition
  # (as regular expressions, separated by | delimiter)
  
  # Then codes exclusions/complications
  # using a similar process.
  
  # Prelims 1: ICD code lists for M3 conditions --------------------------------------------
  
  # ICD_001 Myocardial infarction
  f_ICD_001 <- "I21|I22|I23|I241|I252"
  
  # Congestive heart failure
  f_ICD_002 <- "I099|I110|I130|I132|I255|I420|I425|I426|I427|I428|I429|I43|I50"
  
  # ICD_003 Peripheral vascular
  f_ICD_003 <- "I70|I731|I738|I739|I74|I771|K551|K552|K558|K559"
  
  # ICD_004 Aortic and other aneurysms
  f_ICD_004 <- "I71|I72"
  
  # ICD_005 Venous insufficiency
  f_ICD_005 <- "I830|I832|I872"
  
  # ICD_006 Cerebrovascular disease
  f_ICD_006 <- "I60|I61|I62|I63|I64|I65|I66|I67|I69|G45|G46"
  
  # ICD_007 Dementia
  f_ICD_007 <- "F00|F01|F020|F021|F022|F023|F03|F051|G30|G310|G311"
  
  # ICD_008 Mental and behavioural disorders due to brain damage
  f_ICD_008 <- "F04|F06|F070|F071|F078|F079|F09|G931"
  
  # ICD_009 Chronic pulmonary # Long line!
  f_ICD_009 <- "E84|J40|J41|J42|J43|J44|J45|J46|J47|J60|J61|J62|J63|J64|J65|J66|J67|J684|J701|J703|J84|J961|J980|J982|J983|J984"
  
  # ICD_010 Connective tissue
  f_ICD_010 <- "L93|M05|M06|M08|M120|M123|M30|M31|M32|M33|M34|M350|M351|M352|M353|M354|M355|M356|M358|M359"
  
  # ICD_011 GI ulcer upper GI
  f_ICD_011 <- "K220|K221|K224|K225|K228|K229|K25|K26|K27|K28|K311|K312|K314|K316"
  
  # ICD_012 Diabetes uncomplicated
  f_ICD_012 <- "E100|E101|E109|E110|E111|E119|E120|E121|E129|E130|E131|E139|E140|E141|E149"
  
  # ICD_013 Diabetes complicated
  f_ICD_013 <- "E102|E103|E104|E105|E106|E107|E108|E112|E113|E114|E115|E116|E117|E118|E122|E123|E124|E125|E126|E127|E128|E132|E133|E134|E135|E136|E137|E138|E142|E143|E144|E145|E146|E147|E148"
  
  # ICD_014 Paralysis  
  f_ICD_014 <- "G041|G114|G800|G801|G802|G81|G82|G830|G831|G832|G833|G834|G839"
  
  # ICD_015 Chronic renal
  f_ICD_015 <- "I120|I129|I131|I139|Q60|Q611|Q612|Q613|N032|N033|N034|N035|N036|N037|N038|N039|N042|N043|N044|N045|N046|N047|N048|N049|N052|N053|N054|N055|N056|N057|N058|N059|N11|N18|N19|N250|N258|N259|Z49|Z940|Z992"
  
  # ICD_016 Colorectal cancer
  f_ICD_016 <- "C18|C19|C20|C21"
  
  # ICD_017 Breast cancer
  f_ICD_017 <- "C50"
  
  # ICD_018 Prostate cancer
  f_ICD_018 <- "C61"
  
  # ICD_019 Lung cancer
  f_ICD_019 <- "C33|C34"
  
  # ICD_020 Lymphomas and leukaemias
  f_ICD_020 <- "C81|C82|C83|C84|C85|C91|C92|C93|C94|C95|C96"
  
  # ICD_021 Upper gastrointestinal cancers
  f_ICD_021 <- "C15|C16|C17|C22|C23|C24|C25"
  
  # ICD_022 Malignant melanoma
  f_ICD_022 <- "C43"
  
  # ICD_023 Gynaecological cancers
  f_ICD_023 <- "C51|C52|C53|C54|C55|C56|C57|C58"
  
  # ICD_024 Other cancers
  f_ICD_024 <- "C0|C10|C11|C12|C13|C14|C26|C30|C31|C32|C37|C38|C39|C40|C41|C45|C46|C47|C48|C49|C60|C62|C63|C64|C65|C66|C67|C68|C69|C70|C71|C72|C73|C74|C75|C76|C88|C90"
  
  # ICD_025 Metastatic cancer
  f_ICD_025 <- "C77|C78|C79"
  
  # ICD_026 Liver disease: moderate or severe
  f_ICD_026 <- "I85|I864|I982|K70|K711|K713|K714|K715|K717|K721|K729|K73|K74|K760|K762|K763|K764|K765|K766|K767|K768|K769|Z944"
  
  # ICD_027 AIDS
  f_ICD_027 <- "B20|B21|B22|B23|B24|F024|Z21"
  
  # ICD_028 Angina
  f_ICD_028 <- "I20"
  
  # ICD_029 Hypertension uncomplicated
  f_ICD_029 <- "I10" 
  
  #   ICD_030 Cardiac arrhythmia
  f_ICD_030 <- "I441|I442|I443|I456|I459|I47|I48|I49|T821|Z450|Z950"
  
  # ICD_031 Pulmonary circulation disorder
  f_ICD_031 <- "I26|I27|I280|I281|I288|I289"
  
  # ICD_032 Cardiac valve
  f_ICD_032 <- "I05|I06|I07|I08|I091|I098|I34|I35|I36|I37|I38|T820|Q230|Q231|Q232|Q233|Q238|Q239|Z952|Z953|Z954"
  
  # ICD_033 Bowel disease inflammatory
  f_ICD_033 <- "K50|K51|K522|K528|K529"
  
  # ICD_034 Other neurological disorders exc epilepsy
  f_ICD_034 <- "G10|G110|G111|G112|G113|G118|G119|G12|G13|G20|G21|G23|G255|G312|G318|G319|G35|G36|G37|G90|G934|R470"
  
  # ICD_035 Epilepsy
  f_ICD_035 <- "G400|G401|G402|G403|G404|G406|G407|G408|G409|G41"
  
  # ICD_036 Muscular peripheral nerve disorder
  f_ICD_036 <- "G60|G61|G620|G621|G622|G628|G629|G64|G70|G71|G720|G721|G722|G723|G724|G728|G729|G731"
  
  # ICD_037 Major psychiatric disorder
  f_ICD_037 <- "F20|F22|F25|F28|F29|F302|F31|F321|F322|F323|F328|F329|F33|F39"
  
  # ICD_038 Anxiety and Behavioural disorders
  f_ICD_038 <- "F40|F41|F42|F44|F45|F48|F50|F55|F59|F60|F61|F63|F64|F65|F66|F68|F69"
  
  # ICD_039 Coagulopathy and other blood disorder
  f_ICD_039 <- "D55|D56|D57|D58|D590|D591|D592|D593|D594|D598|D599|D60|D61|D64|D66|D67|D680|D681|D682|D688|D689|D691|D692|D693|D694|D696|D698|D699|D70|D71|D72|D74|D750|D752|D758|D759"
  
  # ICD_040 Anemia deficiency
  f_ICD_040 <- "D50|D51|D52|D53"
  
  # ICD_041 Obesity
  f_ICD_041 <- "E66"
  
  # ICD_042 Alcohol abuse
  f_ICD_042 <- "F101|F102|F103|F104|F105|F106|F107|F108|F109|Z502|Z714"
  
  # ICD_043 Drug abuse
  f_ICD_043 <- "F11|F12|F13|F14|F15|F16|F18|F19|Z503|Z715|Z722"
  
  # ICD_044 Pancreatitis
  f_ICD_044 <- "K85|K860|K861|K868"
  
  # ICD_045 Endocrine disorder
  f_ICD_045 <- "E01|E02|E03|E05|E062|E063|E065|E07|E163|E164|E168|E169|E20|E210|E212|E213|E214|E215|E22|E230|E232|E233|E236|E237|E240|E241|E243|E244|E248|E249|E25|E26|E27|E31|E32|E345|E348|E349"
  
  # ICD_046 Urinary tract problem chronic
  f_ICD_046 <- "N301|N302|N31|N32|N35|N36"
  
  # ICD_047 Tuberculosis
  f_ICD_047 <- "A15|A16|A17|A18|A19|B90"
  
  # ICD_048 Bone disorders
  f_ICD_048 <- "M80|M830|M831|M832|M833|M834|M835|M838|M839|M85|M863|M864|M865|M866|M88"
  
  # ICD_049 Osteoporosis Uncomplicated
  f_ICD_049 <- "M810|M811|M815|M818|M819"
  
  # ICD_050 Immune system disorder
  f_ICD_050 <- "D80|D81|D82|D83|D84|D86|D89"
  
  # ICD_051 Metabolic disorder
  f_ICD_051 <- "E70|E71|E72|E74|E75|E76|E77|E78|E791|E798|E799|E80|E83|E85"
  
  # ICD_052 Mental retardation
  f_ICD_052 <- "F70|F71|F72|F73|F78|F79|F842|F843|F844|E000|E001|E002|E009|Q90"
  
  # ICD_053 Hepatitis Chronic viral
  f_ICD_053 <- "B18|B942|Z225"
  
  # ICD_054 Sleep disorder
  f_ICD_054 <- "F51|G470|G471|G472|G473"
  
  # ICD_055 Inner ear disorder
  f_ICD_055 <- "H80|H81|H83|H90|H910|H911|H913|H918|H919|H930|H931|H932|H933"
  
  # ICD_056 Infection Chronic NOS
  f_ICD_056 <- "A30|A31|A52|B91|B92|B941|B948|B949"
  
  # ICD_057 Malnutrition nutritional
  f_ICD_057 <- "E40|E41|E42|E43|E44|E45|E46|E50|E51|E52|E53|E54|E55|E56|E58|E59|E60|E61|E63|E64"
  
  # ICD_058 Eye problem long term
  f_ICD_058 <- "H16|H181|H184|H185|H186|H201|H212|H301|H311|H312|H313|H314|H330|H332|H333|H334|H335|H34|H35|H43|H46|H47|H49|H50|H51|H530|H531|H532|H533|H534|H536|H538|H539|H54|Q12|Q13|Q14|Q15"
  
  # ICD_059 Cardiac disease other
  f_ICD_059 <- "I119|I248|I249|I250|I251|I253|I254|I256|I258|I259|I310|I311|I421|I422|I424"
  
  # ICD_060 Intestinal disorder
  f_ICD_060 <- "K57|K592|K593|K90"
  
  # ICD_061 Joint spinal disorder
  f_ICD_061 <- "M07|M13|M150|M151|M152|M154|M158|M159|M400|M402|M403|M404|M405|M41|M42|M43|M45|M460|M461|M462|M47|M480|M481|M482|M485|M488|M489|G950|G951"
  
  # Exclusions and complications --------------------------------------------
  
  
  # Diabetes complications (if found, recode uncomplicated diabetes to complicated.)
  c_DIA_Complication <- "I20|I21|I22|I23|I24|I25|I6|I7|N03|N04|N18|G603|G62|G638|H35|H36|L97"
  
  # Osteoporosis exclusions (if found, then remove osteoporosis code.)
  e_OST_Exclusions <- "M80|S220|S320|S52|S720"
  
  # Hypertension exclusions (if found, then remove uncomplicated HT code.)
  e_HYP_Exclusions <- "I11|I12|I13|I20|I21|I22|I23|I24|I25|I6|I70|I71|I72|N03|N04|N18"
  
  # Prelims 2: Make vector of ICD codes for M3 conditions ------------------------------

  f_ICD_codes <- c(f_ICD_001,
                   f_ICD_002,
                   f_ICD_003,
                   f_ICD_004,
                   f_ICD_005,
                   f_ICD_006,
                   f_ICD_007,
                   f_ICD_008,
                   f_ICD_009,
                   f_ICD_010,
                   f_ICD_011,
                   f_ICD_012,
                   f_ICD_013,
                   f_ICD_014,
                   f_ICD_015,
                   f_ICD_016,
                   f_ICD_017,
                   f_ICD_018,
                   f_ICD_019,
                   f_ICD_020,
                   f_ICD_021,
                   f_ICD_022,
                   f_ICD_023,
                   f_ICD_024,
                   f_ICD_025,
                   f_ICD_026,
                   f_ICD_027,
                   f_ICD_028,
                   f_ICD_029,
                   f_ICD_030,
                   f_ICD_031,
                   f_ICD_032,
                   f_ICD_033,
                   f_ICD_034,
                   f_ICD_035,
                   f_ICD_036,
                   f_ICD_037,
                   f_ICD_038,
                   f_ICD_039,
                   f_ICD_040,
                   f_ICD_041,
                   f_ICD_042,
                   f_ICD_043,
                   f_ICD_044,
                   f_ICD_045,
                   f_ICD_046,
                   f_ICD_047,
                   f_ICD_048,
                   f_ICD_049,
                   f_ICD_050,
                   f_ICD_051,
                   f_ICD_052,
                   f_ICD_053,
                   f_ICD_054,
                   f_ICD_055,
                   f_ICD_056,
                   f_ICD_057,
                   f_ICD_058,
                   f_ICD_059,
                   f_ICD_060,
                   f_ICD_061)
  
  # Prelims 3: Make vector of names for M3 conditions ------------------------------
  
  # Labels in same order
  # End with a semi-colon for easy replacement at next step
  # otherwise using str_replace a partial match for J42 as "MY DISEASE" 
  # on a stored data code of J421 will result in "MY DISEASE1"
  f_ICD_label <- c("ICD_001 Myocardial infarction;",
                   "ICD_002 Congestive heart failure;",
                   "ICD_003 Peripheral vascular;",
                   "ICD_004 Aortic and other aneurysms;",
                   "ICD_005 Venous insufficiency;",
                   "ICD_006 Cerebrovascular disease;",
                   "ICD_007 Dementia;",
                   "ICD_008 Mental and behavioural disorders due to brain damage;",
                   "ICD_009 Chronic pulmonary;",
                   "ICD_010 Connective tissue;",
                   "ICD_011 GI ulcer upper GI;",
                   "ICD_012 Diabetes uncomplicated;",
                   "ICD_013 Diabetes complicated;",
                   "ICD_014 Paralysis;",
                   "ICD_015 Chronic renal;",
                   "ICD_016 Colorectal cancer;",
                   "ICD_017 Breast cancer;",
                   "ICD_018 Prostate cancer;",
                   "ICD_019 Lung cancer;",
                   "ICD_020 Lymphomas and leukaemias;",
                   "ICD_021 Upper gastrointestinal cancers;",
                   "ICD_022 Malignant melanoma;",
                   "ICD_023 Gynaecological cancers;",
                   "ICD_024 Other cancers;",
                   "ICD_025 Metastatic cancer;",
                   "ICD_026 Liver disease: moderate or severe;",
                   "ICD_027 AIDS;",
                   "ICD_028 Angina;",
                   "ICD_029 Hypertension uncomplicated;",
                   "ICD_030 Cardiac arrhythmia;",
                   "ICD_031 Pulmonary circulation disorder;",
                   "ICD_032 Cardiac valve;",
                   "ICD_033 Bowel disease inflammatory;",
                   "ICD_034 Other neurological disorders exc epilepsy;",
                   "ICD_035 Epilepsy;",
                   "ICD_036 Muscular peripheral nerve disorder;",
                   "ICD_037 Major psychiatric disorder;",
                   "ICD_038 Anxiety and Behavioural disorders;",
                   "ICD_039 Coagulopathy and other blood disorder;",
                   "ICD_040 Anemia deficiency;",
                   "ICD_041 Obesity;",
                   "ICD_042 Alcohol abuse;",
                   "ICD_043 Drug abuse;",
                   "ICD_044 Pancreatitis;",
                   "ICD_045 Endocrine disorder;",
                   "ICD_046 Urinary tract problem chronic;",
                   "ICD_047 Tuberculosis;",
                   "ICD_048 Bone disorders;",
                   "ICD_049 Osteoporosis Uncomplicated;",
                   "ICD_050 Immune system disorder;",
                   "ICD_051 Metabolic disorder;",
                   "ICD_052 Mental retardation;",
                   "ICD_053 Hepatitis Chronic viral;",
                   "ICD_054 Sleep disorder;",
                   "ICD_055 Inner ear disorder;",
                   "ICD_056 Infection Chronic NOS;",
                   "ICD_057 Malnutrition nutritional;",
                   "ICD_058 Eye problem long term;",
                   "ICD_059 Cardiac disease other;",
                   "ICD_060 Intestinal disorder;",
                   "ICD_061 Joint spinal disorder;")
  
  
  # Prelims 4: Make a named vector of conditions and codes ----------------------------
  
  # Now take the LABELS as a vector
  f_ICD_coding_vector <- f_ICD_label
  # and name these labels with the matching text expressions
  # (counterintuitive... but this is how it works!)
  names(f_ICD_coding_vector) <- f_ICD_codes
  
  # Prelims 5: Make named vectors of complications/exclusions and codes ---------------------------
  
  f_ICD_comp_codes <- "ICD_XXX Diabetes Complications;"
  names(f_ICD_comp_codes) <- c_DIA_Complication
  
  f_ICD_excl_codes  <- c(e_OST_Exclusions, e_HYP_Exclusions)
  f_ICD_excl_labels <- c("ICD_XXY Osteoporosis exclusions;",
                         "ICD_XXZ Hypertension exclusions;")
  
  # Now take the LABELS as a vector
  f_ICD_cod_exc_vector <- f_ICD_excl_labels
  # and name these labels with the matching text expressions
  # (counterintuitive... but this is how it works!)
  names(f_ICD_cod_exc_vector) <- f_ICD_excl_codes
  
  # M3 Step 1: Code for main condition set ----------------------------------

  # Discover and code the ICD codes across to the M3 index conditions
  # This first step is reasonably fully annotated: subs
  
  coded_step1_M3 <- in_df %>% 
    # Two steps here: str_replace_all uses the f_ICD_coding_vector to search and replace M3 codes
    # with their corresponding M3 variable name
    # Wrapped inside str_split_fixed, which was an inelegant way to strip trailing numbers.
    mutate(M3_code = str_split_fixed(str_replace_all(!!sym(clinical_code_col), f_ICD_coding_vector), ";", 2)[, 1]) %>%
    # If not an M3 code: then recode text as z_Exclude
    mutate(M3_code = ifelse(str_sub(M3_code, 1, 4)=="ICD_", M3_code, "z_Exclude")) %>%
    # Filter out exclude rows (can turn off line for testing coding works)
    filter(M3_code != "z_Exclude") %>%
    # Set M3_indi value to 1 for spread function to make a binary variable for each condition
    mutate(M3_indi = 1) %>%
    # Arrange by study ID
    arrange(across(all_of(id_cols))) %>%
    # Just keep study IDs and M3_code variable
    select(all_of(id_cols), M3_code, M3_indi) %>%
    # get rid of duplicate rows (i.e. reduce to one row per condition per patient)
    distinct() %>%
    # UPDATE Nov 2022: to make sure all codes are incorporated as columns
    bind_rows(data.frame(M3_code = str_replace_all(f_ICD_label,  ";", ""), remove_me = 1)) %>% 
    # now pivot_wider (convert LONG table to WIDE)
    # and put a zero in place for not-recorded conditions
    pivot_wider(names_from = M3_code, values_from=M3_indi, values_fill=0, names_sort = TRUE) %>% 
    # remove dummy rows from above
    filter(is.na(remove_me)) %>% select(-remove_me)
  
  # M3 Step 2: Code for diabetes complications ----------------------------------
  coded_step2_Complications <- in_df %>%
    mutate(M3_comp = str_split_fixed(str_replace_all(!!sym(clinical_code_col), f_ICD_comp_codes), ";", 2)[, 1]) %>%
    mutate(M3_comp = ifelse(str_sub(M3_comp, 1, 4)=="ICD_", M3_comp, "z_Exclude")) %>%
    filter(M3_comp != "z_Exclude") %>%
    mutate(M3_indi = 1) %>%
    arrange(across(all_of(id_cols))) %>%
    select(all_of(id_cols), M3_comp, M3_indi) %>%
    distinct() %>%
    bind_rows(data.frame(M3_comp = str_replace_all(f_ICD_comp_codes,  ";", ""), remove_me = 1)) %>% 
    pivot_wider(names_from = M3_comp, values_from=M3_indi, values_fill=0, names_sort = TRUE) %>% 
    filter(is.na(remove_me)) %>% select(-remove_me)


  # M3 Step 3: Code for osteoporosis/hypertension exclusions --------------------
  coded_step3_Exclusions <- in_df %>%
    mutate(M3_excl = str_split_fixed(str_replace_all(!!sym(clinical_code_col), f_ICD_cod_exc_vector), ";", 2)[, 1]) %>%
    mutate(M3_excl = ifelse(str_sub(M3_excl, 1, 4)=="ICD_", M3_excl, "z_Exclude")) %>%
    filter(M3_excl != "z_Exclude") %>%
    mutate(M3_indi = 1) %>%
    arrange(across(all_of(id_cols))) %>%
    select(all_of(id_cols), M3_excl, M3_indi) %>%
    distinct() %>%
    bind_rows(data.frame(M3_excl = str_replace_all(f_ICD_excl_labels,  ";", ""), remove_me = 1)) %>% 
    pivot_wider(names_from = M3_excl, values_from=M3_indi, values_fill=0, names_sort = TRUE) %>% 
    filter(is.na(remove_me)) %>% select(-remove_me)


  # M3 Step 4: bring together codes, complications, exclusions --------------
  coded_step4_ApplyComplicationsExclusions <- coded_step1_M3 %>%
    left_join(coded_step2_Complications) %>%
    left_join(coded_step3_Exclusions) %>%
    # Fill blanks for complications/exclusions (when not found, replace with a zero)
    mutate(`ICD_XXX Diabetes Complications` = ifelse(is.na(`ICD_XXX Diabetes Complications`), 0, `ICD_XXX Diabetes Complications`)) %>%
    mutate(`ICD_XXY Osteoporosis exclusions` = ifelse(is.na(`ICD_XXY Osteoporosis exclusions`), 0, `ICD_XXY Osteoporosis exclusions`)) %>%
    mutate(`ICD_XXZ Hypertension exclusions` = ifelse(is.na(`ICD_XXZ Hypertension exclusions`), 0, `ICD_XXZ Hypertension exclusions`)) %>%
    # now process these complications/exclusions
    # Diabetes complicated: add if has diabetes uncomplicated + one or more complication codes...
    mutate(`ICD_013 Diabetes complicated` = ifelse(`ICD_012 Diabetes uncomplicated`==1 & 
                                                     `ICD_XXX Diabetes Complications`==1, 1, `ICD_013 Diabetes complicated`)) %>%
    # Set Diabetes uncomplicated to zero if 013 diabetes complicated recorded or XXX complications found
    mutate(`ICD_012 Diabetes uncomplicated` = ifelse(`ICD_013 Diabetes complicated`==1 | 
                                                       `ICD_XXX Diabetes Complications`==1, 0, `ICD_012 Diabetes uncomplicated`)) %>%
    mutate(`ICD_049 Osteoporosis Uncomplicated` = ifelse(`ICD_XXY Osteoporosis exclusions`==1, 0, `ICD_049 Osteoporosis Uncomplicated`)) %>%
    mutate(`ICD_029 Hypertension uncomplicated` = ifelse(`ICD_XXZ Hypertension exclusions`==1, 0, `ICD_029 Hypertension uncomplicated`)) %>%
    # Post-testing: drop the complications columns
    select(-starts_with("ICD_XX"))
  
  # M3 Step 5: Deal with metastatic cancer handling --------------
  # Need to remove metastatic cancers from other cancer conditions here
  coded_step5_MetastaticCancers <- coded_step4_ApplyComplicationsExclusions %>%
    # Change all cancer columns EXCEPT for metastatic cancer
    mutate_at(vars(contains("cancer"), contains("leukaemia"), contains("melanoma"), -contains("Metastatic")), 
              # and replace with zeros if metastatic cancer found.
              ~ifelse(`ICD_025 Metastatic cancer`==1, 0, .))
  
  # M3 Step 6: rename some columns, set a maximum length of 32 chars --------------  
  coded_step6_rename <- coded_step5_MetastaticCancers %>% 
    # Change names from ICD_001 to M3_
    rename_with(~paste0("M3_", str_sub(.x, 9, 99)), -any_of(id_cols)) %>% 
    # Replace annoying characters
    rename_with(~str_replace_all(.x, " |:", "_"), -any_of(id_cols)) %>% 
    # Also trim names to 32 characters (or 31 characters if ends in _)
    rename_with(~ifelse(str_sub(.x, 32, 32)=="_",
                        str_sub(.x, 1, 31),
                        str_sub(.x, 1, 32)), -any_of(id_cols))
  
  # M3 Step 7: Make final scored dataset, and prepare for return dataframe --------------
  final_M3_scored <- coded_step6_rename %>% 
    mutate(M3Score = 
             M3_AIDS								          *	0.452647425	+
             M3_Alcohol_abuse					        *	0.576907507	+
             M3_Anemia_deficiency				      *	0.180927466	+
             M3_Angina						          	*	0 			    + #-0.082399267	
             M3_Anxiety_and_Behavioural_disor	*	0.121481351	+
             M3_Aortic_and_other_aneurysms		*	0.260195993	+
             M3_Bone_disorders					      *	0.132827597	+
             M3_Bowel_disease_inflammatory		*	0.086960591	+
             M3_Breast_cancer					        *	0.411891435	+
             M3_Cardiac_arrhythmia				    *	0.173859876	+
             M3_Cardiac_disease_other			    *	0 			    + #-0.104225698
             M3_Cardiac_valve					        *	0.256577208	+
             M3_Cerebrovascular_disease			  *	0.097803808	+
             M3_Chronic_pulmonary				      *	0.6253395	+
             M3_Chronic_renal					        *	0.334155906	+
             M3_Coagulopathy_and_other_blood	*	0.265142145	+
             M3_Colorectal_cancer				      *	0.372878764	+
             M3_Congestive_heart_failure			*	0.539809861	+
             M3_Connective_tissue				      *	0.290446442	+
             M3_Dementia							        *	1.021975368	+
             M3_Diabetes_complicated				  *	0.271607393	+
             M3_Diabetes_uncomplicated			  *	0.299383867	+
             M3_Drug_abuse						        *	0.558979499	+
             M3_Endocrine_disorder				    *	0.112673001	+
             M3_Epilepsy							        *	0.594991823	+
             M3_Eye_problem_long_term			    *	0.179923774	+
             M3_GI_ulcer_upper_GI				      *	0.152986438	+
             M3_Gynaecological_cancers			  *	0.70658858	+
             M3_Hepatitis_Chronic_viral			  *	0.569092852	+
             M3_Hypertension_uncomplicated		*	0.117746303	+
             M3_Immune_system_disorder			  *	0.398529751	+
             M3_Infection_Chronic_NOS			    *	0 			    + #-0.237983891
             M3_Inner_ear_disorder				    *	0.06090681	+
             M3_Intestinal_disorder				    *	0			      + #-0.254089697
             M3_Joint_spinal_disorder			    *	0.095585857	+
             M3_Liver_disease__moderate_or_se	*	0.474321939	+
             M3_Lung_cancer						        *	1.972481401	+
             M3_Lymphomas_and_leukaemias			*	1.190108503	+
             M3_Major_psychiatric_disorder		*	0.212789563	+
             M3_Malignant_melanoma				    *	0.342233292	+
             M3_Malnutrition_nutritional			*	0.331335106	+
             M3_Mental_and_behavioural_disord	*	0.039711074	+
             M3_Mental_retardation				    *	1.405761403	+
             M3_Metabolic_disorder				    *	0.006265195	+
             M3_Metastatic_cancer				      *	2.468586878	+
             M3_Muscular_peripheral_nerve_dis	*	0.208276284	+
             M3_Myocardial_infarction			    *	0.197491908	+
             M3_Obesity							          *	0.248243722	+
             M3_Osteoporosis_Uncomplicated		*	0.083506878	+
             M3_Other_cancers					        *	1.103452294	+
             M3_Other_neurological_disorders	*	0.564391512	+
             M3_Pancreatitis						      *	0 			    + #-0.103132585
             M3_Paralysis						          *	0.281895685	+
             M3_Peripheral_vascular				    *	0.349250005	+
             M3_Prostate_cancer					      *	0.432343447	+
             M3_Pulmonary_circulation_disorde	*	0.398432833	+
             M3_Sleep_disorder					      *	0.245749995	+
             M3_Tuberculosis						      *	0 			    + #-0.104290289
             M3_Upper_gastrointestinal_cancer	*	1.941498638	+
             M3_Urinary_tract_problem_chronic	*	0.046548658	+
             M3_Venous_insufficiency				  *	0.214050369	)
    
  
  # Finally: remove unwanted columns from M3 process as per arguments
  # If didn't want the individual-column-conditions, then drop those
  if(!return_condition_cols){
    final_M3_scored <- select(final_M3_scored, -starts_with("M3_"))
  }
  # If didn't want the M3Score, then drop it
  if(!return_score){
    final_M3_scored <- select(final_M3_scored, -M3Score)
  }
    
  # Return this final object
  final_M3_scored
   
}


# IDI prep dataset helper -------------------------------------------------

# This is a helper function to select the codes of interest from a DIAGS and EVENTS table set
# prior to passing through the code above
idiM3_prep <- function(in_diags, in_events,
                       # Defaults below are for the PUBLIC discharge data tables
                       dia_sys_code_var  = "moh_dia_clinical_sys_code",
                       dia_type_code_var = "moh_dia_diagnosis_type_code",
                       clinical_code_var = "moh_dia_clinical_code",
                       id_vars=c("snz_uid", "snz_moh_uid"), 
                       event_vars=c("moh_evt_event_id_nbr")){
  
  # Select ICD 10 codes, and just the A/B diagnosis codes
  in_diags %>% 
    filter(!!sym(dia_sys_code_var) > 6  & !!sym(dia_type_code_var) %in% c("A","B")) %>% 
    left_join(in_events %>% select(any_of(id_vars), any_of(event_vars))) %>% 
    select(any_of(id_vars), any_of(clinical_code_var)) %>% # moh_dia_condition_onset_code
    distinct()
  
}