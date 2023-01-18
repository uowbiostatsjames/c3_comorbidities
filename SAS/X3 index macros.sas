*--M3 macros for coding ICD records and assigning M3 and C3 index scores--*;
*--Release v. 1.0, 17th Jan 2023 --*;

*-- v0.33: fixed error with misnamed column in convertlong function
			that prevented proper functioning when passing a wide
			dataset for processing.
*-- v0.34: fixed issue with "lung" pre-processing for C3 index
            which mistakenly had an = rather than a - in the ICD
            code specifications

*-- v0.35: fixed C3 index error with selection of OTHER_MALIGNANCY:
             some cancers later in the list (STOMACH, LIVER)
			 were inappropriately missing out RECTAL and ANAL
			 cancers from their pre-selection list.

*-- v0.40: expanded C3 index to include "HEADNECK" and "PROSTATE" primary cancer sites
		(expansion for our use in ELEMENT project)

*-- v0.41: fixed OVARIAN CANCER being missing from selection as an Other Malignancy for some cancers
           and LUNG CANCER being inappropriately included as an Other Malignancy for PROSTATE cancer runs

		   updated some codes for M3 index for EPILEPSY and MAJOR PSYCHIATRIC DISORDER
			where e.g. F32.10 was selected, but F32.1 was not selected

*-- v1.0: 17/01/2023. Set to v. 1.0 now that methods adjusted/set alongside 
*-- v1.01: 18/01/2023 Added option OPTIONS VALIDVARNAME="V7" near start to work around Enterprise Guide naming conventions

--*;


*-- For queries regarding maintenance or errors in the code:
	email james.stanley@otago.ac.nz

	For queries regarding the development of the M3 index
	email james.stanley@otago.ac.nz

	For queries regarding the development of the C3 index
	email diana.sarfati@otago.ac.nz

	This code is based on code for the M3 index orginally written by (alphabetical order)
	James Stanley & Jane Zhang.

	Collated and amended for sharing by James Stanley, Jan/Feb 2018.	
	Partially based on code for the C3 index as written by (alphabetical order)
	Jason Gurney, Clare Salmond, Diana Sarfati, James Stanley & Jane Zhang.
	All authors University of Otago, Wellington, New Zealand.

	An earlier macro code file was available for just the C3 index.

	The C3 and Multimorbidity projects were funded by the New Zealand Health Research Council (HRC)
--*;

/*
All references below to Stanley & Sarfati are to:

Stanley, J & Sarfati, D, 2017, 
The new measuring multimorbidity index predicted mortality better than Charlson and Elixhauser indices 
among the general population. 
Journal of clinical epidemiology, 92, 99-110.
doi: 10.1016/j.jclinepi.2017.08.005
http://www.ncbi.nlm.nih.gov/pubmed/28844785

Please cite Stanley & Sarfati (2017) when using the M3 index.

All references below to Sarfati et al. 2014 are to:

Sarfati, D., Gurney, J., Stanley, J., 
Salmond, C., Crampton, P., Dennett, E., Koea, J. & Pearce, N. (2014). 
Cancer-specific administrative data-based comorbidity indices provided 
valid alternative to Charlson and National Cancer Institute Indices. 
Journal of clinical epidemiology, 67(5), 586-595. 
doi: 10.1016/j.jclinepi.2013.11.012
http://www.ncbi.nlm.nih.gov/pubmed/24582212

Please cite Sarfati et al. (2014) when using the C3 index.

ICD Codes used for Charlson and Elixhauser indices drawn from:
Quan, H., V. Sundararajan, P. Halfon, A. Fong, B. Burnand, 
J. C. Luthi, L. D. Saunders, C. A. Beck, T. E. Feasby and W. A. Ghali (2005). 
Coding algorithms for defining comorbidities in ICD-9-CM and ICD-10 administrative data.
Med Care 43(11): 1130-1139.
http://www.ncbi.nlm.nih.gov/pubmed/16224307

Charlson original paper and weights from:
Charlson, M. E., P. Pompei, K. L. Ales and C. R. MacKenzie (1987). 
A new method of classifying prognostic comorbidity in longitudinal studies: development and validation." 
J Chronic Dis 40(5): 373-383.
http://www.ncbi.nlm.nih.gov/pubmed/3558716

Elixhauser condition selection is in:
Elixhauser A, Steiner C, Harris DR, Coffey RM. 
Comorbidity measures for use with administrative data. 
Med Care. 1998;36(1):8-27.
https://www.ncbi.nlm.nih.gov/pubmed/9431328

and van Walraven developed weights for Elxihauser drawn from:
van Walraven C, Austin PC, Jennings A, Quan H, Forster AJ. 
A modification of the Elixhauser comorbidity measures into a point system for hospital death using administrative data. 
Med Care. 2009;47(6):626-33.
https://www.ncbi.nlm.nih.gov/pubmed/19433995
*/

/*
Structure and order of the file is:
1) top-level macro create_m3_index
2) sub-level macros/code, in the order that they are called
	by create_m3_index.

Users should familiarise themselves with the Create_X3_index macro
options (see below and help documentation). Those interested in
the underlying coding process can review the subsequent macros,
though these are 

As a courtesy, please consider emailing the authors if you use this code 
as the basis for considering comorbidity.
james.stanley@otago.ac.nz
*/

*---------------- create_X3_index macro------------------------------*;
/*
The following macro runs all necessary code to create the 
M3 or C3 index. 

The Charlson and Elixhauser indices can also be added to the output
(note that for the C3 index only the Charlson index is available, and scoring
of Charlson for the C3-index version ignores cancer score)
.*/
%MACRO create_X3_index(indata=, 
						data_format=, outdata= , 
						index_code = M3, 
						addCharlson = 0, addElixhauser = 0,
							
						/*Options for both C3/M3 data*/
						IDvarlist=, 
						pass_ICD_prefix=diag, pass_ICDcolname=ICDcodes,

						/*M3 index options*/
						cancer_dataM3=, 
						CancerSiteName=site, 
						CancerExtentName=extent,
						CancerExtentVal='E',

						/*C3 index options	*/
						subset_precancer=, cancersite=, 

						/*Debugging/testing options*/
						keepconditions = 1, cleanup_tempdata = 1);
/*
Arguments:

Data arguments:

indata: name of SAS data file to process (contains ICD codes and identifying variables)
		(include library: e.g. work.rawdata will use the file rawdata in the work library)
data_format: specifies if 'indata' is in a WIDE or LONG format (case-insensitive)
outdata: name of final dataset to use at end of coding and scoring process.
		(include library: e.g. work.M3_scored_data will use the file M3_scored_data in the WORK library)

index coding arguments:

index_code: C3 or M3. indicates whether to code for C3 index (see below for adding Charlson) or M3 index (one option selected only)
addCharlson:  option to add Charlson index as implemented by Quan et al (1=add, 0=leave out)
addElixauser: option to add Elixhauser index as implemented by Quan et al (1=add, 0=leave out)

specifics about source data:

IDvarlist: list of variables in 'indata' that uniquely identify patients/participants
						Can take one or multiple variables
						If using more than one variable, separate with whitespace e.g.
						IDVarList = HospitalLevelID PatientLevelID 
						if patients are numbered sequentially within hospitals.
pass_ICD_prefix (default prefix to look for is diag)
pass_ICDcolname= (default is to name as ICDcodes)

M3 specific options:

						cancer_dataM3: name of SAS data for in-period Cancer registry data (used for refined detection
						of metastasised cancer)
						CancerSiteName (default is site): name of variable in above data that includes ICD10 code site informaiton;
						CancerExtentName (default is extent): name of variable in above data that gives extent of cancer;
						CancerExtentVal (default is 'E', note that text needs to be wrapped in quotes)
										: Value of CancerExtentName variable that indicates metastasised cancer.

C3 specific options:


subset_precancer= dataset to use to look for diagnoses which are only coded pre-cancer 
cancersite: site of cancer for patients. Options are:

						values allowed: --->  maps to larger group:
						Bladder			--->  Urological
						Breast 			--->  Breast
						Colon			--->  Colorectal
						Kidney			--->  Urological
						Liver			--->  Stomach/Liver
						Ovarian			--->  Gynecological
						Rectal			--->  Colorectal
						Stomach			--->  Stomach/Liver
						Uterine			--->  Gynecological
Added for Sandar, 18/11/2015:
						Lung			--->  Lung
Added for Element project, 20/12/2021 
						HeadNeck        --->  HeadNeck
						Prostate	--->  Prostate

Debugging and testing options:
keepconditions (set this to 1 to keep the columns with the health conditions in the output data file.)
cleanup_tempdata= (set to 1 to remove temporary datasets, which are prefixed with a _C3_) 

*/

/*Firstly force some important macro variables to uppercase*/
%LET index_code = %UPCASE(&index_code); 
%LET cancersite = %UPCASE(&cancersite);
%LET data_format = %UPCASE(&data_format);

/* Just put a reminder up about what's being coded...*/
%PUT coding the &index_code index now.;

/* If site is not in list as specified...then exit with warning message.*/
/* Not using "IN" comparator here as not present in earlier versions of SAS.*/
%IF &index_code EQ C3 AND 
	NOT(&cancersite EQ BLADDER OR &cancersite EQ BREAST  OR &cancersite EQ COLON 
	 OR &cancersite EQ KIDNEY  OR &cancersite EQ LIVER   OR &cancersite EQ OVARIAN 
	 OR &cancersite EQ RECTAL  OR &cancersite EQ STOMACH OR &cancersite EQ UTERINE
	 OR &cancersite EQ LUNG    
         OR &cancersite EQ HEADNECK OR &cancersite EQ PROSTATE)

	%THEN %DO;
		%PUT C3 index requested, BUT... ;
		%PUT Cancer site &cancersite not defined in C3 index.;
		%PUT Please refer to macro for defined conditions.;
		%ABORT;
	%END;

%IF &index_code EQ M3 AND "&cancer_dataM3" EQ "" %THEN %DO;
	%PUT M3 index being coded without cancer registry data;
	%PUT (usually used to code for metastatic cancer);
	%END;

/* Macro: coding for C3 index.*/
/******************************/
%IF &index_code EQ C3 %THEN %DO;

%IF &data_format EQ WIDE %THEN %DO;
	%ConvertLong(inputdata = &indata, outputdata = _C3_converted_long,
					IDvars = &IDvarlist,
					ICD_prefix=&pass_ICD_prefix, ICDcolname=&pass_ICDcolname);
	%LET indata = _C3_converted_long;

	*will need to do same with subset_precancer_diag;
	%ConvertLong(inputdata = &subset_precancer, outputdata = _C3_converted_long_precancer,
					IDvars = &IDvarlist,
					ICD_prefix=&pass_ICD_prefix, ICDcolname=&pass_ICDcolname);
	%LET subset_precancer = _C3_converted_long_precancer;
%END;

/* Macro: Define the ICD codes to use in classifying the C3 conditions.*/
%DefineICDCodes_C3(&cancersite);

/*
CodeConditions macro (which is multi-purpose)
determines whether conditions are present based on ICD10 codes.
*/
%CodeConditions_C3( inputdata=&indata, 
					outputdata=_C3_conditionset1, 
					IDVars=&IDVarList , ICDcolname= &pass_ICDcolname,
					condition_FMT= C3_i10icd);
/* 
Some C3 conditions are defined only if detected *prior* to cancer treatment;
This needs to be supplied as a second dataset.
IT IS UP TO THE USER TO LIMIT THE &subset_cancer DATASET TO THE CORRECT TIME PERIOD.
 */
%CodeConditions_C3( inputdata=&subset_precancer, 
				    outputdata=_C3_conditionset2, 
				    IDVars=&IDVarList, ICDcolname= &pass_ICDcolname,
					condition_FMT= C3_i10icd_Pre);

/*n.b. the C3 version of Charlson skips cancer condition scoring.*/
%IF &AddCharlson EQ 1 %THEN %DO;
%CodeConditions_C3( inputdata=&indata, 
				    outputdata=_Ch_conditionset1, 
				    IDVars=&IDVarList, ICDcolname= &pass_ICDcolname,
					use_prefix = Ch_,
					condition_FMT= Ch_i10icd);

%CodeConditions_C3( inputdata=&subset_precancer, 
				    outputdata=_Ch_conditionset2, 
				    IDVars=&IDVarList, ICDcolname= &pass_ICDcolname,
					use_prefix = Ch_,
					condition_FMT= Ch_i10icd_pre);
%END;

*CodeConditions_C3 processes both sort by IDVarList;
*So no need for additional sorting here;

/*First merge in the above coded files (including Charlson if used)*/
%Merge_C3_Finalset(patient_id=&IDvarlist, passAddCharlson=&addCharlson);

%Define_C3_Weights(site=&cancersite);

%Score_C3_Conditions(_C3_joined_conditions, &outdata, site=&cancersite, retainC3conditions=&keepconditions);

/*If adding Charlson index, then passing and returning same dataset with new columns.*/
%IF &AddCharlson EQ 1 %THEN %DO;
	%Score_C3_Charlson_Conditions(&outdata, &outdata, retainChconditions=&keepconditions);
%END;

/*Now end loop for coding C3 (with or without Charlson.)*/
%END;
/*End of C3 code*/
/******************************/

/******************************/
/*New conditional section for scoring M3 index.*/
/******************************/
%IF &index_code EQ M3 %THEN %DO;

/*
If input data is in wide format (multiple diagnoses form different columns), 
then convert this to long format (multiple diagnoses sit in different rows.) 
*/
%IF &data_format EQ WIDE %THEN %DO;
	%ConvertLong(inputdata = &indata, outputdata = _M3_converted_long,
					IDvars = &IDvarlist,
					ICD_prefix=&pass_ICD_prefix, ICDcolname=&pass_ICDcolname);
	%LET indata = _M3_converted_long;
%END;

/*Load the M3 index formats for diagnostic codes.*/
%DefineICDCodes_M3();

/* First score the main M3 index output. */
%CodeConditions_M3(inputdata=&indata , patient_id =&IDvarlist, 
						clin_cd_name = &pass_ICD_prefix, condition_FMT=ICD_Cat);

/* Now check out the complications/exclusions.*/
/* Diabetes complications*/
%CodeConditions_M3(inputdata=&indata, patient_id =&IDvarlist, 
					clin_cd_name = &pass_ICD_prefix, condition_FMT=DIA_COMP);
/* Osteoporosis exclusions*/
%CodeConditions_M3(inputdata=&indata, patient_id =&IDvarlist, 
					clin_cd_name = &pass_ICD_prefix, condition_FMT=OST_EXCL);
/* Hypertension exclusions*/
%CodeConditions_M3(inputdata=&indata, patient_id =&IDvarlist, 
					clin_cd_name = &pass_ICD_prefix, condition_FMT=HYP_EXCL);

/*N.B. have separate macro for combining and resolving the exclusions. Add here and append macro.*/
/*Now combine and resolve the exclusions/complications.*/
%CombineM3Exclusions(passIDvarlist=&IDvarlist);

/*If adding in cancer registry data,
then do the following.
These macros are currently in a separate holding file.
*/
%IF "&cancer_dataM3" NE "" %THEN %DO;
	%CodeCancerM3(passDataFile=&cancer_dataM3, 
						patient_id = &IDvarlist,
						passCancerSiteName=&CancerSiteName, 
						passCancerExtentName=&CancerExtentName,
						passCancerExtentVal=&CancerExtentVal);
	%CombineCancerM3(patient_id = &IDvarlist);
%END;

/*Finally: score the M3 index.*/
%Score_M3();
PROC SORT DATA=_mm_m3_scored; BY &IDvarlist; RUN;

/*Load formats if coding Charlson or Elixhauser*/
/*n.b. one common file.*/
%IF &addCharlson EQ 1 OR &addElixhauser EQ 1 %THEN %DO;
	%DefineICDCodes_CharlElix();
%END;

/*Coding for Charlson*/
%IF &addCharlson EQ 1 %THEN %DO;
	%CodeConditions_M3(inputdata=&indata, patient_id =&IDvarlist, 
						clin_cd_name = &pass_ICD_prefix, condition_FMT=ICD_Charlson);

	%ScoreCharlson_M3();
	PROC SORT DATA=_mm_char_scored; BY &IDvarlist; RUN;
%END;

%IF &addElixhauser  EQ 1 %THEN %DO;
/* Coding for Elixhauser: note that this needs two passes
   as some ICD codes counted in multiple categories.      */
	%CodeConditions_M3(inputdata=&indata, patient_id =&IDvarlist, 
						clin_cd_name = &pass_ICD_prefix, condition_FMT=ICD_Elixhauser);
	%CodeConditions_M3(inputdata=&indata, patient_id =&IDvarlist, 
						clin_cd_name = &pass_ICD_prefix, condition_FMT=ICD_Elixhauser_B);
	/*The two above files are combined in the following macro:	*/
	%ScoreElixhauser_M3();
	PROC SORT DATA=_mm_elix_scored; BY &IDvarlist; RUN;
%END;

/*Now merge together:
(a) M3 conditions and score;
and if requested:
(b) Charlson conditions and score
(c) Elixhauser conditions and score.
*/

%merge_M3_finalset(passOutData=&outdata, patient_id=&IDvarlist, passkeepconditions=&keepconditions,
					passAddCharlson = &AddCharlson, passAddElixhauser = &AddElixhauser);

/*And the end here of the M3 index scoring code section.*/
%END;
/*End of M3 coding section    */
/******************************/

	/*If value of cleanup_tempdata in macro call is set to 1 
	  -- tidy up the intermediary datasets.*/
		%IF &cleanup_tempdata EQ 1 %THEN %DO;
		PROC DATASETS LIB=work;
    	DELETE _: index_weights;
		QUIT;
		%END;	

%MEND Create_X3_index;

/************************/
/*General macro section */
/*Used by both C3 and M3 indices*/
/************************/

*-----------------ConvertLong macro---------------------*;

/* 
The ConvertLong macro takes a wide-format file,
with one row per admission/event and multiple columns 
(one for each possible recordable diagnosis)
and returns a long format file (with one row per diagnosis)
which is suitable for processing in the next step.

*/
%MACRO ConvertLong(inputdata=, outputdata=, 
					IDvars= , ICD_prefix=diag, ICDcolname=ICDcodes);
/*
inputdata is the name of the datafile (with library specified e.g. work.mydata) to be transformed;
outputdata is the name of the output dataset
IDvars takes a list of all the relevant identifying variables to keep.
	These should be specified with no spaces... 
	Can be used to keep any variables for cross-checking this process.
ICD_prefix is the prefix of the column names for the columns which contain 
	the ICD codes to be searched for presence/absence of comorbidity.
	It is vital that all column names in this datafile start with the same 
	prefix, and that NO OTHER columns in this datafile share this prefix 
	(or they will be treated as though they could contain valid ICD codes.)
ICDcolname is the name given to the new column containing the returned ICD codes.
*/
	DATA &outputdata ;
	SET  &inputdata ;
  
	/*Create array of all variables starting with ICD_prefix (e.g. "diag*")*/
	/*Note use of colon (:) wildcard operator. */
	ARRAY ICDwide(*) &ICD_prefix: ; 
	/*Pre-assign length of ICD column output field to 8 characters.*/
	LENGTH &ICDcolname $ 8;

	/*Loop through all diagnoses.*/
	  DO i = 1 TO dim(ICDwide);
		/*For non-blank diagnostic field, store the ICD code and write out as a dataline.*/
	    IF ICDwide(i) NE "" THEN DO; 
			&ICDcolname = ICDwide(i);
	    	OUTPUT;
		END;
		/*Also keep one blank record for each row with no ICD codes, otherwise 
		we can lose patients (problem in our project where merging NZCR and C3 data.*/
		IF ICDwide(i) EQ "" AND i EQ 1 THEN DO; 
			&ICDcolname = "BLANKICD";
	    	OUTPUT;
		END;

	  END;
	/*Only keep variables listed in the IDvars macro variable*/
	KEEP &IDvars &ICDcolname;
	RUN;

	/*Fix in v 0.33 -- added as a separate step for clarity.*/ /*16/06/2022 ACTUALLY THIS SEEMS TO HAVE BROKEN/CONFUSED THINGS!*/
	DATA &outputdata ;
	SET  &outputdata ;
	&ICD_prefix = &ICDcolname;
	KEEP &IDvars &ICD_prefix ;
	RUN;


%MEND ConvertLong;

/*************************/
/* C3 specific sub-macros,
in order of main macro code
DefineICDCodes_C3
CodeConditions_C3
Define_C3_Weights
Score_C3_Conditions
Score_C3_Charlson_conditions
Merge_C3_finalset
*/
/*************************/

*------C3 Comorbidity coding: FORMAT code-------;

*--Code largely by Clare Salmond and Jason Gurney, 2011-2013--*;
*--JStanley July 2015 -- only cosmetic alterations made to this macro--*;
*-- e.g. renaming of macro variables.-*;

*----macro for the specified cancer site----;
%MACRO DefineICDCodes_C3(site);

%LET site = %UPCASE(&site);

*----The following macro variables define
	 some cancer-site specific ICD codes for the C3 cancer sites ----;

/* PrimaryMalignancy macro variable values are ICD codes pertaining
	to the patient's primary cancer site.
   OtherMalignancy are malignant cancers at other sites in the
	lookback period.
*/

*----make macro variables available to the whole programme----;
%IF &site = BLADDER %THEN %DO; 
	%LET OtherMalignancy = 
	  'C00'-'C21XX','C23'-'C26XX','C30'-'C33XX','C37'-'C39XX',
      'C43'-'C43XX','C45'-'C58XX','C60'-'C66XX','C68'-'C70XX',
      'C72'-'C75XX','C81'-'C85XX','C88'-'C88XX','C90'-'C95XX';
	  *----excludes C67, malignant neoplasm of bladder----;
	%LET PrimaryMalignancy = 'C67'-'C67XX';
  %END; 

%IF &site = BREAST %THEN %DO; 
	%LET OtherMalignancy = 
	  'C00'-'C21XX','C23'-'C26XX','C30'-'C33XX','C37'-'C39XX',
      'C43'-'C43XX','C45'-'C49XX','C51'-'C58XX','C60'-'C70XX',
      'C72'-'C75XX','C81'-'C85XX','C88'-'C88XX','C90'-'C95XX'; 
	  *----excludes C50, malignant neoplasm of breast----;
	%LET PrimaryMalignancy ='C50'-'C50XX';
  %END; 

%IF &site = COLON %THEN %DO; 
	%LET OtherMalignancy = 
	  'C00'-'C17XX','C20'-'C21XX','C23'-'C26XX','C30'-'C33XX',
      'C37'-'C39XX','C43'-'C43XX','C45'-'C58XX','C60'-'C70XX',
      'C72'-'C75XX','C81'-'C85XX','C88'-'C88XX','C90'-'C95XX';
	  *----excludes C18 and C19, malignant neoplasm of colon/rectosigmoid junction----;
	%LET PrimaryMalignancy ='C18'-'C19XX';
  %END; 

%IF &site = KIDNEY %THEN %DO; 
	%LET OtherMalignancy = 
	'C00'-'C21XX','C23'-'C26XX','C30'-'C33XX','C37'-'C39XX',
    'C43'-'C43XX','C45'-'C58XX','C60'-'C63XX','C66'-'C70XX',
    'C72'-'C75XX','C81'-'C85XX','C88'-'C88XX','C90'-'C95XX';
	  *----excludes C64, malignant neoplasm of kidney, except renal pelvis
	  		    and C65, malignant neoplasm of renal pelvis----;
	%LET PrimaryMalignancy ='C64'-'C65XX';
%END; 

%IF &site = OVARIAN %THEN %DO; 
	%LET OtherMalignancy = 
	  'C00'-'C21XX','C23'-'C26XX','C30'-'C33XX','C37'-'C39XX',
      'C43'-'C43XX','C45'-'C55XX','C57'-'C58XX','C60'-'C70XX',
      'C72'-'C75XX','C81'-'C85XX','C88'-'C88XX','C90'-'C95XX';
	  *----excludes C56, malignant neoplasm of ovary----;
	%LET PrimaryMalignancy ='C56'-'C56XX';
%END;

%IF &site = UTERINE %THEN %DO;
	%LET OtherMalignancy=
	  'C00'-'C21XX','C23'-'C26XX','C30'-'C33XX','C37'-'C39XX',
	  'C43'-'C43XX','C45'-'C52XX','C56'-'C58XX','C60'-'C70XX',
	  'C72'-'C75XX','C81'-'C85XX','C88'-'C88XX','C90'-'C95XX'; 
	%LET PrimaryMalignancy = 'C53'-'C55XX'; 
	%END;

%IF &site = ALT_RECTAL %THEN %DO; *remove Alt_ if looking at notes review cohorts;
	%LET OtherMalignancy=
	  'C00'-'C19XX','C21'-'C21XX','C23'-'C26XX','C30'-'C33XX',
	  'C37'-'C39XX','C43'-'C43XX','C45'-'C58XX','C60'-'C70XX',
	  'C72'-'C75XX','C81'-'C85XX','C88'-'C88XX','C90'-'C95XX';
	%LET PrimaryMalignancy = 'C20'-'C20XX'; 
	%END;

%IF &site = RECTAL %THEN %DO; *remove Alt_ if looking at notes review cohorts;
	%LET OtherMalignancy=
	  'C00'-'C19XX','C21'-'C21XX','C23'-'C26XX','C30'-'C33XX',
	  'C37'-'C39XX','C43'-'C43XX','C45'-'C58XX','C60'-'C70XX',
	  'C72'-'C75XX','C81'-'C85XX','C88'-'C88XX','C90'-'C95XX';
	%LET PrimaryMalignancy = 'C20'-'C20XX'; 
	%END;

%IF &site = ALT_STOMACH %THEN %DO;
	%LET OtherMalignancy=
	  'C00'-'C15XX','C17'-'C21XX','C23'-'C26XX','C30'-'C33XX',
	  'C37'-'C39XX','C43'-'C43XX','C45'-'C58XX','C60'-'C70XX',
	  'C72'-'C75XX','C81'-'C85XX','C88'-'C88XX','C90'-'C95XX';
	%LET PrimaryMalignancy = 'C16'-'C16XX'; 
	%END;

%IF &site = STOMACH %THEN %DO;
	%LET OtherMalignancy=
	  'C00'-'C15XX','C17'-'C21XX','C23'-'C26XX','C30'-'C33XX',
	  'C37'-'C39XX','C43'-'C43XX','C45'-'C58XX','C60'-'C70XX',
	  'C72'-'C75XX','C81'-'C85XX','C88'-'C88XX','C90'-'C95XX';
	%LET PrimaryMalignancy = 'C16'-'C16XX'; 
	%END;

%IF &site = ALT_LIVER %THEN %DO;
	%LET OtherMalignancy=
	  'C00'-'C21XX','C23'-'C26XX','C30'-'C33XX',
	  'C37'-'C39XX','C43'-'C43XX','C45'-'C58XX','C60'-'C70XX',
	  'C72'-'C75XX','C81'-'C85XX','C88'-'C88XX','C90'-'C95XX';
	%LET PrimaryMalignancy = 'C22'-'C22XX'; 
	%END;

%IF &site = LIVER %THEN %DO;
	%LET OtherMalignancy=
	  'C00'-'C21XX','C23'-'C26XX','C30'-'C33XX',
	  'C37'-'C39XX','C43'-'C43XX','C45'-'C58XX','C60'-'C70XX',
	  'C72'-'C75XX','C81'-'C85XX','C88'-'C88XX','C90'-'C95XX';
	%LET PrimaryMalignancy = 'C22'-'C22XX'; 
	%END;

/*Cancers added after original C3 index creation*/;

%IF &site = LUNG %THEN %DO;
	%LET OtherMalignancy=
	  'C00'-'C21XX','C23'-'C26XX','C30'-'C33XX',
	  'C37'-'C39XX','C43'-'C43XX','C45'-'C58XX','C60'-'C70XX',
	  'C72'-'C75XX','C81'-'C85XX','C88'-'C88XX','C90'-'C95XX'; 
	%LET PrimaryMalignancy = 'C34'-'C34XX'; 
	%END;

%IF &site = PROSTATE %THEN %DO;
	%LET OtherMalignancy=
	  'C00'-'C21XX','C23'-'C26XX','C30'-'C33XX',
      'C37'-'C39XX','C43'-'C43XX','C45'-'C58XX','C60'-'C60XX', 'C62'-'C70XX',
	  'C72'-'C75XX','C81'-'C85XX','C88'-'C88XX','C90'-'C95XX'; 
	%LET PrimaryMalignancy = 'C61'-'C61XX'; 
	%END;

%IF &site = HEADNECK %THEN %DO;
	%LET OtherMalignancy=
	  'C15'-'C21XX','C23'-'C26XX','C33'-'C33XX',
      'C37'-'C39XX','C43'-'C43XX','C45'-'C58XX','C60'-'C70XX',
	  'C72'-'C75XX','C81'-'C85XX','C88'-'C88XX','C90'-'C95XX';
	%LET PrimaryMalignancy = 'C00'-'C14XX', 'C30'-'C32XX'; 
	%END;

*----Comorbidity list and related ICD-10 dignosis codes for cancer as contained in the
     version of 'Comorbidity list (Sarfati list)' dated 26 March 2012. 

     NOTES: 
	 The comorbidities are in the order originally listed in the Sarfati comorbidity list
     (connected with the literature sources). Comorbidities apparently alcohol-induced  
	 have been separated from those not alcohol-induced, so that the two variable names have 
 	 'AlcInd' or 'not_AlcInd' at the end of the variable names. The variable pairs can be 
	 added together after examination in a further data step, if required. They are coded 
	 with numbers up to 4-digits, the first one or two being the comorbidity code, and the 
	 final two being '00' for 'not AlcInd' and '32' for 'AlcInd' (since the code '32' is for
	 alcohol abuse). These 'AlcInd' and 'not_AlcInd' codes will sort together numerically. 
	 The comorbidity groups are later alphabetised for ease of reading the long output. The 
	 comorbiidity names now have the relevant organ or system name as their first part. The 
     malignancy names (currently analysed, or extra) are '14' with extensions '88' and '00'.
	 The comorbidity that had been coded 36 ('Gall bladder problems') has been dropped.
	 The original comorbidity '6' has been split in two: '61' and '62'----;

PROC FORMAT CNTLOUT=work._C3_lookup_ICD ;
INVALUE C3_i10icd
'I70'-'I71XX','I720'-'I720X','I731'-'I731X','I738'-'I739X','I771'-'I771X',
											 'K551'-'K552X','K558'-'K559X' =3                             
							  'I830'-'I830X','I832'-'I832X','I872'-'I872X' =4
							   'G45'-'G46XX', 'I60'-'I67XX', 'I69'-'I69XX' =5
 'F00'-'F01XX','F020'-'F023X', 'F03'-'F03XX','F051'-'F051X', 'G30'-'G311X' =6
 'F04'-'F04XX', 'F06'-'F06XX', 'F07'-'F071X','F078'-'F079X', 'F09'-'F09XX',
															'G931'-'G931X' =9999 /* was prev 7*/
 'E84'-'E84XX', 'J40'-'J44XX', 'J47'-'J47XX', 'J60'-'J67XX','J684'-'J684X',
'J701'-'J701X','J703'-'J703X', 'J84'-'J84XX','J961'-'J961X','J980'-'J980X',
											  'J982'-'J984X','J45'-'J46XX' =8 /*(7=chronic pulmonary disease combined with 35=asthma)*/															 
 'L93'-'L93XX', 'M05'-'M06XX', 'M08'-'M08XX','M120'-'M120X','M123'-'M123X',
							   'M30'-'M34XX','M350'-'M356X','M358'-'M359X' =9
'K220'-'K221X','K224'-'K225X','K228'-'K229X', 'K25'-'K28XX','K311'-'K312X',
											 'K314'-'K314X','K316'-'K316X' =10
'E100'-'E101X','E109'-'E109X','E110'-'E111X','E119'-'E119X','E120'-'E121X',
'E129'-'E129X','E130'-'E131X','E139'-'E139X','E140'-'E141X','E149'-'E149X' =11														        
'E102'-'E108X','E112'-'E118X','E122'-'E128X','E132'-'E138X','E142'-'E148X' =12  
'G041'-'G041X','G114'-'G114X','G800'-'G802X', 'G81'-'G82XX','G830'-'G834X',
															'G839'-'G839X' =13
'I120'-'I120X','I131'-'I131X','N032'-'N039X','N042'-'N049X','N052'-'N059X',
 'N11'-'N11XX', 'N18'-'N19XX','N250'-'N250X','N258'-'N259X','Z490'-'Z49XX',
											 'Z940'-'Z940X','Z992'-'Z992X' =14
    													  &PrimaryMalignancy = 9999 /* was prev 16*/
'K711'-'K711X','K713'-'K715X','K717'-'K717X','K721'-'K721X','K729'-'K729X', 
 'K73'-'K74XX','K760'-'K760X','K762'-'K769X', 'I85'-'I85XX','I864'-'I864X',
						  	   'I982'-'I982X','Z944'-'Z944X','K70'-'K70XX' =17 /*1500 combined with 1532*/
 							   'B20'-'B24XX','F024'-'F024X', 'Z21'-'Z21XX' =9999 /* was prev 18*/
 															 'I20'-'I20XX' =19 
 'I05'-'I08XX','I091'-'I091X','I098'-'I098X', 'I34'-'I38XX','T820'-'T820X',
			   				  'Q230'-'Q233X','Q238'-'Q239X','Z952'-'Z954X' =23
							   'K50'-'K51XX','K522'-'K522X','K528'-'K529X' =24
 'G10'-'G10XX','G110'-'G113X','G118'-'G119X', 'G12'-'G13XX', 'G20'-'G21XX',
 'G23'-'G23XX','G255'-'G255X','G318'-'G319X', 'G35'-'G37XX', 'G90'-'G90XX',
 							  'G934'-'G934X','R470'-'R470X','G312'-'G312X' =25 /*2300 combined with 2332*/
							  'G400'-'G404X','G406'-'G409X', 'G41'-'G41XX' =26
 'G60'-'G61XX','G620'-'G620X','G622'-'G622X','G628'-'G629X', 'G64'-'G64XX', 
 'G70'-'G71XX','G720'-'G720X','G722'-'G724X','G728'-'G729X','G731'-'G731X',
                                             'G621'-'G621X','G721'-'G721X' =27 /*2500 combined with 2532 */
'F302'-'F302X', 'F31'-'F31XX','F321'-'F323X','F328'-'F329X', 'F33'-'F33XX', 
  'F39'-'F39XX','F20'-'F20XX', 'F22'-'F22XX', 'F25'-'F25XX', 'F28'-'F29XX' =28 /* 'Major psychiatric conditions (bipolar or depressive illness)' combined with 'Major psychiatric conditions (schizophernia and psychosis)'*/
 'D55'-'D58XX','D590'-'D594X','D598'-'D599X', 'D60'-'D61XX', 'D64'-'D64XX',
 'D66'-'D67XX','D680'-'D682X','D688'-'D689X','D691'-'D694X','D696'-'D696X',
'D698'-'D699X', 'D70'-'D72XX', 'D74'-'D74XX','D750'-'D750X','D752'-'D752X',
															'D758'-'D759X' =30
											 				 'E66'-'E66XX' =32
			   'F101'-'F109X','K292'-'K292X','Z502'-'Z502X','Z714'-'Z714X' =33
 'F11'-'F16XX', 'F18'-'F19XX','Z503'-'Z503X','Z715'-'Z715X','Z722'-'Z722X' =9999 /* was prev 34*/
							  'K861'-'K861X','K868'-'K868X','K860'-'K860X' =9999 /* was prev 35*/
							  													 /*3400 combined with 3432 */
 'E01'-'E03XX', 'E05'-'E05XX','E062'-'E063X','E065'-'E065X', 'E07'-'E07XX',
'E163'-'E164X','E168'-'E169X', 'E20'-'E20XX','E210'-'E210X','E212'-'E215X',
 'E22'-'E22XX','E230'-'E230X','E232'-'E233X','E236'-'E237X','E240'-'E241X',
'E243'-'E243X','E248'-'E249X', 'E25'-'E27XX', 'E31'-'E32XX','E345'-'E345X',		   
											 'E244'-'E244X','E348'-'E349X' =36 /*3700 combined with 3732*/
							  'N301'-'N302X', 'N31'-'N32XX', 'N35'-'N36XX' =37
 											  'A15'-'A19XX', 'B90'-'B90XX' =9999 /* was prev 38*/
 'M80'-'M80XX','M810'-'M811X','M815'-'M815X','M818'-'M819X','M831'-'M835X',
			   'M838'-'M839X', 'M85'-'M85XX','M863'-'M866X', 'M88'-'M88XX' =39											 
							   'D80'-'D84XX', 'D86'-'D86XX', 'D89'-'D89XX' =9999 /* was prev 40*/
 'E70'-'E72XX', 'E74'-'E78XX','E791'-'E791X','E798'-'E799X', 'E80'-'E80XX', 
							   'E83'-'E83XX', 'E85'-'E85XX', 'E88'-'E88XX' =41
'E000'-'E002X','E009'-'E009X', 'F70'-'F73XX', 'F78'-'F79XX','F842'-'F844X',
															 'Q90'-'Q90XX' =9999 /* was prev 42*/
							   'B18'-'B18XX','B942'-'B942X','Z225'-'Z225X' =43
											  'F51'-'F51XX','G470'-'G473X' =44
 'H80'-'H81XX', 'H83'-'H83XX', 'H90'-'H90XX','H910'-'H911X','H913'-'H913X',
											 'H918'-'H919X','H930'-'H933X' =45
 'A30'-'A31XX', 'A52'-'A52XX', 'B91'-'B92XX','B941'-'B941X','B948'-'B949X' =9999 /* was prev 46*/
 'H16'-'H16XX','H181'-'H181X','H184'-'H186X','H201'-'H201X','H212'-'H212X',
'H301'-'H301X','H311'-'H314X','H330'-'H330X','H332'-'H335X', 'H34'-'H35XX', 
 'H43'-'H43XX', 'H46'-'H47XX', 'H49'-'H49XX', 'H50'-'H51XX','H530'-'H534X',
			   'H536'-'H536X','H538'-'H539X', 'H54'-'H54XX', 'Q12'-'Q15XX' =48
'I248'-'I249X','I250'-'I251X','I253'-'I254X','I256'-'I256X','I258'-'I259X',
							  'I310'-'I311X','I421'-'I422X','I424'-'I424X' =49											
							   'K57'-'K57XX','K592'-'K593X', 'K90'-'K90XX' =50
 'M07'-'M07XX', 'M13'-'M13XX','M150'-'M152X','M154'-'M154X','M158'-'M159X',
'M400'-'M400X','M402'-'M405X', 'M41'-'M43XX', 'M45'-'M45XX','M460'-'M462X',
 'M47'-'M47XX','M480'-'M482X','M485'-'M485X','M488'-'M489X','G950'-'G951X' =51
																	  ' '  =.
																     Other =9999;

INVALUE C3_i10icd_Pre
                               'I21'-'I23XX','I241'-'I241X','I252'-'I252X' =1
'I099'-'I099X','I110'-'I110X','I130'-'I130X','I132'-'I132X','I255'-'I255X',
'I420'-'I420X','I425'-'I425X','I427'-'I429X', 'I43'-'I43XX', 'I50'-'I50XX',
                                                            'I426'-'I426X' =2 /* 200='Congestive heart failure (not alcohol-induced)' and 232='Congestive heart failure (alcohol-induced)' combined*/                                                
										   &OtherMalignancy =15
 				'I10'-'I10XX','I119'-'I119X','I129'-'I129X','I139'-'I139X' =20
'I441'-'I443X','I456'-'I456X','I459'-'I459X', 'I47'-'I49XX','T821'-'T821X',
											 'Z450'-'Z450X','Z950'-'Z950X' =21
 							   'I26'-'I27XX','I280'-'I281X','I288'-'I289X' =22
 'F40'-'F42XX', 'F44'-'F45XX', 'F48'-'F48XX', 'F50'-'F50XX', 'F55'-'F55XX',
 							   'F59'-'F61XX', 'F63'-'F66XX', 'F68'-'F69XX' =29
 															 'D50'-'D53XX' =31
 'E40'-'E46XX', 'E50'-'E56XX', 'E58'-'E59XX', 'E60'-'E61XX', 'E63'-'E64XX' =47
 																	  ' '  =.
																	  ''  = 9999
																     Other =9999;

*Define a format to apply to the in-formatted ICD codes:
 allows clearer identification of which condition is which;

/* JS July 2015: reconciled names below with names in C3 paper.*/
/* Commented out names are conditions not included in the C3 index.*/

VALUE C3_fcomorb
  3='Peripheral vascular disease'
  4='Venous insufficiency'
  5='Cerebrovascular disease'
  6='Dementia'
/*  7='Mental etc disorders from brain damage'*/
  8='COPD and asthma'
  9='Connective tissue disorders'
  10='Upper GI disorders'
  11='Diabetes no complications'
  12='Diabetes with complications'
  13='Paralysis'
  14='Renal disease'
/*  16='Primary malignancy as in title'*/
  17='Liver disease moderate or severe'
/*  18='AIDS'*/
  19='Angina'
  23='Cardiac valve disorders'
  24='Inflammatory bowel disease'
  25='Neurological excl epilepsy'
  26='Epilepsy'
  27='Peripheral nerve or muscular disorder'
  28='Major psychiatric disorders'
  30='Coagulopathies and other blood disorders'
  32='Obesity'
  33='Alcohol abuse'
/*  34='Drug abuse'*/
/*  35='Pancreatitis (alcohol and not alcohol induced)'*/
  36='Endocrine disorders'
  37='Urinary tract disorder'
/*  38='TB'*/
  39='Osteoporosis and bone disorders'
/*  40='Immune system disorders'*/
  41='Metabolic disorder'
/*  42='Mental retardation'*/
  43='Chronic viral hepatitis'
  44='Sleep disorder' 
  45='Inner ear disorders'
/*  46='Chronic infection NOS'*/
  48='Eye problems'
  49='Other cardiac conditions'
  50='Intestinal disorders'
  51='Joint and spinal disorders'

9999='Other'
;

VALUE C3_fcomorb_Pre
  1='Myocardial infarction'
  2='Congestive heart failure'
  15='Other malignancy'
  20='Hypertension'
  21='Cardiac arrhythmia'
  22='Pulmonary circulation disorder'
  29='Anxiety and behavioral disorders'
  31='Anemia'
  47='Nutritional disorders'
;

/*
Codes here are from Quan et al.'s mapping of 16 Charlson comorbid conditions to ICD-10 codes
Reference:
Quan, H., V. Sundararajan, P. Halfon, A. Fong, B. Burnand, 
J. C. Luthi, L. D. Saunders, C. A. Beck, T. E. Feasby and W. A. Ghali (2005). 
Coding algorithms for defining comorbidities in ICD-9-CM and ICD-10 administrative data.
Med Care 43(11): 1130-1139.
*/

INVALUE Ch_I10ICD
'I70'-'I71XX','I731'-'I731X','I738'-'I739X','I771'-'I771X','I790'-'I790X',
'I792'-'I792X','K551'-'K551X','K558'-'K559X','Z958'-'Z959X'
=3
'G45'-'G46XX','H340'-'H340X','I60'-'I69XX'
=4
'F00'-'F03XX','F051'-'F051X','G30'-'G30XX','G311'-'G311X'
=5
'I278'-'I279X','J40'-'J47XX','J60'-'J67XX','J684'-'J684X','J701'-'J701X','J703'-'J703X'
=6
'M05'-'M06XX','M315'-'M315X','M32'-'M34XX','M351'-'M351X','M353'-'M353X','M360'-'M360X'
=7
'K25'-'K28XX'
=8
'B18'-'B18XX','K700'-'K703X','K709'-'K709X','K713'-'K715X','K717'-'K717X',
'K73'-'K74XX','K760'-'K760X','K762'-'K764X','K768'-'K769X','Z944'-'Z944X'
=9
'E100'-'E101X','E106'-'E106X','E108'-'E111X','E116'-'E116X','E118'-'E121',
'E126'-'E126X','E128'-'E131X','E136'-'E136X','E138'-'E141X','E146'-'E146X',
'E148'-'E149X'
=10
'G041'-'G041X','G114'-'G114X','G801'-'G802X','G81'-'G82XX','G830'-'G834X',
'G839'-'G839X'
=11
'I120'-'I120X','I131'-'I131X','N032'-'N037X','N052'-'N057X','N18'-'N19XX',
'N250'-'N250X','Z490'-'Z492X','Z940'-'Z940X','Z992'-'Z992X'
=12
'E102'-'E105X','E107'-'E107X',
'E112'-'E115X','E117'-'E117X',
'E122'-'E125X','E127'-'E127X',
'E132'-'E135X','E137'-'E137X',
'E142'-'E145X','E147'-'E147X'
=13
&OtherMalignancy
=14
'I850'-'I850X','I859'-'I859X','I864'-'I864X','I982'-'I982X','K704'-'K704X',
'K711'-'K711X','K721'-'K721X','K729'-'K729X','K765'-'K767X'
=15
'B20'-'B22XX','B24'-'B24XX'
=16
' '=.
Other=9999
;

INVALUE Ch_I10ICD_Pre
'I21'-'I22XX','I252'-'I252X'
=1
'I099'-'I099X','I110'-'I110X','I130'-'I130X','I132'-'I132X',
'I255'-'I255X','I420'-'I420X','I425'-'I425X','I426'-'I426X','I427'-'I427X',
'I428'-'I428X','I429'-'I429X','I43'-'I43XX','I50'-'I50XX','P290'-'P290X'
=2
' '=.
Other=9999
;

/*JS: Dropped Charlson_ prefix from the following to work with new code.*/
VALUE Ch_fcomorb
3='PVD'
4='CBVD'
5='Dementia'
6='CPD'
7='RD'
8='PUD'
9='LiverMild'
10='DiabNoComp'
11='HemiPara'
12='Renal'
13='DiabWithComp'
14='Cancer'
15='LiverModSevere'
16='AIDS'
9999='Other'
;

VALUE Ch_fcomorb_Pre
1='MI'
2='CHF'
;

RUN;

%MEND DefineICDCodes_C3;

*-----------------CodeConditions_c3 macro---------------------*;

/*	Adapted by James Stanley, June/July 2015,*/
/*	from Multimorbidity project code developed*/
/*	by JS and Jane Zhang, Jan-May 2015 */

%MACRO CodeConditions_C3(inputdata=, outputdata=, 
						 condition_FMT=C3_i10icd,  use_prefix = C3_,
						 IDvars=, ICDcolname=);

/*
inputdata is name of the datafile (with library included e.g. work.mydata) to be transformed;
outputdata is the name of the output dataset

condition_FMT is the name of the format which is to be used to search
	for the appropriate ICD codes for each condition.

use_prefix is to differentiate between C3 columns (prefix is C3_)
	and Charlson columns (prefix is Ch_ -- must be requested!)

IDvars takes a list of all the relevant identifying variables to keep.
	These should be specified with no spaces... 
	Can be used to keep any variables for cross-checking this process.

ICDcolname is the name of the single column containing all ICD codes.

*/

/* Firstly, load the appropriate formats for this call.*/
/* Time-invariant conditions (any time in study lookback period)*/
%IF "&Condition_FMT" EQ "C3_i10icd" %THEN %DO;
	%LET text_format = C3_fcomorb;

	PROC SORT DATA = work._C3_lookup_ICD 
			   OUT = work._C3_lookup_formats (KEEP=label rename=(label=C3_Cat)) NODUPKEY;
	WHERE FMTName='C3_FCOMORB';
	BY Label ;
	RUN;
%END;

/* Time-dependent conditions (must be prior to cancer treatment.)*/
%IF "&Condition_FMT" EQ "C3_i10icd_Pre" %THEN %DO;
	%LET text_format = C3_fcomorb_Pre;
	PROC SORT DATA = work._C3_lookup_ICD  
			   OUT = work._C3_lookup_formats (KEEP=label rename=(label=C3_Cat)) NODUPKEY;
	WHERE FMTName='C3_FCOMORB_PRE';
	BY Label ;
	RUN;
%END;

/*Charlson: n.b. keep label rename as C3_Cat to make life easier. */
%IF "&Condition_FMT" EQ "Ch_i10icd" %THEN %DO;
	%LET text_format = Ch_fcomorb;
	PROC SORT DATA = work._C3_lookup_ICD  
			   OUT = work._C3_lookup_formats (KEEP=label rename=(label=C3_Cat)) NODUPKEY;
	WHERE FMTName='CH_FCOMORB';
	BY Label ;
	RUN;
%END;

%IF "&Condition_FMT" EQ "Ch_i10icd_pre" %THEN %DO;
	%LET text_format = Ch_fcomorb_Pre;
	PROC SORT DATA = work._C3_lookup_ICD  
			   OUT = work._C3_lookup_formats (KEEP=label rename=(label=C3_Cat)) NODUPKEY;
	WHERE FMTName='CH_FCOMORB_PRE';
	BY Label ;
	RUN;
%END;

/*Make a temporary data file with the called formats.*/
	DATA _C3_lookup_formats;
	SET  _C3_lookup_formats;
	IF C3_Cat NE 'Other';
	RUN;

/* Now to begin coding the data...*/
/* This code works by applying formats to data. */
	DATA _C3_diag_cat;
		SET &inputdata;
		LENGTH ICD_status 3;
		/*Apply input formats to data, also store as text.*/
		C3_Coded = INPUT(&ICDcolname, &condition_FMT..);
		C3_Cat  = PUT(C3_Coded, &text_format..);
		FORMAT C3_Coded &text_format..;
		/*Add a flag to indicate that this condition is present*/
		ICD_status=1;
		/*Drop row if no C3 condition match found in diagnosis formats.*/
		IF C3_Coded EQ 9999 THEN DELETE;
		/*Keep relevant variables only.*/
		KEEP &IDvars C3_Coded C3_Cat ICD_status;
	RUN;

	/* Remove duplicates from data such that the 
	processing file only has one record per comorbidity category per person*/

	PROC SORT DATA=_C3_diag_cat OUT=_C3_diag_cat_uni (KEEP=&IDvars C3_Cat ICD_status ) NODUPKEY;
		BY  &IDvars C3_Coded;
	RUN;

	/* Generate an empty table with all potential comorbidity categories listed*/
	/* with one row per combination of each unique person in original dataset*/
	/* with each C3 condition.*/
	PROC SORT DATA=&inputdata OUT=_C3_ID_uni (KEEP= &IDvars) NODUPKEY;
		BY &IDvars;
	RUN;

	PROC SQL;
		CREATE TABLE _C3_ID_empty AS
			SELECT * FROM _C3_ID_uni
			CROSS JOIN 
			      work._C3_lookup_formats;
	QUIT;
	
	/*Sort the empty table and the diagnoses recorded tables*/
	PROC SORT DATA = _C3_ID_empty;
		BY &IDvars C3_Cat;
	RUN;

	PROC SORT DATA = _C3_diag_cat_uni;
		BY &IDvars C3_Cat;
	RUN;

	/*Merge these two files together to populate table.*/
	DATA _C3_merged_long;
		MERGE _C3_ID_empty _C3_diag_cat_uni;
		BY &IDVars C3_Cat;
		IF ICD_status=. THEN ICD_status=0;
	RUN;

	/*Sort this (possibly superfluous, here as a safety step.)*/
	PROC SORT DATA = _C3_merged_long;
		BY &IDvars C3_Cat;
	RUN;

	/* UPDATE v1.01: Don't allow spaces in variable names! Default in Enterprise Guide 
	is for more flexible name allowances, which breaks the code */
	%LET VALVARSET = %sysfunc(getoption(VALIDVARNAME)); 
	OPTIONS VALIDVARNAME=V7;
	/* Transpose the resulting table: this takes it from a long table
	with one row per person per condition...
	to a wide table, with one row per person, and one column per condition. */
	PROC TRANSPOSE DATA = _C3_merged_long 
		OUT=&outputdata(drop=_name_ compress=binary) PREFIX=&use_prefix;
		ID  C3_Cat; /*names is too long, not good for following processing?*/
		VAR ICD_status;
		BY  &IDvars;
	RUN;
	* Reset the naming option again here;
	OPTIONS VALIDVARNAME=&VALVARSET;
	
	/*Preceding steps may already have handled sorting... just in case...*/
	PROC SORT DATA = &outputdata; BY &IDVars; RUN;

	/*Pop back into data to apply some exclusions here.*/
	%IF &text_format EQ C3_fcomorb OR &text_format EQ Ch_fcomorb %THEN %DO;
		DATA &outputdata;
		SET &outputdata;
		*Make diabetes groups mutually exclusive;
		%IF &text_format EQ C3_fcomorb %THEN %DO;
			IF C3_Diabetes_With_Complications EQ 1 THEN C3_Diabetes_No_Complications = 0;
		%END;

		%IF &text_format EQ Ch_fcomorb %THEN %DO;
			IF Ch_DiabWithComp   EQ 1 THEN Ch_DiabNoComp = 0;
			IF Ch_LiverModSevere EQ 1 THEN Ch_LiverMild = 0;
		%END;
		RUN;
	%END;


%MEND CodeConditions_C3;

*-----------------Define index_weights file ------------------------*;

/*
Define weights to be used in index creation...
*/

%MACRO Define_C3_Weights(site=);

%LET site = %UPCASE(&site);

/*
This macro makes a dataset that contains the C3 index weight components
to be used in creating the index.

The only argument this function takes is the cancer site being processed;
It outputs a dataset called index_weights to the work directory.

site must be specified as some conditions are excluded as 
too closely related to primary cancer e.g. 
Renal_disease for Urological cancer patients (related to primary cancer site)

Index component names are prefixed by
_a_ for all sites
*/

/*Create datafile.*/
DATA index_weights;
/*All sites weights*/
_a_Alcohol_abuse = 1.08;
_a_Anemia = 0.59;
_a_Angina = 0.51;
_a_Anxiety_and_behavioral_disord = 0.57;
_a_Cardiac_arrhythmia = 0.77;
_a_Cardiac_valve_disorders = 1.1;
_a_Cerebrovascular_disease = 1.09;
_a_Chronic_viral_hepatitis = 0.39;
_a_Coagulopathies_and_other_bloo = 0.75;
_a_Congestive_heart_failure = 1.26;
_a_Connective_tissue_disorders = 0.51;
_a_COPD_and_asthma = 1.09;
_a_Dementia = 1.35;
_a_Diabetes_with_complications = 0.88;
_a_Diabetes_no_complications = -0.03;
_a_Endocrine_disorders = 0.77;
_a_Epilepsy = 1.04;
_a_Eye_problems = 0.63;
_a_Upper_GI_disorders = 0.11;
_a_Hypertension = 0.72;
_a_Inflammatory_bowel_disease = 0.52;
_a_Inner_ear_disorders = 0.54;
_a_Intestinal_disorders = 0.11;
_a_Joint_and_spinal_disorders = 0.69;
_a_Liver_disease_moderate_or_sev = 0.92;
_a_Major_psychiatric_disorders = 0.79;
_a_Nutritional_disorders = 1.16;
_a_Metabolic_disorder = 0.61;
_a_Myocardial_infarction = 0.93;
_a_Neurological_excl_epilepsy = 1.06;
_a_Obesity = 0.83;
_a_Osteoporosis_and_bone_disorde = 0.49;
_a_Other_cardiac_conditions = 0.62;
_a_Other_malignancy = 0.17;
_a_Paralysis = 1.03;
_a_Peripheral_nerve_or_muscular = 1.2;
_a_Peripheral_vascular_disease = 0.98;
_a_Pulmonary_circulation_disorde = 0.95;
_a_Renal_disease = 1.38;
_a_Sleep_disorder = 1.41;
_a_Urinary_tract_disorder = 0.12;
_a_Venous_insufficiency = 0.7;

/* 
Over-ruling of all-sites weights for cancer-related condtions:
Note that no such over-ruling applied to gynaecological cancers.
*/

%IF &site EQ COLON OR  &site EQ RECTAL %THEN %DO;
	_a_Anemia = 0;
	_a_Intestinal_disorders = 0;
%END;

%IF &site EQ LIVER OR  &site EQ STOMACH %THEN %DO;
	_a_Anemia = 0; 
	_a_Coagulopathies_and_other_bloo = 0;
	_a_Upper_GI_disorders = 0;
	_a_Liver_disease_moderate_or_sev = 0;
%END;

%IF &site EQ KIDNEY OR  &site EQ BLADDER %THEN %DO;
	_a_Renal_disease = 0; 
	_a_Urinary_tract_disorder = 0;
%END;

OUTPUT;
RUN;

%MEND define_c3_weights;

*-----------------End of definition of index_weights file ------------------------*;

*----------------------score_C3_conditions macro----------------------*;

%MACRO Score_C3_Conditions(inputdata, outputdata, site=NULL, keepweights = 0, retainC3conditions=1);

/*
Takes file with one row per patient/participant 
and one column per C3 condition (as automatically named in the 
earlier processing steps.)

Merges in the C3 study weights (stored in index_weights) as listed in:
Sarfati et al. 2014
*/

/* Arguments:
inputdata: library/file with one row per patient, one column per C3 condition.
			n.b. as produced by CodeConditions_C3 macro.
			default library/file as per call in the create_C3_index macro is 
			work._C3_conditionset1 (or _C3_conditionset2)
outputdata: library/file with output file name. No default.
			e.g. work.final_scored_patients OR finaldat.all_patients_c3_scores
site: name of cancer site for patient group. 
		Usually passed through from call to create_C3_index macro.
keepweights: option to keep weighting columns in output dataset.
				Default is to drop these columns.
retainC3conditions: option to keep binary condition columns in output dataset.
					Default is to keep these columns.
					Note that create_C3_index macro can override this default
					with keepconditions option set to 1 (which is default there!)
*/

	/* Load input data.*/
	DATA &outputdata;
	SET  &inputdata ;

		/*
		Merge in weights to each line.
		Coding tip from http://stackoverflow.com/q/9790868/1834244
		*/
		IF _N_ EQ 1 THEN SET index_weights;

		/* Change site variable to uppercase.*/
		%LET site = %UPCASE(&site);

		/* Set up array of C3 condition names*/
		ARRAY conditions (*) 
			C3_Alcohol_abuse
			C3_Angina
			C3_COPD_and_asthma
			C3_Cardiac_valve_disorders
			C3_Cerebrovascular_disease
			C3_Chronic_viral_hepatitis
			C3_Coagulopathies_and_other_bloo
			C3_Connective_tissue_disorders
			C3_Dementia
			C3_Diabetes_no_complications
			C3_Diabetes_with_complications
			C3_Endocrine_disorders
			C3_Epilepsy
			C3_Eye_problems
			C3_Inflammatory_bowel_disease
			C3_Inner_ear_disorders
			C3_Intestinal_disorders
			C3_Joint_and_spinal_disorders
			C3_Liver_disease_moderate_or_sev
			C3_Major_psychiatric_disorders
			C3_Metabolic_disorder
			C3_Neurological_excl_epilepsy
			C3_Obesity
			C3_Osteoporosis_and_bone_disorde
			C3_Other_cardiac_conditions
			C3_Paralysis
			C3_Peripheral_nerve_or_muscular
			C3_Peripheral_vascular_disease
			C3_Renal_disease
			C3_Sleep_disorder
			C3_Upper_GI_disorders
			C3_Urinary_tract_disorder
			C3_Venous_insufficiency
			C3_Anemia
			C3_Anxiety_and_behavioral_disord
			C3_Cardiac_arrhythmia
			C3_Congestive_heart_failure
			C3_Hypertension
			C3_Myocardial_infarction
			C3_Nutritional_disorders
			C3_Other_malignancy
			C3_Pulmonary_circulation_disorde
			;

		/* 
			Set up array of names for the allsites_weights.
			These were entered manually to ensure column alignment with above
			list -- DO NOT ALTER!
		*/
		ARRAY allsites_weights (*)
			_a_Alcohol_abuse
			_a_Angina
			_a_COPD_and_asthma
			_a_Cardiac_valve_disorders
			_a_Cerebrovascular_disease
			_a_Chronic_viral_hepatitis
			_a_Coagulopathies_and_other_bloo
			_a_Connective_tissue_disorders
			_a_Dementia
			_a_Diabetes_no_complications
			_a_Diabetes_with_complications
			_a_Endocrine_disorders
			_a_Epilepsy
			_a_Eye_problems
			_a_Inflammatory_bowel_disease
			_a_Inner_ear_disorders
			_a_Intestinal_disorders
			_a_Joint_and_spinal_disorders
			_a_Liver_disease_moderate_or_sev
			_a_Major_psychiatric_disorders
			_a_Metabolic_disorder
			_a_Neurological_excl_epilepsy
			_a_Obesity
			_a_Osteoporosis_and_bone_disorde
			_a_Other_cardiac_conditions
			_a_Paralysis
			_a_Peripheral_nerve_or_muscular
			_a_Peripheral_vascular_disease
			_a_Renal_disease
			_a_Sleep_disorder
			_a_Upper_GI_disorders
			_a_Urinary_tract_disorder
			_a_Venous_insufficiency
			_a_Anemia
			_a_Anxiety_and_behavioral_disord
			_a_Cardiac_arrhythmia
			_a_Congestive_heart_failure
			_a_Hypertension
			_a_Myocardial_infarction
			_a_Nutritional_disorders
			_a_Other_malignancy
			_a_Pulmonary_circulation_disorde
		;

		/* Begin calculating weighted index scores. */
		/* Initialise to zero in the first instance.*/
		C3score_allsites = 0;

		DO i = 1 TO DIM(conditions);
		/*
		Error message step: if conditions in lists above don't match by position,
		then abort datastep and macro (otherwise: index weightings will be incorrect.)
		*/
			IF substr(vname(conditions{i}), 4) NE substr(vname(allsites_weights{i}), 4)
			THEN DO;
				ABORT NOLIST;
			END;
			/* Cumulative scores.*/
			/* Multiply condition indicator by appropriate weight.*/
			C3score_allsites     = C3score_allsites     + conditions{i} * allsites_weights{i};
		END;
		DROP i;

		/* Produce categorised versions of these variables...*/
		/* As described in:
		Sarfati, D., Gurney, J., Stanley, J., & Koea, J. (2014). 
		A retrospective cohort study of patients with stomach and liver cancers: 
		the impact of comorbidity and ethnicity on cancer care and outcomes. 
		BMC Cancer, 14, 821. doi: 10.1186/1471-2407-14-821
		*/
		
		/*Allsites categorised...*/
		IF      C3score_allsites LE 0 THEN C3cat_allsites = 0;
		ELSE IF C3score_allsites LE 1 THEN C3cat_allsites = 1;
		ELSE IF C3score_allsites LE 2 THEN C3cat_allsites = 2;
		ELSE IF C3score_allsites GT 2 THEN C3cat_allsites = 3;
		IF C3score_allsites = . THEN C3cat_allsites = . ;

		/* Drop weight variables unless KEEPWEIGHTS macro option == 1*/
		%IF &KEEPWEIGHTS NE 1 %THEN %DO;
			DROP _: ;
		%END;
		/* Drop C3 condition variables unless RETAINC3CONDITIONS macro option == 1*/
		%IF &RETAINC3CONDITIONS NE 1 %THEN %DO;
			DROP C3_: ;
		%END;
	RUN ;

	/*Will keep a transposed version of these weights for review.*/
	PROC TRANSPOSE DATA=index_weights OUT=index_weights_transposed;
	RUN;

	DATA index_weights_transposed;
	SET index_weights_transposed;
		Condition = substr(_NAME_, 4);
		Weight    = COL1;
		DROP _NAME_ COL1;
	RUN;

%MEND score_C3_conditions;

*----------------------score_Charlson_conditions macro----------------------*;
%MACRO Score_C3_Charlson_conditions(inputdata, outputdata, retainChconditions=1);

/* Load input data.*/
DATA &outputdata;
SET  &inputdata ;
	/********************************************/
	/*Largely adapted from Clare Salmond's code.*/
	/********************************************/
	*Calculate Charlson score, using weights from their Charlson et al. (1987);
	CharlsonScore = 
	/* Single weighted conditions */
	Ch_MI  + Ch_CHF + Ch_PVD + Ch_CBVD + Ch_Dementia + Ch_CPD + Ch_RD + 
	Ch_PUD + Ch_LiverMild +	Ch_DiabNoComp + 
	/*Double weighted conditions...*/
	(Ch_HemiPara * 2) + (Ch_Renal * 2) + (Ch_DiabWithComp * 2) +
	/*Higher weighted variables...*/
	(Ch_LiverModSevere * 3) + (Ch_AIDS * 6);
	/* Note that for the C3 index, Cancer is NOT included in scoring.*/

	IF CharlsonScore = . THEN CharlsonScore = 0; 

	LENGTH CharlsonCat $2;
	IF CharlsonScore = 0 THEN CharlsonCat = '0';
	IF CharlsonScore = 1 THEN CharlsonCat = '1';
	IF CharlsonScore = 2 THEN CharlsonCat = '2';
	IF CharlsonScore >= 3 THEN CharlsonCat = '3+';
	/***********************/	

	/* Drop Ch_ condition variables unless RETAINChCONDITIONS macro option == 1*/
	%IF &RETAINChCONDITIONS NE 1 %THEN %DO;
		DROP Ch_: ;
	%END;
RUN ;

%MEND Score_C3_Charlson_conditions;

%MACRO Merge_C3_finalset(patient_id=&IDvarlist, passAddCharlson=&addCharlson);
/*
Merge datasets.
Note that some people may be in the superset (any-time conditions)
but not in the subset (conditions that precede cancer treatment)
or vice versa
and so we need to fill in these blanks with zeros (condition not found).
					*/
DATA _C3_joined_conditions;
MERGE _C3_conditionset1 _C3_conditionset2
	/* Add the Charlson conditions on, if using.*/
	%IF &passAddCharlson EQ 1 %THEN %DO;
		_Ch_conditionset1 _Ch_conditionset2
	%END;
	;
	BY &patient_id;

	ARRAY all_conditions {*} C3_: %IF &AddCharlson EQ 1 %THEN %DO; Ch_: %END; ;
		DO i = 1 to DIM(all_conditions);
			IF all_conditions{i} EQ . THEN all_conditions{i} = 0;
		END;
	DROP i;
RUN;
%MEND Merge_C3_finalset;

/*************************/
/*End of C3 specific code*/
/*************************/


/*************************/
/* M3 specific sub-macros,
in order of call in main code:

CodeConditions_M3
CombineM3Exclusions
CodeCancerM3
DefineICDCodes_M3
Score_M3
DefineICDCodes_CharlElix
ScoreCharlson_M3
ScoreElixhauser_M3
Merge_M3_finalset
*/
/*************************/

%MACRO CodeConditions_M3(inputdata= , patient_id = ,
						 clin_cd_name = ,
						 condition_FMT=ICD_Cat);

	%LET outlib = work;	
	%IF &condition_FMT=ICD_Cat %THEN %DO;
		%LET lookup_table=work._lookup_ICD_CAT;
		%LET outname = _MM_M3;
		%LET var_prefix = M3_;
	%END; 
	%ELSE %IF &condition_FMT=ICD_Cat_Cancer %THEN %DO;
		%LET lookup_table=work._lookup_ICD_CAT_cancer;
		%LET outname = _MM_cancer;
		%LET var_prefix = Canc_;
	%END; 
	
	%ELSE %IF &condition_FMT=DIA_COMP %THEN %DO;
		%LET lookup_table=work._lookup_diabetes_comp;
		%LET outname = _MM_comp_dia;
		%LET var_prefix = DIA_;
	%END; 

	%ELSE %IF &condition_FMT=HYP_EXCL %THEN %DO;
		%LET lookup_table=work._lookup_hypertension_excl;
		%LET outname = _MM_excl_hyp;
		%LET var_prefix = HYP_;
	%END; 

	%ELSE %IF &condition_FMT=OST_EXCL %THEN %DO;
		%LET lookup_table=work._lookup_osteoporosis_excl;
		%LET outname = _MM_excl_ost;
		%LET var_prefix = OST_;
	%END; 

	%ELSE %IF &condition_FMT=ICD_Charlson %THEN %DO;
		%LET lookup_table=work._lookup_icd_charlson;
		%LET outname = _MM_Char;
		%LET var_prefix = CHAR_;
	%END; 
	
	%ELSE %IF &condition_FMT=ICD_Elixhauser %THEN %DO;
		%LET lookup_table=work._lookup_icd_Elixhauser;
		%LET outname = _MM_Elix;
		%LET var_prefix = ELIX_;
	%END; 

		%ELSE %IF &condition_FMT=ICD_Elixhauser_B %THEN %DO;
		%LET lookup_table=work._lookup_icd_Elixhauser_B;
		%LET outname = _MM_ElixB;
		%LET var_prefix = ELIX_;
	%END; 

	PROC CONTENTS DATA=&lookup_table OUT = _vars(KEEP= varnum name) NOPRINT;
	RUN;

%PUT Warning given here about SELECT and ORDER BY is ignorable (confirmed): ;
	PROC SQL NOPRINT;
		SELECT DISTINCT name
		INTO :lookup_variable SEPARATED BY ' '
		FROM _vars
		ORDER BY varnum;
	QUIT; 

/* load diagnosis data*/
	DATA _diag;
	SET &inputdata;
	RUN;
/*	%findfmt($ICD_Cat, F);*/
	DATA _diag_cat;
	SET _diag;
		LENGTH Cond_status 3;
		%PUT &clin_cd_name;
		&lookup_variable = PUT(&clin_cd_name, $&condition_FMT..);
		Cond_status=1;
		IF &lookup_variable='Other' THEN DELETE;
	RUN;

	/* count once each category for each person*/
	PROC SORT DATA=_diag_cat OUT=_diag_cat_uni (KEEP=&patient_id  &lookup_variable Cond_status ) NODUPKEY;
		BY  &patient_id 	&lookup_variable;
	RUN;

	/* get unique people and generate empty table with all potential categories listed*/
	PROC SORT DATA=_diag_cat_uni OUT=_person_uni (KEEP = &patient_id ) NODUPKEY;
		BY   &patient_id;
	RUN;

	PROC SQL;
		CREATE TABLE _person_cond_empty AS
			SELECT * FROM _person_UNI
			CROSS JOIN &lookup_table ;
	QUIT;

	PROC SORT DATA=_person_cond_empty ;
		BY   &patient_id &lookup_variable;
	RUN;

	/* assign NMDS status*/
	PROC SORT DATA=_diag_cat_uni;
		BY   &patient_id &lookup_variable;
	RUN;

	DATA _Cond_file;
		merge _person_cond_empty _diag_cat_uni;
		BY   &patient_id &lookup_variable;
		IF Cond_status=. THEN Cond_status=0;
		/*JS: Following line will make Cond_Cat shorter by removing numbered prefix.
		But note that this doesn't seem to speed up the proc transpose step below.*/
		&lookup_variable = substr(&lookup_variable, 9);
		%IF &condition_FMT EQ DIA_COMP %THEN %DO;
		&lookup_variable = "COMP";
		%END;
		%IF &condition_FMT EQ HYP_EXCL OR &condition_FMT EQ OST_EXCL %THEN %DO;
		&lookup_variable = "EXCL";
		%END;
	RUN;

	PROC SORT DATA=_Cond_file;
		BY   &patient_id &lookup_variable;
	RUN;

	/*Drop prefix option from proc transpose statement
	and add in id variable to give verbose names to columns.*/
	/*For some reason this is REALLY SLOW (takes 13 seconds, rather than 0.5 s)*/

	/* UPDATE v1.01: Don't allow spaces in variable names! Default in Enterprise Guide 
	is for more flexible name allowances, which breaks the code */
	%LET VALVARSET = %sysfunc(getoption(VALIDVARNAME)); 
	OPTIONS VALIDVARNAME=V7;
	PROC TRANSPOSE DATA=_Cond_file OUT=WORK.&outname(DROP=_name_ COMPRESS=binary) PREFIX=&var_prefix;
		ID  &lookup_variable; /*names is too long, not good for following processing?*/
		VAR Cond_status;
		BY  &patient_id;
	RUN;
	* Reset the naming option again here;
	OPTIONS VALIDVARNAME=&VALVARSET;
	
/*Ideas here -- tidy up some intermediary datasets?*/


/*Actually...do this debugging deletion step elsewhere (in master loop?)*/
 
%MEND CodeConditions_M3;

%MACRO CombineM3Exclusions(passIDvarlist=);
/*Code for taking diabetes complications and hypertension/osteoporosis exclusions*/
/*Only coded as a macro to make earlier code simpler to read.*/
/*All should be sorted by &IDvarlist (only dynamic variable)*/
/*And only one row per unique combo of the variables in &IDvarlist*/

DATA _MM_M3_Combined;
MERGE _MM_M3(in=M) _MM_comp_dia _MM_excl_hyp _MM_excl_ost;
		BY &IDvarlist;

		IF M; /*Only include if in MASTER dataset...
				So if no other M3 index conditions, 
				then ignore complications/exclusions data.*/
		
		IF DIA_COMP EQ . THEN DIA_COMP = 0;
		IF HYP_EXCL EQ . THEN HYP_EXCL = 0;
		IF OST_EXCL EQ . THEN OST_EXCL = 0;

		/* update diabetes complication  
		/*	a. have code of diabetes uncomplication, but have additional dia_comp=1 code
			   SET NMDS_Diabetes complicatio=1 and NMDS_Diabetes uncomplicated=0 ?*/
		IF   M3_Diabetes_uncomplicated EQ 1 AND DIA_COMP EQ 1 THEN DO; 
			 M3_Diabetes_complicated=1;
			 M3_Diabetes_uncomplicated=0;
			 /*For troubleshooting...*/
		END;

		/*  Diabetes uncomplicated
			remove diabetes_uncomplicated if NMDS_Diabetes_complication AND NMDS_Diabetes_uncomplicated=1 
			A few cases found from NMDS that had both records in need of updating.*/
		IF  M3_Diabetes_uncomplicated EQ 1 AND M3_Diabetes_complicated EQ 1   THEN DO; 
			 M3_Diabetes_uncomplicated = 0;
			 /*For troubleshooting...*/
		END;

		/* update Hypertension uncomplicated
			remove NMDS_Hypertension_uncomplicated if Hypertension_excl=1*/
		IF   M3_Hypertension_uncomplicated EQ 1 AND HYP_EXCL EQ 1   THEN DO;
			 M3_Hypertension_uncomplicated = 0;
		END;

		/* update Uncomplicated osteoporosis
			remove NMDS_Uncomplicated osteoporosis  if OST_excl=1*/
		IF M3_Osteoporosis_uncomplicated EQ 1 AND OST_EXCL=1   THEN DO;
		   M3_Osteoporosis_uncomplicated = 0;
		END;

		/*update  NMDS_Metastatic cancer
		  remove other cancers if had record of metastatic cancer*/
		ARRAY otherCancer[9] M3_Colorectal_cancer M3_Breast_cancer  M3_Prostate_cancer M3_Lung_cancer
				M3_Lymphomas_and_leukaemias  M3_Upper_gastrointestinal_cancer M3_Malignant_melanoma
			    M3_Gynaecological_cancers M3_Other_cancers;
		
		IF M3_Metastatic_cancer EQ 1 THEN DO i=1 to 9;
			otherCancer[i] = 0;
		END;
		DROP i DIA_COMP HYP_EXCL OST_EXCL;
RUN;
%MEND CombineM3Exclusions;

%MACRO CodeCancerM3(passDataFile=, 
					patient_id = ,
					passCancerSiteName=site, 
					passCancerExtentName=extent,
					passCancerExtentVal='E');

/*Code snippet as macro to code for cancers in additional cancer file, 
then update the main M3 condition file. */


DATA _cancer;
SET &passDataFile;

	/*Code recorded cancer against ICD-10 categories for main M3 index.*/
	Can_cat=put(&passCancerSiteName, $ICD_CAT.);
	/*if extent='E', set metastatic_cancer_flag=1 and can_cat='ICD_025 Metastatic cancer'*/

	IF &passCancerExtentName EQ &passCancerExtentVal THEN DO;
		can_cat='ICD_025 Metastatic cancer';
	END;
	/*Moved this selection step to AFTER selection of ICD_025*/
	/*Ordering shouldn't have any impact as these Other codes shouldn't be 
	possible for them to be metastasised...*/
	IF can_cat EQ 'Other' THEN DELETE;
	status=1;

	/*Check if need to keep status...*/
	KEEP &patient_id can_cat status; 
RUN;


/* transpose to matrix*/
PROC SORT DATA = _cancer OUT = _can_cat(KEEP = can_cat) NODUPKEY;
	BY can_cat;
RUN;

/* get unique people and generate empty table with all potential categories listed*/
PROC SORT DATA = _cancer OUT = _NHI_uni (KEEP = &patient_id) NODUPKEY;
	BY &patient_id ;
RUN;

PROC SQL;
	CREATE TABLE _NHI_cancer_empty AS
		SELECT * FROM _NHI_UNI
		CROSS JOIN _CAN_CAT;
QUIT;
	
PROC SORT DATA =_NHI_cancer_empty;
	BY &patient_id can_cat;
RUN;

PROC SORT DATA = _cancer OUT = _cancer_1 NODUPKEY;
	BY &patient_id can_cat;
RUN;

/*merge in known cancers from NZCR*/
DATA _Cancer_tot;
	MERGE _NHI_cancer_empty _cancer_1;
	BY &patient_id can_cat;
	IF status EQ . THEN status=0;
	Can_cat = SUBSTR(can_cat, 9);
RUN;

/*Added nodupkey below to just select single record per person/cancer combination.*/
PROC SORT DATA =_cancer_tot NODUPKEY;
	BY &patient_id can_cat;
RUN;

PROC TRANSPOSE DATA = _cancer_tot OUT = _Cancer_T(DROP = _name_ ) PREFIX = CR_;
	ID can_cat;  
	VAR status;
	BY &patient_id;
RUN;

PROC SORT DATA = _Cancer_T;
	BY &patient_id;
RUN;

/* link to main M3 data, and generate combine Cancer category*/
PROC SORT DATA = _MM_M3_Combined;
	BY &patient_id;
RUN;

PROC SORT DATA = _Cancer_T;
	BY &patient_id;
RUN;

DATA _MM_M3_Combined;
	MERGE _MM_M3_Combined _Cancer_T ;
	BY &patient_id;
/*
Now to come up with composite cancer categorisation.
Firstly combine hospital and cancer registry records.
Then use metastatic cancer to "trump" site-specific cancer records.
*/

	ARRAY M3_cancers{10}    M3_Metastatic_cancer M3_Breast_cancer M3_Colorectal_cancer M3_Gynaecological_cancers
		  		    M3_Lung_cancer M3_Lymphomas_and_leukaemias M3_Malignant_melanoma
					M3_Other_cancers M3_Prostate_cancer M3_Upper_gastrointestinal_cancer;

	ARRAY CR_cancers{10} 	CR_Metastatic_cancer CR_Breast_cancer CR_Colorectal_cancer CR_Gynaecological_cancers
		  		    CR_Lung_cancer CR_Lymphomas_and_leukaemias CR_Malignant_melanoma
					CR_Other_cancers CR_Prostate_cancer CR_Upper_gastrointestinal_cancer;
	
	/* Use M3 as the main cancer data, so update this.*/
	DO i=1 TO 10;
		/* Now going to set main CANC record to 1 if a corresponding NMDS record is 1. Simple!*/
		IF M3_cancers[i] EQ . THEN M3_cancers[i] = 0;  
		IF CR_cancers[i] EQ 1 THEN M3_cancers[i] = 1;  
	END;

	/* If metastatic cancer detected...	[position 1]*/
	/* Then override cancer site classifications.*/
	/* Which are stored in [position 2 to position 10]*/
	IF M3_cancers[1] EQ 1 THEN DO;  
		DO i=2 TO 10;
			M3_cancers[i] = 0;
		END;
	END;
	/*clean up*/
	DROP i CR_:;

RUN;

%MEND CodeCancerM3;

*------Comorbidity coding for M3: FORMAT code-------;

*--JStanley January 2018 -- ;
*----macro for the specified cancer site----;
%macro DefineICDCodes_M3();
PROC FORMAT library=work CNTLOUT=work._lookup_icd_cat_total;
	VALUE $ICD_CAT  
			'I210'-'I219', 'I220'-'I229', 'I230'-'I239', 'I241', 'I252'
			='ICD_001 Myocardial infarction'

			'I099', 'I110', 'I130', 'I132', 'I255', 'I420','I425', 'I426', 'I427','I428','I429',
			'I430'-'I439', 'I500'-'I509'
			='ICD_002 Congestive heart failure'

			'I700', 'I7000'-'I7099','I731','I738','I739', 'I740'-'I749', 'I771','K551','K5520'-'K5529','K558','K559'
			='ICD_003 Peripheral vascular'

			'I7100'-'I7199','I720'-'I729'
			='ICD_004 Aortic and other aneurysms'

			'I830', 'I832','I872'
			='ICD_005 Venous insufficiency'
		 
			'I600'-'I609','I610'-'I619','I620'-'I629','I630'-'I639', 'I64','I650'-'I659','I660'-'I669',
			'I670'-'I679','I690'-'I699','G450'-'G459','G460'-'G469'
			='ICD_006 Cerebrovascular disease'
 
			'F000'-'F009','F010'-'F019','F020','F021','F022','F023','F03','F051', 'G300'-'G309','G310','G311'
			='ICD_007 Dementia'
 
			'F0400'-'F0499','F060'-'F069', 'F070','F071','F078','F079','F09','G931'
			='ICD_008 Mental and behavioural disorders due to brain damage'

			'E840'-'E849','J40','J410'-'J419','J42','J430'-'J439','J440'-'J449','J450'-'J459','J46','J47',
			'J60','J61','J620'-'J629','J630'-'J639','J64','J65','J660'-'J669','J670'-'J679','J684','J701',
			'J703','J840'-'J849','J961','J980','J982','J983','J984'
			='ICD_009 Chronic pulmonary'

			'L930'-'L939','M0500'-'M0599','M0600'-'M0699','M0800'-'M0899', 'M1200'-'M1209', 'M1230'-'M1239',
			'M300'-'M309', 'M310'-'M319', 'M320'-'M329','M330'-'M339','M340'-'M349','M350','M351','M352',
			'M353','M354','M355','M356','M358','M359'
			='ICD_010 Connective tissue'

			'K220','K221','K224','K225','K228','K229','K250'-'K259','K260'-'K269','K270'-'K279',
			'K280'-'K289','K311','K312','K314','K316'
			='ICD_011 GI ulcer upper GI'

			'E100','E1010'-'E1019','E109','E1100'-'E1109','E1110'-'E1119','E119','E120','E121','E129',
			'E1300'-'E1309','E1310'-'E1319','E139','E1400'-'E1409','E1410'-'E1419','E149'
			='ICD_012 Diabetes uncomplicated'

			'E1020'-'E1029','E1030'-'E1039','E1040'-'E1049','E1050'-'E1059','E1060'-'E1069','E1070'-'E1079',
			'E108','E1120'-'E1129','E1130'-'E1139','E1140'-'E1149','E1150'-'E1159','E1160'-'E1169','E1170'-'E1179',
			'E118','E1220'-'E1229','E1230'-'E1239','E1240'-'E1249','E1250'-'E1259','E1260'-'E1269','E1270'-'E1279',
			'E1280'-'E1289','E1320'-'E1329','E1330'-'E1339','E1340'-'E1349','E1350'-'E1359','E1360'-'E1369',
			'E1370'-'E1379','E138','E1420'-'E1429','E1430'-'E1439','E1440'-'E1449','E1450'-'E1459',
			'E1460'-'E1469','E1470'-'E1479','E148'
			='ICD_013 Diabetes complicated'

			'G041', 'G114','G8000'-'G8009','G801','G802','G810'-'G819','G8200'-'G8299','G830','G831','G832',
			'G833','G834','G839'
			='ICD_014 Paralysis'

			/* JS amended 27/01/2015  */
		  	'I129',  'I139',
			'Q600'-'Q609', 'Q611'-'Q613',
			'N032','N033','N034','N035','N036','N037','N038','N039','N042','N043','N044','N045','N046','N047',
			'N048','N049','N052','N053','N054','N055','N056','N057','N058','N059','N110'-'N119','N180'-'N1899',
			'N19','N250','N258','N259','I120','I131','Z490'-'Z499','Z940','Z992'
			='ICD_015 Chronic renal'

			'C180'-'C189', 'C19','C20', 'C210'-'C219'='ICD_016 Colorectal cancer'

			'C500'-'C509'='ICD_017 Breast cancer'

			'C61'='ICD_018 Prostate cancer'
 
 			'C33', 'C340'-'C3499'= 'ICD_019 Lung cancer'
 
			'C810'-'C859', 'C9100'-'C9599','C960'-'C969'='ICD_020 Lymphomas and leukaemias'

			'C150'-'C179','C220'-'C259'='ICD_021 Upper gastrointestinal cancers'
 
			'C430'-'C439'= 'ICD_022 Malignant melanoma'

			'C510'-'C519','C52','C530'-'C549','C55','C56','C570'-'C579','C58'
			='ICD_023 Gynaecological cancers'

			'C000'-'C009','C01','C020'-'C029','C030'-'C039','C040'-'C049','C050'-'C059','C060'-'C069','C07',
			'C080'-'C089','C090'-'C099','C100'-'C109','C110'-'C119', 'C12','C130'-'C139','C140'-'C149',
		    'C260'-'C269','C300'-'C309','C310'-'C319','C320'-'C329','C37',
			'C380'-'C389', 'C390'-'C399','C40'-'C409', 'C4100'-'C4199', 'C450'-'C459','C460'-'C469','C470'-'C479',
			'C480'-'C489','C490'-'C499', 'C600'-'C609', 'C620'-'C629','C630'-'C639','C64','C65','C66','C670'-'C679',
			'C680'-'C689', 'C690'-'C699','C700'-'C709','C710'-'C719','C720'-'C729','C73','C740'-'C749','C750'-'C759',
			'C760'-'C769', 'C8800'-'C8899','C9000'-'C9099' 
			='ICD_024 Other cancers'
 
			'C770'-'C799'='ICD_025 Metastatic cancer'
 
			'K700'-'K709', 'K711', 'K713', 'K714', 'K715', 'K717', 'K721', 'K729', 'K730'-'K739','K740'-'K749', 
			'K760', 'K762', 'K763', 'K764', 'K765', 'K766', 'K767', 'K768', 'K769','I850'-'I859', 'I864', 
			'I9820'-'I9829', 'Z944'
			='ICD_026 Liver disease: moderate or severe'
 
			'B20', 'B21', 'B22', 'B230'-'B239', 'B24', 'F024', 'Z21'
			='ICD_027 AIDS'

			'I200'-'I209'
			='ICD_028 Angina'

			/*JS amended 27/1/2015*/
			'I10' 
			='ICD_029 Hypertension uncomplicated' /*Hypertension primary*/

			'I441', 'I442', 'I443', 'I456', 'I459', 'I470'-'I479', 'I48', 'I490'-'I499', 'T821', 'Z450', 'Z950'
			='ICD_030 Cardiac arrhythmia'
				
			'I260'-'I269', 'I270'-'I279', 'I280', 'I281', 'I288', 'I289'
			='ICD_031 Pulmonary circulation disorder'

			'I050'-'I059', 'I060'-'I069', 'I070'-'I079', 'I080'-'I089', 'I091', 'I098', 'I340'-'I349', 
			'I350'-'I359', 'I360'-'I369', 'I370'-'I379', 'I38', 'T820', 'Q230', 'Q231', 'Q232', 'Q233', 'Q238',
			'Q239', 'Z952', 'Z953', 'Z954'
			='ICD_032 Cardiac valve'
 
			'K500'-'K509', 'K510'-'K519', 'K522', 'K528', 'K529'
			='ICD_033 Bowel disease inflammatory'

			'G10', 'G110', 'G111', 'G112', 'G113', 'G118', 'G119', 'G120'-'G129', 'G130'-'G139', 'G20', 
			'G210'-'G219', 'G230'-'G239', 'G255', 'G312', 'G318', 'G319', 'G35', 'G360'-'G369','G370'-'G379',
			'G900'-'G909', 'G934', 'R470'
			='ICD_034 Other neurological disorders exc epilepsy'

			'G4000'-'G4009', 'G4010'-'G4019', 'G4020'-'G4029', 'G4030'-'G4039', 'G4040'-'G4049', 'G4060'-'G4069',
			'G4070'-'G4079', 'G4080'-'G4089', 'G4090'-'G4099', 'G410'-'G419',
			/*Added v0.41 -- not quite flexible enough!*/
			'G400', 'G401', 'G402', 'G403', 'G404', 'G406', 'G407', 'G408', 'G409'
			='ICD_035 Epilepsy'

			'G600'-'G609', 'G610'-'G619', 'G620', 'G621', 'G622', 'G628', 'G629', 'G64', 'G700'-'G709', 'G710'-'G719', 
			'G720', 'G721',  'G722', 'G723', 'G724', 'G728', 'G729', 'G731'
			='ICD_036 Muscular peripheral nerve disorder'

			'F200'-'F209', 'F220'-'F229', 'F250'-'F259', 'F28', 'F29', 'F302', 'F310'-'F319', 
            'F3210'-'F3219', 'F3220'-'F3229', 'F3230'-'F3239', 'F3280'-'F3289', 'F3290'-'F3299', 
			'F330'-'F339', 'F39',
			/*Added v0.41 -- not quite flexible enough!*/
            'F321', 'F322','F323','F328',  'F329' 
			='ICD_037 Major psychiatric disorder'
 
			'F4000'-'F4099', 'F410'-'F419', 'F420'-'F429','F440'-'F4499','F450'-'F4599', 'F480'-'F489',
			'F500'-'F509', 'F550'-'F559', 'F59', 'F600'-'F6099', 'F61', 'F630'-'F639', 'F640'-'F649','F650'-'F659',
			'F660'-'F669', 'F680'-'F689', 'F69'
			='ICD_038 Anxiety and Behavioural disorders'

			'D550'-'D559', 'D560'-'D569', 'D570'-'D579', 'D580'-'D589','D590', 'D591','D592','D593','D594', 'D598', 'D599', 
			'D600'-'D609', 'D610'-'D619', 'D640'-'D649', 'D66', 'D67', 'D680', 
			'D681', 'D682', 'D688', 'D689', 'D691' , 'D692', 'D693', 'D694', 'D696', 'D698', 'D699','D70','D71', 
			'D720'-'D729', 'D740'-'D749', 'D750', 'D752', 'D758', 'D759'
			='ICD_039 Coagulopathy and other blood disorder'

			'D500'-'D509', 'D510'-'D519', 'D520'-'D529', 'D530'-'D539'
			='ICD_040 Anemia deficiency'

			'E660'-'E669'
			='ICD_041 Obesity'

			'F101', 'F102', 'F103', 'F104', 'F105', 'F106', 'F107', 'F108', 'F109',  'Z502', 'Z714'
			='ICD_042 Alcohol abuse'
 
			'F110'-'F119', 'F120'-'F129', 'F130'-'F139', 'F140'-'F149', 'F150'-'F159', 'F160'-'F169','F180'-'F189',
			'F190'-'F199', 'Z503', 'Z715', 'Z722'
			='ICD_043 Drug abuse'

			'K85', 'K850'-'K859','K860', 'K861', 'K868'
			='ICD_044 Pancreatitis'

			'E010'-'E019','E02', 'E030'-'E039', 'E050'-'E059','E062', 'E063', 'E065','E070'-'E079','E163','E164',
			'E168','E169', 'E200'-'E209', 'E210', 'E212', 'E213','E214','E215', 'E220'-'E229','E230','E232','E233',
			'E236', 'E237', 'E240', 'E241', 'E243', 'E244', 'E248','E249','E250'-'E259','E260'-'E269','E270'-'E279', 
			'E310'-'E319', 'E320'-'E329', 'E345', 'E348', 'E349'
			='ICD_045 Endocrine disorder'
 
			'N301', 'N302', 'N310'-'N319', 'N320'-'N329', 'N350'-'N359', 'N360'-'N369'
			='ICD_046 Urinary tract problem chronic'

			'A150'-'A159', 'A160'-'A169', 'A170'-'A179', 'A180'-'A189', 'A190'-'A199', 'B900'-'B909'
			='ICD_047 Tuberculosis'

			'M8000'-'M8099', /*JS removed uncomplicated Osteoporosis from here.*/
			'M8300'-'M8309', /*JS added 3/3/2015, altered 17/03/2015*/
			'M8310'-'M8319', 'M8320'-'M8329', 'M8330'-'M8339','M8340'-'M8349','M8350'-'M8359','M8380'-'M8389',
			'M8390'-'M8399', 'M8500'-'M8599','M8630'-'M8639','M8640'-'M8649','M8650'-'M8659', 'M8660'-'M8669',
			'M880'-'M8899'
			='ICD_048 Bone disorders' 

			'M8100'-'M8109', 'M8110'-'M8119', 'M8150'-'M8159', 'M8180'-'M8189', 'M8190'-'M8199'
			='ICD_049 Osteoporosis Uncomplicated' 
 
			'D800'-'D809','D810'-'D819','D820'-'D829','D830'-'D839','D840'-'D849','D860'-'D869','D890'-'D899'
			='ICD_050 Immune system disorder'

			'E700'-'E709','E710'-'E719','E720'-'E729','E740'-'E749','E750'-'E759','E760'-'E769','E770'-'E779',
			'E780'-'E789','E791', 'E798', 'E799', 
			/*Next line was fixed in v 0.31: was E00-E01!*/
			'E800'-'E809', 
			'E830'-'E839', 'E850'-'E859'
			/* Removed at earlier stage.*/
/*			, 'E880'-'E889'*/
			='ICD_051 Metabolic disorder'
 
			'F700'-'F709','F710'-'F719','F720'-'F729', 'F730'-'F739', 'F780'-'F789', 'F790'-'F799','F842',
			'F843', 'F844', 'E000', 'E001', 'E002', 'E009', 'Q900'-'Q909'
			='ICD_052 Mental retardation'

			'B180'-'B189',  'B942', 'Z2250'-'Z2259'
			='ICD_053 Hepatitis Chronic viral'
 
			'F510'-'F519', 'G470', 'G471', 'G472', 'G4730'-'G4739'
			='ICD_054 Sleep disorder'

			'H800'-'H809', 'H810'-'H819', 'H830'-'H839', 'H900'-'H909', 'H910', 'H911', 'H913', 'H918', 
			'H919', 'H930', 'H931', 'H932', 'H933'
			='ICD_055 Inner ear disorder'

			'A300'-'A309', 'A310'-'A319', 'A520'-'A529', 'B91', 'B92', 'B941', 'B948', 'B949'
			='ICD_056 Infection Chronic NOS'

			'E40', 'E41', 'E42', 'E43','E440'-'E449','E45','E46','E500'-'E509','E510'-'E519','E52','E530'-'E539', 
			'E54','E550'-'E559','E560'-'E569','E58','E59','E60','E610'-'E619', 'E630'-'E639', 'E640'-'E649'
			='ICD_057 Malnutrition nutritional'

			'H160'-'H169','H181','H184', 'H185', 'H186', 'H201','H212','H301','H311','H312','H313','H314','H330',
			'H332','H333','H334','H335','H340'-'H349','H350'-'H359','H430'-'H439','H46','H470'-'H479','H490'-'H499',
			'H500'-'H509','H510'-'H519', 'H530', 'H531', 'H532','H533','H534','H536','H538', 'H539', 'H540'-'H549',
			'Q120'-'Q129','Q130'-'Q139', 'Q140'-'Q149','Q150'-'Q159'
			='ICD_058 Eye problem long term'
 
			'I119', 'I248', 'I249', 'I250', 'I2510'-'I2519', 'I253', 'I254', 'I256', 'I258', 'I259', 'I310', 'I311', 'I421',
			'I422', 'I424'
			='ICD_059 Cardiac disease other'

			'K5700'-'K5799', 'K592', 'K593', 'K900'-'K909'
			='ICD_060 Intestinal disorder'
 
			'M0700'-'M0799', 'M130'-'M1399', 'M150', 'M151', 'M152','M154','M158','M159','M4000'-'M4009',
			'M4020'-'M4029', 'M4030'-'M4039', 'M4040'-'M4049','M4050'-'M4059','M4100'-'M4199','M4200'-'M4299', 
			'M4300'-'M4399','M4500'-'M4599', 'M4600'-'M4609','M461','M4620'-'M4629','M4700'-'M4799','M4800'-'M4809',
			'M4810'-'M4819','M4820'-'M4829','M4850'-'M4859','M4880'-'M4889', 'M4890'-'M4899', 'G950', 'G951'
			='ICD_061 Joint spinal disorder'

			Other='Other';

RUN;	  
  
PROC SORT DATA=work._lookup_icd_cat_total  
			OUT = work._lookup_icd_cat (KEEP=label RENAME=(label=NMDS_Cat)) NODUPKEY;
	WHERE FMTName='ICD_CAT';
	BY Label ;
RUN;

DATA work._lookup_ICD_CAT;
	SET work._lookup_ICD_CAT;
	IF NMDS_CAT NE 'Other';
RUN;

PROC FORMAT CNTLOUT=work._lookup_complications;
	VALUE $DIA_COMP  
			'I200'-'I259', 'I600'-'I699', 'I700'-'I799', 
			'N030'-'N039', 'N040'-'N049', 'N180'-'N1891',
			'G603', 'G620'-'G629', 
			/* Fixed as at v 0.32 was G630 to G638*/
			'G638',
			'H350'-'H359', 'H36',  'H368','L97' 
			='DIA_Complication'
			Other='Other';

	VALUE $OST_EXCL
			'M8000'-'M8099', 
			'S2200'-'S2209', 
			'S3200'-'S3209', 
			'S5200'-'S5299', 
			'S7200'-'S7209'
			='OST_Exclusions'
			Other='Other';
 
	VALUE $HYP_EXCL
			'I110'-'I119',
			'I120'-'I129',
			'I130'-'I139',
			'I200'-'I259',
			'I600'-'I699',
			'I700'-'I709',
			'I710'-'I719',
			'I720'-'I729',
			'N030'-'N039', 
			'N040'-'N049', 
			'N180'-'N1891'
			='HYP_Exclusions'
			Other='Other';
	;
RUN;
  
PROC SORT DATA=work._lookup_complications
			OUT = work._lookup_diabetes_comp(KEEP=label RENAME=(label=DIA_Comp)) NODUPKEY;
	WHERE FMTName='DIA_COMP';
	BY Label ;
RUN;

DATA work._lookup_diabetes_comp;
	SET work._lookup_diabetes_comp;
	IF DIA_COMP NE 'Other';
RUN;

PROC SORT DATA=work._lookup_complications
			OUT = work._lookup_osteoporosis_excl(KEEP=label RENAME=(label=OST_Excl)) NODUPKEY;
	WHERE FMTName='OST_EXCL';
	BY Label ;
RUN;

DATA work._lookup_osteoporosis_excl;
	SET work._lookup_osteoporosis_excl;
	IF OST_EXCL NE 'Other';
RUN;

PROC SORT DATA=work._lookup_complications
			OUT = work._lookup_hypertension_excl(KEEP=label RENAME=(label=HYP_Excl)) NODUPKEY;
	WHERE FMTName='HYP_EXCL';
	BY Label ;
RUN;

DATA work._lookup_hypertension_excl;
	SET work._lookup_hypertension_excl;
	IF HYP_EXCL NE 'Other';
RUN;

%MEND DefineICDCodes_M3;


%MACRO Score_M3();

DATA _mm_m3_scored;
SET  _mm_m3_combined;

/*Make dynamic array of all morbidity conditions (using M3 wildcard with : colon).*/
ARRAY all_M3{*} M3: ;

/*Step through all conditions in array above...*/
DO i = 1 to dim(all_M3);
	/*If a record is blank (only happens when no hospital contacts for any condition in list) 	*/
	/*set the value to zero instead.															*/
	IF all_M3{i} EQ . THEN all_M3{i} = 0;
END;
DROP i; /*Tidy up loop counter variable.*/


/*Manual scoring of M3 index... will be needed for macro.*/
M3Score = 
		M3_AIDS								*	0.452647425	+
		M3_Alcohol_abuse					*	0.576907507	+
		M3_Anemia_deficiency				*	0.180927466	+
		M3_Angina							*	0 			+ /*-0.082399267	+*/
		M3_Anxiety_and_Behavioural_disor	*	0.121481351	+
		M3_Aortic_and_other_aneurysms		*	0.260195993	+
		M3_Bone_disorders					*	0.132827597	+
		M3_Bowel_disease_inflammatory		*	0.086960591	+
		M3_Breast_cancer					*	0.411891435	+
		M3_Cardiac_arrhythmia				*	0.173859876	+
		M3_Cardiac_disease_other			*	0 			+ /*-0.104225698	+*/
		M3_Cardiac_valve					*	0.256577208	+
		M3_Cerebrovascular_disease			*	0.097803808	+
		M3_Chronic_pulmonary				*	0.6253395	+
		M3_Chronic_renal					*	0.334155906	+
		M3_Coagulopathy_and_other_blood		*	0.265142145	+
		M3_Colorectal_cancer				*	0.372878764	+
		M3_Congestive_heart_failure			*	0.539809861	+
		M3_Connective_tissue				*	0.290446442	+
		M3_Dementia							*	1.021975368	+
		M3_Diabetes_complicated				*	0.271607393	+
		M3_Diabetes_uncomplicated			*	0.299383867	+
		M3_Drug_abuse						*	0.558979499	+
		M3_Endocrine_disorder				*	0.112673001	+
		M3_Epilepsy							*	0.594991823	+
		M3_Eye_problem_long_term			*	0.179923774	+
		M3_GI_ulcer_upper_GI				*	0.152986438	+
		M3_Gynaecological_cancers			*	0.70658858	+
		M3_Hepatitis_Chronic_viral			*	0.569092852	+
		M3_Hypertension_uncomplicated		*	0.117746303	+
		M3_Immune_system_disorder			*	0.398529751	+
		M3_Infection_Chronic_NOS			*	0 			+ /*-0.237983891	+*/
		M3_Inner_ear_disorder				*	0.06090681	+
		M3_Intestinal_disorder				*	0			+ /*-0.254089697	+*/
		M3_Joint_spinal_disorder			*	0.095585857	+
		M3_Liver_disease__moderate_or_se	*	0.474321939	+
		M3_Lung_cancer						*	1.972481401	+
		M3_Lymphomas_and_leukaemias			*	1.190108503	+
		M3_Major_psychiatric_disorder		*	0.212789563	+
		M3_Malignant_melanoma				*	0.342233292	+
		M3_Malnutrition_nutritional			*	0.331335106	+
		M3_Mental_and_behavioural_disord	*	0.039711074	+
		M3_Mental_retardation				*	1.405761403	+
		M3_Metabolic_disorder				*	0.006265195	+
		M3_Metastatic_cancer				*	2.468586878	+
		M3_Muscular_peripheral_nerve_dis	*	0.208276284	+
		M3_Myocardial_infarction			*	0.197491908	+
		M3_Obesity							*	0.248243722	+
		M3_Osteoporosis_Uncomplicated		*	0.083506878	+
		M3_Other_cancers					*	1.103452294	+
		M3_Other_neurological_disorders		*	0.564391512	+
		M3_Pancreatitis						*	0 			+ /*-0.103132585	+*/
		M3_Paralysis						*	0.281895685	+
		M3_Peripheral_vascular				*	0.349250005	+
		M3_Prostate_cancer					*	0.432343447	+
		M3_Pulmonary_circulation_disorde	*	0.398432833	+
		M3_Sleep_disorder					*	0.245749995	+
		M3_Tuberculosis						*	0 			+ /*-0.104290289	+*/
		M3_Upper_gastrointestinal_cancer	*	1.941498638	+
		M3_Urinary_tract_problem_chronic	*	0.046548658	+
		M3_Venous_insufficiency				*	0.214050369	;

RUN;

%MEND Score_M3;

/*Format code for Charlson and Elixhauser index: M3 implementation (i.e. includes cancer)*/
%MACRO DefineICDCodes_CharlElix();
PROC FORMAT library=work CNTLOUT=work._lookup_icd_charlson;

/*
Codes here are from Quan et al.'s mapping of 16 Charlson comorbid conditions to ICD-10 codes
Reference:
Quan, H., V. Sundararajan, P. Halfon, A. Fong, B. Burnand, 
J. C. Luthi, L. D. Saunders, C. A. Beck, T. E. Feasby and W. A. Ghali (2005). 
Coding algorithms for defining comorbidities in ICD-9-CM and ICD-10 administrative data.
Med Care 43(11): 1130-1139.

NOTE THAT CONTRARY TO C3, MI AND CHF ARE CODED IRRESPECTIVE TO TIMING.
SINCE NO "EVENT" AT INDEX.
*/

VALUE $ICD_Charlson
'I21'-'I22XX','I252'-'I252X'
="CHA_001 MI"

'I099'-'I099X','I110'-'I110X','I130'-'I130X','I132'-'I132X',
'I255'-'I255X','I420'-'I420X','I425'-'I425X','I426'-'I426X','I427'-'I427X',
'I428'-'I428X','I429'-'I429X','I43'-'I43XX','I50'-'I50XX','P290'-'P290X'
="CHA_002 CHF"

'I70'-'I71XX','I731'-'I731X','I738'-'I739X','I771'-'I771X','I790'-'I790X',
'I792'-'I792X','K551'-'K551X','K558'-'K559X','Z958'-'Z959X'
="CHA_003 PVD"

'G45'-'G46XX','H340'-'H340X','I60'-'I69XX'
="CHA_004 CBVD"

'F00'-'F03XX','F051'-'F051X','G30'-'G30XX','G311'-'G311X'
="CHA_005 Dementia"

'I278'-'I279X','J40'-'J47XX','J60'-'J67XX','J684'-'J684X','J701'-'J701X','J703'-'J703X'
="CHA_006 CPD"

'M05'-'M06XX','M315'-'M315X','M32'-'M34XX','M351'-'M351X','M353'-'M353X','M360'-'M360X'
="CHA_007 RD"

'K25'-'K28XX'
="CHA_008 PUD"

'B18'-'B18XX','K700'-'K703X','K709'-'K709X','K713'-'K715X','K717'-'K717X',
'K73'-'K74XX','K760'-'K760X','K762'-'K764X','K768'-'K769X','Z944'-'Z944X'
="CHA_009 LiverMild"

'E100'-'E101X','E106'-'E106X','E108'-'E111X','E116'-'E116X','E118'-'E121',
'E126'-'E126X','E128'-'E131X','E136'-'E136X','E138'-'E141X','E146'-'E146X',
'E148'-'E149X'
="CHA_010 DiabNoComp"

'G041'-'G041X','G114'-'G114X','G801'-'G802X','G81'-'G82XX','G830'-'G834X',
'G839'-'G839X'
="CHA_011 HemiPara"

'I120'-'I120X','I131'-'I131X','N032'-'N037X','N052'-'N057X','N18'-'N19XX',
'N250'-'N250X','Z490'-'Z492X','Z940'-'Z940X','Z992'-'Z992X'
="CHA_012 Renal"

'E102'-'E105X','E107'-'E107X',
'E112'-'E115X','E117'-'E117X',
'E122'-'E125X','E127'-'E127X',
'E132'-'E135X','E137'-'E137X',
'E142'-'E145X','E147'-'E147X'
="CHA_013 DiabWithComp"

'C00'-'C21XX','C23'-'C26XX','C30'-'C33XX','C37'-'C39XX',
'C43'-'C43XX','C45'-'C49XX','C50'-'C50XX',
'C51'-'C58XX','C60'-'C70XX', 'C72'-'C75XX','C81'-'C85XX',
'C88'-'C88XX','C90'-'C95XX'
="CHA_014 Cancer"

'I850'-'I850X','I859'-'I859X','I864'-'I864X','I982'-'I982X','K704'-'K704X',
'K711'-'K711X','K721'-'K721X','K729'-'K729X','K765'-'K767X'
="CHA_015 LiverModSevere"

'B20'-'B22XX','B24'-'B24XX'
="CHA_016 AIDS"

'C770'-'C80X'
="CHA_017 CancerMetastatic"

Other="Other"
;

RUN;

PROC SORT DATA=work._lookup_icd_charlson  
			OUT = work._lookup_icd_charlson (KEEP=label RENAME=(label=CHAR_Cat)) NODUPKEY;
	WHERE FMTName='ICD_CHARLSON';
	BY Label ;
RUN;

DATA work._lookup_ICD_Charlson;
	SET work._lookup_ICD_Charlson;
	IF CHAR_CAT NE 'Other';
RUN;

/*Now Elixhauser*/
PROC FORMAT library=work CNTLOUT=work._lookup_icd_elixhauser;

/*
Codes here are from Quan et al.'s mapping of 31 Elixhauser comorbid conditions to ICD-10 codes
Reference:
Quan, H., V. Sundararajan, P. Halfon, A. Fong, B. Burnand, 
J. C. Luthi, L. D. Saunders, C. A. Beck, T. E. Feasby and W. A. Ghali (2005). 
Coding algorithms for defining comorbidities in ICD-9-CM and ICD-10 administrative data.
Med Care 43(11): 1130-1139.
*/

VALUE $ICD_Elixhauser
'I099'-'I099X','I110'-'I110X','I130'-'I130X','I132'-'I132X',
'I255'-'I255X','I420'-'I420X','I425'-'I425X','I426'-'I426X','I427'-'I427X',
'I428'-'I428X','I429'-'I429X','I43'-'I43XX','I50'-'I50XX','P290'-'P290X'
="ELX_001 CHF"

'I441' - 'I443X',
'I456'-'I456X', 'I459'-'I459X', 'I47'-'I49XX', 
'R000' - 'R000X', 'R001' - 'R001X', 'R008' - 'R008X', 'T821' - 'T821X',
'Z450' - 'Z450X', 'Z950' - 'Z950X'
="ELX_002 Cardiac arrythmia"

'A520'-'A520X', 'I05'-'I08XX', 'I091'-'I091X', 'I098'-'I098X', 'I34'-'I39XX', 'Q230'-'Q233X', 'Z952X'-'Z954X'
="ELX_003 Valvular disease"

'I26'-'I26XX', 'I27'-'I27XX', 'I280'-'I280X', 'I288'-'I288X', 'I289'-'I289X'
="ELX_004 Pulmonary circulation disorders"

'I70'-'I71XX','I731'-'I731X','I738'-'I738X','I739'-'I739X', 'I771'-'I771X','I790'-'I790X',
'I792'-'I792X','K551'-'K551X','K558'-'K559X','Z958'-'Z959X'
="ELX_005 PVD"

'I10'-'I10XX'
="ELX_006 Hypertension uncomplicated"

'G041'-'G041X','G114'-'G114X','G801'-'G801X', 'G802'-'G802X','G81'-'G82XX','G830'-'G834X',
'G839'-'G839X'
="ELX_008 Paralysis"

'E100'-'E101X',
'E109'-'E111X',
'E119'-'E121X',
'E129'-'E131X',
'E139'-'E141X',
'E149'-'E149X'
="ELX_011 DiabNoComp"

'E102'-'E108X',
'E112'-'E118X',
'E122'-'E128X',
'E132'-'E138X',
'E142'-'E148X'
="ELX_012 DiabWithComp"

'E00'-'E03X', 'E890'-'E890X'
="ELX_013 Hypothyroidism"

'I120'-'I120X','I131'-'I131X',
'N18'-'N19XX', 'N250'-'N250X',
'Z490'-'Z492X','Z940'-'Z940X','Z992'-'Z992X'
="ELX_014 Renal"

'B18'-'B18XX','I85'-'I85XX',
'I864'-'I864X','I982'-'I982X',
'K70'-'K70XX', 
'K711'-'K711X',
'K713'-'K715X','K717'-'K717X',
'K72'-'K74XX',
'K760'-'K760X',
'K762'-'K769X',
'Z944'-'Z944X'
="ELX_015 Liver"

'K257'-'K257X',
'K259'-'K259X',
'K267'-'K267X',
'K269'-'K269X',
'K277'-'K277X',
'K279'-'K279X',
'K287'-'K287X',
'K289'-'K289X'
="ELX_016 PUD"

'B20'-'B22XX','B24'-'B24XX'
="ELX_017 AIDS"

'C81'-'C85X','C88'-'C88X',
'C96'-'C96X', 'C900'-'C900X',
'C902'-'C902X'
="ELX_018 Lymphoma"

'C77'-'C80X'
="ELX_019 CancerMetastatic"

'C00'-'C26XX',
'C30'-'C34XX',
'C37'-'C41XX',
'C43'-'C43XX',
'C45'-'C58XX',
'C60'-'C76XX', 
'C97'-'C97XX'
="ELX_020 Cancer solid no metastasis"

'L940'-'L940X','L941'-'L941X','L943'-'L943X',
'M05'-'M05XX','M06'-'M06XX',
'M08'-'M08XX','M120'-'M120X','M123'-'M123X',
'M30'-'M30XX','M31'-'M313X','M32'-'M35XX',
'M45'-'M45XX','M461'-'M461X','M468'-'M468X','M469'-'M469X'
="ELX_021 Rheumatoid arthritis"

'D65'-'D68XX',
'D691'-'D691X',
'D693'-'D696X'
="ELX_022 Coagulopathy"

'E66'-'E66X'
="ELX_023 Obesity"

'E40'-'E46X', 'R634'-'R634X', 'R64'-'R64X'
="ELX_024 Weight loss"

'E222'-'E222X', 'E86'-'E86XX', 'E87'-'E87XX'
="ELX_025 Fluid and electrolyte disorders"

'D500' - 'D500X'
="ELX_026 Blood loss anemia"

'D508' - 'D508X', 'D509' - 'D509X', 
'D51' - 'D53X' 
="ELX_026 Deficiency anemia"


'F11'-'F16X', 'F18'-'F18X', 'F19'-'F19X', 'Z715', 'Z722'
="ELX_028 Drug abuse"

'F20'-'F20X', 'F22'-'F25X', 'F28'-'F28X', 'F29'-'F29X', 
'F302', 'F312','F315'
="ELX_029 Psychoses"

Other="Other"
;
RUN;

PROC SORT DATA=work._lookup_icd_Elixhauser
			OUT = work._lookup_icd_Elixhauser (KEEP=label RENAME=(label=ELIX_Cat)) NODUPKEY;
	WHERE FMTName='ICD_ELIXHAUSER';
	BY Label ;
RUN;

DATA work._lookup_ICD_Elixhauser;
	SET work._lookup_ICD_Elixhauser;
	IF ELIX_CAT NE 'Other';
RUN;


/*REPEAT THIS FOR CONDITIONS THAT OVERLAP WITH FIRST LIST...*/


PROC FORMAT library=work CNTLOUT=work._lookup_icd_elixhauser_B;
/*Second version of formats for overlapping conditions*/
VALUE $ICD_Elixhauser_B
/*N.B. NON-EXCLUSIVE CODING WITH CHF*/
'I11'-'I13XX', 'I15'-'I15XX'
="ELX_007 Hypertension complicated"

/*N.B. NON-EXCLUSIVE CODING WITH PARALYSIS*/
'G10'-'G13X', 'G20'-'G22X', 'G254'-'G254X', 'G255'-'G255X', 'G312'-'G312X', 'G318'-'G318X', 'G319'-'G319X', 'G32'-'G32X', 
'G35'-'G37X', 'G40'-'G40XX', 'G41'-'G41X', 'G931'-'G931X', 'G934'-'G934X', 'R470'-'R470X', 'R56'-'R56X'
="ELX_009 Other Neurological"

/*Overlaps with Pulmonary Circulation diesease*/
'I278'-'I278X','I279'-'I279X','J40'-'J47XX','J60'-'J67XX','J684'-'J684X','J701'-'J701X','J703'-'J703X'
="ELX_010 CPD"


/*Overlaps with LIVER*/
'F10' - 'F10X', 'E52' - 'E52X', 
'G621'-'G621X', 'I426'-'I426X', 'K292'-'K292X', 'K700'-'K700X', 'K703'-'K703X', 'K709'-'K709X', 
'T51'-'T51X', 'Z502'-'Z502X', 'Z714'-'Z714X', 'Z721'-'Z721X'
="ELX_027 Alcohol abuse"

/*Overlaps with Psychoses*/
'F204'-'F204X',
'F313'-'F315X', 'F32'-'F32X', 'F33'-'F33X', 
'F341'-'F341X', 'F412'-'F412X', 'F432'-'F432X'
="ELX_030 Depression"

Other="Other"
;

RUN;


PROC SORT DATA=work._lookup_icd_Elixhauser_B
			OUT = work._lookup_icd_Elixhauser_B (KEEP=label RENAME=(label=ELIX_Cat_B)) NODUPKEY;
	WHERE FMTName='ICD_ELIXHAUSER_B';
	BY Label ;
RUN;

DATA work._lookup_ICD_Elixhauser_B;
	SET work._lookup_ICD_Elixhauser_B;
	IF ELIX_CAT_B NE 'Other';
RUN;

%MEND DefineICDCodes_CharlElix;

%MACRO ScoreCharlson_M3();
/*This macro implements the scoring of Charlson 
as per the M3 coding structure, which
includes scoring for cancer.*/

/*First merge the Charlson conditions 
onto the main M3 conditions/score file.*/

/*Check these are sorted first?*/

DATA _mm_char_scored;
SET _mm_char;

ARRAY all_CHAR{*} CHAR: ;

/*Step through all conditions in array above...*/
DO i = 1 TO dim(all_CHAR);
	/*If a record is blank (only happens when no hospital contacts for any condition for that patient)	*/
	/*set the value to zero instead.																	*/
	IF all_CHAR{i} EQ . THEN all_CHAR{i} = 0;
END;
DROP i;	

/*Now score Charlson...*/
		/********************************************/
		/*Largely adapted from Clare Salmond's code.*/
		/********************************************/
		*Calculate Charlson score, using weights from their Charlson et al. (1987);
		CharlsonScore = 
		/* Single weighted conditions */
		CHAR_MI  + CHAR_CHF + CHAR_PVD + CHAR_CBVD + 
		CHAR_Dementia + CHAR_CPD + CHAR_RD +  CHAR_PUD + 	 
		/*Higher weighted conditions...*/
		(CHAR_HemiPara * 2) + (CHAR_Renal * 2) + (CHAR_AIDS * 6) +
		/*Either/or option for diabetes, liver disease, cancer	*/
		MAX(CHAR_DiabNoComp, CHAR_DiabWithComp * 2)   +
		MAX(CHAR_LiverMild,  CHAR_LiverModSevere * 3) +
		MAX(CHAR_Cancer * 2 , CHAR_CancerMetastatic * 6)
		;
RUN;

%MEND ScoreCharlson_M3;

%MACRO ScoreElixhauser_M3();
/*This macro implements the scoring of Elixhauser
following the van Walraven weights pre condition.*/

/*Need to join these two files!*/

DATA _mm_elix_scored;
MERGE _mm_elix _mm_elixB;

BY &IDvarlist;/* Add in Elixhauser scoring...*/

/*This next replace-missing-with-zero step is needed 
for Elixhauser (but not Charlson) as Elixhauser is compiled in
		two steps, and a person *may* be missing in either set 
		(I don't think this can happen, but this is a safeguard)
	*/
ARRAY all_ELIX{*} ELIX: ;

/*Step through all conditions in array above...*/
DO i = 1 TO dim(all_ELIX);
	/*If a record is blank (only happens when no hospital contacts for any condition in list) 	*/
	/*set the value to zero instead.															*/
	IF all_ELIX{i} EQ . THEN all_ELIX{i} = 0;
END;
DROP i; /*Tidy up loop counter variable.*/

/*Now score Elixhauser.*/
ElixhauserScore = 	  ELIX_CHF * 7 +
					+ ELIX_Cardiac_arrythmia * 5					
					+ ELIX_Valvular_disease * -1
					+ ELIX_Pulmonary_circulation_disor * 4
					+ ELIX_PVD * 2
					+ ELIX_Hypertension_uncomplicated * 0
					+ ELIX_Hypertension_complicated * 0 
					+ ELIX_Paralysis *7 
					+ ELIX_Other_Neurological * 6
					+ ELIX_CPD * 3 
					+ ELIX_DiabNoComp  * 0
					+ ELIX_DiabWithComp * 0
					+ ELIX_Hypothyroidism  * 0 
					+ ELIX_Renal  * 5
					+ ELIX_Liver * 11
					+ ELIX_PUD * 0
					+ ELIX_AIDS * 0 
					+ ELIX_Lymphoma * 9
					+ ELIX_CancerMetastatic * 12 
					+ ELIX_Cancer_solid_no_metastasis * 4
					+ ELIX_Rheumatoid_arthritis * 0
					+ ELIX_Coagulopathy * 3 
					+ ELIX_Obesity * -4
					+ ELIX_Weight_loss * 6
					+ ELIX_Fluid_and_electrolyte_disor * 5
					+ ELIX_Blood_loss_anemia * -2
					+ ELIX_Deficiency_anemia * -2
					+ ELIX_Alcohol_abuse * 0
					+ ELIX_Drug_abuse * -7
					+ ELIX_Psychoses * 0
					+ ELIX_Depression  * -3 
					;
RUN;

%MEND ScoreElixhauser_M3;


%MACRO Merge_M3_finalset(passoutdata=, patient_id=, passkeepconditions=,
						 passAddCharlson=, passAddElixhauser=);

DATA &passoutdata;
%IF &passAddCharlson EQ 0 AND &passAddElixhauser EQ 0 %THEN %DO;
	SET 										  %END;
%IF &passAddCharlson EQ 1 OR  &passAddElixhauser EQ 1 %THEN %DO;
	MERGE                                         %END;
	/*Add base of M3 for everyone*/
	_mm_m3_scored
%IF &passAddCharlson EQ 1   %THEN %DO;
	_mm_char_scored		    %END;
%IF &passAddElixhauser EQ 1 %THEN %DO;
	_mm_elix_scored		    %END;
; *Close line off for set/merge statement here;
/*Following not strictly needed if going with just the M3 index*/
BY &patient_id;


ARRAY all_conds {*} 
				M3_: M3Score /* Always include M3 conditions*/
				%IF &passAddCharlson   EQ 1   %THEN %DO; CHAR_: CharlsonScore   %END;
				%IF &passAddElixhauser EQ 1   %THEN %DO; ELIX_: ElixhauserScore %END;
				;
/*Step through all conditions in array above...*/
DO i = 1 TO dim(all_conds );
	/*If a record is blank (only happens when no hospital contacts for any condition for that patient)	*/
	/*set the value to zero instead.																	*/
	IF all_conds{i} EQ . THEN all_conds{i} = 0;
END;
DROP i;

%IF &passkeepconditions EQ 0 %THEN %DO;
	DROP M3_ : CHAR_: ELIX_: ;
%END;
RUN;
%MEND Merge_M3_finalset;

/*************************/
/*End of M3 specific code*/
/*************************/

/* End of macro file*/
