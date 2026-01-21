------------------ 
TABLE OF CONTENTS

* GENERAL INFORMATION
* FILE ORGANIZATION
* METHODOLOGICAL INFORMATION 
* DATA DICTIONARY
* SHARING AND ACCESSING INFORMATION
* ADDITIONAL NOTES/COMMENTS

--------------------
GENERAL INFORMATION

1- DATASET TITLE: 
Individual participant data (IPD) corresponding to the article:
“Combinatorial benefit without synergy in recent clinical trials of immune checkpoint inhibitors” 
by Adam C. Palmer, Benjamin Izar , Haeun Hwangbo and Peter K. Sorger. 
https://www.medrxiv.org/content/10.1101/2020.01.31.20019604v2

2- AUTHORS AND AFFILIATIONS:
Name: Adam C. Palmer 
Organization/institution: University of North Carolina at Chapel Hill
Email: palmer@unc.edu

Name: Peter K. Sorger 
Organization/institution: Harvard Medical School
Email: peter_sorger@hms.harvard.edu (cc: sorger_admin@hms.harvard.edu)

3- DATE OF DATA COLLECTION:
2021-10-19

-------------------
FILE ORGANIZATION

1- FILE NAMING:
Files are found in folders labeled by disease setting and drug combination.
The number preceding the folder name corresponds to the identifier in the metadata file. 
Each file follows the following naming convention:
CombinationName.csv

2- FILE SUMMARY COUNT: 
15 IPD files and 1 metadata file

3- FILE LIST :
Trial information summarized in metadata file:
IPD metadata.xlsx
    
4- RELATIONSHIP BETWEEN FILES:
Each .csv file contains IPD from a different figure from a published clinical trial. 

---------------------------
METHODOLOGICAL INFORMATION

1- DESCRIPTION OF METHODS USED FOR COLLECTION/GENERATION OF DATA:
For a description on the individual participant imputation, please see the following publications: 
Guyot, P., Ades, A., Ouwens, M.J., and Welton, N.J. (2012). Enhanced secondary analysis of survival data: reconstructing the data from published Kaplan-Meier survival curves. BMC Med Res Methodol 12, 9.
Fell, G., Redd, R.A., Vanderbeek, A.M., Rahman, R., Louv, B., McDunn, J., Arfè, A., Alexander,B.M., Ventz, S., Trippa, L. (2021). KMDATA: a curated database of reconstructed individual patient-level data from 153 oncology clinical trials, Database, Volume 2021, 2021, baab037.

2-  PEOPLE INVOLVED WITH SAMPLE COLLECTION, PROCESSING, ANALYSIS AND/OR SUBMISSION:
Adam C. Palmer, Benjamin Izar , Haeun Hwangbo and Peter K. Sorger- data collection, analysis, article writing, and editing
Deborah Plana- data organization for release

----------------
DATA DICTIONARY

1- VARIABLE LIST: 
IPD files:
Time	Timing of progression or censoring events. Time presented in the same units as reported in the original trial publication.
Event	Censoring (0) or progression (1).
Arm	Name of treatment arm. 

Metadata file IPD metadata.xlsx

--------------------------------
SHARING & ACCESSING INFORMATION

1- LICENSES/RESTRICTIONS PLACED ON THE DATA:
CC BY-SA 4.0

2- LINKS TO OTHER PUBLICLY ACCESSIBLE LOCATIONS OF THE DATA:
cancertrials.io

3- RECOMMENDED CITATION FOR THIS DATASET:
-cancertrials.io
-“Combinatorial benefit without synergy in recent clinical trials of immune checkpoint inhibitors” 
by Adam C. Palmer, Benjamin Izar , Haeun Hwangbo and Peter K. Sorger. 
https://www.medrxiv.org/content/10.1101/2020.01.31.20019604v2

--------------------------
ADDITIONAL NOTES/COMMENTS

Please let Adam Palmer (palmer@unc.edu) know if any errors are found in this data.  