#Example Input File Read Me

## Input Files: 
1. SkyMat_3isotopes_test_neg.csv
* Quantification table containing areas for light and heavy isotopes and transitions exported from Skyline for the negative polarity mode LC-MS/MS data aquisition.  
2. SkyMat_3isotopes_test_pos.csv
* Quantification table containing areas for light and heavy isotopes and transitions exported from Skyline for the positive polarity mode LC-MS/MS data aquisition.
3. SkyMat_3isotopes_test_pos_and_neg.csv
* Sequence sheet exported from the LC-MS instrument with added columns: ionMode (pos/neg), goodData (0/1), sType (blank,rep,std,pool), runOrder (numeric)
* Column headers include: Sample Type, File Name, Path, Instrument Method, Position, Inj Vol, Sample Name, ionMode, goodData, sType, runOrder
4. TransitionList_SkyMat_Example.xlsx - Transition List 
* Transition list containing metabolites and their precursor masses, fragments, molecular weights, etc. This file should be used in Skyline for determine your integration peak list and for SkyMat convertMoles.m function

## Note: Skyline Quantification Tables require specific export columns in order for SkyMat to work. 
* Required columns exported from Skyline: 
  * Replicate Name
  * File Name - name of file(s) imported into Skyline for peak integration.
  * Sample Type - Double Blank, Blank, Quality Control, Standard, or Unknown
  * Isotope Label Type - light,heavyD5,heavyC13
  * Fragment Ion Type - custom (fragment), precursor
  * Transition Note - light/heavy
  * Analyte Concentration - concentration of standards
  * Molecule List Name - name of group specific molecules are grouped into (e.g., molecules1) 
  * Molecule Name - name of molecule
  * Precursor Mz - mass:charge ratio of the precursor ion
  * Product Mz - mass:charge ratio of the fragment ion
  * Points Across Peak - number of points across each integrated peak
  * Background - background noise level calculated by Skyline
  * Area - Area of integrated peak in Skyline
  * Fwhm - full width half max of each integrated peak in Skyline

## Where to use input files: 
* Input files are used in both riSkyline_C13.m and riSkyline_D5.m depending on which isotope(s) you have. In this example, we have both isotopes so both codes will be used. 
* riSkyline_C13: Section 2 (line 16)
  * fName = SkyMat_3isotopes_test_pos_and_neg.xlsx
* riSkyline_C13: Section 3 (line 23)
  * dfile_pos = SkyMat_3isotopes_test_pos.csv
  *dfile_neg = SkyMat_3isotopes_test_neg.csv
* riSkyline_C13: Section 9 (line 270)
  * tFile = string([tDir filesep 'TransitionList_SkyMat_Example.xlsx']);
