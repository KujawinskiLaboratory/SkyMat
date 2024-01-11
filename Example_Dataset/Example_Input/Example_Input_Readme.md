# Example Input File Read Me

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

## Skyline Quantification Tables require specific export columns in order for SkyMat to work. 
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
 
## Sequence Sheet format is important in order for SkyMat to work. Columns where value input is essential are described below:
* `Sample Type` column should be used accordingly:
   * Blank - any blank that you don't want in your final output table, blanks are solely used for assessing contamination and carryover in Skyline
   * Std Bracket - this should be used to denote any of your standard curve samples (including your "zero standard")
   * Unknown - any pooled QC sample, replicate samples, or blank that you want exported in your final quantification table.
* `Sample Name` column
   * For Skymat to conduct polarity merging sample names must be identical for both positive and negative mode, excludiing their suffix which would be ` pos` or ` neg` respectively. This syntax is looked for in SkyMat when conducting the polarity merging, if there are not equal numbers of paired samples in positive and negative mode the code will throw an error.
     * Example from `SkyMat_3isotopes_test_pos_and_neg.xlsx`: `Syn_10-30_filtrate_B1 neg` (row #20) and `Syn_10-30_filtrate_B1 pos` (row #42)
     * Note: The code will not break if pooled samples are not named in this convention, however the word "pool" must be included in the `Sample Name` for pooled samples which riSkyline_${SILISType}.m will then rename in a separate column to ensure they match downstream (Example if sample is solely labeled "pool QC" for all injections, riSkyline_${SILISType}.m will rename the first instance of "pool QC" to "pool01" in positive and negative mode, respectively and use as a pair.
* ionMode - `pos` or `neg`
* goodData - Boolean argument where samples listed as `1` are included in downstream SkyMat analyses and those listed as `0` are excluded.
   * Example of files that should be exluded and listed as `0` incuded: column conditioning samples, samples where an injection failure occured, or a system suitability standard.  
* sType - sample type inputs include:
   * rep - unknown/study samples or "replicates"
   * blank - MQ or other 'blank' samples excluding the 0 concentration sample included in the standard curve
   * std - standard curve samples
   * pool - pooled QCs
* runOrder - column with a numerical indicator of the order in which the samples were injected onto the LC-MS instrument

## Where to use input files: 
* Input files are used in both riSkyline_C13.m and riSkyline_D5.m depending on which isotope(s) you have. In this example, we have both isotopes so both codes will be used. 
* riSkyline_${SILISType}: Section 2 (line 16)
  * fName = SkyMat_3isotopes_test_pos_and_neg.xlsx
* riSkyline_${SILISType}: Section 3 (line 23)
   * dfile_pos = SkyMat_3isotopes_test_pos.csv
   * dfile_neg = SkyMat_3isotopes_test_neg.csv
* riSkyline_${SILISType}: Section 9 (line 270)
  * tFile = string([tDir filesep 'TransitionList_SkyMat_Example.xlsx']);
