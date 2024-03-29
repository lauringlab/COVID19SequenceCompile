# README for Compilation Pipeline:

### Set-Up:

<b>Folder Structure:</b>

Folder structure lives within the MED-LauringLab/SEQUENCING DropBox folder.

* INFLUENZA_A
  * 0_ProcessingScripts
  * 2_Plate_Maps
  * 3_ProcessedGenomes
  * 4_SequenceSampleMetadata
  * 5_GISAID_Uploads
  * 6_SampleSheets
  * 7_BarcodeMaps
* SARSCOV2
  * 00_ARCHIVE_SarsCoV2_2020: Contains sample analysis from early 2020 sequencing work
  * 0_ProcessingScripts: Contains processing code used in original sequence processing
  * 1_OriginalManifests_CBR: Original manifest documents provided by CBR
  * 1_OriginalManifests_Lynx: Original manifest documents provided by LynxDx
  * 1_OriginalManifests_SPH: Original manifest documents provided by SPH
  * 1_OriginalManifests_UHS: Original manifest documents provided by UHS
  * 2_PlateMaps: Original plate map files are first placed here by laboratory team
  * 3_ProcessedGenomes: Original sequence files are placed here, organized by plate.
    * Within each plate folder = .all.consensus.final.fasta (fasta file of sequences with >90% completeness), all.consensus.renamed.full.fasta (fasta file of all sequences)
  * 4_SequenceSampleMetadata
      * FinalSummary: Contains the final compiled data full_compiled_data.csv
          * ReportNotifications: Contains negative_control_warnings.tsv
          * secret: Contains a duplicate final file (full_compiled_data.csv) for comparison in case of changes made manually to the compiled file
          * CDC_IVY_UPLOADS: Contains the files created in order to upload to the CDC IVY RedCap database
          * RVTN_UPLOADS: Contains the files created in order to upload to the RVTN RedCap database
      * Manifests: contains all manifest files, organized by received source
          * ASCENSION: Contains manifests from Ascension Hospital system, one of the MI-SAPPHIRE sites
          * CBR: Contains manifest files from the University of Michigan Campus BioRepository
          * CDCIVY: Contains manifest files from the CDC IVY project
              * Full_IVY_Set: Contains IVY_sample_full_manifest_list.csv, the compiled CDC IVY manifest lists with all variables sent (a subset of these are kept for the final compiled sequence list)
              * Keys: Contains CDC_SiteCodebook.csv, the list of CDC IVY sites, their locations, and their coded name
          * CSTP: Contains manifest files from the COVID-19 Sampling & Tracking Program (ultimately from LynxDx)
          * EDIDNOW: Contains manifest files from the ED IDnow testing program at Michigan Medicine
          * HENRYFORD: Contains manifests from the Henry Ford Health System, one of the MI-SAPPHIRE sites
          * IVYIC: Contains manifests from and IVY4 sub study on immunocompromised individuals
          * Martin: Contains manifest files from the Martin Lab at the School of Public Health
          * MDHHS: Contains manifests from the Michigan Department of Health and Human Services
          * PUIMISC: Contains manifests for one-off projects, PUIs, etc.
          * RVTN: Contains manifest files from the RCTN project
          * RVTN_TEST: Contains manifests from initial RVTN test round
          * TRINITY: Contains manifests ffrom the Trinity Health System, one of the MI-SAPPHIRE sites
          * UHS: Contains manfests from University Health Services samples
          * VIEW: Contains manifest files from the VIEW study
          * ManifestsComplete: Contains the compiled manifest list sample_full_manifest_list.csv, the archived manifest list sample_full_manifest_list_prev.csv, and the list of archived manifest files as current_manifest_list.RDS
      * PlateMaps: Contains all plate map files as YYYYMMDD_SC2_Illumina_Plate_##.xlsx or YYYYMMDD_SC2_Nanopore_Run_##.xlsx
          * PlateMapsComplete: Contains the compiled plate map list sample_full_plate_list.csv
      * PreviousLists: Contains the original processed sample list ProcessedSampleCumulativeList_20210326.csv maintained by Andrew Valesano prior to code implementation
      * SequenceOutcomes
          * gisaid: Contains the GISAID metadata file metadata_YYYY-MM-DD_MM-DD.tsv; also contains the GISAID upload template file (GISAID_UPLOAD_TEMPLATE.xls)
          * nextclade: Contains the NextClade output files as either YYYYMMDD_Plate_##\_##\_nextclade.tsv or YYYYMMDD\_Nanopore\_Run_##_nextclade.tsv
          * pangolin: Contains the Pangolin output files as either YYYYMMDD_Plate_##\_##\_pangolin.csv or YYYYMMDD\_Nanopore\_Run_##_pangolin.csv; also contains the code and file(s) for the updated current pangolin status system
          * SequenceOutcomeComplete: Contains the sequence outcome final compiled files of sample_full_gisaid_list.csv, sample_full_nextclade_list.csv, and sample_full_pangolin_list.csv
    * 5_GISAID_Uploads: contains each plate's necessary files to submit for GISAID submission
    * 6_Sequencing_ReportsAndCharts
    * 7_QAQC_CDCIVY
    * 8_QPCR_IVY
    * 9_Coverage
    * 10_transfer
    * 11_GridION_RunReports
    * 12_SampleSheets

<b>R Libraries Necessary:</b>

* tidyverse: https://www.tidyverse.org/
* lubridate: https://lubridate.tidyverse.org/
* janitor: https://cran.r-project.org/web/packages/janitor/index.html
* withr: https://cran.r-project.org/web/packages/withr/index.html
* openxlsx: https://cran.r-project.org/web/packages/openxlsx/index.html
* reshape2: https://cran.r-project.org/web/packages/reshape2/index.html

<b>Python Libraries Necessary:</b>

* import argparse: https://docs.python.org/3/library/argparse.html
* import glob: https://docs.python.org/3/library/glob.html
* import pandas as pd: https://pandas.pydata.org/
* from Bio import SeqIO: https://biopython.org/wiki/SeqIO
* from Bio.Seq import Seq: https://biopython.org/docs/1.75/api/Bio.Seq.html
* from Bio.SeqRecord import SeqRecord https://biopython.org/docs/1.75/api/Bio.SeqRecord.html
* import os: https://docs.python.org/3/library/os.html

<b>Code:</b>

In order to run the compilation code pipeline, download the following code sets:

* <span style="background-color: #8AA29E">pipeline_functions.R</span>
* manifest_code.R
* plate_map_code.R
* pangolin_code.R
* nextclade_code.R
* gisaid_code.R
* determine_current_manifest_state.R
* compile_components_code.R
* <span style="background-color: #F0B7B3">full_run_code.R</span>
* <b>OutsidePipeline</b>
  * cdc_ivy_upload_code.R
  * checking_compiled_files.R
  * gisaid_upload_file_creation.R
  * moving_nextclade_output.R
  * moving_pangolin_output.R
  * moving_plate_map_files.R
  * rvtn_upload_code.R
  * subset_compiled_for_fasta.R
  * <span style="background-color: #8AA29E">summary_report_samples.Rmd</span>
  * <span style="background-color: #8AA29E">who_variant_crosswalk.csv</span>
    * <b>ProcessingFASTA</b>
      * checking_sampleids.py
      * prep_fasta_NumberTwo.py
      * prep_fasta_for_gisaid.py
      * prep_fasta_gisaid_flu.py
      * prep_fasta_NumberOne_2.py
      * replace_barcodes_fasta.py
      * SelectSequences.py
    * <b>UpToDatePangolin</b>
      * compile_all_fastas.py

And put them all in the same folder structure on your computer. Use full_run_code.R to run everything in the correct order. Ensure that all file paths and folder names mentioned in each code set are correct.

<b>Additional Notes:</b>

* At the beginning of full_run_code.R, there is a section of code to set the ```starting_path```, ```code_path```, ```batch_path```, and ```influenza_path``` for your particular computer. This only works for adjusting the portions of the necessary code that are called from the full_run_code.R file. In any instance, always check all file paths before running code files.
* Ensure that you have access to either the LauringLab DropBox folder, or your own version of the pipeline folder structure, from your computer.

---

### Data Dictionary for SARS-CoV-2 Compiled File:

The final created file is called <b>full_compiled_data.csv</b>.

| Column | Source | Data Style | Join Variable | Definition |
| --- | --- | --- | --- | --- |
| sample_id | sample_full_manifest_list.csv | numeric | Yes | Unique identification number for a particular sample; One individual (subject_id) can have multiple samples, each with their own sample_id number |
| subject_id | sample_full_manifest_list.csv | numeric | | Unique identifcation number for an individual; received source determines the type. CSTP = umid, Martin = other, UHS = MRN, CBR = MRN |
| coll_date | sample_full_manifest_list.csv | date | | Date the sample was collected from the individual |
| flag | sample_full_manifest_list.csv | character | | Notes for each sample, some warnings introduced by the pipeline |
| received_source | sample_full_manifest_list.csv | character | | Place the sample was received from; options include = UHS, CBR, CSTP, Martin, etc. |
| SampleBarcode | sample_full_plate_list.csv | character | | Barcode type identifier for a particular sample; Not all samples have a barcode style identifier |
| PlateDate | sample_full_plate_list.csv | date | | Date the plate was created on |
| PlatePlatform | sample_full_plate_list.csv | character | | Testing system the plate was run through; Options include = Nanopore, Illumina, etc. |
| PlateNumber | sample_full_plate_list.csv | character | | ##_## format, indicates what distinct plate the sample was run with |
| pangolin_lineage | sample_full_pangolin_list.csv | character | | Strain type of the sample; Ex. B.1.1.7 |
| pangolin_probability | sample_full_pangolin_list.csv | numeric | | Probability of correct lineage assignment (pangolin_lineage) by PANGOLIN. Archived as of 3 May 2021 update |
| pangolin_status | sample_full_pangolin_list.csv | character | | Pangolin classifier assigns either “passed_qc” or “fail” to the sequenced SARS-CoV-2 strain based on sequence completeness and quality. If the consensus SARS-CoV-2 sequence has 50 percent or more ambiguous bases (50%N's) in the consensus fasta or if the consensus length is less than or equal to 10,000bp only, the status will change to fail. |
| pangolin_note | sample_full_pangolin_list.csv | character | | Additional information provided by the pangolin system |
| nextclade_clade | sample_full_nextclade_list.csv | character | | The global clade identified for the sequenced sample |
| nextclade_totalMissing | sample_full_nextclade_list.csv | numeric | | The number of bases missing from the nextclade analysis |
| nextclade_completeness | sample_full_nextclade_list.csv | numeric | | Percentage; how complete the nextclade coverage was of the sample genome |
| gisaid_strain | sample_full_gisaid_list.csv | character | | Virus name, as listed on GISAID. Format is USA/MI-UM-sample_id/YYYY. |
| gisaid_epi_isl | sample_full_gisaid_list.csv | character | | GISAID database accession number. General format is EPI_ISL_NNNNNN. |
| gisaid_clade | sample_full_gisaid_list.csv | character | | Clade, as found in the GISAID metadata download |
| gisaid_pango_lineage | sample_full_gisaid_list.csv | character | | Pangolin lineage, as found in the GISAID metadata download |
| received_date | sample_full_manifest_list.csv | date | Yes | Date the sample was received at the lab from the received_source |
| position | sample_full_manifest_list.csv | character | | Position the sample was in, in the box received from the received source |
| SiteName | sample_full_manifest_list.csv | character | | Research site where the sample was collected from; Only applicable for the CDC IVY study |
| subject_id_length | sample_full_manifest_list.csv | numeric | | Number of characters in the subject id; Calculated before a leading zero is added, so in some cases may be one short from the current length |
| PlateName | sample_full_plate_list.csv | character | | Full plate name of the sample test run; generally corresponds to plate file name |
| PlatePosition | sample_full_plate_list.csv | character | | Position the sample was in, in the plate that it was tested on |
| SampleSourceLocation | sample_full_plate_list.csv | character | | Where the sample came from (should correlate with received source) |
| pangoLEARN_version | sample_full_pangolin_list.csv | date | | The version date of the multinomial logistic regression model that was used for global lineage assignment. |
| pangolin_conflict | sample_full_pangolin_list.csv | numeric | | New as of pangolin software update (notified on 3 May 2021);  the number of conflicts in the algorithm's decision tree |
| pango_version | sample_full_pangolin_list.csv | character | | New as of pangolin software update (notified on 3 May 2021); additional software version information |
| pangolin_version | sample_full_pangolin_list.csv | character | | New as of pangolin software update (notified on 3 May 2021); additional software version information |
| pangolin_runDate | sample_full_pangolin_list.csv | date | | Date of pangolin output receipt; pulled from the pangolin filename the sample_id was reported on |
| PlateToPangolin_days | {calculated in compile_components_code.R} | numeric | | Number of days between pangolin_runDate and PlateDate; Put in place to add more context for samples when joining multiple sources, as some sample aliquots are received more than once from the same source |
| nextclade_qcOverallScore | sample_full_nextclade_list.csv | numeric | | Quality control overall score (smaller is better) |
| nextclade_qcOverallStatus | sample_full_nextclade_list.csv | character | | Quality control overall status |
| nextclade_totalMutations| sample_full_nextclade_list.csv| numeric | | Total number of nucleotide substitutions relative to Wuhan/Hu-1 reference. Does not include ambiguous nucleotides (see below). |
| nextclade_totalNonACGTNs | sample_full_nextclade_list.csv | numeric | | Number of ambiguous nucleotide characters, i.e. not A, C, T, G, or N. |
| nextclade_runDate | sample_full_nextclade_list.csv | date | | Date of nextclade output receipt; pulled from the nextclade filename the sample_id was reported on |
| PlateToNextclade_days | {calculated in compile_components_code.R} | numeric | | Number of days between nextclade_runDate and PlateDate; Put in place to add more context for samples when joining multiple sources, as some sample aliquots are received more than once from the same source |
| sample_per_subject | {calculated in compile_components_code.R} | numeric | | Numbering of sample_ids per subject_id of the listed sample; samples are ordered by collection date and assigned a number based on 1st, 2nd, 3rd, etc. sampling |
| multiSamples | {calculated in compile_components_code.R} | numeric | | 1,0 binary; has a value of 1 if the max sample_per_subject of a given sample's subject_id is greater than 1, has a value of 0 if the subject_id has only one sample_id associated with it. |
| daysFromPrevious | {calculated in compile_components_code.R} | numeric | | Number of days since the previous sample coll_date to the current sample coll_date |
| ninetyDayFromPrevious | {calculated in compile_components_code.R} | numeric | |  binary 1,0; if daysFromPrevious is greater than 90, then is equal to 1, otherwise it is 0 |
| previousLineageDifferentThanCurrent | {calculated in compile_components_code.R} | numeric | | binary 1,0; If the previous sample's pangolin_lineage value is different than the current pangolin_lineage, then has a value of 1, otherwise is 0 |
| previousCladeDifferentThanCurrent | {calculated in compile_components_code.R} | numeric | | binary 1,0; If the previous sample's nextclade_clade value is different than the current nextclade_clade, then has a value of 1, otherwise is 0 |
| sample_id_lauring | {calculated in compile_components_code.R} | character | | "L" followed by 13 digits, used to hide sample ids in public sequence information  when necessary |
| data_quality_rule | {calculated in compile_components_code.R} | character | | "pass" or "not passed"; in order to pass, the sequence must have (pangolin_status %in% c("pass", "passed_qc")) & (nextclade_qcOverallStatus %in% c("good", "mediocre")) & (nextclade_completeness > 80) |
| newest_pangolin_lineage | {compiled in compile_components_code.R} | character | | updated pangolin lineage based on newest available pangolin program |
| newest_pangolin_date | {calculated in compile_components_code.R} | numeric | | date in YYYYMMDD format of when the sample's fasta file was re-run through pangolin |


---

<b>Additional Information</b>
* Latest epidemiological lineages of SARS-CoV-2: https://cov-lineages.org/index.html
* Additional Definition Source: https://perkinelmer-appliedgenomics.com/wp-content/uploads/marketing/Coronavirus/NEXTFLEX_Variant-Seq_SARS-CoV-2_Software-from-CosmosID.pdf
* Nextclade GitHub: https://github.com/nextstrain/nextclade
* CDC's Monitoring of Variant Proportions : https://covid.cdc.gov/covid-data-tracker/#variant-proportions

---

### SARS-CoV-2 Processing Code Order:

Use full_run_code.R to run all the pieces of the pipeline in order. This order is outlined within the full_run_code.R file. These steps are also outlined within the [Google Tracking Document for SARS-CoV-2](<https://docs.google.com/spreadsheets/d/1GuPIPou3Y15_TH2cZbNJ1Y6BLTHNllD-2yvwmuPhfEM/edit#gid=744361978>)

---

### Manifests

Manifests are received from the following sources:

* Martin Lab at the University of Michigan School of Public Health (Martin)
* University of Michigan Central Biorepository (CBR)
* CDC IVY Project (CDCIVY) - samples from 21 sites sent to Vanderbilt, then to University of Michigan (Currently IVY4 as of Sept. 2022)
* CDC IVY Sub-study of Viral Evolution in Immunocompromised Individuals (IVYIC)
* University Health Services (UHS)
* CDC Respiratory Virus Transmission Network (RVTN)
* Michigan Sequencing Academic Partnership for Public Health Innovation and Response (MI-SAPPHIRE) - samples from Henry Ford, Trinity, and Ascension health systems
* VIEW Project

Manifests are reviewed, checked and renamed as necessary, and placed in the appropriate Manifests folder [within DropBox/MED-LauringLab/4_SequenceSampleMetadata/Manifests]. Manifest documentation to exchange with new testing sites, go to Dropbox (University of Michigan)/MED-LauringLab/LabManagement/Documentation_Notes/ManifestFileFormatting.

##### Basic Manifest Column Format

| Columns | Data Type	| Variable Description |
| --- | --- | --- |
| position | character | Where the sample is located in the box sent |
| sample_id | numeric | Identification number for the sample; unique |
| subject_id | numeric | Identification number for the subject (individual) - there may be multiple samples per subject; Should not be identifiable; Should be available with barcode scan |
| coll_date | date | Date the sample was collected on; M/D/YY format |
| flag | character | Notes |

##### Manifest Column Format (for IVY & RVTN)

Note: these manifests are generally send as Excel files.

| Columns | Data Type	| Variable Description |
| --- | --- | --- |
| Box # | Numeric | The identifying number on the box sent to the laboratory |
| Position # | Numeric | Location of the individual sample within the shipping box |
| Site Name | Character | Hidden; Can be left blank, as site can be determined from study ID, filled with site name the sample came from |
| Study ID | Numeric | Identifying subject number |
| Collection Date | Date | Date the sample was collected on |
| Aliquot ID | Alpha-numeric | Sample ID, non identifying |
| Specimen Type | Character | Type of sample collected (e.g. "Upper Respiratory") |
| FluSC2 RNP Ct | NA | Hidden, blank |
| SC2 Ct | NA | Hidden, blank |
| RNP Ct | Numeric | Ct value for the sample |
| Covid-19 N1 | Numeric	| Ct value for the sample |
| Covid-19 N2 | Numeric | Ct value for the sample | 	
| Comments | Character | Any notes on the sample, if provided; Can be blank | 


##### Manifest File Name Format

nameOfSource_YYYYMMDD_#.csv or nameOf Source_YYYYMMDD_#.xlsx

nameOfSource = Where the samples came from, corresponds to the name of the project folders inside the Manifests folder (CBR, CSTP, Martin, CDCIVY, EDIDNOW, etc.)

YYYYMMDD = Year, Month, and Day of when the samples arrived/the manifest file was received

\# = Number of the manifest; Will usually be a "1", but if two batches of samples arrive on the same day from the same source, with two separate associated manifest files, then these would be numbered accordingly

---

##### Samples

Samples are received, and sequenced on plates using Nanopore or Illumina systems. Plate Map files are generated, matching samples to their location on those plates. Original plate map files are placed in [DropBox/MED-LauringLab/Plate Maps] as .xlsx files.

---

##### Plate Map File Column Format

| Columns | Data Type	| Variable Description |
| --- | --- | --- |
| Processing Plate | character | Name of the plate. Should match the file name |
| Slot | numeric (int) | Well position on the plate that the sample was placed in |
| Sample ACCN | character/numeric | Corresponds to sample_id |
| Sample MRN | character/numeric | If filled in, corresponds to subject_id |
| Sample Order# | blank | blank |
| Barcode | character | Barcode information, will be output in .fastq file to identify samples |
| Source | character | Where the sample came from; In general, should correspond to the Manifest source names |

These columns will occupy columns A-G of the Excel file. Columns H-T contain the sample_id information laid out in the format of the plate itself. This information is ignored in data processing, but can be a good point of reference for Quality Assurance checking.

---

##### Plate Map File Name Format

YYYYMMDD_sequenceSystem_#.xlsx

YYYYMMDD = Year, Month, and Day of when the Plate Map file was created and the samples began the sequencing process.

sequenceSystem = Nanopore_Run or Illumina_Plate, depending on whichever system was used. If a new system is introduced, the main system information should be first. The character string after the "_" character is occasionally dropped.

\# = The number assigned to the testing plate.

---

These Nanopore or Illumina systems generate .fastq formatted files, containing the raw sequence data for each sample. There is a python script on the lab computer that converts these files to .fasta files, with a ">" symbol marking the beginning of each new sequence. After that symbol will be a string denoting the barcode value for the corresponding sample, and the sequence will proceed on a new line following that information.

~~The .fasta files generated are called <b>plateMapName.all.consensus.fasta & plateMapName.all.consensus.final.fasta</b> and are placed in the corresponding plateMapName folder within [DropBox/MED-LauringLab/ProcessedGenomes/].~~

~~plateMapName.all.consensus.fasta = the sequence information for all samples in a given plate run. Some samples may not appear if there was essentially no matching genetic material in the sample.~~

~~plateMapName.all.consensus.final.fasta = filtered version of plateMapName.all.consensus.fasta (deprecated 26 May 2021)~~

As of August 19, 2021, the .fasta files that are generated are called <b>plateMapName.all.consensus.final.fasta & plateMapName.all.consensus.renamed.full.fasta</b>. These now have the ">" separators use the sample_id as the identifier, rather than the barcode identifier.

These files are placed in the corresponding plateMapName folder within [DropBox/MED-LauringLab/ProcessedGenomes/].

plateMapName.all.consensus.renamed.full.fasta = the sequence information for all samples in a given plate run. Some samples may not appear if there was essentially no matching genetic material in the sample.

plateMapName.all.consensus.final.fasta = filtered version of plateMapName.all.consensus.fasta (not used in any subsequent analyses)

---

~~There is another python script that takes the generated .fasta file and replaces the barcode string with the matching sample ID string. Once this .fasta file has been created, it will be used in the Pangolin and NextClade systems [prep_fasta_NumberOne.py].~~

~~<b>Steps:</b>~~
1. ~~Compile full data set using pipeline code.~~
2. ~~Ensure there is a folder in [DropBox/MED-LauringLab/ProcessedGenomes] named the same as the originating PlateMap file~~
3. ~~In that folder, there will be two files: plateMapName.all.consensus.fasta & plateMapName.all.consensus.final.fasta; only plateMapName.all.consensus.fasta is necessary.~~
4. ~~Using the subset_compiled_for_fasta.R code, subset the newly generated full_compiled_data.csv file to only contain the relevant plate samples, copy it here [DropBox/MED-LauringLab/ProcessedGenomes/plateMapName] and rename it to plateMapName.meta.csv (This is because the barcode locations repeat on sample runs, so this is necessary to ensure proper barcode to sample_id matching)~~
5. ~~Run prep_fasta_NumberOne.py~~

~~This creates the following file:~~
* ~~plateMapName.all.consensus.renamed.full.fasta~~

~~That should be placed into [DropBox/MED-LauringLab/ProcessedGenomes/plateMapName]~~

<b>New Steps (now that prep_fasta_NumberOne.py no longer needs to be run):</b>

1. Compile full data set using pipeline code [full_run_code.R]
2. Ensure there is a folder in [DropBox/MED-LauringLab/ProcessedGenomes] named the same as the originating PlateMap file
3. In that folder, there will be two files: plateMapName.all.consensus.final.fasta & plateMapName.all.consensus.renamed.full.fasta; only plateMapName.all.consensus.renamed.full.fasta is necessary
4. Using the subset_compiled_for_fasta.R code, subset the newly generated full_compiled_data.csv file to only contain the relevant plate samples. Check for missing manifest entries [by ensuring that there are no missing coll_date entries in seq_list2]. This code will write out a copy of that information to [DropBox/MED-LauringLab/ProcessedGenomes/plateMapName] as plateMapName.meta.csv.

_At this point, the [DropBox/MED-LauringLab/ProcessedGenomes/plateMapName] folder should contain plateMapName.all.consensus.final.fasta, plateMapName.all.consensus.renamed.full.fasta, and plateMapName.meta.csv._

---

##### Pangolin

Instructions for installing the command line version of Pangolin can be found here: https://cov-lineages.org/pangolin_docs/installation.html

Note: This can only be completed using Linux or OS operating systems. (The environment.yml file to install pangolin requires minimap2 and gofasta libraries, both of which are only built for those systems, not for Windows) If you only have a Windows operating system, you can install a Windows Subsystem for Linux (directions here: https://docs.microsoft.com/en-us/windows/wsl/install-win10) which will allow you to proceed with the pangolin installation process.

Pangolin can also be used at the following site: https://pangolin.cog-uk.io/

Steps:
1. In command line, navigate to [DropBox/MED-LauringLab/ProcessedGenomes/plateMapName]
2. Activate the pangolin environment ```conda activate pangolin```
3. Run the command ``` pangolin plateMapName.all.consensus.renamed.full.fasta```
4. This will output a file called lineage_report.csv in [DropBox/MED-LauringLab/ProcessedGenomes/plateMapName]
5. Rename this file to plateMapName_pangolin.csv
6. Copy this renamed file to [DropBox/MED-LauringLab/SequenceSampleMetadata/SequenceOutcomes/pangolin]

---

##### Next Clade

Next Clade can be accessed at this site: https://clades.nextstrain.org/

After the process is complete, the output should be downloaded as a .tsv file.

Steps:
1. Navigate to https://clades.nextstrain.org/
2. In the SARS-CoV-2 box in the lower right-hand corner, ensure the "From File" tab is selected and either drag & drop or select the <b>plateMapName.all.consensus.renamed.full.fasta</b> from the [DropBox/MED-LauringLab/ProcessedGenomes/plateMapName] folder you've been working from
3. Once complete, in the top bar of the processing screen, click the arrow next to "Export to CSV" and select "Export to TSV"
4. The file will go to the Downloads file of your computer as nextclade.tsv - Copy the file to [DropBox/MED-LauringLab/ProcessedGenomes/plateMapName] and rename it to plateMapName_nextclade.tsv
5. Copy that re-named file to [DropBox/MED-LauringLab/SequenceSampleMetadata/SequenceOutcomes/nextclade]

---

##### GISAID

GISAID can be accessed at this site: https://www.gisaid.org/

Note: You will need to register - https://www.gisaid.org/registration/register/

1. Compile full data set using pipeline code.
2. Run gisaid_upload_file_creation.R; Ensure that file names/folder paths are all accurate; This will create a YYYYMMDD_Lauring_gisaid_upload_metadata_run_#.xlsx file
3. Run prep_fasta_NumberTwo.py; Ensure that file names/folder paths are all accurate; This will create a .all.consensus.final.gisaid.fasta file
4. Copy that gisaid.fasta file to the proper run folder in [DropBox/MED-LauringLab/GISAID_Uploads/]
5. Double check the fasta file and the .xlsx file have the same number of records, no duplicates, proper rows.
6. Re-save the .xslx file as an .xls file for submission to GISAID
7. Go to https://www.gisaid.org/ and log-in
8. Click on "Upload" then "Batch Upload"
9. Upload the .xls and .fasta file in their corresponding sections, and press "Check and Submit" in the lower right-hand corner

Email correspondence with the GISAID curators will determine the remainder of the GISAID submission process.

Once email confirms sequence release, those will show up within the download portion of GISAID after 24 hours. When that time period has passed:

~~1. Go to https://www.gisaid.org/ and log-in~~
~~2. With the EpiCoV tab highlighted, click the "Downloads" tab~~
~~3. Under the header "Download packages", click on the "metadata" box~~
~~4. Accept the "Terms & Conditions" in the pop-up, and click "Download" - a file called metadata_tsv_YYYY_MM_DD.tar.xz will go to your Downloads folder~~
~~5. Unzip the tar.xz file - This can be done by right-clicking on the file and using 7-Zip (on a Windows computer); You'll need to do it twice - once to get out of .xz, and once to get out of .tar~~
~~6. Inside of these, will be a file called metadata.tsv. Copy that file to [DropBox/MED-LauringLab/SequenceSampleMetadata/SequenceOutcomes/gisaid] and rename it as metadata_YYYY-MM-DD.tsv~~
~~7. Move the previous metadata_YYYY-MM-DD.tsv to the ARCHIVE folder (These files are large, so only keeping the 2-3 most recent in the ARCHIVE folder is sufficient)~~

1. Go to https://www.gisaid.org/ and log-in
2. With the EpiCoV tab highlighted, click the "Search" tab
3. Under "Virus Name", search for "IVY", "RVTN", "MI-UM", select all the sequences, and click "Download" (Separately for each search term)

Note: For "MI-UM", you will need to include a date range, as there is a max download limit of 10,000 sequences.

4. Select "Patient Status Metadata" in the pop-up window, and click "Download".
5. Move the file that is downloaded to [Dropbox (University of Michigan)/MED-LauringLab/SEQUENCING/SARSCOV2/4_SequenceSampleMetadata/SequenceOutcomes/gisaid/new_gisaid_try] and add "_um", "_ivy", or "_rvtn" to the end of the file name.
