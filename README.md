# README for Compilation Pipeline:

### Set-Up:

<b>Folder Structure:</b>

Folder structure lives within the MED-LauringLab DropBox folder.

* SequenceSampleMetadata
    * FinalSummary: Contains the final compiled data full_compiled_data.csv
        * ReportNotifications: Contains out_of_range_alert.csv
        * secret: Contains a duplicate final file (full_compiled_data.csv) for comparison in case of changes made manually to the compiled file
        * CDC_IVY_UPLOADS: Contains the files created in order to upload to the CDC IVY RedCap database
    * Manifests
        * CBR: Contains manifest files from the University of Michigan Campus BioRepository
        * CDCIVY: Contains manifest files from the CDC IVY project
            * Full_IVY_Set: Contains IVY_sample_full_manifest_list.csv, the compiled CDC IVY manifest lists with all variables sent (a subset of these are kept for the final compiled sequence list)
            * Keys: Contains CDC_SiteCodebook.csv, the list of CDC IVY sites, their locations, and their coded name
        * CSTP: Contains manifest files from the COVID-19 Sampling & Tracking Program (ultimately from LynxDx)
        * EDIDNOW: Contains manifest files from the ED IDnow testing program at Michigan Medicine
        * Martin: Contains manifest files from the Martin Lab at the School of Public Health
        * ManifestsComplete: Contains the manifest_output_report_YYYYMMDD.xlsx reports, as well as the compiled manifest list sample_full_manifest_list.csv
    * PlateMaps: Contains all plate map files as YYYYMMDD_Illumina_Plate_##.xlsx or YYYYMMDD_Nanopore_Run_##.xlsx
        * PlateMapsComplete: Contains the compiled plate map list sample_full_plate_list.csv
    * PreviousLists: Contains the original processed sample list ProcessedSampleCumulativeList_20210326.csv maintained by Andrew Valesano prior to code implementation
    * SequenceOutcomes
        * gisaid: Contains the GISAID metadata file metadata_YYYY-MM-DD_MM-DD.tsv
        * nextclade: Contains the NextClade output files as either YYYYMMDD_Plate_##\_##\_nextclade.tsv or YYYYMMDD\_Nanopore\_Run_##_nextclade.tsv
        * pangolin: Contains the Pangolin output files as either YYYYMMDD_Plate_##\_##\_pangolin.csv or YYYYMMDD\_Nanopore\_Run_##_pangolin.csv
        * SequenceOutcomeComplete: Contains the sequence outcome final compiled files of sample_full_gisaid_list.csv, sample_full_nextclade_list.csv, and sample_full_pangolin_list.csv

<b>R Libraries:</b>

* tidyverse: https://www.tidyverse.org/
* lubridate: https://lubridate.tidyverse.org/
* janitor: https://cran.r-project.org/web/packages/janitor/index.html
* withr: https://cran.r-project.org/web/packages/withr/index.html
* openxlsx: https://cran.r-project.org/web/packages/openxlsx/index.html

<b>Code:</b>

In order to run the compilation code pipeline, download the following code sets:

* <span style="background-color: #8AA29E">pipeline_functions.R</span>
* manifest_code.R
* plate_map_code.R
* pangolin_code.R
* nextclade_code.R
* gisaid_code.R
* compile_components_code.R
* <span style="background-color: #F0B7B3">full_run_code.R</span>

And put them all in the same folder on your computer. Use full_run_code.R to run everything in the correct order. Ensure that all file paths and folder names mentioned in each code set are correct.

<b>Additional Notes:</b>

* All code sets have a variable called ```starting_path``` that will need to be changed to your individual path to the level just above the  SampleMetaDataOrganization folder.
* Ensure that you have access to either the LauringLab DropBox folder, or your own version of the pipeline folder structure, on your computer.


### Code Order:

Use full_run_code.R to run all the pieces of the pipeline in order. This order is:

* manifest_code.R
* plate_map_code.R
* pangolin_code.R
* nextclade_code.R
* gisaid_code.R
* compile_components_code.R

### Data Dictionary for Compiled File:

The final created file is called <b>full_compiled_data.csv</b>.

| Column | Source | Data Style | Join Variable | Definition |
| --- | --- | --- | --- | --- |
| sample_id | sample_full_manifest_list.csv | numeric | Yes | Unique identification number for a particular sample; One individual (subject_id) can have multiple samples, each with their own sample_id number |
| subject_id | sample_full_manifest_list.csv | numeric | | Unique identifcation number for an individual; received source determines the type. CSTP = umid, Martin = other, UHS = MRN, CBR = MRN |
| coll_date | sample_full_manifest_list.csv | date | | Date the sample was collected from the individual |
| flag | sample_full_manifest_list.csv | character | | Notes for each sample, some warnings introduced by the pipeline |
| received_source | sample_full_manifest_list.csv | character | | Place the sample was received from; options include = UHS, CBR, CSTP, Martin, etc. |
| SiteName | sample_full_manifest_list.csv | character | | Research site where the sample was collected from; Only applicable for the CDC IVY study |
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
| received_date | sample_full_manifest_list.csv | date | Yes | Date the sample was received at the lab from the received_source |
| position | sample_full_manifest_list.csv | character | | Position the sample was in, in the box received from the received source |
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
| IlluminaPangolin_OutOfRange | {calculated in compile_components_code.R} | numeric | | 1,0 binary; If PlatePlatform is Illumina and PlateToPangolin is more than 8, then the column is marked (1) as potentially being out of range/an incorrect sample to data match; Marked rows are output in SampleMetadataOrganization/FinalSummary/ReportNotifications/out_of_range_alert.csv |
| NanoporePangolin_OutOfRange | {calculated in compile_components_code.R} | numeric | | 1,0 binary; If PlatePlatform is Nanopore and PlateToPangolin is more than 4, then the column is marked (1) as potentially being out of range/an incorrect sample to data match; Marked rows are output in SampleMetadataOrganization/FinalSummary/ReportNotifications/out_of_range_alert.csv |
| IlluminaNextclade_OutOfRange | {calculated in compile_components_code.R} | numeric | | 1,0 binary; If PlatePlatform is Illumina and PlateToNextclade is more than 8, then the column is marked (1) as potentially being out of range/an incorrect sample to data match; Marked rows are output in SampleMetadataOrganization/FinalSummary/ReportNotifications/out_of_range_alert.csv |
| NanoporeNextclade_OutOfRange | {calculated in compile_components_code.R} | numeric | | 1,0 binary; If PlatePlatform is Nanopore and PlateToNextclade is more than 4, then the column is marked (1) as potentially being out of range/an incorrect sample to data match; Marked rows are output in SampleMetadataOrganization/FinalSummary/ReportNotifications/out_of_range_alert.csv |

---

<b>Additional Information</b>
* https://cov-lineages.org/index.html
* https://perkinelmer-appliedgenomics.com/wp-content/uploads/marketing/Coronavirus/NEXTFLEX_Variant-Seq_SARS-CoV-2_Software-from-CosmosID.pdf
* https://github.com/nextstrain/nextclade

---

## Outside the Pipeline

#### Checking Compiled Files

The checking_compiled_files.R code file can be used to see if the "main" version of full_compiled_data.csv matches the "secret" version. The two files should always be the same, but it is possible that the "main" version could differ if individuals manually change the data.

#### Generating .csv File for RedCap Upload

The cdc_ivy_upload_code.R code file is used to generate the new rows of data that need to be manually uploaded to the CDC IVY RedCap database.

### Folder: ProcessingFASTA

#### Converting Names in .FASTA files

The prep_fasta_for_gisaid.py code converts the genome names within .fasta files for use in pangolin/nextclade/gisaid.

#### Subset .FASTA files by ID

The SelectSequences.py code is for getting a subset of a .fasta file by ID.

---

### Process

Manifests are received from the following sources:

* COVID-19 Sampling & Tracking Program (CSTP) - samples from this source are processed first by LynxDx
* Martin Lab at the University of Michigan School of Public Health (Martin)
* University of Michigan Central Biorepository (CBR)
* CDC IVY Project (CDCIVY) - samples from 21 sites sent to Vanderbilt, then to University of Michigan
* Michigan Medicine ED ID Now project (EDIDNOW)

Dr. Adam Lauring reviews these manifests, checks and renames columns as necessary, renames the file, and places them in the appropriate Manifests folder [within DropBox/MED-LauringLab/SequenceSampleMetadata/Manifests].

---

##### Manifest Column Format (for all except CDC IVY)

| Columns | Data Type	| Variable Description |
| --- | --- | --- |
| position | character | Where the sample is located in the box sent |
| sample_id | numeric | Identification number for the sample; unique |
| subject_id | numeric | Identification number for the subject (individual) - there may be multiple samples per subject |
| coll_date | date | Date the sample was collected on; M/D/YY format |
| flag | character | Notes |

##### Manifest File Name Format

nameOfSource_YYYYMMDD_#.csv

nameOfSource = Where the samples came from, corresponds to the name of the project folders inside the Manifests folder (CBR, CSTP, Martin, CDCIVY, EDIDNOW)

YYYYMMDD = Year, Month, and Day of when the samples arrived/the manifest file was received

\# = Number of the manifest; Will usually be a "1", but if two batches of samples arrive on the same day from the same source, with two separate associated manifest files, then these would be numbered accordingly

---

##### Samples

Samples are received, and sequenced on plates using Nanopore or Illumina systems. Plate Map files are generated, matching samples to their location on those plates. Original plate map files are placed in [DropBox/MED-LauringLab/Plate Maps] as .xlsx files. Before running the data compilation pipleline, these files will need to be copied to [DropBox/MED-LauringLab/SequenceSampleMetadata/PlateMaps] and checked for format in column order/naming and file naming.

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

---

There is another python script that takes the generated .fasta file and replaces the barcode string with the matching sample ID string. Once this .fasta file has been created, it will be used in the Pangolin, Next Clade, and GISAID systems [prep_fasta_for_gisaid.py].

---

##### Pangolin

Instructions for installing the command line version of Pangolin can be found here: https://cov-lineages.org/pangolin_docs/installation.html

Note: This can only be completed using Linux or OS operating systems.

Pangolin can also be used at the following site: https://pangolin.cog-uk.io/

##### Next Clade

Next Clade can be accessed at this site: https://clades.nextstrain.org/

After the process is complete, the output should be downloaded as a .tsv file.

##### GISAID

GISAID can be accessed at this site: https://www.gisaid.org/

Note: You will need to register - https://www.gisaid.org/registration/register/
