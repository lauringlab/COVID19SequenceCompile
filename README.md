# README for Compilation Pipeline:

### Set-Up:

<b>Folder Structure:</b>

* SampleMetaDataOrganization
    * FinalSummary: Contains the final compiled data full_compiled_data.csv
        * ReportNotifications: Contains out_of_range_alert.csv
        * secret: Contains a duplicate final file (full_compiled_data.csv) for comparison in case of changes made manually to the compiled file
    * Manifests
        * CBR: Contains manifest files from the University of Michigan Campus BioRepository
        * CDCIVY: Contains manifest files from the CDC IVY project
            * Full_IVY_Set: Contains IVY_sample_full_manifest_list.csv, the compiled CDC IVY manifest lists with all variables sent (a subset of these are kept for the final compiled sequence list)
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

<b>Additional Notes:</b>

All code sets have a variable called ```starting_path``` that will need to be changed to your individual path to the level just above the  SampleMetaDataOrganization folder.

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

### Outside the Pipeline

#### Checking Compiled Files

The checking_compiled_files.R code file can be used to see if the "main" version of full_compiled_data.csv matches the "secret" version. The two files should always be the same, but it is possible that the "main" version could differ if individuals manually change the data.
