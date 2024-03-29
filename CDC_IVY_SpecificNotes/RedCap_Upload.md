# README for CDC IVY RedCap Dataset:

### Data Dictionary:

| Original from Manifest | Column | Data Type | Definition |
| --- | --- | --- | --- |
| aliquot.id | sample_id | character | CUID alpha numeric – printed in an alpha numeric cycling sequence that should be unique to the sample; One individual (subject_id) can have multiple samples, each with their own sample_id number |
| study.id | subject_id | character | Unique identification string for an individual; Identifies the year collected, the site, participant # and the specimen type (nomenclature decided by the CDC/IVY Study Leads) |
| collection.date | coll_date | date | Date the sample was collected from the individual |
| | flag | character | Notes for each sample, some warnings introduced by the pipeline |
| | received_source | character | Place the sample was received from; all should be CDCIVY |
| site.name | sitename | character | Research site where the sample was collected from; coded with site name and state |
| | samplebarcode | character |  |
| | platedate | date | Date the plate was created on |
| | plateplatform | character | Testing system the plate was run through; Options include = Nanopore, Illumina, etc. |
| | platenumber | character | ##_## format, indicates what distinct plate the sample was run with |
| | pangolin_lineage | character | Strain type of the sample; Ex. B.1.1.7 |
| | pangolin_probability | numeric | Probability of correct lineage assignment (pangolin_lineage) by PANGOLIN. Archived as of 3 May 2021 update; Blank in all IVY records |
| | pangolin_status | character | Pangolin classifier assigns either “passed_qc” or “fail” to the sequenced SARS-CoV-2 strain based on sequence completeness and quality. If the consensus SARS-CoV-2 sequence has 50 percent or more ambiguous bases (50%N's) in the consensus fasta or if the consensus length is less than or equal to 10,000bp only, the status will change to fail. |
| | pangolin_note | character | Additional information provided by the pangolin system |
| | nextclade_clade | character | The global clade identified for the sequenced sample |
| | nextclade_totalmissing | numeric | The number of bases missing from the nextclade analysis |
| | nextclade_completeness | numeric | Percentage; how complete the nextclade coverage was of the sample genome |
| | gisaid_strain | character | Virus name, as listed on GISAID. Format is USA/MI-UM-sample_id/YYYY. |
| | gisaid_epi_isl | character | GISAID database accession number. General format is EPI_ISL_NNNNNN. |
| | received_date | date | Date the sample was received at the lab from the received_source |
| position.# | position | character | Location of the sample in the box received at the University of Michigan |
| | platename | character | Full plate name of the sample test run; generally corresponds to plate file name |
| | plateposition | character | Position the sample was in, in the plate that it was tested on |
| | samplesourcelocation | character | Where the sample came from (should correlate with received source) |
| | pangolearn_version | date | The version date of the multinomial logistic regression model that was used for global lineage assignment. |
| | pangolin_conflict | numeric | New as of pangolin software update (notified on 3 May 2021);  the number of conflicts in the algorithm's decision tree |
| | pango_version | character | New as of pangolin software update (notified on 3 May 2021); additional software version information |
| | pangolin_version | character | New as of pangolin software update (notified on 3 May 2021); additional software version information |
| | pangolin_rundate | date | Date of pangolin output receipt; pulled from the pangolin filename the sample_id was reported on |
| | platetopangolin_days | numeric | Number of days between pangolin_runDate and PlateDate; Put in place to add more context for samples when joining multiple sources, as some sample aliquots are received more than once from the same source |
| | nextclade_qcoverallscore | numeric | Quality control overall score (smaller is better) |
| | nextclade_qcoverallstatus | character | Quality control overall status |
| | nextclade_totalmutations | numeric | Total number of nucleotide substitutions relative to Wuhan/Hu-1 reference. Does not include ambiguous nucleotides (see below). as of 6/18/2021 new nextclade records will combine totalmuatations and totalsubstitutions columns to accomodate nextclade changes |
| | nextclade_totalnonacgtns | numeric | Number of ambiguous nucleotide characters, i.e. not A, C, T, G, or N. |
| | nextclade_runDate | date | Date of nextclade output receipt; pulled from the nextclade filename the sample_id was reported on |
| | platetonextclade_days | numeric | Number of days between nextclade_runDate and PlateDate; Put in place to add more context for samples when joining multiple sources, as some sample aliquots are received more than once from the same source |

---

<b>Additional Information</b>
* https://cov-lineages.org/index.html
* https://perkinelmer-appliedgenomics.com/wp-content/uploads/marketing/Coronavirus/NEXTFLEX_Variant-Seq_SARS-CoV-2_Software-from-CosmosID.pdf
* https://github.com/nextstrain/nextclade

---

Completed data rows are updated as information becomes available - If you need a new copy of the data as it was on a particular date, please reach out to Kim Hart (kim.hart@vumc.org) and Julie (Jules) Gilbert (juliegil@umich.edu). Samples are received from Vanderbilt University Medical Center. Please direct questions regarding the dataset to Julie (Jules) Gilbert (juliegil@umich.edu).

---

Specific Notes:

1.  The results from 21020015 and 2102016 should be merged/assigned to the IVY 3 counterparts: 2102015 -> merged with 2102007 and 2102016 -> merged with 2102008
2. GISAID results are uploaded upon availability
3. sample_id ZZX9KW9M resides in RedCap as blank row (withdrawn participant sample)
