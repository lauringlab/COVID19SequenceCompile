# CDC IVY Information

## Data Dictionaries

#### IVY_sample_full_manifest_list.csv
"position.#"      "site.name"       "study.id"        "specimen.type"   "collection.date"
 [6] "aliquot.id"      "rnp.ct"          "covid-19.n1"     "covid-19.n2"     "freezer.box"
[11] "cell.#"
| Column | Data Style | Equivalent To | Definition |
| --- | --- | --- | --- |
| position.# | numeric (integer) | position | Location of the sample in the box received at the University of Michigan |
| site.name | character | | Name of the site that the sample originated from; format = <location>_<state abbreviation>; used in GISAID identifier |
| study.id | character | subject_id | Identifies the year collected, the site, participant # and the specimen type (nomenclature decided by the CDC/IVY Study Leads) |
| specimen.type | character | | Type of specimen collected from the study participant; ex. Nasal/Throat Swab |
| collection.date | date | coll_date | Date the sample was collected on |
| aliquot.id | character | sample_id | CUID alpha numeric – printed in an alpha numeric cycling sequence that should be unique to the sample |
| rnp.ct | numeric | | |
| covid-19.n1 | numeric | | |
| covid-19.n2 | numeric | | |
| freezer.box | numeric (integer) | | Freezer location from previous site |
| cell.# | numeric (integer) | | Cell location from previous site |
