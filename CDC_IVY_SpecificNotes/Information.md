# CDC IVY Information

## Data Dictionaries

#### IVY_sample_full_manifest_list.csv

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

## CDC IVY Sites

[Centers for Disease Control and Prevention in Atlanta, GA]

| Number | Institution | Location | Site Name Code |
| --- | --- | --- | --- |
| 1 | Baylor, Scott & White | Temple, TX | BSWH_TX |
| 2 | Baystate Medical Center | Springfield, MA | Baystate_MA |
| 3 | Beth Israel Deaconess | Boston, MA | |
| 4 | Cleveland Clinic | Cleveland, OH | Cleveland_OH |
| 5 | Emory University | Atlanta, GA | |
| 6 | Hennepin County Medical Center | Minneapolis, MN | Hennepin_MN |
| 7 | Intermountain Medical Center | Murray, UT | |
| 8 | Johns Hopkins | Baltimore, MD | |
| 9 | Montefiore Medical Center | Bronx, NY | |
| 10 | Ohio State Medical Center | Columbus, OH | |
| 11 | Oregon Health & Sciences University | Portland, OR | |
| 12 | Stanford University | Stanford, CA | |
| 13 | University of California, Los Angeles | Los Angeles, CA | |
| 14 | University of Colorado | Aurora, CO | Univ. Colorado |
| 15 | University of Iowa | Iowa City, IA | Univ. of Iowa_IA |
| 16 | University of Miami | Miami, FL | |
| 17 | University of Michigan | Ann Arbor, MI | |
| 18 | University of Washington | Seattle, WA | Univ of Wash_WA |
| 19 | Vanderbilt University Medical Center | Nashville, TN | Vanderbilt_TN |
| 20 | Wake Forest University | Winston-Salem, NC | |
| 21 | Washington University | St. Louis, MO | |

Coordinating Center: Vanderbilt University Medical Center in Nashville, TN
