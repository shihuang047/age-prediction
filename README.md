Human skin, oral, and gut microbiomes predict chronological age
-----------------------
This study performed Random Forest regression analyses of human microbiota from multiple body sites (gut, mouth and skin).
This repository included source codes for generation of all results in the manuscript.


## The data resource and cohort characteristics
### Qiita study IDs involved in the meta-analysis: 
* Gut microbiota:
[10317](https://qiita.ucsd.edu/study/description/10317),
[11757](https://qiita.ucsd.edu/study/description/11757)
* Oral microbiota:
550, 1841, 1774, 2010, 2024, 2136, 10317, 11052, 10052
* Skin microbiota:
1841, 2010, 10317, 11052
### The age distribution of all samples in gut, oral and skin datasets: 


## R scripts
There are some other R scripts and files in this repository that were used in the process of preparing the manuscript, also. Here I'll try to explain some of these.

#### Usage requirements and dependencies
This meta-analysis depends on the self-developed R package [`crossRanger`](https://github.com/shihuang047/crossRanger) that can be downloaded as following.
``` r 
## install.packages('devtools') # if devtools not installed
devtools::install_github('shihuang047/crossRanger')
```
#### What analysis has been done by R script `Age.crossRF_reg.ranger.R`?
The R script `Age.crossRF_reg.ranger.R` performs the meta-analysis of microbiota data from gut, mouth and skin. For each dataset, the main analyzing processes here include 
* Data trimming (such as sample filtering by NA values in the metadata). 
* RF modeling and performance evaluation for the whole dataset. 
* RF modeling and performance evaluation for the sub-datasets. To test if confounders (such as sex) affected the modeling, we first trained the age model within a sub-dataset stratified by a confounder, then applied it on all the other sub-datasets. For both model training and testing, we evaluated regression performance using mean absolute error (MAE). 
* Cross-application of RF models built on the sub-datasets and evaluated the performance using MAE

#### What inputs are neccessary for this R script? 

| Input | gut_data |oral_data | skin_data | Description |
| ------------------ | ------------------ |------------------ |------------------ |-------------------------- |
| datafile  |  gut_data/gut_4575_rare.biom | oral_data/oral_4014.biom | skin_data/skin_4168.biom | Biom-table path |
| sample_metadata  | gut_data/gut_4575_rare_map.txt | oral_data/oral_4014_map.txt | skin_data/skin_4168_map.txt | Metadata path |
| feature_metadata |  gut_data/gut_taxonomy.txt | oral_data/oral_taxonomy.txt | skin_data/skin_taxonomy.txt | Feature metadata path |
| prefix_name  | gut_4575 | oral_4014 | skin_4168 | The prefix of datasets |
| s_category  |  `c("cohort", "sex")` | `"qiita_host_sex"` | `c("body_site","qiita_host_sex")` | The metadata category for dividing datasets |
| c_category  |  `"age"` | `"qiita_host_age"` | `"qiita_host_age"` | The targeted metadata category for RF modeling |


## About the `Input/` folder
This folder includes all the input files (biom table, sample metadata and feature metadata files) necessary for the RF regression analysis. 

## About the `Output/` folder
This folder contains all of the output from the main R script `Age.crossRF_reg.ranger.R`.

## About the `Figures/` folder
This folder contains selected output figures from the `Output` folder for our manuscript.

## Acknowledgements

 This work is supported by IBM Research AI through the AI Horizons Network. For
 more information visit the [IBM AI Horizons Network website](https://www.research.ibm.com/artificial-intelligence/horizons-network/).
