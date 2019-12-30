Human skin, oral, and gut microbiomes predict chronological age
-----------------------
This study performed Random Forest regression analyses of human microbiota from multiple body sites (gut, mouth and skin).

## What is this repository?

This repository included source codes for generation of all results in this meta-analysis study.
Qiita study IDs involved in the meta-analysis: 
* Gut microbiota:
[10317](https://qiita.ucsd.edu/study/description/10317),
[11757](https://qiita.ucsd.edu/study/description/11757)
* Oral microbiota:
550, 1841, 1774, 2010, 2024, 2136, 10317, 11052, 10052
* Skin microbiota:
1841, 2010, 10317, 11052

There are some other R scripts and files in this repository that were used in
the process of preparing the manuscript, also. Here I'll try to explain some of
these.

## R scripts
The meta-analysis will dependent on the self-developed R package [`crossRanger`] (https://github.com/shihuang047/crossRanger).
``` r 
## install.packages('devtools') # if devtools not installed
devtools::install_github('shihuang047/crossRanger')
```

## About the `input/` folder


## About the `output/` folder


## About the `figures/` folder
