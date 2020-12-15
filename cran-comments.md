## Test environments
* local R installation, R 4.1.0
* ubuntu 16.04 (on travis-ci), R 4.1.0
* win-builder (devel)

## R CMD check results

###  Mon Dec 14 19:08:13 2020
── R CMD check results ────────────────────────────────────── gac 0.0.90014 ────
Duration: 6m 25.4s

❯ checking whether package ‘gac’ can be installed ... WARNING
W  checking whether package ‘gac’ can be installed (16.6s)
   Found the following significant warnings:
     Warning: replacing previous import ‘GenomicRanges::union’ by ‘dplyr::union’ when loading ‘gac’
     Warning: replacing previous import ‘GenomicRanges::intersect’ by ‘dplyr::intersect’ when loading ‘gac’
     Warning: replacing previous import ‘GenomicRanges::setdiff’ by ‘dplyr::setdiff’ when loading ‘gac’
   See ‘/private/var/folders/f7/v_vt_5bn07z9tljnh56n9crnrzdbg9/T/Rtmpuo5hhr/gac.Rcheck/00install.out’ for details.
❯ checking installed package size ... NOTE
    installed size is  7.5Mb
    sub-directories of 1Mb or more:
      data      1.7Mb
      extdata   4.8Mb

❯ checking R code for possible problems ... NOTE
  HeatmapCNR: no visible global function definition for ‘gpar’
  doKSpectral: no visible binding for global variable ‘topEigenValues’
  doKSpectral: no visible binding for global variable ‘k’
  doKSpectral: no visible binding for global variable ‘dLambdaMax’
  Undefined global functions or variables:
    dLambdaMax gpar k topEigenValues

❯ checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ‘cnr.rds’

0 errors ✔ | 1 warning ✖ | 3 notes ✖



### Mon Oct  5 21:02:11 2020
── R CMD check results ─────────────────────────────────────── gac 0.0.9009 ────
Duration: 3m 18.8s

❯ checking installed package size ... NOTE
    installed size is  8.4Mb
    sub-directories of 1Mb or more:
      data      2.7Mb
      extdata   4.8Mb

❯ checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ‘cnr.rds’

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

### Tue Sep 15 18:53:38 2020
 
* checking package dependencies ... ERROR
Namespace dependency not required: ‘SCclust’

See section ‘The DESCRIPTION file’ in the ‘Writing R Extensions’
manual.
* DONE
Status: 1 ERROR


Duration: 3m 15.3s


### 
❯ checking installed package size ... NOTE
    installed size is  8.9Mb
    sub-directories of 1Mb or more:
      data      2.7Mb
      doc       1.2Mb
      extdata   4.8Mb

❯ checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ‘cnr.rds’ ‘dna.rds’

0 errors ✔ | 0 warnings ✔ | 2 notes ✖


* This work is still on alpha-release.
