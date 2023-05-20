## Test environments
* local R installation, R 4.1.0
* ubuntu 16.04 (on travis-ci), R 4.1.0
* CentOS7 Institutional HPC Cluster


## R CMD check results
## Fri May 19 21:19:39 EDT 2023
── R CMD check results ─────────────────────────────────────────────── gac 0.0.9032 ────
Duration: 21m 30.7s

❯ checking installed package size ... NOTE
    installed size is  5.9Mb
    sub-directories of 1Mb or more:
      data   3.5Mb
      doc    2.0Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖

## Fri Mar 10 16:36:02 EST 2023
── R CMD check results ────────────────────────────────────────────────────── gac 0.0.9031 ────
Duration: 9m 53.2s

❯ checking installed package size ... NOTE
    installed size is  5.9Mb
    sub-directories of 1Mb or more:
      data   3.5Mb
      doc    2.0Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖

## Thu Mar  2 13:28:46 EST 2023
── R CMD check results ───────────────────────────────────────────────────── gac 0.0.90030 ────
Duration: 9m 42.8s

❯ checking installed package size ... NOTE
    installed size is  5.9Mb
    sub-directories of 1Mb or more:
      data   3.5Mb
      doc    2.0Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖


### Tue Oct 11 12:13:35 EDT 2022
── R CMD check results ───────────────────────────────────────────────────── gac 0.0.90029 ────
Duration: 9m 44.9s

❯ checking installed package size ... NOTE
    installed size is  6.0Mb
    sub-directories of 1Mb or more:
      data   3.5Mb
      doc    2.1Mb


### Thu Sep 22 12:09:26 EDT 2022
── R CMD check results ───────────────────────────────────────────────────── gac 0.0.90029 ────
Duration: 31m 45.2s

❯ checking installed package size ... NOTE
    installed size is  6.2Mb
    sub-directories of 1Mb or more:
      data   3.5Mb
      doc    2.2Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖

### Fri Sep  9 08:23:10 EDT 2022
── R CMD check results ───────────────────────────────────────────────────── gac 0.0.90028 ────
Duration: 15m 5.3s

❯ checking installed package size ... NOTE
    installed size is  6.0Mb
    sub-directories of 1Mb or more:
      data   3.5Mb
      doc    2.1Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖

── R CMD check results ────────────────────────────────────────────────────── gac 0.0.90028 ────
Duration: 15m 52.7s

❯ checking installed package size ... NOTE
    installed size is  6.0Mb
    sub-directories of 1Mb or more:
      data   3.5Mb
      doc    2.1Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖

### Thu Sep  8 11:32:33 EDT 2022
── R CMD check results ────────────────────────────────────────────────────── gac 0.0.90027 ────
Duration: 12m 30.7s

❯ checking PDF version of manual ... WARNING
  LaTeX errors when creating PDF version.
  This typically indicates Rd problems.
  LaTeX errors found:
  ! LaTeX Error: File `inconsolata.sty' not found.
  
  Type X to quit or <RETURN> to proceed,
  or enter new name. (Default extension: sty)
  
  ! Emergency stop.
  <read *> 
           
  l.297 ^^M
           
  !  ==> Fatal error occurred, no output PDF file produced!

❯ checking installed package size ... NOTE
    installed size is  6.0Mb
    sub-directories of 1Mb or more:
      data   3.5Mb
      doc    2.1Mb

❯ checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ‘gac-manual.tex’

0 errors ✔ | 1 warning ✖ | 2 notes ✖


### Wed Sep  7 16:51:29 EDT 2022
── R CMD check results ────────────────────────────────────────────────────── gac 0.0.90027 ────
Duration: 9m 44.6s

❯ checking installed package size ... NOTE
    installed size is  6.0Mb
    sub-directories of 1Mb or more:
      data   3.5Mb
      doc    2.1Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖

### Wed Aug 31 23:02:21 EDT 2022
── R CMD check results ─────────────────────────────────────────────────────── gac 0.0.90026 ────
Duration: 1h 17m 19.9s

❯ checking whether package ‘gac’ can be installed ... WARNING
  See below...
     Warning: package ‘vegan’ was built under R version 4.1.1
	 ...

0 errors ✔ | 1 warning ✖ | 0 notes ✔

### Mon May 24 23:26:20 EDT 2021
── R CMD check results ────────────────────────────────────── gac 0.0.90022 ────
Duration: 8m 54.7s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔


### Mon Apr 19 16:18:56 2021
── R CMD check results ────────────────────────────────────── gac 0.0.90021 ────
Duration: 9m 1.4s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

### Sun Apr 11 00:16:00 2021 
── R CMD check results ────────────────────────────────────── gac 0.0.90020 ────
Duration: 9m 10.8s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

### Sun Mar  7 12:40:35 2021
── R CMD check results ────────────────────────────────────── gac 0.0.90019 ────
Duration: 9m 4.9s

❯ checking R code for possible problems ... NOTE
  doKSpectral: no visible binding for global variable ‘topEigenValues’
  doKSpectral: no visible binding for global variable ‘k’
  doKSpectral: no visible binding for global variable ‘dLambdaMax’
  Undefined global functions or variables:
    dLambdaMax k topEigenValues

0 errors ✔ | 0 warnings ✔ | 1 note ✖


### Sun Jan 31 10:42:04 2021 
── R CMD check results ────────────────────────────────────── gac 0.0.90018 ────
Duration: 6m 31.2s

❯ checking R code for possible problems ... NOTE
  doKSpectral: no visible binding for global variable ‘topEigenValues’
  doKSpectral: no visible binding for global variable ‘k’
  doKSpectral: no visible binding for global variable ‘dLambdaMax’
  Undefined global functions or variables:
    dLambdaMax k topEigenValues

### Sun Jan 17 12:04:21 2021
── R CMD check results ────────────────────────────────────── gac 0.0.90017 ────
Duration: 6m 40.1s

❯ checking installed package size ... NOTE
    installed size is  7.7Mb
    sub-directories of 1Mb or more:
      data      1.7Mb
      extdata   4.8Mb

❯ checking dependencies in R code ... NOTE
  Package in Depends field not imported from: ‘circlize’
    These packages need to be imported from (in the NAMESPACE file)
    for when this namespace is loaded but not attached.

❯ checking R code for possible problems ... NOTE
  HeatmapCNR: no visible global function definition for ‘gpar’
  doKSpectral: no visible binding for global variable ‘topEigenValues’
  doKSpectral: no visible binding for global variable ‘k’
  doKSpectral: no visible binding for global variable ‘dLambdaMax’
  Undefined global functions or variables:
    dLambdaMax gpar k topEigenValues

0 errors ✔ | 0 warnings ✔ | 3 notes ✖

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
