## Test environments
* local R installation, R 4.1.0
* ubuntu 16.04 (on travis-ci), R 4.1.0
* win-builder (devel)

## R CMD check results
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
