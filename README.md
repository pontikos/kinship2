# kinship2

This is a fork of  a read-only mirror of the CRAN R package repository.  kinship2 — Pedigree Functions. Homepage: http://r-forge.r-project.org
extended with additional pedigree functions.

```
R CMD build --no-build-vignettes kinship2
R CMD INSTALL kinship2
```
 
# Additional functions
 
 Are in file ```pedigree-functions.R```
 
 ```R
 read.pedigree
 ```
 expects a CSV file with at least the following headers:
 
 ```
 "ID","Family","Father","Mother","Gender","Affection"
 ```
