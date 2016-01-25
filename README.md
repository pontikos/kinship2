# kinship2

Fork of  a read-only mirror of the CRAN R package repository.  kinship2 — Pedigree Functions. Homepage: http://r-forge.r-project.org  

This is a fork of the cran kinship2 package extended with additional pedigree functions.

```
R CMD build --no-build-vignettes kinship2
R CMD INSTALL kinship2
```
 
 ```R
 read.pedigree
 ```
 expects
 
 ```
 "ID","Family","Father","Mother","Gender","Affection"
 ```
