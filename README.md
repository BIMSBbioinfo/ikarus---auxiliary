# Ikarus 
This repository contains code to generate figures for [Ikarus](https://github.com/BIMSBbioinfo/ikarus) paper. 

## Source code for figures
Repository is divided into two parts (part 1 and part 2) that separate code to Python and R.
### Part 1
Python part contains a collection of jupyter notebooks. Running those with a provision of the respective *Data* will render figures in pdf.
### Part 2
R part contains a collection of R scripts. Each script takes as an argument a location of respective *Data* folder and, when executed, will render the corresponding figure in pdf.
For example:
```R
Rscript F1_C.R ./Data
```
All plots were generated with R versiob 4.0.5.


## Data
The tar balls with data to recapitulate the analysis can be found:
| Link | Description |
| --- | --- |
| [part 1](https://bimsbstatic.mdc-berlin.de/akalin/Ikarus/part_1.tar.gz) | Python part |
| [part 2](https://bimsbstatic.mdc-berlin.de/akalin/Ikarus/part_2.tar.gz) | R part |
