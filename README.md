# PDAC-Profiler

This repository contains a framework for predicting subtypes of Pancreatic Ductal Adenocarcinoma (PDAC) using transcriptome data (RNA-seq or microarray) and pre-defined PDAC subtype signatures (see [References](#References)).

## Citation

If you use **PDAC-Profiler** in your work, please cite the original publication:
Hafner et al., Manuscript in prep.

## Getting Started

The provided scripts constitute the official bioinformatic pipeline for processing normalized transcriptome data to predict PDAC subtypes.

**Composing of Results:**
TO DO



### Prerequisites
Requiered software and databases

1. Software and dependencies
	* R (4.4.0 or later)

2. R packages
	* openxlsx
  * tidyverse
    
3. Additional Files (provided in data folder)
	* PDAC_signatures.tsv

## Directory structure

```
PDAC-Profiler
│
├── data
|   ├── toy_data.tsv
│   └── PDAC_signatures.tsv
|
└── script   
    ├── annotateGenes.R
    ├── bed2sequence.R
    └── bedTools_fct.R


```
**Remark:** Before using the **PDAC-Profiler**, ensure that your transcriptome data are properly **normalized**. We recommend **state-of-the-art normalization methods** such as:
- **Library size normalization** and **TMM** (Trimmed Mean of M-values) for RNA-seq data. <br/>
- **RMA** (Robust Multi-array Average) for microarray data.
  
**Batch correction** may be required depending on the dataset. Users are responsible for determining whether batch effects are present and, if so, for applying appropriate correction methods.

**Important:** Use **HUGO Gene Symbols** (Gene Symbol) as unique identifiers for your genes.


## Running PDAC-Profiler
Make sure all required R packages are installed. The **PDAC-Profiler** can then be executed using:

```
Rscript ./PDAC-Profiler.R --input_file "./data/toy_data.tsv"\
                          --signature_file "./data/PDAC_signatures.tsv"\
                          --output_dir "./output"
```

**--input_file** input file with normalized mRNA intensity (e.g. log2 CPM or log2 RMA), expected genes in rows, samples in columns<br/>
**--signature_file** <br/>
**--output_dir** <br/>

### Parameters



## Authors

* Geoffroy Andrieux
* Tonmoy Das


## License
This software is under AGPL3 license.

## Acknowledgments
We thank all members of our laboratories for constructive discussions and suggestions.


## References




