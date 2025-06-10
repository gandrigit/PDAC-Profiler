# PDAC-Profiler

This repository contains a framework for predicting subtypes of Pancreatic Ductal Adenocarcinoma (PDAC) using transcriptome data (RNA-seq or microarray) and pre-defined PDAC subtype signatures (see [References](#References)).

## Citation

If you use **PDAC-Profiler** in your work, please cite the original publication:
Hafner et al., Manuscript in prep.

[![DOI](https://zenodo.org/badge/999408057.svg)](https://doi.org/10.5281/zenodo.15632386)

## Getting Started

The provided scripts constitute the official bioinformatic pipeline for processing normalized transcriptome data to predict PDAC subtypes.

### Prerequisites
Requiered software and databases

1. Software and dependencies
	* R (4.4.0 or later)

2. R packages
	* openxlsx
  * tidyverse
  * optparse
  * org.Hs.eg.db
  * org.Mm.eg.db
  * homologene
  * clusterProfiler
  * gridExtra
  * parallel
  * circlize
  * ComplexHeatmap
  * limma

3. Additional Files (provided in data folder)
	* PDAC_signatures.tsv

## Directory structure

```
PDAC-Profiler
│
├── data
│  ├── toy_data.tsv
│  └── PDAC_signatures.tsv
│
└── R   
    ├── enrichment_plot_helper_fct.R
    └── PDAC-Profiler.R
```

**Remark:** Before using the **PDAC-Profiler**, ensure that your transcriptome data are properly **normalized**. We recommend **state-of-the-art normalization methods** such as:
- **Library size normalization** and **TMM** (Trimmed Mean of M-values) for RNA-seq data. <br/>
- **RMA** (Robust Multi-array Average) for microarray data.
  
**Batch correction** may be required depending on the dataset. Users are responsible for determining whether batch effects are present and, if so, for applying appropriate correction methods.

**Important:** Use **HUGO Gene Symbols** (Gene Symbol) as unique identifiers for your genes.


## Running PDAC-Profiler
Make sure all required R packages are installed. The **PDAC-Profiler** can then be executed using:

```
Rscript ./R/PDAC-Profiler.R --input_file "./data/toy_data.tsv"\
                          --signature_file "./data/PDAC_signatures.tsv"\
                          --output_dir "./output"
                          --species_name "human"
```

### Parameters

**--input_file** input file with normalized mRNA intensity (e.g. log2 CPM or log2 RMA), expected genes in rows, samples in columns<br/>
**--signature_file** PDAC subtype signatures .tsv file. One column per subtype. Use **HUGO Gene Symbols** (Human Gene Symbol) as unique identifiers. <br/>
**--output_dir** output directory<br/>
**--species_name** name of the current species (human or mouse).

### Output directory structure

```
output_dir
└── fgsea
     └── PDAC_subtype
          ├── sample_A_fgsea.xlsx
          ├── sample_B_fgsea.xlsx
          ├── ...
          ├── heatmaps
          └── PDAC_subtype.xlsx

```


## Authors

* Geoffroy Andrieux
* Tonmoy Das


## License
This software is under AGPL3 license.

## Acknowledgments
We thank all members of our laboratories for constructive discussions and suggestions.


## References
Moffitt RA, Marayati R, Flate EL, et al. Virtual microdissection identifies distinct tumor- and stroma-specific subtypes of pancreatic ductal adenocarcinoma. Nat Genet. 2015;47(10):1168-1178. doi:10.1038/ng.3398 [Moffitt et al.](https://pubmed.ncbi.nlm.nih.gov/26343385/)

Bailey P, Chang DK, Nones K, et al. Genomic analyses identify molecular subtypes of pancreatic cancer. Nature. 2016;531(7592):47-52. doi:10.1038/nature16965 [Bailey et al.](https://pubmed.ncbi.nlm.nih.gov/26909576/)

Collisson EA, Sadanandam A, Olson P, et al. Subtypes of pancreatic ductal adenocarcinoma and their differing responses to therapy. Nat Med. 2011;17(4):500-503. doi:10.1038/nm.2344 [Collisson et al.](https://pubmed.ncbi.nlm.nih.gov/21460848/)

Puleo F, Nicolle R, Blum Y, et al. Stratification of Pancreatic Ductal Adenocarcinomas Based on Tumor and Microenvironment Features. Gastroenterology. 2018;155(6):1999-2013.e3. doi:10.1053/j.gastro.2018.08.033 [Puleo et al.](https://pubmed.ncbi.nlm.nih.gov/30165049/)

Chan-Seng-Yue M, Kim JC, Wilson GW, et al. Transcription phenotypes of pancreatic cancer are driven by genomic events during tumor evolution [published correction appears in Nat Genet. 2020 Apr;52(4):463. doi: 10.1038/s41588-020-0588-3.]. Nat Genet. 2020;52(2):231-240. doi:10.1038/s41588-019-0566-9 [Chan-Seng-Yue et al.](https://pubmed.ncbi.nlm.nih.gov/31932696/)


