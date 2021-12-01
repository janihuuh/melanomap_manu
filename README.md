# Evolution and modulation of antigen-specific T cell responses in melanoma patients

Scripts to reproduce figures and analyses in the manuscript "Evolution and modulation of antigen-specific T cell responses in melanoma patients" Huuhtanen et al., submitted 

# Installation

Clone this github repository, e.g., in terminal by

$ git clonehttps://github.com/janihuuh/melanomap_manu

All the other dependencies (R-packages, Python modules) are mentioned in the code and their installation guides can be found in their own respective guides. The typical installation time for all the software on a standard laptop (e.g., Apple M1 16Gb) should not exceed 1 hour. 

# Pseudocode for MelanoMAP

## scTCRab-seq preprocessing



## scRNAseq preprocessing

<pre>
<b>for</b> Huuhtanen et al., Durante et al data, <i><b>do</b></i>  
  <b>read</b> 10X CellRanger output 
  <b>filter</b> based on the quality of cells 

<b>for</b> Jerby-Arnon et al, Li et al., Sade-Feldman et al, <i><b>do</b></i>  
  <b>read</b> count matrix
  <b>filter</b> pre-recieved QC

### No batch correction (Jerby-Arnon et al., Li et al., Sade-Feldman et al)
<pre>
<b>for</b> scTCRseq data, <i><b>do</b></i>  
  <b>read</b> count matrix
  <b>filter</b> ## based on the quality of cells (incomple TCRab information, low confidence data)  
  <b>Normalize and scale</b> 
  <b>Find HVGs</b> ## select top 2000, exclude V(D)J genes
  <b>PCA on HVGs</b> ## retain PCs >2 SD
  <b>Graph-based clustering on PCs</b>
  <b>annotate</b> 
</pre>

### Batch correction (Huuhtanen et al., Durante et al.)
<pre>
<b>for</b> scTCRseq data, <i><b>do</b></i>  
  <b>read</b> count matrix
  <b>filter</b> ## based on the quality of cells (incomple TCRab information, low confidence data)  
  <b>scVI on all cells </b> ## calculate 30 latent embeddings
  <b>Graph-based clustering on latent embeddings</b>
  <b>annotate</b> 
</pre>

	
  
  


## To reproduce the results:

### 1) Clone this repository

```
git clone https://github.com/janihuuh/melanomap_manu
cd path/to/melanomap_manu/
```

### 2) Obtain the data

* The processed scRNA+TCRab-seq data can be received from EGA (accession number EGAS00001005580) and links listed in Supplementary Data 1. 
* The TCRb-seq data can be received from immuneAccess from links listed in Supplementary Data 1

### 3) Create seurat-objects

```
Rscript R/main.R ## init the helper-functions, coloring, etc.
Rscript R/helper/run_preprocessTCRab.R ## preprocess the TCRab-seq
python python/run_scvi_example.py ## obtain the latent embeddings to use in clustering, UMAPs. This is just an example script
Rscript R/rnaseq/run_createSeurat.R ## read the CellRanger output, filter, merge TCRab-data, merge scVI latent embeddings, run singleR, DEs, pathways, etc.

```

### 4) Run additional RNAseq analyses

```
Rscript R/rnaseq/run_cellphone.R ## init data for CellPhoneDb
bash bash/run_cellphonedb.sh ## Run CellPhoneDb
Rscript R/rnaseq/run_scenic.R ## run Scenic-analysis
bash bash/run_vartrix.sh ## Run Vartrix
Rscript R/rnaseq/run_edgeR.R ## run bulk-RNAseq analysis
```

### 5) Run additional TCRseq analyses 

```
bash bash/run_vdjtools.sh ## preprocess the bulk-TCRb-data, subsample, calcualte diverisities, etc.
bash bash/run_gliph2_expample.sh ## run GLIPH2 collectively and on individual samples; this is just an example script
python python/run_tcrgp.py ## run TCRGP collectively on samples
Rscript R/tcrseq/run_antigen_drive.R ## run antigen-drive based analyses
```
