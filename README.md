# A transcriptional and regulatory map of mouse somitogenesis

This repository contains the code related to the publication **A transcriptional and regulatory map of mouse somitogenesis** (Ibarra-Soria, Thierion et al., 2023).


Data
----

The **raw and processed data** from this study are available in the ArrayExpress repository and can be accessed through the BioStudies database (https://www.ebi.ac.uk/biostudies/):

- RNA-seq dataset: [E-MTAB-12511](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12511?accession=E-MTAB-12511).
- ATAC-seq dataset: [E-MTAB-12539](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12539?accession=E-MTAB-12539).

These include the raw FASTQ files as well as gene/peak count tables provided as *processed data*.

Code
----
All code used to pre-process and analyse the data is provided, in the `scripts` directories within each of the data modality folders.

- `RNA-seq` and `ATAC-seq` directories contain the primary analysis of each data type.
- `RNA+ATAC` directory contains analyses that integrate both information sources. 


Results
-------

`RNA-seq/results` contains all **differential expression** results (between somite trios and across developmental stages).

`ATAC-seq/results` contains all **differential accessibility** results (between somite trios and across developmental stages).

