

## Title: Investigating Changes in Stress Hallmarks During *C. elegans* Aging

As stress is one of the major contributing factors to aging and its related diseases, the goal of our project is to investigate changes in stress hallmarks in the nervous system of *C. elegans*. We will use single-cell RNA sequencing data from the 2023 paper by Antoine Emile Roux and Han Yuan titled "Individual cell types in *C. elegans* age differently and activate distinct cell-protective responses."  Additionally, studying stress in the *C. elegans* nervous system can serve as a bridge to understanding the mammalian nervous system and enable comparative studies between worms and mice.

## Example published figure:[Figure 5 Global aging characterization reveals difference in magnitude and transcriptional noise during aging across cell types]<img width="3404" height="2470" alt="image" src="https://github.com/user-attachments/assets/26d21109-be23-417a-b638-02fc04452740" />



## Datasets:

(Raw single-cell RNA sequencing data (Day 1-15)):https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208154

(Spatial Atlas of Molecular Cell Types and AAV Accessibility across the Whole Mouse Brain): [https://singlecell.broadinstitute.org/single_cell/study/SCP1830/spatial-atlas-of-molecular-[…]pes-and-aav-accessibility-across-the-whole-mouse-brain](https://singlecell.broadinstitute.org/single_cell/study/SCP1830/spatial-atlas-of-molecular-cell-types-and-aav-accessibility-across-the-whole-mouse-brain#study-download)

## Software: 
scverse(https://scverse.org), 

scanpy(https://scanpy.readthedocs.io/en/latest/), 

bioconductor(https://bioconductor.org/books/3.21/OSCA.intro/getting-scrna-seq-datasets.html#from-hdf5-based-formats)

## Proposed goals:

First goal: Visualize stress hallmark genes in the aging nervous system of C. elegans. To do this, we will use Scanpy (part of the scverse) to perform the initial loading and preprocessing of the dataset from GEO. Once we have the processed dataset, we will identify 20+ genes that are either upregulated or downregulated with aging by comparing day 1 and day 15 samples using Scanpy. Finally, we will use the Zellkonverter package to visualize the expression of these stress-related genes across the two time points.

First "strech" goal: Using the processed data from Scanpy, we will calculate the Maximum Mean Discrepancy (MMD) between young and old nervous system cells using NumPy, and visualize the results using Pandas. 

Second "strech" goal: We will perform a comparative analysis of stress hallmark genes between the C. elegans and mouse nervous systems using the two datasets we have identified.









