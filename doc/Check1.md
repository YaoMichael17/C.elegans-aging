## How we have addressed prior feedback 

We have fixed and placed the title at the very top of the README.md, replacing "C.elegans-aging." We updated the source of the C. elegans data to "ad_worm_aging.h5ad" from the Calico Research website downloads section, including a description of the sequencing data. We removed Zellkonverter, as we were able to complete the extraction and sorting of the data using Scanpy. In our goals, we highlighted specific types of stress genes and our predictions for their expression across the aging timepoints of C. elegans.


## New progression

In regards to the python script, we parsed the pertinent information from the h5ad file. Leveraging the scanpy and pandas library in particular, we were able to isolate the relevant stress hallmark genes of relevance along with their expression data for the some odd 47.5k cells in the dataset. Along with expression data, we correlated the timepoint upon which the cell was taken from (i.e. day1, day8, etc.).

Using Rstudio, we load in the CSV file that contains the gene expression of our selected stress genes across the C.elegans timepoint. Where we then converted the data into a long format for plotting with ggplot. To generate the mean expression of the genes at each timepoint, we used group_by of the timepoint and genes and then used summarize. To create a line plot to visualize the increase and decrease of the stress genes we used geom_line and geom_point. In the end, we generated a line graph that showcases the changes in expression of stress genes across the lifespan of C.elegans.

<img width="1562" height="817" alt="Stress genes expression" src="https://github.com/user-attachments/assets/686f158d-6aaa-4476-ae55-9e13df263816" />

## Project Organization

Michael Yao selected the C. elegans dataset, organized the GitHub repository, and created the figures in RStudio. Anirudh Seshadri selected the mouse dataset, chose the stress genes, and extracted and sorted the data using Python code.


## Struggles we are encountering and questions we would like advice on 

Some struggles we ran into were being able to view the .h5ad file and having python write all the data from all 47.5k cells to a .csv file.
The h5ad is a large 4.13 GB file that is not easily viewable. Trying to read the file in UNIX provided little valuable information. Luckily we got around this by being able to view parts o the dataset such as the main variables and column headers through scanpy in Python. 
Python would also not write the full data table of expression data to the .csv file, and the work around took a very long time to figure out. It required us to use the pandas library to ask python to effectively provide the maximum columns and rows
Within RStudio, we had to convert the data into long format, as we previously ran into errors due to misaligned data. Since each gene has multiple expression measurements at each timepoint, we had to summarize the data points into a mean expression. For the next checkpoint, we would like to perform a comparative analysis of stress hallmark genes between the C. elegans and mouse nervous systems. We are currently struggling with how to properly compare the two datasets and highlight significant differences or similarities between them.

