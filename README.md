# Cellular Proteomic Analysis

## Objective
The primary objective is to understand the typical cellular homeostasis of healthy and deleterious cells and observe how cellular proteomic homeostasis transforms over time in response to experimental conditions. Ultimately, the project aims to discover ways to direct deleterious cellular states towards non-deleterious states using data-driven identification and control of high-dimensional dynamical systems. The overarching goal is to contribute to the fight against cancer and improve the quality of life for patients and their loved ones.

## Dataset
Our dataset come from the article <[Link](https://www.biorxiv.org/content/10.1101/2021.12.06.471514v1.full)> finding that the “AP-1 transcription factor network” (i.e., the relative distributions and dependency relationships of transcription factors) are predictive of “cellular plasticity in melanoma” (i.e., how easily changeable the phenotype are melanoma cell lines) 

Our project dataset is available for download here: <[Link](https://drive.google.com/uc?id=1m-bc56NfKErzkxdlHXBLWQg14W2R2vd8&export=download)>

This dataset consists of measurements of 22 AP-1 transcription factors and 4 phenotype proteins simultaneously. These measurements are taken on individual cells under different experimental conditions, allowing the observation of protein levels and their relationships over time. To measure the progression over time, cells are split into groups and measurements are made at different time points, providing insights into the temporal dynamics of the experimental conditions.

## Statistical Methods
The project involves several statistical analyses and data science techniques, including hypothesis testing, correlation estimation, regression, and classification. These analyses aim to answer questions about changes in protein levels over time, differences in protein levels between experimental conditions, relationships between proteins at specific time points, and the predictability of cellular phenotypes based on transcription factors. Additionally, the project explores the patterns and meta-analyses of the obtained results.

## Tools
To support our analysis, we will utilize R programming language and various packages such as Ggplot2, Knitr, Dplyr, Tidyverse, and Rpart. These tools will facilitate data visualization, data manipulation, and the generation of informative reports. 

## References
- NCI Dictionaries, National Cancer Institute
  - https://www.cancer.gov/publications/dictionaries
- The Neural Crest and Cancer: A Developmental Spin on Melanoma, National Library of Medicine
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3809092
- What Is Melanoma Skin Cancer?
  - https://www.cancer.org/cancer/melanoma-skin-cancer/about/what-is-melanoma.html
- AP-1 transcription factor network explains diverse patterns of cellular plasticity in melanoma
  - https://www.biorxiv.org/content/10.1101/2021.12.06.471514v1.full
- MITF gene
  - https://medlineplus.gov/genetics/gene/mitf
- SOX10 gene
  - https://medlineplus.gov/genetics/gene/sox10/
- Nerve Growth Factor Receptor
  - https://www.sciencedirect.com/topics/medicine-and-dentistry/nerve-growth-factor-receptor
- AXL receptor tyrosine kinase as a promising anti-cancer approach
  - https://molecular-cancer.biomedcentral.com/articles/10.1186/s12943-019-1090-3

## Contributors
- Matthew Qiankun Yu [![Foo](https://img.shields.io/badge/github-%23121011.svg?style=for-the-badge&logo=github&logoColor=white)](https://github.com/MatthewQiankunYu)
