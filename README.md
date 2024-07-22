# mscell
# Cell Lineage and Trajectory Analysis of Single-cell RNA-seq Data of Multiple Sclerosis Patients
# Content
# Project Description 
The aim of this project is to employ single-cell RNA sequencing (scRNA-Seq) to investigate the lineage relationships and developmental trajectories of rare immune cell populations in multiple sclerosis (MS). By utilizing advanced bioinformatics tools for trajectory analysis, we seek to understand how these rare populations evolve in response to inflammatory signals and their potential roles in disease progression. The ultimate goal is to enhance our understanding of cellular diversity and functional states within the MS microenvironment, and to identify novel therapeutic targets that could inform future treatment strategies.
# Introduction 
Multiple sclerosis (MS) is a debilitating neurodegenerative disease affecting millions worldwide, primarily young adults and women. It involves an immune-mediated attack on the central nervous system, leading to the destruction of the myelin sheath that protects nerve fibers. While significant progress has been made in understanding MS, the roles of rare immune cell populations remain largely unexplored. This project aims to leverage single-cell RNA sequencing (scRNA-Seq) to uncover the lineage relationships and developmental trajectories of these rare cell populations in MS. By doing so, we hope to provide deeper insights into the disease mechanisms and identify potential therapeutic targets.
# Workflow
The scRNA-Seq data analysis workflow utilizes several key R tools. 
1. Data Download and Preparation is supported by GEOquery.
2. Data Preprocessing: This step involves FastQC for quality control and Seurat for normalization. In Quality Control, Seurat is again used to filter low-quality cells.
3. Dimensionality Reduction and Clustering, Seurat facilitates PCA, UMAP/t-SNE, and clustering.
4. Cell Type Annotation employs SingleR and CellMarker.
5. Differential Expression Analysis is conducted with DESeq2 and edgeR.
6. Trajectory and Lineage Analysis utilizes Monocle 3.
7. Functional and Pathway Analysis is performed using ClusterProfiler.
8. Integration and Visualization are achieved through Seurat and ComplexHeatmap for effective data representation.
