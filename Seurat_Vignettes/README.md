Ben Mescher

Dataflow and Availability
-------------------------

#### 10x Genomics Datasets, Cell ranger, [Mouse brain Serial Section 1 (sagittal-Anterior)](https://www.10xgenomics.com/resources/datasets/)

-   FASTQ input, 29GB
-   Genome-Aligned BAM & BAM index, 21GB
-   Feature / Cell (spot) matrix "raw", 62MB \[22MB HDF5\]
-   Feature / Cell (spot) matrix "filtered", 52MB \[19MB HDF5\] ^ 10X
    also provides a clustering analysis (27MB) and per-molecule read
    information (309MB), and additional images (9.66MB)

#### Seurat Vignette for [Spatial Analysis](https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html)

-   Can load 10X data from spaceranger pipeline
    -   Or can use SeuratData package to download filtered feature/cell
        matrix from within the R notebook
-   Normalize and scale using scTransform
    -   Use NormalizeData() on scRNA-seq (log normalizes), but use
        scTransform on spatial data (different samples have different
        cell densities & read volumes!)
    -   scTransform normalizes, detecs high-variance features (the top
        ~2000 genes that will be used in downstream steps), and saves as
        a new Assay in the Seurat object
-   Dim reduction and clustering
    -   RunPCA, FindNeighbors, FindClusters, RunUMAP... results are
        saved in the Sueurat object itself
-   Identify spatially variable genes
    -   Can use the cluster borders from the previous step as labeling
        different sample conditions for DE
    -   Or can use markvariogram (Trendsceek) mark point process
    -   Others not implemented in Seurat: SpatialDE, Splotch
-   Integration with Single-Cell Data
    -   Deconvolute the spatial samples into component cell types
    -   Or can use anchor-based integration methods ("label trransfer")
        for matching to reference single cells' annotations
        -   Can then use markvariogram again to identify **cell types**
            that are spatially variable

Read in the brain data from the Seurat (v3.2) [vignette on Spatial
analysis](https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html):

    library(Seurat)
    brain <- readRDS("brain_final.rds")

The data was already normalized and scaled using Seurat's `sctransform`.
The original data is saved into the `Spatial` assay, while the
transformed data is saved into the `SCT` assay:

    brain@assays

    ## $Spatial
    ## Assay data with 31053 features for 2696 cells
    ## First 10 features:
    ##  Xkr4, Gm1992, Gm37381, Rp1, Sox17, Gm37323, Mrpl15, Lypla1, Gm37988,
    ## Tcea1 
    ## 
    ## $SCT
    ## Assay data with 17668 features for 2696 cells
    ## Top 10 variable features:
    ##  Ttr, Plp1, Hbb-bs, Hba-a1, Mbp, Penk, Ptgds, Hba-a2, S100a5, Ppp1r1b

The data contains gene counts for each sample spot (n=2,696). Samples
are named AAACAAGTATC... and are referred to as 'cells' in the data
structure (since these data structures were originally designed for
scRNA-seq data).

    as.matrix(brain@assays$SCT@data[1:5,1:5])

    ##        AAACAAGTATCTCCCA-1 AAACACCAATAACTGC-1 AAACAGAGCGACTCCT-1
    ## Xkr4             0.000000          0.0000000          0.0000000
    ## Sox17            0.000000          0.6931472          0.0000000
    ## Mrpl15           1.386294          0.0000000          1.6094379
    ## Lypla1           0.000000          0.6931472          0.6931472
    ## Tcea1            1.386294          0.6931472          0.6931472
    ##        AAACAGCTTTCAGAAG-1 AAACAGGGTCTATATT-1
    ## Xkr4            0.0000000          0.0000000
    ## Sox17           0.0000000          0.0000000
    ## Mrpl15          0.0000000          0.6931472
    ## Lypla1          0.6931472          0.0000000
    ## Tcea1           1.0986123          1.3862944

Compare the normalized data above with the original count data for the
same samples and genes:

    target.genes.from.SCT <- brain@assays$Spatial@data@Dimnames[[1]] %in% row.names(brain@assays$SCT@data)[1:5]
    as.matrix(brain@assays$Spatial@data[target.genes.from.SCT,1:5])

    ##        AAACAAGTATCTCCCA-1 AAACACCAATAACTGC-1 AAACAGAGCGACTCCT-1
    ## Xkr4                    0                  0                  0
    ## Sox17                   0                  1                  0
    ## Mrpl15                  2                  1                  4
    ## Lypla1                  0                  1                  1
    ## Tcea1                   2                  2                  1
    ##        AAACAGCTTTCAGAAG-1 AAACAGGGTCTATATT-1
    ## Xkr4                    0                  0
    ## Sox17                   0                  0
    ## Mrpl15                  0                  1
    ## Lypla1                  2                  0
    ## Tcea1                   3                  4

### Upstream Data

The data from this experiment is also available [directly from
10xGenomics](https://www.10xgenomics.com/resources/datasets/). I think I
found the counts matrix (downloaded the filtered version), which is
saved as an ".mtx" file. I think this is a 10X file format produced by
their CellRanger program. Seurat can probably load mtx:

    target.dir <- paste(getwd(), '/Data files From 10x/filtered_feature_bc_matrix/', sep="")

Seurat vignette 1 actually describes the read10X function for reading in
data from 10X's cellranger pipeline. The data from cellranger is
organized into: 1. genes.tsv.gz - list of features 2. barcodes.tsv.gz -
list of samples 3. matrix.mtx.gz - matrix-like of counts

The filtered\_feature\_bc\_matrix from 10X follows this format and
should be readable by Seruat's Read10X()

    library(dplyr)

    ## Warning: package 'dplyr' was built under R version 3.6.3

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    library(Seurat)
    library(patchwork)

    ## Warning: package 'patchwork' was built under R version 3.6.3

    from.10x <- Read10X(data.dir=target.dir)
    str(from.10x)

    ## Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    ##   ..@ i       : int [1:15319170] 6 9 12 18 19 28 32 36 38 61 ...
    ##   ..@ p       : int [1:2697] 0 4242 12102 18434 26391 34182 40473 47960 55003 59879 ...
    ##   ..@ Dim     : int [1:2] 31053 2696
    ##   ..@ Dimnames:List of 2
    ##   .. ..$ : chr [1:31053] "Xkr4" "Gm1992" "Gm37381" "Rp1" ...
    ##   .. ..$ : chr [1:2696] "AAACAAGTATCTCCCA" "AAACACCAATAACTGC" "AAACAGAGCGACTCCT" "AAACAGCTTTCAGAAG" ...
    ##   ..@ x       : num [1:15319170] 2 2 5 1 1 1 2 2 1 6 ...
    ##   ..@ factors : list()

Confirm that the data from 10x (filtered version) matches the Seurat
spatial vignette dataset:

    genes <- from.10x@Dimnames[[1]] %in% row.names(brain@assays$SCT@data)[1:5]
    samples <- paste(from.10x@Dimnames[[2]],"-1",sep="") %in% colnames(brain@assays$SCT@data)[1:5]
    from.10x[genes,samples]

    ## 5 x 5 sparse Matrix of class "dgCMatrix"
    ##        AAACAAGTATCTCCCA AAACACCAATAACTGC AAACAGAGCGACTCCT AAACAGCTTTCAGAAG
    ## Xkr4                  .                .                .                .
    ## Sox17                 .                1                .                .
    ## Mrpl15                2                1                4                .
    ## Lypla1                .                1                1                2
    ## Tcea1                 2                2                1                3
    ##        AAACAGGGTCTATATT
    ## Xkr4                  .
    ## Sox17                 .
    ## Mrpl15                1
    ## Lypla1                .
    ## Tcea1                 4

From the Seurat vignette:

    as.matrix(brain@assays$Spatial@data[target.genes.from.SCT,1:5])

    ##        AAACAAGTATCTCCCA-1 AAACACCAATAACTGC-1 AAACAGAGCGACTCCT-1
    ## Xkr4                    0                  0                  0
    ## Sox17                   0                  1                  0
    ## Mrpl15                  2                  1                  4
    ## Lypla1                  0                  1                  1
    ## Tcea1                   2                  2                  1
    ##        AAACAGCTTTCAGAAG-1 AAACAGGGTCTATATT-1
    ## Xkr4                    0                  0
    ## Sox17                   0                  0
    ## Mrpl15                  0                  1
    ## Lypla1                  2                  0
    ## Tcea1                   3                  4

    dim(from.10x)

    ## [1] 31053  2696

    dim(brain@assays$Spatial)

    ## [1] 31053  2696

In the Seurat tutorial vignette (\#1), the 10X cellranger data is then
converted into a Large Seurat object. Can try to do that here. They
require genes to be counted in at least 3 different cells (sample spots)
and require each sample to have at least 200 different counted genes.

*The Large Seurat object holds the data and also the results from
downstream analysis (like PCA).* The Seurat object's `@meta.data` holds
QC metrics like number of reads (nCount\_RNA) and number of genes in
that sample (nFeature\_RNA).

For spatial data, should use `Load10X_Spatial()` to load from the
spaceranger pipeline and return a Seurat object with the expression
levels and an image of the tissue slide.

    # Initialize the Seurat object with the raw (non-normalized data).
    ben.brain <- CreateSeuratObject(counts = from.10x, project = "ben.brain", min.cells = 3, min.features = 200)
    ben.brain

    ## An object of class Seurat 
    ## 18449 features across 2696 samples within 1 assay 
    ## Active assay: RNA (18449 features, 0 variable features)

The min.cells (samples) filter lowered the feature list length from
31,053 down to 18,449. The actual spatial vignette did not apply
min.cells like the tutorial vignette did.

    ben.genes <- ben.brain[["RNA"]]@data@Dimnames[[1]] %in% row.names(brain@assays$SCT@data)[1:5]
    ben.samples <- paste(ben.brain[["RNA"]]@data@Dimnames[[2]],"-1",sep="") %in% colnames(brain@assays$SCT@data)[1:5]
    ben.brain@assays$RNA@data[ben.genes,ben.samples]

    ## 5 x 5 sparse Matrix of class "dgCMatrix"
    ##        AAACAAGTATCTCCCA AAACACCAATAACTGC AAACAGAGCGACTCCT AAACAGCTTTCAGAAG
    ## Xkr4                  .                .                .                .
    ## Sox17                 .                1                .                .
    ## Mrpl15                2                1                4                .
    ## Lypla1                .                1                1                2
    ## Tcea1                 2                2                1                3
    ##        AAACAGGGTCTATATT
    ## Xkr4                  .
    ## Sox17                 .
    ## Mrpl15                1
    ## Lypla1                .
    ## Tcea1                 4
