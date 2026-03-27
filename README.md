# QeITH

QeITH (Quantifying Ecosystem Intratumor Heterogeneity), a universally applicable computational framework designed to quantify the ecosystem heterogeneity index across single-cell, bulk, and spatial resolutions. 

<img width="5322" height="4948" alt="QeITH" src="https://github.com/user-attachments/assets/236705b4-2552-45e0-a04b-71612c9f0dc5" />

&nbsp;

# Description

**QeITH (Quantifies  Ecosystem Intratumor Heterogeneity)** is an entropy-based computational framework designed to quantify ecosystem heterogeneity by measuring the diversity and distribution of cellular populations and associated functional states within the tumor ecosystem.

**Key features:**

- **Multi-modal support**: Works with bulk transcriptomics, single-cell annotations, and spatial transcriptomics data
- **Simple unified interface**: Single function `QeITH()` for all data types

&nbsp;

# Details

The function `QeITH()` is used to calculate ITH score. It supports **three common data modalities** in cancer research, containing two parameters:

**Parameters**

The **input dataset** format depends on the **data modality** specified by the `type` parameter:

+ | Parameter | Type              | Description       | Requirements                                              |
  | :-------- | :---------------- | :---------------- | :-------------------------------------------------------- |
  | `input`   | matrix/data.frame | The input dataset | Format depends on `type` (see below)                      |
  | `type`    | character         | Analysis mode     | Must be one of: `"bulk"`, `"single-cell"`, or `"spatial"` |

+ **Input** Format Requirements by `type`

  | `type`          | Required Format            | Column/Row Requirements                                      | Example           |
  | :-------------- | :------------------------- | :----------------------------------------------------------- | :---------------- |
  | `"bulk"`        | Gene expression matrix     | **Rows**: Gene symbols **Columns**: Sample IDs               | `bulk_example`    |
  | `"single-cell"` | Cell annotation data.frame | **Must contain two columns**: - `Samples`: Sample identifiers - `Celltype`: Cell type annotations | `sc_example`      |
  | `"spatial"`     | Cell proportion matrix     | **Rows**: Samples/spots **Columns**: Cell types **Values**: Proportions (0-1), each row sums to 1 | `spatial_example` |

+ **Output**

  Returns a data.frame with two columns:

  - `Sample`: Sample identifier
  - `ITH_Score`: Calculated intratumor heterogeneity score


  ## Interpreting QeITH Scores

  - **Higher scores** indicate greater cellular diversity and more evenly distributed cell type proportions or functional states within a sample;
  - **Lower scores** suggest dominance by a few cell types or functional states, indicating reduced heterogeneity;
  - **Note**: Scores from different datasets should be interpreted with caution, as they may not be directly comparable due to differences in underlying calculations.

&nbsp;

# Installation

- Users can install the released version of **QeITH** with:
  &nbsp;

```R
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# Install QeITH
devtools::install_github("XS-Wang-Lab/QeITH")

# Load the package
library(QeITH)
```

&nbsp;



# Quick Start Examples

Below are working examples for each data type.



### Load Example Data&nbsp;

```R
library(QeITH)

# Load built-in example datasets
example_file_path <- system.file("extdata", "example.RData", package = "QeITH")
load(example_file_path)
#ls()
#"bulk_example"  "sc_example" "spatial_example"
# The workspace now contains three example datasets:
# - bulk_example     : Bulk transcriptomics matrix
# - sc_example       : Single-cell annotations
# - spatial_example  : Spatial cell type proportions
```

&nbsp;



### **Apply QeITH in bulk transcriptome**

**Input format**: Gene expression matrix with genes as rows and samples as columns.

```R
# Examine the data
bulk_example[1:5, 1:5]

# Calculate ITH scores
QeITH_bulk <- QeITH(bulk_example, type = "bulk")

# View results
head(QeITH_bulk)
```

&nbsp;

**bulk_example**

```R
bulk_example[1:5,1:5]
```

| row.names  | MB-0362  | MB-0346  | MB-0386  | MB-0574  | MB-0185  |
| ---------- | -------- | -------- | -------- | -------- | -------- |
| **RERE**   | 8.676978 | 9.653589 | 9.033589 | 8.814855 | 8.736406 |
| **RNF165** | 6.075331 | 6.687887 | 5.910885 | 5.62874  | 6.392422 |
| **PHF7**   | 5.83827  | 5.600876 | 6.030718 | 5.849428 | 5.542133 |
| **CIDEA**  | 6.397503 | 5.246319 | 10.11182 | 6.116868 | 5.184098 |
| **TENT2**  | 7.906217 | 8.267256 | 7.959291 | 9.206376 | 8.162845 |

&nbsp;

**QeITH_bulk**

```R
QeITH_bulk[1:5,]
```

| **Sample**  | ITH_Score |
| :---------- | :-------: |
| **MB-0362** | 3.717366  |
| **MB-0346** | 3.669606  |
| **MB-0386** | 3.728047  |
| **MB-0574** | 3.679429  |
| **MB-0185** | 3.716163  |

&nbsp;



## Apply QeITH in a single-cell dateset

**Input format**: Data frame with **exactly two columns** named `Samples` and `Celltype`. Each row represents one cell.

```R
# Examine the data
sc_example[1:5,]

# Calculate ITH scores
QeITH_sc <- QeITH(sc_example, type = "single-cell")

# View results
head(QeITH_sc)
```

&nbsp;

**sc_example**

```R
sc_example[1:5,]
```

| **row.names**              | Samples | Celltype  |
| -------------------------- | ------- | --------- |
| **P01_AAACCTGAGGACACCA-1** | P01     | Malignant |
| **P01_AAACCTGCAAGTCTGT-1** | P01     | NK_cell   |
| **P01_AAACCTGCACCGAAAG-1** | P01     | B_cell    |
| **P01_AAACCTGCAGTATAAG-1** | P01     | NK_cell   |
| **P01_AAACCTGCATGCATGT-1** | P01     | B_cell    |

&nbsp;

**QeITH_sc**

```R
QeITH_sc[1:5,]
```

| **Sample** | ITH_Score |
| :--------- | :-------: |
| **P01**    | 2.750233  |
| **P02**    | 1.719219  |
| **P03**    | 2.227526  |
| **P04**    | 2.260293  |
| **P05**    | 2.524838  |

&nbsp;

## **Apply QeITH to a spatial transcriptome data**set

**Input format**: Cell type proportion matrix with samples/spots as rows, cell types as columns. **Each row must sum to 1**.

```R
# This is a conceptual example showing how to obtain a proportion matrix using the CARD algorithm.
# For full usage details, please refer to the CARD package documentation: https://github.com/YMa-lab/CARD
# Note: Users can also choose other deconvolution algorithms.

#library(CARD) 
#CARD_obj <- createCARDObject(
#	sc_count = sc_count,
#	sc_meta = sc_meta,
#	spatial_count = spatial_count,
#	spatial_location = spatial_location,
#	ct.varname = "Celltype",
#	ct.select = unique(sc_meta$Celltype),
#	sample.varname = "Samples",
#	minCountGene = 100,
#	minCountSpot = 5) 
# CARD_obj <- CARD_deconvolution(CARD_object = CARD_obj)
# spatial_example <- as.data.frame(CARD_obj@Proportion_CARD) 

# Examine the data
spatial_example[1:5, 1:5]

# Calculate ITH scores
QeITH_ST <- QeITH(spatial_example, type = "spatial")

# View results
head(QeITH_ST)
```

&nbsp;

**spatial_example**

```R
spatial_example [1:5,1:5]
```

| **row.names**            | NK         | CAF         | Treg        | Plasmablasts |
| ------------------------ | ---------- | ----------- | ----------- | ------------ |
| **AAACAGCTTTCAGAAG-1_1** | 0.4873921  | 0.07230232  | 0.1175443   | 0.03521669   |
| **AAACAGTGTTCCTGGG-1_1** | 0.3972463  | 0.000000000 | 0.000000000 | 0.1248904    |
| **AAACCTCATGAAGTTG-1_1** | 0.01487764 | 0.000000000 | 0.000000000 | 0.000000000  |
| **AAACTCGGTTCGCAAT-1_1** | 0.1186421  | 0.006878719 | 0.0724607   | 0.04485068   |
| **AAAGGGATGTAGCAAG-1_1** | 0.3091865  | 0.000000000 | 0.000000000 | 0.09319035   |

&nbsp;

**QeITH_ST**

```R
QeITH_ST[1:5,]
```

| **Sample**               | ITH_Score |
| :----------------------- | :-------: |
| **AAACAGCTTTCAGAAG-1_1** | 1.3530028 |
| **AAACAGTGTTCCTGGG-1_1** | 0.7936923 |
| **AAACCTCATGAAGTTG-1_1** | 0.8356857 |
| **AAACTCGGTTCGCAAT-1_1** | 1.6211930 |
| **AAAGGGATGTAGCAAG-1_1** | 0.7807895 |

# Contact

E-mail any questions to Xiaosheng Wang (xiaosheng.wang@hotmail.com)
