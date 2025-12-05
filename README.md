# QeITH

Here, we introduce QeITH to quantify Intratumor heterogeneity (ITH) within the tumor ecosystem.

&nbsp
&nbsp;

# Description

QeITH, a novel entropy-based computational framework designed to quantify ITH by measuring the diversity and distribution of cellular populations and associated functional states within the tumor ecosystem. QeITH applies to single-cell, bulk, and spatial transcriptomic data, which revealed that elevated QeITH scores are intrinsic markers of malignant transformation and consistently correlate with increased tumor aggressiveness.



&nbsp;

# Details

The function `QeITH()` is used to calculate ITH score. It supports **three common data modalities** in cancer research, containing two parameters:

**Parameters**

+ **`data`** : The input dataset, with format **strictly dependent** on the `type`:

  * **`type = "bulk"`** : A gene expression **matrix** or **data.frame**. Rows are gene symbols, columns are sample IDs.

  * **`type = "single-cell"`** : A cell annotation **data.frame** that **must contain** the two columns: `Samples` and `Celltype`.

  * **`type = "spatial"`** : A cell proportion **matrix** or **data.frame**. Rows are samples/spots, ^columns are cell types. Each element is the proportion of a cell type in a sample, and **each row must sum to 1**.

+ **`type`** : A character string specifying the analysis mode. **Must be one of:** `"bulk"`, `"single-cell"`, or `"spatial"`

  + 

    

# Installation

- Users can install the released version of **QeITH** with:
  &nbsp;

```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("XS-Wang-Lab/QeITH")
```

&nbsp;
&nbsp;

# Examples

&nbsp
&nbsp;

### **Apply QeITH in bulk transcriptome**

```R
library(QeITH)
example_file_path <- system.file("extdata", "example.RData", package = "QeITH")
load(example_file_path)
#ls()
#"bulk_example"  "sc_example" "spatial_example"
QeITH_bulk <- QeITH(bulk_example, type = "bulk")
```



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




## Apply QeITH to single-cell 

```R
QeITH_sc <- QeITH(sc_example, type = "single-cell")
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



## **Apply QeITH to spatial transcriptome data**

```R
#Using the deconvolution algorithm to infer cell proportion. Here we used CARD algorithm (https://github.com/YMa-lab/CARD). Users can choose other algorithm to instead.
#Example
#library(CARD) 
#CARD_obj = createCARDObject(
#	sc_count = sc_count,
#	sc_meta = sc_meta,
#	spatial_count = spatial_count,
#	spatial_location = spatial_location,
#	ct.varname = "Celltype",
#	ct.select = unique(sc_meta$Celltype),
#	sample.varname = "Samples",
#	minCountGene = 100,
#	minCountSpot = 5) 
# CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
# spatial_example = as.data.frame(CARD_obj@Proportion_CARD) 

QeITH_ST <- QeITH(spatial_example, type = "spatial")
```



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
