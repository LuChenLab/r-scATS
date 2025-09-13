## scATS
![GitHub release (latest SemVer)](https://img.shields.io/badge/Version-v0.5.4-yellowgreen) ![GitHub release (latest SemVer)](https://img.shields.io/badge/Language-R-yellowgreen)


**An R package for quantifying alternative transcription start site (ATS) in 5' single-cell RNA-seq data**

Here, we developed scATS, which introduces two novel metrics (**α** and **β**) to precisely quantify RNA degradation, and provides both raw (**ψ**) and corrected (**θ**) metrics to correct ATS expression profiles from resulting distortions. 

overall, scATS is a stepwise computational method, scATS, to de novo infer the relative usage of TSSs and correct the RNA degradation in the paired-end 5’ scRNA-seq.

## Contents

- [Overview](#Overview)
- [Requirements](#Requirements)
- [Installation](#Installation)
- [Tutorial](#Tutorial)
- [LRS model](#LRS_model)
- [Citing](#Citing)


## <span id="Overview">Overview</span>

`scATS` contains three main steps (**Fig. a**): 
1. identifing TSS in each gene according to the eCDF of the start sites of all read 1 (R1) reads.
2. quantifying an EM model based RNA degradation rate (**α**) and a degradation metric (**β**), representing the relative area under the cumulative distribution curve of each TSS.
3. providing raw percent spliced in (PSI, **ψ**) and corrected PSI values (**θ**) of each TSS account for RNA degradation.

We further applied scATS on three different 5’ scRNA-seq datasets, including mouse hematopoietic stem and progenitor cells (HSPCs), and human clinical samples from COVID-19 and lung cancer patients (**Fig. b**). Benchmark was performed in mouse HSPCs datasets. And the accuracy of our model was assessed by promoter motifs, epigenetic information, conservation and long-read sequencing (**Fig. c**). We also developed a machine-learning method to screen for the functional TSSs in lung cancer and guide subsequent experimental validation and mechanism study (**Fig. d**).

![scATS](./scATS.png)


## <span id="Requirements">Requirements</span>

scATS has been developed with `R 4.0.0` and the following packages are needed to be installed in `R` (R scripts in `scATS` automatically check to see if an R library is already installed and then install those that are needed. So no need for manual preinstallation):

|        module        |    version   |
| -------------------- | ------------ |
| data.table           | >= 1.15.4    |
| fitdistrplus         | >= 1.2-1     |
| GenomicAlignments    | >= 1.22.1    |
| GenomicFeatures      | >= 1.38.2    |
| GenomicRanges        | >= 1.38.0    |
| ggbio                | >= 1.42.0    |
| mclust               | >= 6.1.1     |
| mixtools             | >= 2.0.0     |
| parallel             | >= 3.6.0     |
| rtracklayer          | >= 1.46.0    |
| Seurat               | >= 4.4.0     |
| SummarizedExperiment | >= 1.16.1    |

## <span id="Installation">Installation</span>

To install `scATS`, you have two options: either install directly from GitHub or use the compressed source file:

```r
# Install from GitHub if remotes package is not installed
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("LuChenLab/r-scATS/scATS/")
```

Alternatively, you can install `scATS` using the source file downloaded from the [repository](https://github.com/LuChenLab/r-scATS/blob/main/scATS_0.5.4.tar.gz) :
```r
# Install scATS from a downloaded source file
R CMD INSTALL scATS_0.5.4.tar.gz
```


## <span id="Tutorial">Tutorial</span>


To use scATS, you can follow the instructions below. The detailed installation and usage of scATS are available in our [document](https://r-scats.readthedocs.io/en/latest/).

Demo data can be obtained from [repository](https://github.com/LuChenLab/r-scATS/tree/main/demo).


## <span id="LRS_model">Tutorial</span>

The training dataset contains 33 features. Detailed feature categories can be found in the column headers of the data file.

For model development:

The dataset was split into training, validation, and test sets using **StratifiedKFold cross-validation** with `n_splits=10`, `shuffle=True`, and `random_state=42`.

Features were standardized using **StandardScaler**.

A **Random Forest** Classifier (**RandomForestClassifier**) was used as the base model.

Both preprocessing and the classifier were wrapped into a **Pipeline**, ensuring they can be saved and reloaded together.

Resources in this repository:

Training data: `dataset.xlsx`

Trained model file: `rf_pipeline_model.pkl`

You can directly load the trained model from [here](https://github.com/LuChenLab/r-scATS/tree/main/LRS) for predictions without applying feature scaling again. Example:

```python
import pickle

# Load the trained pipeline
with open("rf_pipeline_model.pkl", "rb") as f:
    model = pickle.load(f)
```

## <span id="Citing">Citing</span>

A pre-print is going to be uploaded soon.



