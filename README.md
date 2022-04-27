# ShallowChrome

## About
This repository contains the implementation of the `ShallowChrome` modeling pipeline presented in the paper:
```Frasca F., Matteucci M., Leone M., Morelli M. J. and Masseroli M. "Accurate and highly interpretable prediction of gene expression from histone modifications", 2022; 23: 151``` available [here] (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04687-x).
`ShallowChrome` is a novel computational pipeline for accurate and fully interpretable modeling of epigenetic gene transcriptional regulation operated by Histone Mark (HM) modifications. `ShallowChrome` leverages on the procedure of 'peak calling' to retrieve gene-wise, significant and dynamically located HM features that can strongly predict the transcriptional state of genes. In our modeling pipeline we:
- Fit logistic regression models on these extracted features to solve the task of binary classification of gene transcriptional state over 56 cell-types from the REMC database;
- Analyse and rigorously interpret the obtained models by extracting insightful gene-specific regulative patterns;
- Compare the extracted patterns with the characteristic chromatin state emissions from ChromHMM (Ernst et al., 2012), showing that `ShallowChrome` is able to coherently rank groups of chromatin states w.r.t. their transcriptional activity. 

More on how to replicate paper results is in the following.

## Structure
```
ShallowChrome/
|-- README.md
|-- LICENSE
|-- .gitignore
|-- notebooks/
|   |-- utils.py
|   |-- model fitting.ipynb
|   |-- model inspection.ipynb
|   |-- model validation.ipynb
|   |-- model fitting - valley thresholding.ipynb
|   |-- data extraction.ipynb
|-- scores/
|   |-- DeepChrome_scores.txt
|-- data/
|   |-- - splits/
|   |   |-- iteration_0/
|   |   |-- iteration_1/
|   |   |-- iteration_2/
|   |   ...
|   |-- - targets/
|   |   |-- E003/
|   |   |-- E004/
|   |   ...
|   |-- cells.csv
|   |-- gene_list.txt
|   |-- GeneFile.txt
|   |-- names.csv
```

- `README.md` this file
- `LICENSE` MIT license file
- `.gitignore` standard .gitignore file for Python projects
- `notebooks/` folder containing Python notebooks to run the modeling pipeline
- `notebooks/utils.py` core Python routines called from within the notebooks to perform modeling and analyses
- `notebooks/model fitting.ipynb` notebook where ShallowChrome models are fitted to solve binary classification of gene transcriptional state; reproduces Tables 2 and S1 and Figure 2 of the paper
- `notebooks/model inspection.ipynb` notebook to inspect ShallowChrome models and to extract and plot gene-wise regulative patterns; reproduces Figure 3 of the paper
- `notebooks/model validation.ipynb` notebook to compare ShallowChrome regulative patterns with ChromHMM chromatin state emissions; reproduces Figure 4 of the paper
- `notebooks/model fitting - valley thresholding.ipynb` here the classification task is solved with an alternative approach to define target classes; reproduces Table S4 and Figures S2 and S3 of the paper
- `notebooks/data extraction.ipynb` notebook to perform data extraction with pygmql and pandas
- `scores/` default folder where numerical results from the modeling pipeline are stored
- `scores/DeepChrome_scores.txt` test scores from the DeepChrome model (Singh et al., 2016)
- `data/` default folder where data and reference files are stored
- `data/- splits/` folder containing random split indices for model fitting
- `data/- targets/` folder containing RPKM target values for each epigenome
- `data/cells.csv` csv file enumerating the 56 epigenomes object object of the present study
- `data/gene_list.txt/` txt file containing the ordered list of genes considered in the present study
- `data/GeneFile.txt` txt file containing promoter window coordinates for each of the considered genes
- `data/names.csv` csv file enumerating the Histone Mark modifications considered in the present study

## Requirements
In order to run the `ShallowChrome` model fitting and analyses (notebooks `model fitting.ipynb`, `model inspection.ipynb`, `model validation.ipynb`, `model fitting - valley thresholding.ipynb`), the following libraries are required:
- matplotlib
- numpy
- scikit
- scipy
- jupyter
We suggest installing them within a Python virtual environment via pip. Paper results can be reproduced with the following versions on `Python 2.7.15`:
```
matplotlib==2.2.4      
numpy==1.16.6     
scikit-learn==0.20.4     
scipy==1.2.2      
```
In order to run the _de novo_ data extraction notebook (`data extraction.ipynb`) the following libraries are required:
- numpy
- gmql
- pandas
- jupyter 
We employed the following versions on `Python 3.6.12`
```
numpy==1.16.0
gmql==0.1.1
pandas==1.1.5
```
NB: `gmql` additionally requires Java. Please follow the installation procedure [here](https://pygmql.readthedocs.io/en/latest/installation.html).

## Reproducing paper results
1. Download the pre-processed data from [here](https://zenodo.org/record/4445287);
2. Uncompress the downloaded ".zip" file in folder "data";
3. Run the `notebooks/model fitting.ipynb` notebook to reproduce Tables 2 and S1 and Figure 2 of the paper;
4. Run the `notebooks/model inspection.ipynb` notebook to reproduce Figure 3 of the paper;
5. Run the `notebooks/model inspection.ipynb` notebook with variable `target_only` set to `False`: this will perform model selection and fitting for _all_ epigenomes over the 'standard' DeepChrome split, dumping all fitted models to disk;
6. Run the `notebooks/model validation.ipynb` notebook to reproduce Figure 4 of the paper;
7. Run the `notebooks/model fitting - valley thresholding.ipynb` notebook to reproduce Table S4 and Figures S2 and S3 of the paper.
Alternatively:
1. Run the `notebooks/data extraction.ipynb` to prepare all necessary data for model fitting and analyses;
2. Go to step 3. above.
