# space_for_time_analysis
file created by Lena and Noelle March 2022, last edit : 19.04.22

This folder contains parts of the analysis code for the space for time analysis. In this folder specifically : the matrix permutation in order to calculate p-values, variance importance and deviance explained.

The data files are stored in a separate folder "../data/" (not uploaded to GitHub) since some were too big to be uploaded, and not finally decided where to deposit. 

Note : this readme was written for local use, and might contain several files not under version control (file format not suitable).

## Overview over folder content

**Functions and function development**
- Edited matrix permutation : split into 2 functions, as in the `gdm` package :
    - `gdm.varIMP_edit.R` contains the function `gdm.varImp.synthesis.edit()`, which is called by the user.
    - `matrix_perm_permutateSitePair.R` contains the function `permutateSitePair_edit()` and `permutateVarSitePair_edit()`, which are called within `gdm.varImp.synthesis.edit()`
- `22-02_test_varimp.edit_function.R` was used while developing the function, running some simpler tests.
- old versions : for nostalgic and documentation reasons, some selected depreciated material is kept :
    - `21-07_DEPRECIATED_old_gdm.varimp.R` : the edited gdm.varimp function, as used in July 2021.
    -  `22-02_DEPRECIATED_permutatesitepair_without_respecting_nonindependence_of_rows.R` : depreciated function permutatesitepair, before the non-independence of rows was introduced : This function permutes rows without going back to the site-predictor matrix (stays in dissimilarity-matrix format). This function could probably be used in another project (BetaDivMultifun)


**Scripts for running the functions**
- 6 Scripts `GDM_varIMP_<...>.R` contain code to run matrix permutation on all GDM models in the manuscript.
    - `<...>` is one of : 
        - plant_alpha_scaled : plant alphadiversity, space and time
        - plants             : plant betadiveresity, space and time
        - herb_alpha_scaled  : insect herbivores alphadiversity, space and time
        - herb_scaled        : insect herbivores betadiversity, space and time
        - pred_alpha_scaled  : insect secondary consumers alpha diversity, space and time
        - pred_scaled        : insect secondary consumers beta diversity, space and time
    - Each script calls one organism and either alpha or betadiversity GDMs, and all Chao Numbers (plus beta sim for betadiversity), as well as both LUI and it's components (Mowing, Grazing, Fertilisation)
        - within script : results from temporal models contain "time_ayrs" in the name, if not the variable refers to a spatial variable.

**Results**
- output from the permutations are stored in the folder
    - `Mar-01_results/` and reported in the Table `GDM_varEX_pValues.xlsx`

**Documentation**
- At the moment in .docx format locally
