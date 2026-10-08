# Gene-family copy numbers and fungal growth form

This workflow uses OrthoFinder gene-family copy numbers to classify fungal growth forms with XGBoost random forests (XGBRF). The main analysis compares strict unicellular yeasts with strict multicellular fungi. Three sensitivity analyses include dimorphic species under different label definitions.

## Inputs and labels

The gene-count matrix comes from the OrthoFinder analysis in directory 3: `orthofinder_result/Orthogroups.GeneCount.tsv`. The two label tables are in `related_files/`. These analysis-specific labels distinguish strict yeast, dimorphic and strict multicellular species.

| Analysis | Species | Label definitions |
| --- | ---: | --- |
| Main | 55 | 0: strict multicellular fungi (33); 1: strict unicellular yeasts (22). Dimorphic species are excluded. |
| Scenario A | 76 | 0: non-yeast-forming fungi (33); 1: yeast-forming fungi, including strict yeasts and dimorphic species (43). |
| Scenario B | 76 | 0: strict multicellular and dimorphic fungi (54); 1: strict yeasts (22). |
| Scenario C | 76 | 0: strict multicellular fungi (33); 1: strict yeasts (22); 2: dimorphic fungi (21). |

## 1. Prepare the copy-number matrices

Run the scripts and notebooks from this directory.

1. **`main_combine_label.py`** matches the 55 main-analysis species to the OrthoFinder matrix and writes `related_files/ml_55_matrix_main_labeled.tsv`.
2. **`sens_combine_label.py`** builds the corresponding 76-species matrix with `sens_A`, `sens_B` and `sens_C` labels. It also requires a unique match for every species.
3. **`main_preprocess.ipynb`** and **`sens_preprocess.ipynb`** remove constant OGs and retain OGs present in at least 5% of species. All retained copy numbers are transformed with `log1p` and written to the corresponding `*_log1p_cleaned.tsv` files. The skewness summary is diagnostic; it does not decide whether to apply the transformation. The optional PCA preview uses the untransformed, filtered counts.

Both starting matrices contain 38,091 OGs. Filtering retains 15,950 OGs for Main and 15,487 for the sensitivity analyses. These matrices are in `related_files/`.

## 2. Train and evaluate the classifiers

Run **`main_xgb.ipynb`** for Main, then **`sens_xgb.ipynb`** for Scenarios A–C. Results are written to `results/ml_output_main_xgb/` and `results/ml_output_sens_xgb/sens_A/`, `sens_B/` and `sens_C/`.

All analyses use five-fold stratified cross-validation repeated 20 times, giving 100 outer test folds. Within each outer training set, three-fold grid search selects the model settings. Binary models optimize average precision and use a probability threshold of 0.5. Scenario C optimizes macro F1 and predicts one of three classes. Class weights are calculated from each outer training set.

The search covers 300 or 500 trees, depth 2 or 3, minimum child weight 1 or 2, subsampling 0.6 or 0.8, and per-node feature sampling 0.3 or 0.6. The splits are random and stratified by class; they are not phylogenetic blocks. Constant/sparse feature filtering is performed before cross-validation.

Each notebook saves:

- Fold metrics and their mean, SD, median, minimum and maximum.
- Out-of-fold predictions, pooled confusion counts and class-level performance. Each species contributes 20 test predictions, so pooled counts are not counts of independent species.
- Gain and absolute SHAP importance, summarized by mean, SD, median and the fraction of outer folds with nonzero importance. SHAP values are calculated on the outer test folds; Scenario C also averages absolute values across classes.
- ROC/precision–recall curves, confusion matrices and top-30 importance plots. Curve shading represents variation across folds, not a confidence interval.

The field `pr_auc` is calculated with `average_precision_score`. Main reports positive-class precision, recall and F1; the sensitivity analyses report macro-averaged versions. Scenario C uses macro one-versus-rest ROC-AUC and macro average precision.

## 3. Combine the importance rankings

**`RRA_rank_calcu.py`** combines four lists: Main gain, Main SHAP, Scenario A gain and Scenario A SHAP. It reranks the 13,800 shared OGs within this common feature set, using the minimum rank for ties.

For each OG, normalized ranks are combined using beta order-statistic probabilities. The minimum probability is multiplied by the number of lists, followed by BH correction across the shared OGs. The four lists are related analyses, so these scores summarize rank consistency.

Core candidates must rank within the top 30 in all four common-feature lists, have nonzero importance in at least 90% of outer folds in every list, and have positive median importance in every list. The supplied results contain eight core OGs. RRA P or Q values are used for ranking and reporting, not as an additional selection threshold.

**`RRA_rank_pic.py`** plots the core ranking and the four component ranks. Its `TARGET_OG` setting can highlight a selected OG; an empty value leaves all candidates in the same style. RRA tables and plots are saved to `results/ml_output_rra/`.

## Software

Python requires NumPy, pandas, scikit-learn, XGBoost, SHAP, matplotlib, seaborn and Jupyter/IPython.
