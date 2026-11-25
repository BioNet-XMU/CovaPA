# CovaPA

Covariance-aware pathway analysis (CovaPA) compares metabolite variability and covariation between groups at the pathway level. This repository contains the CovaPA reference implementation, two abundance-based comparators, and prepared inputs for the three clinical cohorts analyzed in *Metabolite Variability and Covariation Reveal Pathway-Level Changes beyond Mean Abundance*.

## Run the cohort analyses

Use Python 3.10 or later and R 4.0 or later. The analysis uses NumPy and base R. OpenPyXL is used by the workbook preparation script.

```bash
python -m pip install -r requirements.txt
python run_analysis.py --permutations 49
python run_analysis.py --permutations 999 --output results
```

The 49-permutation run is a quick check. Use `--cohort CRC`, `--cohort METH`, or `--cohort T2DM` to run one cohort. If Rscript is outside your PATH, pass `--rscript "/path/to/Rscript"`. Files are written by contrast as `{cohort}_{case}_vs_{reference}.csv`. Each includes pathway IDs, member counts, the CovaPA squared Bures distance, the Global Test statistic, ORA hit and background counts, and method-specific p and BH q values.

## Cohorts

| CSV | Specimen and groups | Feature records | Tested pathways | Sources |
| --- | --- | ---: | ---: | --- |
| [CRC.csv](data/CRC.csv) | Serum. 66 CRC cases (C), 76 participants with polyps (P), and 92 healthy controls (H) | 113 | 42 per contrast | [ST000284](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST000284) · [Data DOI](https://doi.org/10.21228/M8FG61) · [Zhu et al.](https://doi.org/10.1021/pr500494u) |
| [METH.csv](data/METH.csv) | Saliva. 170 methamphetamine users (M) and 50 controls (H) | 229 | 6 | [MTBLS7807](https://www.ebi.ac.uk/metabolights/MTBLS7807) · [Wang et al., Supplementary Table 2](https://doi.org/10.1128/msystems.00991-23) |
| [T2DM.csv](data/T2DM.csv) | Serum. 100 participants with T2DM (D) and 200 healthy controls (H) | 146 | 40 | [ST003390](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST003390) · [Data DOI](https://doi.org/10.21228/M81V7G) · [Ma et al.](https://doi.org/10.1038/s41597-025-04652-2) |

Rows are samples, and the first column, `Group`, gives the cohort label. The other columns are feature records in the units of the source analysis matrices. Original sample identifiers and clinical covariates are omitted. [feature_map.csv](data/feature_map.csv) maps features to KEGG Compound IDs. Empty mappings can participate in feature screening but do not enter pathway sets. [pathway_memberships.csv](data/pathway_memberships.csv) records the measured pathway sets and full-family degree weights.

To regenerate the CSVs from the prepared study workbooks, run `python code/prepare_data.py PATH_TO_WORKBOOKS --output PATH_TO_OUTPUT`. The input directory should contain `CRC_data.xlsx` (Sheet1), `MTBLS7807-Supplementary table.xlsx` (table2, published as `msystems.00991-23-s0003.xlsx`), and `YLab_T2DM_FinalData_KEGG_pathways_metabolism_only.xlsx` (Data). The script uses the bundled METH feature mapping and preserves the workbook sample order.

These matrices are derived from the cited datasets. Cite the original data and studies when reusing the matrices. Data permissions are separate from permissions for the analysis code. Consult the [Metabolomics Workbench terms](https://www.metabolomicsworkbench.org/about/termsofuse.php) and [EMBL-EBI terms](https://www.ebi.ac.uk/about/terms-of-use/) before publicly redistributing sample-level data.

## Core methods

- [CovaPA](code/covapa_reference.py) centers the abundance columns within each group, weights them by reciprocal pathway degree, and compares the two regularized sample covariance matrices with the squared Bures distance. The cohort runner uses natural logs, including a CRC floor of 1e-12. Duplicate T2DM records are averaged on the log scale for pathway analysis. Whole sample profiles are permuted across groups with fixed group sizes. Valid inference requires exchangeability under the permutation null.
- [Global Test](code/global_test.R) uses the intercept-only binomial score from [Goeman et al.](https://doi.org/10.1093/bioinformatics/btg382) and the shared-permutation implementation checked against the numerical score in Bioconductor `globaltest` 5.67.0. It uses pooled sample-SD autoscaling of the log intensities in CRC and METH and of the reported concentrations in T2DM. Duplicate T2DM records are represented by their geometric mean. The 999-permutation METH Global Test uses the archived seed and workbook order.
- [ORA](code/ora.R) screens feature records using two-sided Welch tests with BH q < 0.05, maps the selected features to unique compounds, and evaluates enrichment with the upper hypergeometric tail, equivalent to a one-sided Fisher exact test. It screens log abundances in CRC and METH and reported concentrations in T2DM. All 146 T2DM records are screened before mapping hits to unique compounds. The CRC and T2DM backgrounds contain 94 and 141 mapped compounds in the tested pathways, respectively. For METH, screening all 229 records recovers the archived 28 mapped hits. The historical METH output uses a 91-compound effective detected background. The bundled map verifies 90 identities, and the documented total of 91 is used for the enrichment denominator. The ORA function accepts an explicit background for other applications.

[Archived cohort results](reference/cohort_pathway_tests.csv) are preserved for comparison. The cohort runner performs fresh permutations with a specified seed. It provides a reproducible reference calculation, with Monte Carlo p and q values that can differ from the archived runs.

The bundled workflow was checked across all five contrasts and 172 pathway tests using 999 permutations. ORA p and q values in all five contrasts and Global Test p and q values for METH agree with the archived values at their reported precision.

The primary T2DM input keeps the reported fumarate concentrations. The manuscript's chromatogram and injection-order audit addresses sensitivity to quantification. This cohort workflow covers the primary pathway tests. The manuscript and supplementary materials describe the clinical adjustments, raw chromatogram audit, simulation designs and external-cohort comparisons.
