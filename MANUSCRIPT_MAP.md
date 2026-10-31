# Manuscript correspondence

Checked against the author-provided manuscript, Supplementary Table 1 and the completed server analysis on 2026-10-31. Source names below are relative to the study directory; source files were not modified.

## Included scripts

| Packaged file | Original source |
| --- | --- |
| `01_protein_bmi_mr.R` | `result/clump_radial_20261019_095058/code/pro_bmi_clump_radial_20261019.R` |
| `02_protein_t2dm_mr.R` | `result/clump_radial_20261019_095058/code/step2_clump_radial_20261019.R` |
| `03_select_candidates.R` | `result/clump_radial_20261019_095058/code/prepare_coloc_55_proteins.R` |
| `04_coloc_bmi.R` | `result/clump_radial_20261019_095058/code/coloc_bmi_clump_radial_20261019.R` |
| `05_coloc_t2dm.R` | `result/clump_radial_20261019_095058/code/coloc_t2dm_clump_radial_20261019.R` |
| `06_multitrait_coloc.R` | `result/clump_radial_20261019_095058/code/moloc_clump_radial_20261019.R` |
| `moloc_typed.R` | `result/clump_radial_20261019_095058/code/moloc_test_typed_20261020.R` |
| `07_finemap_adh1b.R` | `fine_mapping_strict_cap_20261019_113316/code/ADH1B_strict_numcausal5.R` |
| `08_finemap_entpd6.R` | `fine_mapping_strict_cap_20261019_113316/code/ENTPD6_strict_numcausal5.R` |
| `finemap_helpers.R` | `fine_mapping_strict_cap_20261019_113316/code/helpers.R` |
| `CARMA_strict.R` | `fine_mapping_strict_cap_20261019_113316/code/CARMA_strict.R` |
| `09_export_cs95.R` | `fine_mapping_strict_cap_20261019_113316/code/export_supplementary_cs95.R` |

## Verified correspondence and retained details

- The completed MR summary reports 1,944 BMI analyses, 387 selected BMI-associated proteins, 373 T2DM analyses, 81 selected T2DM-associated proteins and 55 concordant candidates, matching the manuscript. This was checked from existing results, without rerunning MR.
- Both MR stages preserve P < 5e-8 selection, P=0 tie ordering by absolute beta/SE, harmonization, Steiger calculation, and one attempted Radial IVW screen only when the IVW Q-test P < 0.05. Outlier alpha is 0.05; failure/no usable outliers preserves the original instruments. The saved `mr_keep` flag is used; the scripts do not discard rows based on `steiger_dir`.
- `ld_clump()` retains its original implicit defaults. The inspected ieugwasr 1.1.0 defaults are r²=0.001 and 10,000 kb, matching the manuscript; other versions must be checked before reuse.
- F statistics are not an additional explicit filter in these two source MR scripts. No new F-statistic calculation/filter has been inserted to match the wording of the manuscript. Instrument preparation upstream and the reported F statistics must be considered separately.
- Candidate selection preserves stage-specific BH/FDR <0.05, residual IVW Q and Egger-intercept filters where available, and matching directions of effect. Its original expected-set validation uses the accompanying 55 protein identifiers.
- The pairwise result files each contain 55 proteins. PP.H4 >0.80 identifies 9 BMI and 11 T2DM proteins, with 7 overlapping, matching the manuscript.
- Multi-trait priors remain 1e-4, 1e-5 and 1e-6; explicit trait types remain quantitative, quantitative and case-control. The source additionally requires at least 50 aligned shared variants.
- The multi-trait script's auxiliary significant-results export uses PP >0.70, while the manuscript uses PP >0.80. This difference is retained. Across all 55 proteins, these thresholds yield 6 and 5 proteins, respectively; intersecting with the 7 pairwise-supported proteins yields the same 5 at either threshold: ADH1B, AP3B1, APOE, ENTPD6 and RBP1. Use the complete posterior table when assessing the manuscript's threshold.
- Fine-mapping retains the original strict five-variant implementation, outlier detection, allele/frequency alignment, seeds and ±250-kb windows. The native CARMA fits use 99% credible sets. The separate original export script extracts 95% sets from those saved fits for the supplementary tables; the native 99% sets must not be relabelled as 95% sets.

## Release scope

Fenland replication, mediation and single-cell analyses are excluded at the author's request. Linear/Cox regression code is pending. The inspected FinnGen script reads an obesity endpoint despite the author confirming continuous BMI_IRN, and the available network/PPI/enrichment scripts include legacy or incomplete workflows; they are not substituted into this release. PheWAS portal queries are not represented by a verified final script.

No participant-level records, association datasets, private machine paths, analysis logs or manuscript files are included. The packaged code has been parsed and compared with the original executable expressions after path/library normalization; this is not an end-to-end rerun or a guarantee that every manuscript method is covered.
