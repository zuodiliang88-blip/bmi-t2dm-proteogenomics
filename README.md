# Analysis templates

| Code | Purpose |
| --- | --- |
| `01_protein_bmi_mr.R` | Protein → BMI MR with LD clumping and conditional Radial IVW. |
| `02_protein_t2dm_mr.R` | Protein → T2DM MR for BMI-associated proteins. |
| `04_coloc_bmi.R` | Protein/BMI colocalization; quantitative/quantitative traits. |
| `05_coloc_t2dm.R` | Protein/T2DM colocalization; quantitative/case-control traits. |
| `06_multitrait_coloc.R` | Protein/BMI/T2DM colocalization; helper included. |
| `07_finemap_adh1b.R` | ADH1B fine-mapping using the original CARMA package. |
| `08_cox.R` | Three-model Cox regression using prepared two-year landmark data. |
| `09_linear.R` | BMI linear regression; Models 1–3 match the manuscript, Model4 is additional. |

Set BMI_T2DM_ROOT to the absolute project directory and adjust the input/output paths. Prepare the processed data and output directories before running. PLINK must be on PATH. Candidate-selection instructions are in the colocalization templates; no separate selection script is included. See the English comments in each file for input requirements. These are analysis templates, not an end-to-end pipeline.

