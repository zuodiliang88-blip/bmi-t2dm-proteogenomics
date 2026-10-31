# BMI and type 2 diabetes proteogenomics

Analysis code accompanying *Proteogenomic analysis identifies plasma proteins associated with body mass index and type 2 diabetes risk*.

This release contains the manuscript-matched two-stage protein MR, candidate selection, colocalization and strict five-causal-variant fine-mapping analyses. Computation and statistical parameters are retained from the study scripts. Packaging changes are limited to paths, script names, environment-specific library setup, comments and whitespace.

## Analysis order

| Script | Analysis |
| --- | --- |
| `01_protein_bmi_mr.R` | UKB-PPP cis-pQTL instruments → MVP BMI; LD clumping and conditional Radial IVW |
| `02_protein_t2dm_mr.R` | BMI-associated proteins → European-ancestry T2DGGI T2DM |
| `03_select_candidates.R` | FDR, heterogeneity, Egger-intercept and concordant-direction filters |
| `04_coloc_bmi.R`, `05_coloc_t2dm.R` | Pairwise colocalization in ±500-kb windows |
| `06_multitrait_coloc.R`, `moloc_typed.R` | Protein/BMI/T2DM colocalization; explicit quantitative/quantitative/case-control types |
| `07_finemap_adh1b.R`, `08_finemap_entpd6.R` | Fine-mapping at rs1229984 and rs6050446, ±250 kb, maximum five causal variants |
| `CARMA_strict.R`, `finemap_helpers.R` | Original study-specific CARMA modification, LD alignment and model-space checks |
| `09_export_cs95.R` | 95% credible-set extraction from saved fits without rerunning CARMA |

Fenland replication, mediation and single-cell scripts are excluded from this release. Observational linear/Cox regression code is pending. FinnGen replication, PheWAS and PPI/enrichment are not included because a final manuscript-matched computational workflow was not established in the inspected scripts. This repository therefore does not yet reproduce every analysis in the manuscript.

## Inputs

Place the following **processed inputs** below the repository root, or symlink them there. Original GWAS downloads require column/variant/build preprocessing before these analysis scripts. No raw or participant-level data are distributed.

```text
data/
  gwas/
    BMI_GCST90475156.tsv.gz
    t2dm_EUR.txt.gz
  pqtl/
    protein_file_list/      # original cis1mb_after_clump filenames (*.txt.gz)
    ukb_cis_1mb/            # corresponding uncompressed *.txt pQTL files
  ld/
    EUR.bed
    EUR.bim
    EUR.fam
```

- GWAS columns: `SNP`, `chr`, `bp`, `A1`, `A2`, `freq`, `beta`, `se`, `p`, `n`.
- pQTL columns: `SNP`, `Chr`, `Pos`, `Effect_allele`, `Other_allele`, `Eaf`, `Beta`, `Se`, `P`, `N`, `Protein`.
- GWAS `A1` and pQTL `Effect_allele` identify the effect allele. Preserve each source's genome build and use the same preprocessed coordinate conventions as the study; do not simply relabel genomic coordinates.
- The MR starting protein universe is determined by the filenames in `protein_file_list/`. Each `.gz` filename must correspond to a file without the final `.gz` in `ukb_cis_1mb/`. The list is not regenerated from every available pQTL file.
- `expected_dual_role_candidates.csv` contains only the 55 protein file identifiers from the completed study run. It is used by the original selection script for a set-equality check, not for effect estimation or selection.

Original data: [MVP BMI, GCST90475156](https://www.ebi.ac.uk/gwas/studies/GCST90475156), [T2DGGI European-ancestry GWAS](https://www.diagram-consortium.org/downloads.html), [UKB-PPP](https://www.synapse.org/Synapse:syn51365301), and [1000 Genomes](https://www.internationalgenome.org/data/). Access and reuse follow the providers' terms.

## Requirements and execution

Use R with `TwoSampleMR`, `MRInstruments`, `ieugwasr`, `plinkbinr`, `dplyr`, `data.table`, `RadialMR`, `coloc`, `moloc`, `CARMA`, `Matrix` and `jsonlite`; install PLINK 1.9 on `PATH`. CARMA and its compiled dependencies require a compatible compiler toolchain. Packages must already be installed; these scripts do not install software.

Versions observed on the source server during packaging: R 4.5.3, TwoSampleMR 0.6.8, ieugwasr 1.1.0, RadialMR 1.1, coloc 5.2.3, moloc 0.1.0 and CARMA 1.0. This is an environment observation, not a historical lockfile. In particular, `moloc_typed.R` intentionally patches exactly one expression in the installed moloc function and stops if the expected expression is absent. The manuscript currently reports R 4.5.1 for single-cell analyses, which are outside this release.

Run from the repository root after preparing the inputs:

```bash
bash run_analysis.sh
```

For an individual stage, first set `export BMI_T2DM_ROOT="$PWD"` and create the output directories used in `run_analysis.sh`. This absolute root is required because the fine-mapping scripts change working directory. PLINK is resolved from `PATH`.

The complete analyses were **not rerun during packaging**. All 12 packaged R scripts passed syntax and parsed-expression equivalence checks against the source snapshot, after normalizing only the documented path and library-setup changes. See [MANUSCRIPT_MAP.md](MANUSCRIPT_MAP.md) for result correspondence and retained differences.

## Third-party code

The study's modified CARMA source is retained in `CARMA_strict.R`; replacing it with the unmodified package would change the analysis. Its original authors, licensing and modification notice are recorded in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) and [LICENSE-CARMA](LICENSE-CARMA). A blanket licence for the other author scripts has not been added.
