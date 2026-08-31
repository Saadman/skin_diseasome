# A Multilayer Human Skin Diseasome

Data and analysis code for the manuscript *A Multilayer Human Skin Diseasome Links Molecular
Architecture, Clinical Co-occurrence, and Drug-Reversal Signatures*.

Six evidence layers — differential expression, GO enrichment, interactome separation (S_AB),
DisGeNET disease-gene enrichment, real-world co-occurrence in TriNetX, and LINCS L1000
drug-signature reversal — are compared across psoriasis (PS), atopic dermatitis (AD),
hidradenitis suppurativa (HS), vitiligo (VL), chronic spontaneous urticaria (CSU) and
sodium lauryl sulfate–induced irritant contact dermatitis (SLS-ICD).

Preprint: *[bioRxiv DOI to be added]*

## Repository layout

| Folder | Contents |
|---|---|
| `DEG/` | Per-disease differential-expression tables and the intersections derived from them |
| `go_enrichment/` | GO biological-process enrichment output |
| `network/` | Interactome, `separation.py`, and the gene lists used to compute S_AB |
| `networksim/` | S_AB values for all fifteen disease pairs, plus the raw run logs |
| `sensitivity/` | Harmonized equal-N gene sets and the threshold sensitivity analysis |
| `disgenet/` | DisGeNET enrichment and the diseasome edge list |
| `trinetx/` | Aggregate TriNetX co-occurrence counts |
| `signaturesearch/` | LINCS reversal rankings and the cross-disease comparisons |
| `scripts/` | All analysis and figure code |
| `images/high_res/` | Publication figures |
| `supplementary/` | Supplementary Data S1 and Tables S1–S3 |

## Figures

| Figure | Script | Main input |
|---|---|---|
| 1. DEG overlap (UpSet) | `Common_DEG_list.R` | `DEG/GEX_datasets/` |
| 2. GO biological process | `figure_go_dotplot.R` | `go_enrichment/` |
| 3. Interactome separation | `figure_sab_heatmap.R` | `networksim/sab_values.csv` |
| 4. Diseasome network | `figure_network_static.R` | `disgenet/`, `networksim/` |
| 5. Drug-repurposing overlap | `figure_drug_repurposing.R` | `signaturesearch/` |

All scripts are in `scripts/`. Current figures are `upset_plot`, `go_dotplot`, `sab_heatmap`,
`diseasome_network_static` and `drug_repurposing_tiers`; `images/high_res/` also holds
earlier versions of some panels.

## Reported numbers and where they come from

| Result | Source |
|---|---|
| DEG intersections (2,060 / 490 / 464 / 402 / 379 / 14) | `DEG/GEX_datasets/`, via `Common_DEG_list.R` |
| S_AB for all fifteen pairs | `networksim/sab_values.csv`, reproducible with `network/sab_fast.py` |
| DisGeNET counts at FDR < 0.25 (33 / 145 / 136 / 12 / 749 / 259) | `disgenet/skin_disease_enrichment_outerjoin_FDR_1030disease_131122.xlsx` |
| TriNetX co-occurrence (55 of 220; 42 diseases) | `trinetx/genetic_rwe_fdr_relation.xlsx` |
| LINCS reversal candidates per disease | `signaturesearch/*_lincs_result.csv` |
| Compound overlap (PS–CSU 19, PS–HS 2; rho = −0.22) | `signaturesearch/tier2_*.csv` |
| Comorbidity drug overlap (AD–liver carcinoma 25, CSU–HIV 11) | `signaturesearch/tier3_overlap_results.csv` |
| Threshold and set-size sensitivity | `sensitivity/`, summarized in Supplementary Data S1 |

## Running the code

Scripts run from the repository root and each reads files already present here, so any
stage can be rerun on its own:

```r
setwd("/path/to/skin_diseasome")
source("scripts/signaturesearch_tier2.R")
```

The one exception is `diseasome_cyto_filePrep.R`, which reads a table that `disgenet.R`
writes and so must follow it.

**Requirements.** R 4.4 with `readxl`, `dplyr`, `ggplot2`, `clusterProfiler`, `org.Hs.eg.db`,
`igraph`, `signatureSearch` and `signatureSearchData`, plus `UpSetR` and `ComplexUpset` for
Figure 1. Exact versions are in `SESSION_INFO.txt`. The network layer runs under Python 3
with `networkx` and `numpy`; see `network/README.md`.

**Two external resources** are needed only to regenerate their own stage, since the outputs
are checked in:

- **LINCS** — `signaturesearch_tier1.R` and `tier3.R` load the `lincs` dataset
  (ExperimentHub EH3226, ~2.3 GB), downloaded on first use.
- **DisGeNET** — `disgenet.R` needs an API key in `DISGENET_API_KEY`.

## Data sources

| Layer | Source |
|---|---|
| Transcriptomes | NCBI GEO: GSE30999 (PS), GSE32924 (AD), GSE148027 (HS), GSE75819 (VL), GSE57178 (CSU), GSE168735 (SLS-ICD) |
| Interactome | Menche et al., *Science* 2015 — redistributed in `network/`; please cite that paper |
| Perturbation signatures | LINCS L1000, Subramanian et al., *Cell* 2017 |
| Disease-gene associations | DisGeNET, queried November 2022 with `disgenet2r` (curated source, HGNC vocabulary) |
| Clinical co-occurrence | TriNetX — see the note below |

**TriNetX.** Individual-level records cannot be redistributed under the platform licence.
`trinetx/genetic_rwe_fdr_relation.xlsx` holds the aggregate counts reported in the
manuscript: one row per queried association, with the DisGeNET FDR that selected it and the
number of patients in the index cohort carrying the candidate diagnosis. Counts are rounded
by the platform and the smallest reportable non-zero value is ten, so an absent count cannot
be distinguished from a count below that threshold.

The 220 queried associations were taken from the DisGeNET enrichment at FDR < 0.05 for PS
(16), AD (23), HS (45), CSU (94) and SLS-ICD (29), and at FDR < 0.25 for VL (12), no VL
association having reached FDR < 0.05. Two composite DisGeNET concepts in the SLS-ICD set
were split into their five component terms so that each could be mapped to an ICD-10-CM
code, and two VL concepts were not carried forward.

## Licence

MIT — see `LICENSE`. Third-party inputs are not covered by it and remain under the terms of
their original sources: the interactome and `separation.py` (Menche et al. 2015), the GEO
series listed above, LINCS L1000 (Subramanian et al. 2017), and DisGeNET, whose API requires
a licence from the provider.

## Citation

*[Full citation and bioRxiv DOI to be added on posting.]*
