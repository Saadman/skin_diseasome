# Mixed-threshold sensitivity analysis

Response to the editorial comment that mixed statistical thresholds threaten validity, because
DEG-set size drives overlap, enrichment, module size and drug hits.

Generated 2026-08-15. Scripts: `scripts/harmonize_build_sets.R`, `scripts/harmonize_go.R`,
`scripts/sensitivity_equalN.R`, `scripts/sensitivity_equalN_v2.R`. S_AB computed with the
published Menche et al. 2015 `separation.py` (Python-3 converted; interactome unchanged,
13,460 nodes / 141,296 links).

## 1. Can a single threshold be used for all six diseases?

No. At |log2FC| > 0.585, unique gene symbols surviving each rule:

| Disease | Published list | Nominal p < 0.05 | FDR < 0.05 |
|---|---|---|---|
| PS | 6,522 | 6,423 | 6,410 |
| AD | 1,527 | 1,507 | **0** |
| HS | 6,801 | 6,875 | 6,783 |
| VL | 3,040 | 3,741 | **3,040** |
| CSU | 976 | 983 | 819 |
| SLS | 33 | 32 | **0** |

A uniform FDR threshold is arithmetically impossible: it eliminates every gene in AD and SLS.
VL matches FDR exactly (3,040), so VL cannot have used nominal p. The published per-disease
criteria are therefore the only ones consistent with the reported gene lists. The fold-change
criterion was uniform across all six.

## 2. Harmonization procedure

Identical rule applied to every disease, removing both the threshold and the size confound:

1. restrict to genes measurable on all platforms (10,340 genes for six diseases; 11,037 for five)
2. nominal p < 0.05
3. rank surviving genes by |log2FC| descending, take top N

Ranking uses effect size rather than p-value because p-value depends on sample size (PS n = 170,
CSU n = 18), which would systematically favour the large-cohort diseases. A significance filter is
applied first, because ranking on effect size alone selects noise — under that variant SLS's
separation disappeared entirely, confirming the filter is load-bearing.

**N is capped at 582** by SLS, the most it can supply at p < 0.05. Excluding SLS, the five clinical
diseases support N up to 2,720. Both are reported.

## 3. Results

### Pairwise overlap (`harmonized_pairwise_overlap.csv`)

Six diseases, genuinely equal N:

| N | Rank 1 | PS–HS rank | SLS mean overlap | Non-SLS mean |
|---|---|---|---|---|
| 250 | PS–AD (71) | 2nd (26) | 10.0 | 21.3 |
| 400 | PS–AD (130) | 2nd (53) | 20.8 | 44.0 |
| 582 | PS–AD (202) | 2nd (99) | 35.0 | 74.3 |

Five clinical diseases, larger N — PS–HS converges on PS–AD as depth increases:

| N | PS–AD | PS–HS | gap |
|---|---|---|---|
| 1,000 | 361 | 242 | 33% |
| 2,000 | 797 | 672 | 16% |
| 2,720 | 1,101 | 1,086 | **1.4%** |

On the published (unequal) sets, size-normalised overlap still ranks PS–HS clearly first:
Jaccard 0.329 vs 0.176 for PS–AD.

### Interactome separation (`harmonized_sab.csv`, `harmonized5_sab.csv`)

Five clinical diseases:

| Pair | Published | N=1,000 | N=2,720 |
|---|---|---|---|
| PS–AD | −0.267 | **−0.573** | **−0.538** |
| PS–HS | **−0.594** | −0.356 | −0.519 |
| AD–HS | −0.148 | −0.260 | −0.401 |
| all 10 pairs negative | yes | **yes** | **yes** |

Six diseases at N ≤ 582: SLS pairs are **not** consistently positive — 4 of 5 are negative at
N = 582 (PS–SLS −0.129, AD–SLS −0.063, HS–SLS −0.043, CSU–SLS −0.034; only VL–SLS positive).

### GO biological process (`harmonized_GO_N582.csv`, `harmonized5_GO_N2720.csv`)

| Disease | Top harmonized term | Neutrophil terms |
|---|---|---|
| PS | keratinization | none |
| AD | monocarboxylic acid metabolic process | none |
| HS | regulation of cell activation | migration, rank 72 |
| VL | ribonucleoprotein complex biogenesis | none |
| CSU | ribonucleoprotein complex biogenesis | activation, rank 95 |
| SLS | *no significant terms* | none |

## 4. What survives, what does not

**Robust to harmonization**

- All ten pairings among PS, AD, HS, VL and CSU have negative S_AB — reproduced at both N.
- PS–HS and PS–AD are the two closest pairs by both overlap and S_AB, at every N tested.
- VL's distinctive ribosome-biogenesis / RNA-metabolism signature.

**Not robust**

- PS–HS as uniquely the closest pair. PS–AD ties or leads under every harmonization; the two are
  within 1.4% (overlap) and 0.019 (S_AB) at the deepest level tested.
- SLS as separated from all other diseases. Under equal-sized sets its S_AB values go negative.
- Neutrophil activation as a leading shared pathway.

**Interpretation**

The SLS result cannot be separated from the weakness of its dataset: SLS yields zero genes at
FDR < 0.05 and only 582 at p < 0.05. Its published 33-gene module is small enough that d_AA is
inflated, which is what produces the positive S_AB. Whether SLS is biologically distinct or simply
underpowered is not decidable from this dataset.

The GO shift is substantially explained by the ranking statistic. Ranking on |log2FC| selects the
largest-fold-change genes, which in psoriasis are keratins and S100 proteins, whereas
neutrophil/chemokine genes have moderate fold changes but high significance. This is a property of
the harmonization, not evidence that the original enrichment was wrong. A rank-based method
without cutoffs (GSEA preranked) is the appropriate test for the pathway claims and has not yet
been run.

## 5. Files

| File | Contents |
|---|---|
| `harmonized_pairwise_overlap.csv` | overlap + Jaccard, 15 pairs × N ∈ {250, 400, 582} |
| `harmonized_sab.csv` | S_AB, six diseases, 15 pairs × N ∈ {250, 400, 582} |
| `harmonized5_sab.csv` | S_AB, five clinical diseases, 10 pairs × N ∈ {1000, 2720} |
| `harmonized_GO_N582.csv` | GO BP enrichment, six diseases |
| `harmonized5_GO_N2720.csv` | GO BP enrichment, five clinical diseases |
| `harmonized_symbols_N*.rds`, `harmonized5_symbols_N*.rds` | the harmonized gene sets |
| `*_N*.txt`, `*_5d_N*.txt` | Entrez ID lists used as separation.py input |
