# Interactome and network-separation code

Everything required to recompute the network-separation statistic S_AB. The resulting
values live in `../networksim/`.

## Contents

| File | What it is |
|---|---|
| `interactome.tsv` | Human protein–protein interaction network from Menche et al., *Science* 2015 (13,460 nodes, 141,296 links), used unmodified |
| `separation.py` | The authors' original implementation, as published with that paper (Python 2) |
| `separation_py3.py` | The same code converted to run under Python 3 (`print` statements, `dict.has_key`, `Graph.selfloop_edges`); the algorithm is unchanged |
| `sab_fast.py` | A faster reimplementation of the same statistic, computing all fifteen pairs in one pass |
| `gene_lists/*.txt` | The per-disease Entrez ID lists used as input, as archived from the 2022 analysis |
| `sab_recomputed.csv` | Output of `sab_fast.py` on the files in this folder |
| `menche_original_readme.txt` | The original readme distributed with `separation.py` |

`gene_lists/` holds the differential-expression gene sets described in the Methods, mapped
from symbols to Entrez IDs; they were written by
`../scripts/go_enrichment.R` from
`Venn_allDisease_final.xlsx`.

## Usage

```
python3 sab_fast.py interactome.tsv \
  PS gene_lists/PS.txt AD gene_lists/AD.txt HS gene_lists/HS.txt \
  VL gene_lists/VL.txt CSU gene_lists/CSU.txt SLS gene_lists/SLS.txt
```

or, for a single pair with the original code:

```
python3 separation_py3.py -n interactome.tsv --g1 gene_lists/PS.txt --g2 gene_lists/HS.txt -o ps_hs.txt
```

`separation_py3.py` requires `networkx` and `numpy` (tested with networkx 3.6.1,
numpy 2.5.2, Python 3.13); `sab_fast.py` needs only the standard library.

The harmonized, equal-sized gene sets used for the sensitivity analysis are in
`../sensitivity/` (`*_N*.txt` for the six-disease analysis, `*_5d_N*.txt` for the five
clinical diseases) and take the same input format.

## Verification

`separation_py3.py` loads the network as 13,460 nodes and 141,296 links, as reported.
`sab_fast.py` reproduces all fifteen values in `../networksim/sab_values.csv` exactly at
three decimal places, and those values match the archived `separation.py` run logs in
`../networksim/separation_outputs_2022-08/`.

## Citation

Menche J, Sharma A, Kitsak M, et al. Uncovering disease-disease relationships through the
incomplete interactome. *Science*. 2015;347(6224):1257601.

`interactome.tsv` and `separation.py` are redistributed here as the supplementary material
of that publication so that the analysis can be rerun as performed; please cite the original
paper for both.
