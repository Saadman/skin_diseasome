# Interactome separation (S_AB)

`sab_values.csv` is the canonical table of network-separation values for all fifteen
disease pairs, and is what Figure 3, Figure 4 and the Tier 2 drug-overlap comparison read.

| Column | Meaning |
|---|---|
| `from`, `to` | the disease pair |
| `Network_separation_s_AB` | S_AB; negative = modules closer than random expectation, positive = more separated |
| `Mean_shortest_distance_d_AB` | mean shortest distance between the two modules |
| `network_diameter_d_A`, `_d_B` | mean shortest distance within each module |

## Provenance

The values are taken from the per-pair `separation.py` outputs in
`separation_outputs_2022-08/`, which are the raw run logs and record the full-precision
numbers. Those runs used the differential-expression gene sets described in the Methods,
mapped to Entrez IDs; the exact input lists are archived in `../network/gene_lists/`.

An earlier pass computed the same statistic from a different set of inputs and is
superseded; its values are not used anywhere in this repository or in the manuscript.

## Reproducing it

```
python3 ../network/sab_fast.py ../network/interactome.tsv \
  PS ../network/gene_lists/PS.txt AD ../network/gene_lists/AD.txt \
  HS ../network/gene_lists/HS.txt VL ../network/gene_lists/VL.txt \
  CSU ../network/gene_lists/CSU.txt SLS ../network/gene_lists/SLS.txt
```

reproduces every value in `sab_values.csv` exactly at three decimal places. See
`../network/README.md` for the original Menche et al. implementation.
