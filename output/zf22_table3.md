# ZF22 Table 3 -- D2DGA reproduction (M6-T1)

Mesh `n_phi = 20`, `n_xi = 200`, CFL = 0.5.

`t_br` is reported at three outlet-concentration thresholds; see the
module docstring and assumptions.md NUM-21 for why one number is not
enough. `1-D shock` is the arrival time of the main shock of BF25
(3.7)'s planar reduction, which isolates the gap-scale closures from
the 2-D solve.

| case | e | m | b (ZF22) | b (BF25) | t_br @0.01 | @0.1 | @0.5 | ZF22 t_br | 1-D shock | eta_E | ZF22 eta_E |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.8 | 0.2 | -50 | -22.4 | 0.053 | 0.054 | 0.059 | 0.44 | 7.332 | 0.375 | 0.66 |
| 2 | 0.8 | 0.2 | 100 | 44.7 | 0.877 | 0.924 | 0.929 | 0.95 | 0.950 | 0.952 | 0.95 |
| 3 | 0.6 | 0.2 | 10 | 4.5 | 0.425 | 0.736 | 0.841 | 0.79 | 0.859 | 0.897 | 0.92 |
| 4 | 0.6 | 0.2 | 1000 | 447.2 | 0.948 | 0.958 | 0.958 | 0.99 | 0.982 | 0.983 | 1.0 |
| 5 | 0.4 | 0.2 | 100 | 44.7 | 0.804 | 0.928 | 0.936 | 0.95 | 0.950 | 0.957 | 0.97 |
| 6 | 0.4 | 5.0 | 100 | 223.6 | 0.817 | 0.913 | 0.917 | 0.93 | 0.935 | 0.940 | 0.93 |
| 7 | 0.2 | 0.5 | 10 | 7.1 | 0.500 | 0.737 | 0.859 | 0.78 | 0.856 | 0.898 | 0.9 |
| 8 | 0.2 | 2.0 | 10 | 14.1 | 0.499 | 0.746 | 0.820 | 0.78 | 0.823 | 0.858 | 0.84 |
| 9 | 0.1 | 0.2 | 100 | 44.7 | 0.813 | 0.930 | 0.938 | 0.97 | 0.950 | 0.959 | 0.97 |
| 10 | 0.1 | 0.2 | 100 | 44.7 | 0.813 | 0.930 | 0.938 | 0.97 | 0.950 | 0.959 | 0.97 |

Volume conservation (max relative drift over the run):

* case 1: 6.62e-15
* case 2: 2.18e-15
* case 3: 1.50e-15
* case 4: 7.06e-15
* case 5: 2.71e-15
* case 6: 6.66e-15
* case 7: 8.74e-16
* case 8: 1.09e-15
* case 9: 1.58e-15
* case 10: 1.58e-15
