# ZF22 Table 3 -- D2DGA reproduction (M6-T1)

Mesh `n_phi = 20`, `n_xi = 200`, CFL = 0.5.

`t_br` is reported at three outlet-concentration thresholds; see the
module docstring and assumptions.md NUM-21 for why one number is not
enough. `1-D shock` is the arrival time of the main shock of BF25
(3.7)'s planar reduction, which isolates the gap-scale closures from
the 2-D solve. BOTH the leading wave and the main shock are quoted:
for an adverse density difference (case 1) the largest shock travels
BACKWARDS, so its arrival time says nothing about breakthrough and only
the leading wave does.

| case | e | m | b (ZF22) | b (BF25) | t_br @0.01 | @0.1 | @0.5 | ZF22 t_br | 1-D lead | 1-D main shock | eta_E | ZF22 eta_E |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.8 | 0.2 | -50 | -22.4 | 0.054 | 0.054 | 0.060 | 0.44 | 0.391 | 7.332 | 0.366 | 0.66 |
| 2 | 0.8 | 0.2 | 100 | 44.7 | 0.880 | 0.927 | 0.932 | 0.95 | 0.667 | 0.950 | 0.951 | 0.95 |
| 3 | 0.6 | 0.2 | 10 | 4.5 | 0.426 | 0.737 | 0.842 | 0.79 | 0.667 | 0.859 | 0.897 | 0.92 |
| 4 | 0.6 | 0.2 | 1000 | 447.2 | 0.953 | 0.963 | 0.963 | 0.99 | 0.669 | 0.982 | 0.983 | 1.0 |
| 5 | 0.4 | 0.2 | 100 | 44.7 | 0.806 | 0.931 | 0.938 | 0.95 | 0.667 | 0.950 | 0.957 | 0.97 |
| 6 | 0.4 | 5.0 | 100 | 223.6 | 0.822 | 0.918 | 0.923 | 0.93 | 0.667 | 0.935 | 0.940 | 0.93 |
| 7 | 0.2 | 0.5 | 10 | 7.1 | 0.501 | 0.738 | 0.860 | 0.78 | 0.667 | 0.856 | 0.897 | 0.9 |
| 8 | 0.2 | 2.0 | 10 | 14.1 | 0.500 | 0.748 | 0.822 | 0.78 | 0.667 | 0.823 | 0.857 | 0.84 |
| 9 | 0.1 | 0.2 | 100 | 44.7 | 0.815 | 0.932 | 0.941 | 0.97 | 0.667 | 0.950 | 0.959 | 0.97 |
| 10 | 0.1 | 0.2 | 100 | 44.7 | 0.815 | 0.932 | 0.941 | 0.97 | 0.667 | 0.950 | 0.959 | 0.97 |

Volume conservation (max relative drift over the run):

* case 1: 1.07e-14
* case 2: 3.03e-15
* case 3: 1.28e-15
* case 4: 4.15e-15
* case 5: 2.03e-15
* case 6: 4.23e-15
* case 7: 9.26e-16
* case 8: 7.57e-16
* case 9: 1.21e-15
* case 10: 1.21e-15
