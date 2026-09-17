# Phase 1 -- K-GEP-1 synthetic-wall A-L sweep

Newtonian pair, `m = 0.5`, `b = 10.0`; wall mode `enlargement`; mesh `20 x 400`; CFL 0.5; run to 1.2 pumped volumes.

`delta/pi` and `|dr_o/dxi|` are the Hele-Shaw validity diagnostics and
come FIRST on purpose: a case that displaces well while failing them is
a warning, not a result. `t_br` is quoted at three outlet thresholds
(assumptions.md NUM-21).  The ZF23 dispersion metrics are evaluated while
the tip is still at 80% of the domain: after breakthrough the axial
profile is truncated and sigma_w+r collapses to zero, and sampling right
at the outlet clips the w_r+ tail exactly where it is largest.

| A (in) | L (m) | delta/pi | \|dr/dxi\| | e range | eta_E | t_br@0.01 | @0.1 | @0.5 | narrow min | residual | sigma_w+r | dispersive |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.0 | 5.0 | 0.0622 | 0.0000 | 0.30–0.30 | 0.9045 | 0.491 | 0.827 | 0.871 | 0.642 | 0.0000 | 0.2047 | True |
| 0.0 | 20.0 | 0.0622 | 0.0000 | 0.30–0.30 | 0.9045 | 0.491 | 0.827 | 0.871 | 0.642 | 0.0000 | 0.2047 | True |
| 0.0 | 60.0 | 0.0622 | 0.0000 | 0.30–0.30 | 0.9045 | 0.491 | 0.827 | 0.871 | 0.642 | 0.0000 | 0.2047 | True |
| 1.5 | 5.0 | 0.0825 | 0.0120 | 0.21–0.30 | 0.9067 | 0.558 | 0.851 | 0.880 | 0.740 | 0.0000 | 0.1886 | True |
| 1.5 | 20.0 | 0.0825 | 0.0030 | 0.21–0.30 | 0.9071 | 0.533 | 0.856 | 0.881 | 0.696 | 0.0000 | 0.1902 | True |
| 1.5 | 60.0 | 0.0825 | 0.0010 | 0.21–0.30 | 0.9074 | 0.526 | 0.831 | 0.880 | 0.683 | 0.0000 | 0.1875 | True |
| 3.0 | 5.0 | 0.0999 | 0.0239 | 0.16–0.30 | 0.9108 | 0.616 | 0.878 | 0.901 | 0.777 | 0.0000 | 0.1581 | True |
| 3.0 | 20.0 | 0.0999 | 0.0060 | 0.16–0.30 | 0.9104 | 0.583 | 0.882 | 0.901 | 0.725 | 0.0000 | 0.1735 | True |
| 3.0 | 60.0 | 0.0999 | 0.0020 | 0.16–0.30 | 0.9109 | 0.573 | 0.855 | 0.899 | 0.701 | 0.0000 | 0.1429 | True |
| 6.0 | 5.0 | 0.1279 | 0.0479 | 0.11–0.30 | 0.9165 | 0.704 | 0.932 | 0.951 | 0.804 | 0.0000 | 0.0975 | True |
| 6.0 | 20.0 | 0.1279 | 0.0120 | 0.11–0.30 | 0.9158 | 0.678 | 0.934 | 0.950 | 0.758 | 0.0000 | 0.1131 | True |
| 6.0 | 60.0 | 0.1279 | 0.0040 | 0.11–0.30 | 0.9167 | 0.671 | 0.910 | 0.947 | 0.728 | 0.0000 | 0.0630 | False |
