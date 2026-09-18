"""
The figures, drawn from a run's own arrays.

No physics here.  Every derived quantity comes from `d2dga.postprocess` or
`d2dga.geometry`; this module only arranges them.
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402

from d2dga.postprocess import axial_profile, narrow_side_profile   # noqa: E402
from ui import theme                     # noqa: E402


def _depth_edges(geometry):
    """Depth (m below surface) at the xi CELL EDGES, deepest first.

    xi is measured UP from bottom hole, so depth = total_depth - xi_hat and the
    array runs from TD down to the casing shoe.
    """
    g = geometry.grid
    return geometry.well.total_depth_m - geometry.xi_hat(g.xi_edges)


def field_with_envelope(geometry, c, figsize=(5.4, 6.6), stretch=False):
    """
    The mandated pairing (build spec section 3).

    Left: the unwrapped annulus -- azimuth across, WIDE SIDE LEFT, narrow side
    right, depth down -- which is how every source paper plots it.

    Right: delta/pi on the SAME depth axis, aligned row for row, coloured
    against ZF23's validated 0.038.  Shared y, one figure, so the alignment is
    structural rather than eyeballed: the reader must be able to see where the
    front is and where the model is untrustworthy in one glance.
    """
    with plt.rc_context(theme.mpl_rc()):
        fig, (ax, axe) = plt.subplots(
            1, 2, figsize=figsize, sharey=True,
            gridspec_kw={"width_ratios": [3.1, 1.0], "wspace": 0.06})

        depth = _depth_edges(geometry)
        d_top, d_bot = float(np.min(depth)), float(np.max(depth))

        # c is (n_phi, n_xi); transpose so rows are depth and columns azimuth.
        # The absolute 0-1 ramp is the default and the honest one.  A nearly
        # complete displacement is then a flat block -- true, and useless as a
        # figure -- so `stretch` rescales to the data range.  The colour bar
        # always states which range is in use, so a stretched field cannot be
        # mistaken for a poor displacement.
        lo, hi = (float(c.min()), float(c.max())) if stretch else (0.0, 1.0)
        if hi - lo < 1e-9:
            lo, hi = 0.0, 1.0
        im = ax.pcolormesh(np.linspace(0.0, 1.0, c.shape[0] + 1), depth, c.T,
                           cmap=theme.CMAP, vmin=lo, vmax=hi,
                           shading="flat", rasterized=True)
        ax.set_xlabel("azimuth  $\\varphi$   wide $\\rightarrow$ narrow")
        ax.set_ylabel("depth below surface (m)")
        ax.set_ylim(d_bot, d_top)          # depth increases downward
        ax.invert_yaxis()
        ax.set_xticks([0.0, 0.5, 1.0])
        ax.set_title("gap-averaged concentration $\\bar c$", loc="left")

        cb = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.09,
                          fraction=0.045, ticks=[lo, 0.5 * (lo + hi), hi])
        cb.outline.set_edgecolor(theme.RULE)
        cb.ax.set_xticklabels([f"{lo:.4g}", f"{0.5 * (lo + hi):.4g}", f"{hi:.4g}"])
        cb.set_label("mud $\\leftarrow$  $\\bar c$  $\\rightarrow$ cement"
                     + ("   (scale stretched to this run)" if stretch else
                        "   (absolute scale)"), fontsize=8)

        xi_c = geometry.grid.xi_centres
        dpi = np.asarray(geometry.narrow_gap_parameter(xi_c), dtype=float)
        depth_c = geometry.well.total_depth_m - geometry.xi_hat(xi_c)
        colours = np.where(dpi > theme.DELTA_PI_VALIDATED, theme.AMBER, theme.TEAL)
        height = abs(depth[1] - depth[0]) * 0.92
        axe.barh(depth_c, dpi, height=height, color=colours, linewidth=0)
        axe.axvline(theme.DELTA_PI_VALIDATED, color=theme.TEAL_DARK, lw=1.0)
        axe.text(theme.DELTA_PI_VALIDATED, d_top, " 0.038", fontsize=7,
                 color=theme.TEAL_DARK, va="top", ha="left")
        axe.set_xlim(0.0, max(0.05, float(dpi.max()) * 1.18))
        axe.set_xlabel("$\\delta/\\pi$")
        frac = float(np.mean(dpi > theme.DELTA_PI_VALIDATED))
        axe.set_title(f"{frac:.0%} outside", loc="left",
                      color=theme.AMBER_INK if frac else theme.TEAL_DARK)
        for s in ("top", "right"):
            axe.spines[s].set_visible(False)
    return fig


def history(times, efficiency, outlet, Z, t_br=None, figsize=(6.4, 2.5)):
    """Efficiency and outlet concentration against pumped volumes."""
    with plt.rc_context(theme.mpl_rc()):
        fig, ax = plt.subplots(figsize=figsize)
        vol = np.asarray(times, dtype=float) / Z
        ax.plot(vol, efficiency, color=theme.TEAL, lw=1.6, label="$\\eta_E$")
        ax.plot(vol, outlet, color=theme.AMBER, lw=1.4, ls="--",
                label="outlet $\\bar c$")
        if t_br is not None and np.isfinite(t_br):
            ax.axvline(t_br / Z, color=theme.MUTED_2, lw=0.9, ls=":")
            ax.text(t_br / Z, 1.02, f" $t_{{br}}$ {t_br / Z:.3f}", fontsize=7,
                    color=theme.MUTED, va="bottom")
        ax.set_xlabel("pumped volumes  $t/Z$")
        ax.set_ylim(0.0, 1.05)
        ax.set_xlim(0.0, float(vol.max()))
        ax.legend(loc="center right")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    return fig


def profiles(geometry, c, figsize=(6.4, 2.5)):
    """Volume-weighted azimuthal mean, and the narrow-side column, against depth."""
    with plt.rc_context(theme.mpl_rc()):
        fig, ax = plt.subplots(figsize=figsize)
        xi_c = geometry.grid.xi_centres
        depth = geometry.well.total_depth_m - geometry.xi_hat(xi_c)
        ax.plot(depth, axial_profile(geometry, c), color=theme.TEAL, lw=1.6,
                label="volume-weighted mean")
        ax.plot(depth, narrow_side_profile(geometry, c), color=theme.AMBER,
                lw=1.4, ls="--", label="narrow side  $\\varphi = 1$")
        ax.set_xlabel("depth below surface (m)")
        ax.set_ylabel("$\\bar c$")
        ax.set_ylim(-0.02, 1.05)
        ax.invert_xaxis()
        ax.legend(loc="lower right")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    return fig
