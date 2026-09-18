"""
The visual system, defined once.

Build spec section 6: warm neutral ground, single teal accent, monospace for all
numerics, serif for headings, and a concentration ramp running tan (mud,
c_bar = 0) to teal (cement, c_bar = 1) that is **defined here and shared by
every plot and legend**.  Nothing in `ui/` may declare a second ramp or a
second accent.
"""

from matplotlib.colors import LinearSegmentedColormap

# -- ground and ink --------------------------------------------------------
GROUND = "#F5F4F0"
PANEL = "#FFFFFF"
BORDER = "#DDDCD5"
RULE = "#E7E5DD"
INK = "#1A1A17"
INK_2 = "#3B3B35"
MUTED = "#5C5C55"
MUTED_2 = "#6B6B63"

# -- the single accent -----------------------------------------------------
TEAL = "#0E6E66"
TEAL_DARK = "#0A5049"
TEAL_WASH = "#EAF2F0"
TEAL_EDGE = "#C3DCD7"

# -- warning, for things outside their validated range ---------------------
AMBER = "#B5742A"
AMBER_INK = "#8A5310"
AMBER_WASH = "#FDF8F0"
AMBER_EDGE = "#E8D3AC"

# -- concentration ---------------------------------------------------------
MUD = "#D9C9A3"          # c_bar = 0, displaced fluid
CEMENT = TEAL            # c_bar = 1, displacing fluid
CMAP = LinearSegmentedColormap.from_list("d2dga_c", [MUD, "#7EA694", CEMENT])

# -- ZF23's validated narrow-gap parameter ---------------------------------
# Frigaard's validated experiments sit at delta/pi ~ 0.038 (ZF23 section 2.1).
# Everything above it is a prediction outside the demonstrated envelope, which
# is why the envelope is plotted beside the field rather than in a settings tab.
DELTA_PI_VALIDATED = 0.038


def mpl_rc() -> dict:
    """Matplotlib settings that make a figure match the interface."""
    return {
        "figure.facecolor": PANEL,
        "axes.facecolor": "#FBFAF7",
        "axes.edgecolor": RULE,
        "axes.labelcolor": INK_2,
        "axes.titlecolor": INK,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.color": MUTED_2,
        "ytick.color": MUTED_2,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "grid.color": RULE,
        "font.family": "monospace",
        "font.monospace": ["IBM Plex Mono", "DejaVu Sans Mono", "monospace"],
        "legend.fontsize": 8,
        "legend.frameon": False,
        "savefig.facecolor": PANEL,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
    }


CSS = f"""
<style>
  .stApp {{ background: {GROUND}; }}
  h1, h2, h3 {{ font-family: 'IBM Plex Serif', Georgia, serif !important;
                letter-spacing: -0.3px; color: {INK}; }}
  code, kbd, .mono {{ font-family: 'IBM Plex Mono', monospace !important; }}
  .panel {{ background: {PANEL}; border: 1px solid {BORDER}; border-radius: 8px;
            padding: 14px 16px; margin-bottom: 12px; }}
  .kv {{ display: flex; justify-content: space-between; align-items: baseline;
         padding: 3px 0; }}
  .kv span:first-child {{ font-size: 12px; color: {INK_2}; }}
  .kv span:last-child  {{ font-family: 'IBM Plex Mono', monospace;
                          font-size: 13px; color: {INK}; }}
  .big {{ font-family: 'IBM Plex Mono', monospace; font-size: 26px;
          font-weight: 500; color: {INK}; line-height: 1.1; }}
  .lbl {{ font-size: 11px; color: {MUTED_2}; text-transform: uppercase;
          letter-spacing: 0.7px; }}
  .note {{ font-size: 11px; color: {MUTED}; line-height: 1.5; }}
  .warn {{ background: {AMBER_WASH}; border: 1px solid {AMBER_EDGE};
           border-radius: 6px; padding: 10px 13px; font-size: 12px;
           color: {MUTED}; line-height: 1.55; }}
  .warn strong {{ color: {AMBER_INK}; }}
  .good {{ background: {TEAL_WASH}; border: 1px solid {TEAL_EDGE};
           border-radius: 6px; padding: 10px 13px; font-size: 12px;
           color: {TEAL_DARK}; line-height: 1.55; }}
  .unavail {{ border: 1px dashed {BORDER}; border-radius: 6px;
              padding: 16px; font-size: 12px; color: {MUTED_2};
              line-height: 1.55; background: transparent; }}
</style>
<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Serif:wght@500;600&display=swap" rel="stylesheet">
"""
