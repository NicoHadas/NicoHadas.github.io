"""Figure for the eigendecomposition post.

Left:  unit circle with the eigenvectors of A.
Right: the ellipse A maps the circle to; eigenvectors are only scaled,
       ordinary vectors change direction.

Run:  _figures/.venv/bin/python _figures/eigendecomposition_figure.py
Needs: numpy, matplotlib (see requirements.txt)
Output: assets/img/eigendecomposition-geometry.{png,svg}
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Polygon

# Palette taken from assets/css/style.css
TEXT, MUTED, RULE = "#111111", "#666666", "#e5e5e5"
C1, C2 = "#0b57d0", "#d9480f"  # v1 (blue = site link colour), v2 (orange)
GHOST = "#b5b5b5"

plt.rcParams.update({
    "font.family": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "mathtext.fontset": "stixsans",
    "text.color": TEXT,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})

A = np.array([[2.0, 1.0], [1.0, 2.0]])
evals, evecs = np.linalg.eigh(A)
order = np.argsort(evals)[::-1]  # lambda_1 = 3, lambda_2 = 1
evals, evecs = evals[order], evecs[:, order]
evecs[:, 0] *= np.sign(evecs[0, 0])  # v1 = (1, 1)/sqrt(2)
evecs[:, 1] *= np.sign(evecs[0, 1])  # v2 = (1, -1)/sqrt(2)
v1, v2 = evecs[:, 0], evecs[:, 1]

t = np.linspace(0, 2 * np.pi, 400)
circle = np.vstack([np.cos(t), np.sin(t)])
ellipse = A @ circle
L = 3.35  # axis half-range


def arrow(ax, vec, color, lw=2.4, ls="-", z=5):
    ax.annotate(
        "", xy=vec, xytext=(0, 0), zorder=z,
        arrowprops=dict(arrowstyle="-|>,head_width=0.28,head_length=0.55",
                        color=color, lw=lw, ls=ls, shrinkA=0, shrinkB=0,
                        capstyle="round"),
    )


def label(ax, xy, text, color, dx=0.0, dy=0.0, size=18):
    ax.text(xy[0] + dx, xy[1] + dy, text, color=color, fontsize=size,
            ha="center", va="center", zorder=6)


def style(ax, title):
    ax.set_aspect("equal")
    ax.set_xlim(-L, L)
    ax.set_ylim(-L, L)
    ax.axis("off")
    ax.set_title(title, fontsize=17, color=MUTED, pad=10)
    # eigen-axes as faint guide lines
    for v, c in ((v1, C1), (v2, C2)):
        ax.plot(*np.outer(v, [-L * 1.2, L * 1.2]), color=c, lw=0.8,
                alpha=0.25, zorder=1)
    # minimal coordinate cross
    ax.plot([-L, L], [0, 0], color=RULE, lw=0.8, zorder=0)
    ax.plot([0, 0], [-L, L], color=RULE, lw=0.8, zorder=0)


fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10, 5.2))

# ---- Left: before ----------------------------------------------------------
style(ax0, "Before")
ax0.add_patch(Circle((0, 0), 1, fc="#0b57d0", ec="none", alpha=0.06, zorder=1))
ax0.plot(*circle, color=TEXT, lw=1.4, zorder=2)
arrow(ax0, v1, C1)
arrow(ax0, v2, C2)
label(ax0, v1, r"$v_1$", C1, dx=-0.38, dy=0.5)
label(ax0, v2, r"$v_2$", C2, dx=0.4, dy=-0.05)

# ---- Right: after ----------------------------------------------------------
style(ax1, r"After applying $A$")
ax1.plot(*circle, color=GHOST, lw=1.0, ls=(0, (3, 3)), zorder=2)  # ghost circle
ax1.add_patch(Polygon(ellipse.T, closed=True, fc="#0b57d0", ec="none",
                      alpha=0.06, zorder=1))
ax1.plot(*ellipse, color=TEXT, lw=1.4, zorder=3)
arrow(ax1, A @ v1, C1)
arrow(ax1, A @ v2, C2)
label(ax1, A @ v1, rf"$Av_1 = {evals[0]:.0f}\,v_1$", C1, dx=-0.45, dy=0.75)
label(ax1, A @ v2, rf"$Av_2 = {evals[1]:.0f}\,v_2$", C2, dx=0.95, dy=-0.3)

# ordinary (non-eigen) vectors, before and after
for ax, M in ((ax0, np.eye(2)), (ax1, A)):
    for e in (np.array([1.0, 0.0]), np.array([0.0, 1.0])):
        arrow(ax, M @ e, GHOST, lw=1.4, ls=(0, (3, 2)), z=4)

fig.tight_layout(w_pad=2)

out_dir = Path(__file__).resolve().parent.parent / "assets" / "img"
out_dir.mkdir(parents=True, exist_ok=True)
for ext in ("png", "svg"):
    path = out_dir / f"eigendecomposition-geometry.{ext}"
    fig.savefig(path, dpi=220, bbox_inches="tight", pad_inches=0.15)
    print("saved", path)
