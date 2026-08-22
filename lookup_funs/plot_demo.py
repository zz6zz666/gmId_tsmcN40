"""
Comprehensive validation script for the Python lookup interface.
Plots a 3x3 grid of the most common gm/ID design charts.

Features:
- Pan / Zoom via the native matplotlib toolbar
- Click on any curve to see exact coordinates (like Matlab Data Cursor)
"""

import matplotlib

matplotlib.use("qtagg")

import matplotlib.pyplot as plt
import numpy as np
import mplcursors

from lookup_table import loadmat, lookup

# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------
nch = loadmat(r"D:\tsmcN40_lookup\nch_tt.h5")

VDS = 0.55
VSB = 0
lengths = [0.04, 0.1, 0.18, 0.35, 1.0]  # µm
colors = plt.cm.viridis(np.linspace(0, 0.9, len(lengths)))

# ------------------------------------------------------------------
# Create a 3x3 figure
# ------------------------------------------------------------------
fig, axes = plt.subplots(3, 3, figsize=(14, 12))
fig.suptitle(
    f"TSMC 40 nm nch  ({nch.CORNER.upper()} corner,  VDS={VDS} V,  VSB={VSB} V)",
    fontsize=14,
    y=0.98,
)

for ax in axes.flat:
    ax.grid(True, ls="--", alpha=0.4)

all_lines = []

# Helper to collect lines and add legend
def plot_sweep(ax, ydata, xdata, label, ylabel, logy=False, ylim=None):
    for L, c in zip(lengths, colors):
        if logy:
            (line,) = ax.semilogy(xdata, ydata[L], color=c, lw=1.5, label=f"L={L:.2f} µm")
        else:
            (line,) = ax.plot(xdata, ydata[L], color=c, lw=1.5, label=f"L={L:.2f} µm")
        all_lines.append(line)
    ax.set_xlabel(r"$V_{GS}$ (V)")
    ax.set_ylabel(ylabel)
    ax.set_title(label)
    ax.legend(loc="best", fontsize=7)
    if ylim is not None:
        ax.set_ylim(ylim)


def plot_parametric(ax, xdata_dict, ydata_dict, label, xlabel, ylabel, logy=False, ylim=None):
    for L, c in zip(lengths, colors):
        x = np.asarray(xdata_dict[L]).flatten()
        y = np.asarray(ydata_dict[L]).flatten()
        valid = np.isfinite(x) & np.isfinite(y) & (y > 0 if logy else True)
        if logy:
            (line,) = ax.semilogy(x[valid], y[valid], color=c, lw=1.5, label=f"L={L:.2f} µm")
        else:
            (line,) = ax.plot(x[valid], y[valid], color=c, lw=1.5, label=f"L={L:.2f} µm")
        all_lines.append(line)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(label)
    ax.legend(loc="best", fontsize=7)
    if ylim is not None:
        ax.set_ylim(ylim)


# ------------------------------------------------------------------
# Pre-compute all sweeps (VGS sweep, fixed VDS/VSB/L)
# ------------------------------------------------------------------
vgs = nch.VGS

id_w = {}
gm_id = {}
vt = {}
vd_sat = {}
cgg_w = {}
cgd_w = {}
ft = {}
gain = {}

for L in lengths:
    id_w[L] = lookup(nch, "ID_W", "VGS", vgs, "VDS", VDS, "VSB", VSB, "L", L, WARNING="off") * 1e6
    gm_id[L] = lookup(nch, "GM_ID", "VGS", vgs, "VDS", VDS, "VSB", VSB, "L", L, WARNING="off")
    vt[L] = lookup(nch, "VT", "VGS", vgs, "VDS", VDS, "VSB", VSB, "L", L, WARNING="off")
    vd_sat[L] = lookup(nch, "VDSAT", "VGS", vgs, "VDS", VDS, "VSB", VSB, "L", L, WARNING="off")
    cgg_w[L] = lookup(nch, "CGG", "VGS", vgs, "VDS", VDS, "VSB", VSB, "L", L, WARNING="off") / nch.W * 1e15
    cgd_w[L] = lookup(nch, "CGD", "VGS", vgs, "VDS", VDS, "VSB", VSB, "L", L, WARNING="off") / nch.W * 1e15
    ft[L] = lookup(nch, "FT", "VGS", vgs, "VDS", VDS, "VSB", VSB, "L", L, WARNING="off") / 1e9
    gain[L] = lookup(nch, "GM_GDS", "VGS", vgs, "VDS", VDS, "VSB", VSB, "L", L, WARNING="off")

# ------------------------------------------------------------------
# Row 1: vs VGS
# ------------------------------------------------------------------
plot_sweep(axes[0, 0], id_w, vgs, "(a) Drain Current Density", r"$I_D/W$ ($\mu$A/$\mu$m)", logy=True)
plot_sweep(axes[0, 1], gm_id, vgs, "(b) Transconductance Efficiency", r"$g_m/I_D$ (S/A)")
plot_sweep(axes[0, 2], vt, vgs, "(c) Threshold Voltage", r"$V_T$ (V)")

# ------------------------------------------------------------------
# Row 2: parametric vs gm/ID
# ------------------------------------------------------------------
plot_parametric(axes[1, 0], gm_id, ft, "(d) Transit Frequency", r"$g_m/I_D$ (S/A)", r"$f_T$ (GHz)", logy=True)
plot_parametric(axes[1, 1], gm_id, gain, "(e) Intrinsic Gain", r"$g_m/I_D$ (S/A)", r"$g_m/g_{ds}$ (V/V)", logy=True)
plot_parametric(axes[1, 2], gm_id, vd_sat, "(f) Saturation Voltage", r"$g_m/I_D$ (S/A)", r"$V_{DSAT}$ (V)")

# ------------------------------------------------------------------
# Row 3: more parametric / mixed
# ------------------------------------------------------------------
plot_parametric(axes[2, 0], gm_id, id_w, "(g) Current Density", r"$g_m/I_D$ (S/A)", r"$I_D/W$ ($\mu$A/$\mu$m)", logy=True)
plot_parametric(axes[2, 1], gm_id, cgg_w, "(h) Gate Capacitance", r"$g_m/I_D$ (S/A)", r"$C_{gg}/W$ (fF/$\mu$m)")

# (i) Capacitance ratio CGD/CGG vs gm/ID
cgd_cgg = {}
for L in lengths:
    cgg = lookup(nch, "CGG", "VGS", vgs, "VDS", VDS, "VSB", VSB, "L", L, WARNING="off")
    cgd = lookup(nch, "CGD", "VGS", vgs, "VDS", VDS, "VSB", VSB, "L", L, WARNING="off")
    cgd_cgg[L] = np.asarray(cgd).flatten() / np.asarray(cgg).flatten()

plot_parametric(axes[2, 2], gm_id, cgd_cgg, "(i) Capacitance Ratio", r"$g_m/I_D$ (S/A)", r"$C_{gd}/C_{gg}$")

# ------------------------------------------------------------------
# Data cursor: click a curve to show a floating label
# ------------------------------------------------------------------
cursor = mplcursors.cursor(all_lines, hover=False)


@cursor.connect("add")
def on_cursor_add(sel):
    x, y = sel.target
    sel.annotation.set_text(f"x = {x:.4f}\ny = {y:.4f}")
    sel.annotation.get_bbox_patch().set_alpha(0.9)


plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()
