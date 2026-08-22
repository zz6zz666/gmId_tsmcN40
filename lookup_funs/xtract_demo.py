"""
EKV parameter extraction demo (Appendix A.1).
Reproduces the key figures from Section A.1.4:
    - (a) Drain current density  : original vs. reconstructed
    - (b) Transconductance efficiency : original vs. reconstructed
    - (c) Percentage difference between the two

Usage
-----
    python xtract_demo.py

The script assumes that a TSMC N40 lookup table (.h5) is available;
adjust the ``loadmat`` path below.
"""

import numpy as np
import matplotlib.pyplot as plt
from lookup_table import loadmat, lookup
from ekv_extract import XTRACT

# ------------------------------------------------------------------
# Load lookup table (adjust path to your installation if necessary)
# ------------------------------------------------------------------
nch = loadmat(r"D:\tsmcN40_lookup\nch_tt.h5")

# ------------------------------------------------------------------
# Example: N40 core device  (L = 0.10 µm, VDS = 0.55 V, VSB = 0 V)
# ------------------------------------------------------------------
# ------------------------------------------------------------------
# IMPORTANT: pick an L that actually exists in your lookup table.
# For the N40 grid the shortest channel is L = 0.04 µm.
# ------------------------------------------------------------------
L = 0.1           # µm
VDS = 0.55        # V
VSB = 0.0         # V
rho = 0.6
TEMP = 300.0

y = XTRACT(nch, L, VDS, VSB, rho=rho, TEMP=TEMP)

print("=" * 50)
print("XTRACT output (L = 0.10 µm, VDS = 0.55 V)")
print("=" * 50)
print(f"  VDS           = {y[0]:.4f} V")
print(f"  n             = {y[1]:.4f}")
print(f"  VT            = {y[2]:.4f} V")
print(f"  JS            = {y[3] * 1e6:.4f} µA/µm")
print(f"  dn/dVDS       = {y[4]:.4f}")
print(f"  dVT/dVDS      = {y[5]:.4f}")
print(f"  dlogJS/dVDS   = {y[6]:.4f}")
print(f"  d²n/dVDS²     = {y[7]:.4f}")
print(f"  d²VT/dVDS²    = {y[8]:.4f}")
print(f"  d²logJS/dVDS² = {y[9]:.4f}")
print("=" * 50)

# ------------------------------------------------------------------
# Reconstruct JD and gm/ID from the extracted basic EKV model
# ------------------------------------------------------------------
k = 1.380649e-23
qe = 1.602176634e-19
UT = k * TEMP / qe

n = y[1]
VT = y[2]
JS = y[3]

# Sweep normalized charge density q  (A.1.1 – A.1.2)
q = np.logspace(-3, 1, 200)
i = q ** 2 + q
VP = UT * (2.0 * (q - 1.0) + np.log(q))

VGS_ekv = n * VP + VT          # (A.1.3)
JD_ekv = i * JS                # J_D = i · J_S   [A/µm]
gm_ID_ekv = 1.0 / (n * UT * (1.0 + q))   # analytical gm/ID for the basic EKV

# ------------------------------------------------------------------
# Original lookup-table data
# ------------------------------------------------------------------
vgs = nch.VGS
JD_lut = lookup(nch, 'ID_W', 'VGS', vgs, 'VDS', VDS, 'VSB', VSB, 'L', L,
                WARNING='off')
gm_ID_lut = lookup(nch, 'GM_ID', 'VGS', vgs, 'VDS', VDS, 'VSB', VSB, 'L', L,
                   WARNING='off')

JD_lut = np.asarray(JD_lut).flatten()
gm_ID_lut = np.asarray(gm_ID_lut).flatten()

# Locate the reference point on the reconstructed curves for plotting
# (interpolate VGS from gm_ID_ekv)
gm_id_inc = gm_ID_ekv[::-1]
vgs_inc = VGS_ekv[::-1]
gm_id_inc, uniq = np.unique(gm_id_inc, return_index=True)
vgs_inc = vgs_inc[uniq]
VGS_o = float(np.interp(rho * np.max(gm_ID_lut), gm_id_inc, vgs_inc))
JD_o = float(np.interp(VGS_o, VGS_ekv, JD_ekv))

# ------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

# (a) Drain current density
ax = axes[0]
ax.semilogy(vgs, np.abs(JD_lut) * 1e6, 'b-', lw=1.5, label='Lookup data')
ax.semilogy(VGS_ekv, JD_ekv * 1e6, 'r--', lw=1.5, label='Basic EKV')
ax.plot(VGS_o, JD_o * 1e6, 'ko', ms=7, mfc='none', mew=1.2,
        label=r'Ref. ($\rho=0.6$)')
ax.set_xlabel(r'$V_{GS}$ (V)')
ax.set_ylabel(r'$J_D$ ($\mu$A/$\mu$m)')
ax.set_title('(a) Drain Current Density')
ax.legend(loc='best', fontsize=9)
ax.grid(True, ls='--', alpha=0.4)

# (b) Transconductance efficiency
ax = axes[1]
ax.plot(vgs, gm_ID_lut, 'b-', lw=1.5, label='Lookup data')
ax.plot(VGS_ekv, gm_ID_ekv, 'r--', lw=1.5, label='Basic EKV')
ax.plot(VGS_o, rho * np.max(gm_ID_lut), 'ko', ms=7, mfc='none', mew=1.2,
        label=r'Ref. ($\rho=0.6$)')
ax.set_xlabel(r'$V_{GS}$ (V)')
ax.set_ylabel(r'$g_m/I_D$ (S/A)')
ax.set_title('(b) Transconductance Efficiency')
ax.legend(loc='best', fontsize=9)
ax.grid(True, ls='--', alpha=0.4)

# (c) Percentage difference
ax = axes[2]
JD_ekv_on_lut = np.interp(vgs, VGS_ekv, JD_ekv)
with np.errstate(divide='ignore', invalid='ignore'):
    D = (JD_ekv_on_lut - JD_lut) / JD_lut * 100.0

# Only plot where the original current is physically meaningful
valid = JD_lut > 0
ax.plot(vgs[valid], D[valid], 'g-', lw=1.5)
ax.set_xlabel(r'$V_{GS}$ (V)')
ax.set_ylabel(r'$D$ (%)')
ax.set_title('(c) Percentage Difference')
ax.set_ylim([-20, 20])
ax.axhline(0, color='k', ls='-', lw=0.5)
ax.grid(True, ls='--', alpha=0.4)

plt.tight_layout()
plt.show()
