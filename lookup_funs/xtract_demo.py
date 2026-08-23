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
from scipy.optimize import brentq


def invq(x):
    x = np.atleast_1d(np.asarray(x, dtype=float))
    out = np.empty_like(x)
    for i, xx in enumerate(x):
        out[i] = brentq(
            lambda q: 2.0 * (q - 1.0) + np.log(q) - float(xx), 1e-12, 80.0
        )
    return float(out[0]) if out.size == 1 else out


def EKV_params(dev, L, VDS, VSB=0.0, rho=0.6, TEMP=300.0):
    return XTRACT(dev, L, VDS, VSB, rho=rho, TEMP=TEMP)


def _param_matrix(dev, L, VDS, VSB, rho, TEMP):
    return np.atleast_2d(XTRACT(dev, L, VDS, VSB, rho=rho, TEMP=TEMP))


def ekv_idvds(dev, L, VDS, VGS, VSB=0.0, rho=0.6, TEMP=300.0):
    UT = 1.380649e-23 * TEMP / 1.602176634e-19
    y = _param_matrix(dev, L, VDS, VSB, rho, TEMP)
    n, VT, JS = y[:, 1], y[:, 2], y[:, 3]
    qS = invq(((float(VGS) - VT) / n) / UT)
    ID = dev.W * JS * (qS**2 + qS)
    return ID, qS, y


def ekv_gds(dev, L, VDS, VGS, VSB=0.0, rho=0.6, TEMP=300.0):
    UT = 1.380649e-23 * TEMP / 1.602176634e-19
    y = _param_matrix(dev, L, VDS, VSB, rho, TEMP)
    if y.shape[0] < 2:
        raise ValueError("ekv_gds requires a VDS vector of at least two points.")
    SVT, SIS = y[:, 5], y[:, 6]
    ID, qS, _ = ekv_idvds(dev, L, VDS, VGS, VSB=VSB, rho=rho, TEMP=TEMP)
    gm = dev.W * y[:, 3] / (y[:, 1] * UT) * qS
    x = (float(VGS) - y[:, 2]) / y[:, 1] / UT
    return -gm * (SVT + x * y[:, 4]) + ID * SIS


def ekv_aint(dev, L, VDS, VGS, VSB=0.0, rho=0.6, TEMP=300.0):
    UT = 1.380649e-23 * TEMP / 1.602176634e-19
    y = _param_matrix(dev, L, VDS, VSB, rho, TEMP)
    if y.shape[0] < 2:
        raise ValueError("ekv_aint requires a VDS vector of at least two points.")
    SVT, SIS = y[:, 5], y[:, 6]
    n = y[:, 1]
    qS = invq(((float(VGS) - y[:, 2]) / n) / UT)
    gmid = 1.0 / (n * UT * (1.0 + qS))
    x = ((float(VGS) - y[:, 2]) / n) / UT
    return 1.0 / (SIS / gmid - SVT - x * y[:, 4])

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

vds_vec = nch.VDS[nch.VDS >= 0.3]
id_ekv, q_s, params = ekv_idvds(nch, L, vds_vec, 0.4, VSB=VSB, rho=rho, TEMP=TEMP)
id_lut = np.array([lookup(nch, 'ID', 'VGS', 0.4, 'VDS', vd, 'VSB', VSB,
                          'L', L, WARNING='off') for vd in vds_vec])
gds_ekv = ekv_gds(nch, L, vds_vec, 0.4, VSB=VSB, rho=rho, TEMP=TEMP)
gds_lut = np.array([lookup(nch, 'GDS', 'VGS', 0.4, 'VDS', vd, 'VSB', VSB,
                            'L', L, WARNING='off') for vd in vds_vec])

fig2, axes2 = plt.subplots(1, 2, figsize=(11, 4.5))
axes2[0].plot(vds_vec, id_lut * 1e6, 'b-', label='Lookup')
axes2[0].plot(vds_vec, id_ekv * 1e6, 'r--', label='EKV')
axes2[0].set(xlabel='VDS (V)', ylabel='ID (uA)', title='ID vs VDS')
axes2[0].legend(); axes2[0].grid(True, ls='--', alpha=0.4)
axes2[1].plot(vds_vec, gds_lut * 1e6, 'b-', label='Lookup')
axes2[1].plot(vds_vec, gds_ekv * 1e6, 'r--', label='EKV')
axes2[1].set(xlabel='VDS (V)', ylabel='gds (uS)', title='gds vs VDS')
axes2[1].legend(); axes2[1].grid(True, ls='--', alpha=0.4)
fig2.tight_layout()

vds_der = nch.VDS[(nch.VDS >= 0.5) & (nch.VDS <= 0.7)]
der = EKV_params(nch, L, vds_der, VSB, rho=rho, TEMP=TEMP)
i_der = np.argmin(np.abs(vds_der - 0.6))
n_der, vt_der = der[i_der, 1], der[i_der, 2]
dn_der, svt, sis = der[i_der, 4], der[i_der, 5], der[i_der, 6]
gmid_lut = np.asarray(lookup(nch, 'GM_ID', 'VGS', nch.VGS, 'VDS', 0.6,
                             'VSB', VSB, 'L', L, WARNING='off')).flatten()
gds_id_lut = np.asarray(lookup(nch, 'GDS', 'VGS', nch.VGS, 'VDS', 0.6,
                               'VSB', VSB, 'L', L, WARNING='off')).flatten() / np.asarray(
    lookup(nch, 'ID', 'VGS', nch.VGS, 'VDS', 0.6, 'VSB', VSB, 'L', L,
           WARNING='off')).flatten()
gain_lut = np.asarray(lookup(nch, 'GAIN', 'VGS', nch.VGS, 'VDS', 0.6,
                             'VSB', VSB, 'L', L, WARNING='off')).flatten()
q_gain = np.maximum(1.0 / (n_der * UT * gmid_lut) - 1.0, np.finfo(float).tiny)
x_gain = 2.0 * (q_gain - 1.0) + np.log(q_gain)
gds_id_ekv = -gmid_lut * (svt + x_gain * dn_der) + sis
gain_ekv = 1.0 / (sis / gmid_lut - svt - x_gain * dn_der)

fig3, axes3 = plt.subplots(1, 2, figsize=(11, 4.5))
axes3[0].plot(gmid_lut, gds_id_lut, 'b-', label='Lookup')
axes3[0].plot(gmid_lut, gds_id_ekv, 'r--', label='EKV')
axes3[0].set(xlabel='gm/ID (S/A)', ylabel='gds/ID (1/V)', title='gds/ID vs gm/ID')
axes3[0].legend(); axes3[0].grid(True, ls='--', alpha=0.4)
axes3[1].semilogy(gmid_lut, gain_lut, 'b-', label='Lookup')
axes3[1].semilogy(gmid_lut, gain_ekv, 'r--', label='EKV')
axes3[1].set(xlabel='gm/ID (S/A)', ylabel='Intrinsic gain', title='Intrinsic gain')
axes3[1].legend(); axes3[1].grid(True, ls='--', alpha=0.4)
fig3.tight_layout()
plt.show()
