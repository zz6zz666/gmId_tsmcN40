"""
Python equivalents of Murmann's Matlab XTRACT.m and XTRACT2.m
(Appendix A.1 — The EKV Parameter Extraction Algorithm).
"""

import numpy as np
from lookup_table import lookup
from scipy.interpolate import PchipInterpolator


def XTRACT(dev, L, VDS, VSB, rho=0.6, TEMP=300.0):
    """
    Extract basic EKV parameters from lookup-table data.

    Syntax (Matlab-compatible) ::

        y = XTRACT(dev, L, VDS, VSB, rho=0.6, TEMP=300.0)

    Parameters
    ----------
    dev : LookupTable
        Lookup-table structure (e.g. ``nch``).
    L : float
        Channel length (scalar, µm).
    VDS : float or 1-D array_like
        Drain-to-source voltage (scalar or column vector, V).
    VSB : float
        Source-to-bulk voltage (scalar, V).
    rho : float, optional
        Normalized transconductance efficiency that defines the
        reference point.  Default is 0.6 (moderate inversion).
    TEMP : float, optional
        Absolute temperature (K).  Default is 300 K.

    Returns
    -------
    np.ndarray
        * If ``VDS`` is a scalar — shape ``(10,)`` ::

            [VDS, n, VT, JS,
             dn/dVDS, dVT/dVDS, dlogJS/dVDS,
             d²n/dVDS², d²VT/dVDS², d²logJS/dVDS²]

        * If ``VDS`` is a vector — shape ``(len(VDS), 10)`` with one row
          per VDS value.
    """
    k = 1.380649e-23       # Boltzmann constant  [J/K]
    qe = 1.602176634e-19   # Elementary charge   [C]
    UT = k * TEMP / qe     # Thermal voltage     [V]

    vds_query = np.atleast_1d(VDS).flatten().astype(float)
    L = float(L)
    VSB = float(VSB)

    # ------------------------------------------------------------------
    # 1) Extract n, VT, JS on an independent grid.  Derivatives are
    #    evaluated on half-grid points below, matching pygmid's XTRACT.
    # ------------------------------------------------------------------
    vds_grid = np.asarray(dev.VDS, dtype=float).flatten()
    if len(vds_grid) < 3:
        raise ValueError("XTRACT requires at least three VDS grid points.")
    vds_grid = vds_grid[1:]
    results = []
    for vds in vds_grid:
        vgs = dev.VGS.copy()
        # Ensure monotonically increasing for safe interpolation
        sort_idx = np.argsort(vgs)
        vgs = vgs[sort_idx]

        # Lookup drain-current density and gm/ID
        try:
            jd = lookup(dev, 'ID_W', 'VGS', vgs, 'VDS', vds,
                        'VSB', VSB, 'L', L, WARNING='off')
            gm_id = lookup(dev, 'GM_ID', 'VGS', vgs, 'VDS', vds,
                           'VSB', VSB, 'L', L, WARNING='off')
        except Exception:
            results.append([vds, np.nan, np.nan, np.nan])
            continue

        jd = np.asarray(jd).flatten()[sort_idx]
        gm_id = np.asarray(gm_id).flatten()[sort_idx]

        # Discard non-finite points
        valid = np.isfinite(jd) & np.isfinite(gm_id)
        vgs = vgs[valid]
        jd = jd[valid]
        gm_id = gm_id[valid]

        if len(vgs) == 0:
            results.append([vds, np.nan, np.nan, np.nan])
            continue

        # Maximum of gm/ID
        idx_max = np.argmax(gm_id)
        M = float(gm_id[idx_max])

        if M <= 0 or not np.isfinite(M):
            results.append([vds, np.nan, np.nan, np.nan])
            continue

        # Subthreshold slope factor  (A.1.5)
        n = 1.0 / (M * UT)

        # Reference gm/ID
        gm_id_ref = rho * M

        # Interpolate VGS_o from the right-hand side of the peak
        gm_id_side = gm_id[idx_max:]
        vgs_side = vgs[idx_max:]

        valid_side = np.isfinite(gm_id_side)
        gm_id_side = gm_id_side[valid_side]
        vgs_side = vgs_side[valid_side]

        if len(gm_id_side) < 2:
            results.append([vds, np.nan, np.nan, np.nan])
            continue

        # Reverse so that gm_id is monotonically *increasing* for np.interp
        gm_id_inc = gm_id_side[::-1]
        vgs_inc = vgs_side[::-1]

        # Remove possible duplicates in the independent variable
        gm_id_inc, uniq_idx = np.unique(gm_id_inc, return_index=True)
        vgs_inc = vgs_inc[uniq_idx]

        if gm_id_ref < gm_id_inc[0] or gm_id_ref > gm_id_inc[-1]:
            VGS_o = np.nan
        else:
            # Smooth interpolation avoids slope jumps when VGS_o crosses a LUT node.
            VGS_o = float(PchipInterpolator(gm_id_inc, vgs_inc)(gm_id_ref))

        if not np.isfinite(VGS_o):
            results.append([vds, np.nan, np.nan, np.nan])
            continue

        # Drain-current density at the reference point
        JD_o = float(PchipInterpolator(vgs, jd)(VGS_o))
        if not np.isfinite(JD_o) or JD_o <= 0:
            results.append([vds, np.nan, np.nan, np.nan])
            continue

        # Normalized mobile charge density  (A.1.6)
        q_o = 1.0 / rho - 1.0

        # Pinch-off voltage  (A.1.2)
        VP_o = UT * (2.0 * (q_o - 1.0) + np.log(q_o))

        # Threshold voltage  (A.1.7)
        VT = VGS_o - n * VP_o

        # Normalized drain current  (A.1.1)
        i_o = q_o ** 2 + q_o

        # Specific current density  (A.1.8)
        JS = JD_o / i_o

        results.append([vds, n, VT, JS])

    results = np.array(results, dtype=float)

    valid_rows = np.isfinite(results[:, 1:4]).all(axis=1) & (results[:, 3] > 0)
    if np.sum(valid_rows) < 4:
        out = np.full((len(vds_query), 10), np.nan)
        out[:, 0] = vds_query
        return out[0] if out.shape[0] == 1 else out

    fit_vds = results[valid_rows, 0]
    n_vec, VT_vec, JS_vec = results[valid_rows, 1:].T
    fit_vds1 = 0.5 * (fit_vds[:-1] + fit_vds[1:])
    fit_vds2 = 0.5 * (fit_vds1[:-1] + fit_vds1[1:])
    d1n = np.diff(n_vec) / np.diff(fit_vds)
    d1VT = np.diff(VT_vec) / np.diff(fit_vds)
    d1logJS = np.diff(np.log(JS_vec)) / np.diff(fit_vds)
    d2n = np.diff(d1n) / np.diff(fit_vds1)
    d2VT = np.diff(d1VT) / np.diff(fit_vds1)
    d2logJS = np.diff(d1logJS) / np.diff(fit_vds1)

    # PCHIP on the parameter grid and on the half-grid derivatives avoids
    # endpoint artifacts when the caller requests a VDS subrange.
    interp = lambda x, y, z: PchipInterpolator(x, y)(
        np.clip(z, x[0], x[-1])
    )
    n_q = interp(fit_vds, n_vec, vds_query)
    VT_q = interp(fit_vds, VT_vec, vds_query)
    JS_q = interp(fit_vds, JS_vec, vds_query)
    derivs = np.column_stack([
        interp(fit_vds1, d1n, vds_query),
        interp(fit_vds1, d1VT, vds_query),
        interp(fit_vds1, d1logJS, vds_query),
        interp(fit_vds2, d2n, vds_query),
        interp(fit_vds2, d2VT, vds_query),
        interp(fit_vds2, d2logJS, vds_query),
    ])
    out = np.column_stack([vds_query, n_q, VT_q, JS_q, derivs])
    if out.shape[0] == 1:
        return out.flatten()
    return out
def XTRACT2(VGS, ID, rho=0.6, TEMP=300.0):
    """
    Extract basic EKV parameters from directly-supplied I_D(V_GS) data.

    Syntax (Matlab-compatible) ::

        p = XTRACT2(VGS, ID, rho=0.6, TEMP=300.0)

    Parameters
    ----------
    VGS : array_like
        Gate-to-source voltage (column vector, V).
    ID : array_like
        Drain current (column vector or matrix, A).
        If a matrix, each *column* corresponds to a different transistor
        curve (same VGS sweep).
    rho : float, optional
        Normalized transconductance efficiency.  Default is 0.6.

    Returns
    -------
    np.ndarray
        * If ``ID`` is a 1-D vector — shape ``(3,)`` with ``[n, VT, IS]``.
        * If ``ID`` is a 2-D matrix — shape ``(M, 3)`` where ``M`` is the
          number of curves.  Each row holds ``[n, VT, IS]``.
    """
    VGS = np.asarray(VGS).flatten()
    ID_arr = np.asarray(ID)
    squeeze_output = ID_arr.ndim == 1

    # Ensure ID is 2-D with shape (N, M) — N VGS points, M curves
    if ID_arr.ndim == 1:
        ID_arr = ID_arr.reshape(-1, 1)
    elif ID_arr.ndim == 2:
        if ID_arr.shape[0] != len(VGS) and ID_arr.shape[1] == len(VGS):
            ID_arr = ID_arr.T
    else:
        raise ValueError("ID must be a vector or a 2-D matrix.")

    if ID_arr.shape[0] != len(VGS):
        raise ValueError(
            f"ID length along first axis ({ID_arr.shape[0]}) must match "
            f"VGS length ({len(VGS)})."
        )

    N, M = ID_arr.shape

    k = 1.380649e-23
    qe = 1.602176634e-19
    UT = k * TEMP / qe

    out = []
    for m in range(M):
        id_vec = ID_arr[:, m].copy()
        vgs = VGS.copy()

        valid = np.isfinite(id_vec) & (id_vec > 0)
        vgs = vgs[valid]
        id_vec = id_vec[valid]

        if len(vgs) < 3:
            out.append([np.nan, np.nan, np.nan])
            continue

        # Transconductance efficiency
        # Use log-gradient for numerical stability (weak-inversion ID is exponential)
        gm_id = np.gradient(np.log(id_vec), vgs)

        # Maximum of gm/ID
        idx_max = np.argmax(gm_id)
        M_val = float(gm_id[idx_max])

        if M_val <= 0 or not np.isfinite(M_val):
            out.append([np.nan, np.nan, np.nan])
            continue

        # Subthreshold slope factor  (A.1.5)
        n = 1.0 / (M_val * UT)

        # Reference point
        gm_id_ref = rho * M_val

        # Interpolate VGS_o on the right-hand side of the peak
        gm_id_side = gm_id[idx_max:]
        vgs_side = vgs[idx_max:]

        valid_side = np.isfinite(gm_id_side)
        gm_id_side = gm_id_side[valid_side]
        vgs_side = vgs_side[valid_side]

        if len(gm_id_side) < 2:
            out.append([np.nan, np.nan, np.nan])
            continue

        gm_id_inc = gm_id_side[::-1]
        vgs_inc = vgs_side[::-1]

        gm_id_inc, uniq_idx = np.unique(gm_id_inc, return_index=True)
        vgs_inc = vgs_inc[uniq_idx]

        if gm_id_ref < gm_id_inc[0] or gm_id_ref > gm_id_inc[-1]:
            VGS_o = np.nan
        else:
            VGS_o = float(PchipInterpolator(gm_id_inc, vgs_inc)(gm_id_ref))

        if not np.isfinite(VGS_o):
            out.append([np.nan, np.nan, np.nan])
            continue

        ID_o = float(PchipInterpolator(vgs, id_vec)(VGS_o))
        if not np.isfinite(ID_o) or ID_o <= 0:
            out.append([np.nan, np.nan, np.nan])
            continue

        q_o = 1.0 / rho - 1.0
        VP_o = UT * (2.0 * (q_o - 1.0) + np.log(q_o))
        VT = VGS_o - n * VP_o
        i_o = q_o ** 2 + q_o
        IS = ID_o / i_o

        out.append([n, VT, IS])

    out = np.array(out, dtype=float)
    if squeeze_output:
        return out.flatten()
    return out
