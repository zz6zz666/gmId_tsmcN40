"""
Python equivalent of Murmann's Matlab lookup.m and lookupVGS.m
for gm/ID lookup tables (HDF5 .h5 or MAT .mat files).

The 4-D variable arrays follow the official Murmann lookup-table layout
(L, VGS, VDS, VSB) — the same order used by techsweep_spectre_run.m.
"""

import os
import numpy as np
from scipy.interpolate import interpn, PchipInterpolator, interp1d
import h5py


class LookupTable:
    """
    Container for a single MOSFET lookup table loaded from a MAT (.mat)
    or HDF5 (.h5) file.
    """

    def __init__(self, filepath):
        ext = os.path.splitext(filepath)[1].lower()
        if ext in (".h5", ".hdf5"):
            data = self._read_h5(filepath)
        elif ext == ".mat":
            data = self._read_mat(filepath)
        else:
            raise ValueError(
                f"Unsupported lookup-table format '{ext}' (expected .mat or .h5)."
            )

        # Normalise shape: vectors -> 1-D, scalars -> Python scalar
        for attr in ["VSB", "VGS", "VDS", "L"]:
            if hasattr(self, attr):
                setattr(self, attr, np.atleast_1d(getattr(self, attr)).flatten())

        for attr in ["W", "NFING", "TEMP"]:
            if hasattr(self, attr):
                v = getattr(self, attr)
                if isinstance(v, np.ndarray):
                    v = v.item()
                setattr(self, attr, v)

        # Collect all 4-D sweep variables
        self._vars = {}
        for key in [k for k in dir(self) if not k.startswith("_")]:
            val = getattr(self, key)
            if isinstance(val, np.ndarray) and val.ndim == 4:
                self._vars[key] = val

        # Axes must match the dimension order of the 4-D arrays.
        # The official Murmann lookup-table layout is (L, VGS, VDS, VSB),
        # matching the data files produced by techsweep_spectre_run.m.
        self._axes = (self.L, self.VGS, self.VDS, self.VSB)

    # ------------------------------------------------------------------
    # file readers
    # ------------------------------------------------------------------
    def _read_h5(self, filepath):
        with h5py.File(filepath, "r") as f:
            for key in f.keys():
                val = f[key]
                if isinstance(val, h5py.Dataset):
                    arr = np.array(val)
                    if arr.dtype == object:
                        s = arr[()]
                        if isinstance(s, bytes):
                            s = s.decode("utf-8")
                        setattr(self, key, s)
                    else:
                        arr = np.atleast_1d(arr).squeeze()
                        setattr(self, key, arr)

    def _read_mat(self, filepath):
        """Read a MATLAB v5 (.mat) lookup table.

        Variable arrays are stored in the official (L, VGS, VDS, VSB) order.
        """
        import scipy.io

        m = scipy.io.loadmat(filepath, squeeze_me=True, struct_as_record=False)
        key = next(k for k in m if not k.startswith("__"))
        s = m[key]
        for name in dir(s):
            if name.startswith("_"):
                continue
            val = getattr(s, name)
            if isinstance(val, str):
                setattr(self, name, val)
            elif isinstance(val, np.ndarray):
                if val.dtype.kind in "SUO":
                    setattr(self, name, str(np.asarray(val).item()))
                else:
                    setattr(self, name, np.atleast_1d(val).squeeze())
            elif np.isscalar(val):
                setattr(self, name, val)

    def _get_single(self, name):
        """Return a raw 4-D array by name, or None."""
        if name in self._vars:
            return self._vars[name]
        if hasattr(self, name):
            val = getattr(self, name)
            if isinstance(val, np.ndarray) and val.ndim == 4:
                return val
        return None

    def get_array(self, name):
        """
        Resolve a variable name or a ratio string (e.g. 'GM_ID', 'GM_CGG').

        Stored primitives are returned directly. Ratios are recomputed on the
        fly by the generic ``A/B`` rule below, using the textbook naming:

            GM_ID  = GM/ID
            GM_CGG = GM/CGG      (=> fT = GM_CGG/(2 pi))
            GM_GDS = GM/GDS      (=> intrinsic gain)
            ID_W   = ID/W
        """
        arr = self._get_single(name)
        if arr is not None:
            return arr

        if "_" in name:
            parts = name.split("_")
            for i in range(1, len(parts)):
                a = "_".join(parts[:i])
                b = "_".join(parts[i:])
                num = self._get_single(a)
                if num is None:
                    continue
                if b == "W":
                    den = self.W
                else:
                    den = self._get_single(b)
                if den is None:
                    continue
                with np.errstate(divide="ignore", invalid="ignore"):
                    return num / den

        raise KeyError(f"Variable or ratio '{name}' not found in lookup table.")


def _parse_params(args, kwargs):
    """Convert Matlab-style key-value pairs into a dict."""
    if len(args) % 2:
        raise ValueError("Arguments must be supplied as name-value pairs.")
    params = {}
    i = 0
    while i < len(args):
        key = str(args[i]).upper()
        val = np.atleast_1d(args[i + 1])
        params[key] = val.flatten()
        i += 2
    for k, v in kwargs.items():
        params[k.upper()] = np.atleast_1d(v).flatten()
    return params


def lookup(data, outvar, *args, **kwargs):
    """
    Python equivalent of Murmann's Matlab ``lookup.m``.

    Usage modes
    -----------
    1. Basic lookup ::

        id = lookup(data, 'ID', 'VGS', 0.6, 'VDS', 0.55, 'L', 0.1, 'VSB', 0)

    2. Ratio lookup ::

        gm_id = lookup(data, 'GM_ID', 'VGS', 0.6, 'VDS', 0.55, 'L', 0.1)

    3. Cross-lookup (one ratio versus another) ::

        wt = lookup(data, 'GM_CGG', 'GM_ID', np.array([5, 10, 15]),
                    'VDS', 0.55, 'L', 0.1, 'VSB', 0)

    Parameters
    ----------
    data : LookupTable
        Loaded lookup table.
    outvar : str
        Output variable or ratio (e.g. ``'ID'``, ``'GM_ID'``, ``'GM_CGG'``).
    *args : pairs of (str, array_like)
        Additional arguments as key-value pairs.
    **kwargs : dict
        Additional arguments as keyword parameters.

    Returns
    -------
    np.ndarray
        Interpolated result (singleton dimensions are squeezed out).
    """
    params = _parse_params(args, kwargs)

    method = params.pop("METHOD", np.array(["pchip"])).flatten()
    method = str(method[0]).lower() if len(method) > 0 else "pchip"

    warning_flag = params.pop("WARNING", np.array(["on"])).flatten()
    warning_flag = str(warning_flag[0]).lower() if len(warning_flag) > 0 else "on"

    default_params = {
        "L": np.array([np.min(data.L)]),
        "VGS": data.VGS,
        "VDS": np.array([np.max(data.VDS) / 2.0]),
        "VSB": np.array([0.0]),
    }

    coord_keys = {"L", "VGS", "VDS", "VSB", "METHOD", "WARNING"}
    non_coord_keys = [k for k in params if k not in coord_keys]

    # ------------------------------------------------------------------
    # Mode 3: cross-lookup of one ratio against another
    # ------------------------------------------------------------------
    if len(non_coord_keys) == 1:
        in_ratio = non_coord_keys[0]
        if "_" not in outvar or "_" not in in_ratio:
            raise ValueError(
                "Invalid syntax or usage mode: cross-lookup requires both "
                "the output and input variables to be ratios."
            )
        in_values = params.pop(in_ratio)
        out_ratio = outvar

        vsb = np.atleast_1d(params.pop("VSB", default_params["VSB"])).flatten()
        vds = np.atleast_1d(params.pop("VDS", default_params["VDS"])).flatten()
        l = np.atleast_1d(params.pop("L", default_params["L"])).flatten()
        vgs_full = np.atleast_1d(params.pop("VGS", data.VGS)).flatten()

        if params:
            raise ValueError(
                f"Invalid syntax or usage mode! Unexpected keys: {list(params.keys())}"
            )

        # Build broadcast grid so that vector L/VDS/VSB are all supported
        vsb_grid, vds_grid, l_grid = np.meshgrid(vsb, vds, l, indexing="ij")
        coord_shape = vsb_grid.shape

        in_values_arr = np.atleast_1d(in_values).flatten()
        result = np.empty((*coord_shape, *in_values_arr.shape), dtype=float)

        for idx in np.ndindex(coord_shape):
            vsb_i = float(vsb_grid[idx])
            vds_i = float(vds_grid[idx])
            l_i = float(l_grid[idx])

            # Suppress noisy warnings from the internal sweeps
            x = lookup(
                data, in_ratio, "VGS", vgs_full, "VDS", vds_i, "VSB", vsb_i, "L", l_i,
                WARNING="off",
            )
            y = lookup(
                data, out_ratio, "VGS", vgs_full, "VDS", vds_i, "VSB", vsb_i, "L", l_i,
                WARNING="off",
            )
            x = np.asarray(x).flatten()
            y = np.asarray(y).flatten()

            # Remove NaN / inf (e.g. 0/0 at VGS = 0)
            valid = np.isfinite(x) & np.isfinite(y)
            x = x[valid]
            y = y[valid]

            if len(x) < 2:
                result[idx] = np.nan
                continue

            # Handle non-monotonicity for known problematic ratios
            if in_ratio == "GM_ID":
                idx_max = np.argmax(x)
                x = x[idx_max:]
                y = y[idx_max:]
            elif in_ratio in ("GM_CGG", "GM_CGS"):
                idx_max = np.argmax(x)
                x = x[: idx_max + 1]
                y = y[: idx_max + 1]
            else:
                dx = np.diff(x)
                if np.any(dx > 0) and np.any(dx < 0):
                    raise ValueError(
                        f"lookup: Error! There are multiple curve intersections for '{in_ratio}'. "
                        f"Try to reduce the search range by specifying the VGS vector explicitly. "
                        f"Example: lookup(data, '{out_ratio}', '{in_ratio}', val, "
                        f"'VGS', data.VGS[10:])"
                    )

            # Remove duplicate x values (required for interpolators)
            x, uniq_idx = np.unique(x, return_index=True)
            y = y[uniq_idx]

            if len(x) < 2:
                result[idx] = np.nan
                continue

            if method == "pchip":
                interpolator = PchipInterpolator(x, y, extrapolate=False)
            else:
                interpolator = interp1d(
                    x, y, kind=method, bounds_error=False, fill_value=np.nan
                )

            result[idx] = interpolator(in_values_arr)

        # Squeeze out singleton dimensions (mimics original scalar behaviour)
        result = np.squeeze(result)

        # Auto-extract scalar so callers don't need .item()
        if isinstance(result, np.ndarray) and result.ndim == 0:
            result = result.item()

        if warning_flag == "on" and np.any(np.isnan(result)):
            print(f"lookup warning: {in_ratio} input out of range! (output is NaN)")

        return result

    # ------------------------------------------------------------------
    # Mode 1 & 2: direct multidimensional interpolation
    # ------------------------------------------------------------------
    elif len(non_coord_keys) == 0:
        axes_names = ["L", "VGS", "VDS", "VSB"]
        query_points = []
        for ax_name in axes_names:
            if ax_name in params:
                vec = params[ax_name]
            else:
                vec = default_params[ax_name]
            query_points.append(vec)

        values = data.get_array(outvar)

        grids = np.meshgrid(*query_points, indexing="ij")
        xi = np.stack(grids, axis=-1)

        result = interpn(
            data._axes, values, xi, method="linear", bounds_error=False, fill_value=np.nan
        )
        result = result.reshape([len(q) for q in query_points])

        # Squeeze out singleton dimensions (mimics Matlab scalar/vector output)
        result = np.squeeze(result)

        if warning_flag == "on" and np.any(np.isnan(result)):
            # Distinguish between true out-of-bounds and ratio artefacts
            nan_frac = np.sum(np.isnan(result)) / result.size
            if nan_frac < 1.0:
                print(
                    "lookup warning: some requested points are out of bounds (NaN returned)."
                )

        return result

    else:
        raise ValueError(
            f"Invalid syntax or usage mode! Unexpected keys: {non_coord_keys}"
        )


def lookupVGS(data, *args, **kwargs):
    """
    Python equivalent of Murmann's Matlab ``lookupVGS.m``.

    Finds VGS for a given inversion level (``GM_ID``) or current density
    (``ID_W``) and terminal voltages.

    Mode 1 – known source terminal ::

        VGS = lookupVGS(data, 'GM_ID', 10, 'VDS', 0.55, 'VSB', 0.1, 'L', 0.1)

    Mode 2 – unknown source voltage ::

        VGS = lookupVGS(data, 'GM_ID', 10, 'VDB', 0.9, 'VGB', 1.0, 'L', 0.1)

    Parameters
    ----------
    data : LookupTable
        Loaded lookup table.
    *args : pairs of (str, array_like)
        Additional arguments as key-value pairs.
    **kwargs : dict
        Additional arguments as keyword parameters.

    Returns
    -------
    np.ndarray
        VGS values (1-D array).
    """
    params = _parse_params(args, kwargs)

    method = params.pop("METHOD", np.array(["pchip"])).flatten()
    method = str(method[0]).lower() if len(method) > 0 else "pchip"

    l_values = np.atleast_1d(
        params.pop("L", np.array([np.min(data.L)]))
    ).astype(float).flatten()

    # ------------------------------------------------------------------
    # Mode 2: unknown source (VDB / VGB supplied)
    # ------------------------------------------------------------------
    if "VDB" in params or "VGB" in params:
        if "VDB" not in params or "VGB" not in params:
            raise ValueError("lookupVGS mode 2 requires both 'VDB' and 'VGB'")
        vdb = float(params.pop("VDB")[0])
        vgb = float(params.pop("VGB")[0])
        if len(l_values) != 1:
            raise ValueError("lookupVGS mode 2 requires scalar 'L'")
        l = float(l_values[0])

        target_key = None
        target_val = None
        for k in list(params.keys()):
            if k in ("GM_ID", "ID_W"):
                target_key = k
                target_val = params.pop(k)
                break

        if target_key is None:
            raise ValueError("lookupVGS mode 2 requires 'GM_ID' or 'ID_W'")
        if params:
            raise ValueError(
                f"Invalid syntax or usage mode! Unexpected keys: {list(params.keys())}"
            )

        target_val = np.atleast_1d(target_val).flatten()

        # Constraints:
        #   VSB = VGB - VGS
        #   VDS = VDB - VSB = VDB - VGB + VGS
        # We search over the intersection of data.VGS with the feasible region.
        vgs_feasible = data.VGS.copy()
        vsb_search = vgb - vgs_feasible
        vds_search = vdb - vsb_search

        mask = (
            (vsb_search >= np.min(data.VSB))
            & (vsb_search <= np.max(data.VSB))
            & (vds_search >= np.min(data.VDS))
            & (vds_search <= np.max(data.VDS))
        )
        vgs_search = vgs_feasible[mask]
        vsb_search = vsb_search[mask]
        vds_search = vds_search[mask]

        if len(vgs_search) < 2:
            return np.full_like(target_val, np.nan)

        # Vectorised interpolation over the search sweep
        # (L, VGS, VDS, VSB) axis order
        xi = np.column_stack(
            [np.full_like(vgs_search, l), vgs_search, vds_search, vsb_search]
        )
        y_search = interpn(
            data._axes,
            data.get_array(target_key),
            xi,
            method="linear",
            bounds_error=False,
            fill_value=np.nan,
        )

        # Remove NaNs (out-of-bounds artefacts)
        valid = ~np.isnan(y_search)
        if np.sum(valid) < 2:
            return np.full_like(target_val, np.nan)
        vgs_search = vgs_search[valid]
        y_search = y_search[valid]

        # Non-monotonicity handling (same logic as lookup mode 3)
        if target_key == "GM_ID":
            idx_max = np.argmax(y_search)
            y_search = y_search[idx_max:]
            vgs_search = vgs_search[idx_max:]
        elif target_key in ("GM_CGG", "GM_CGS"):
            idx_max = np.argmax(y_search)
            y_search = y_search[: idx_max + 1]
            vgs_search = vgs_search[: idx_max + 1]

        if len(y_search) < 2:
            return np.full_like(target_val, np.nan)

        # Remove duplicate y values
        y_search, uniq_idx = np.unique(y_search, return_index=True)
        vgs_search = vgs_search[uniq_idx]

        if method == "pchip":
            interp = PchipInterpolator(y_search, vgs_search, extrapolate=False)
        else:
            interp = interp1d(
                y_search, vgs_search, kind=method, bounds_error=False, fill_value=np.nan
            )

        result = np.squeeze(interp(target_val))
        if target_key == "GM_ID" and np.any(np.isnan(result)):
            print("lookupVGS: GM_ID input larger than maximum!")
        if isinstance(result, np.ndarray) and result.ndim == 0:
            result = result.item()
        return result

    # ------------------------------------------------------------------
    # Mode 1: known source terminal
    # ------------------------------------------------------------------
    else:
        vds_values = np.atleast_1d(
            params.pop("VDS", np.array([np.max(data.VDS) / 2.0]))
        ).astype(float).flatten()
        vsb_values = np.atleast_1d(
            params.pop("VSB", np.array([0.0]))
        ).astype(float).flatten()

        target_key = None
        target_val = None
        for k in list(params.keys()):
            if k in ("GM_ID", "ID_W"):
                target_key = k
                target_val = params.pop(k)
                break

        if target_key is None:
            raise ValueError("lookupVGS requires 'GM_ID' or 'ID_W'")
        if params:
            raise ValueError(
                f"Invalid syntax or usage mode! Unexpected keys: {list(params.keys())}"
            )

        target_val = np.atleast_1d(target_val).flatten()
        vector_count = sum(
            len(v) > 1 for v in (target_val, l_values, vds_values, vsb_values)
        )
        if vector_count > 1:
            raise ValueError("lookupVGS accepts at most one vector input.")

        biases = np.array(
            np.meshgrid(l_values, vds_values, vsb_values, indexing="ij")
        ).reshape(3, -1).T
        rows = []
        for l, vds, vsb in biases:
            vgs_full = data.VGS.copy()
            y_search = np.asarray(
                lookup(
                    data, target_key, "VGS", vgs_full, "VDS", vds,
                    "VSB", vsb, "L", l, WARNING="off",
                )
            ).flatten()
            valid = np.isfinite(y_search)
            x = y_search[valid]
            y = vgs_full[valid]
            if len(x) < 2:
                rows.append(np.full_like(target_val, np.nan, dtype=float))
                continue
            if target_key == "GM_ID":
                idx_max = np.argmax(x)
                x, y = x[idx_max:], y[idx_max:]
            x, uniq_idx = np.unique(x, return_index=True)
            y = y[uniq_idx]
            if len(x) < 2:
                rows.append(np.full_like(target_val, np.nan, dtype=float))
                continue
            if method == "pchip":
                interp = PchipInterpolator(x, y, extrapolate=False)
            else:
                interp = interp1d(
                    x, y, kind=method, bounds_error=False, fill_value=np.nan
                )
            rows.append(np.asarray(interp(target_val), dtype=float))

        result = np.squeeze(np.asarray(rows))
        if target_key == "GM_ID" and np.any(np.isnan(result)):
            print("lookupVGS: GM_ID input larger than maximum!")
        if isinstance(result, np.ndarray) and result.ndim == 0:
            result = result.item()
        return result


def loadmat(filepath):
    """Convenience wrapper: load a MAT (.mat) or HDF5 (.h5) lookup table."""
    return LookupTable(filepath)
