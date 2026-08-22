#!/usr/bin/env python3
"""
Unified extraction script for TSMC N40 (CRN40LP) 1.1V core devices.
Produces .h5 lookup tables.

Storage format:
  - HDF5 container, one file per device+corner (e.g. nch_tt.h5)
  - 4-D variable arrays (VSB, VDS, VGS, L) stored as float32 with gzip
    compression + chunking (per-VSB chunks, one L per L-slice)
  - metadata as scalar/string datasets: CORNER, DEVICE, INFO, TEMP, W, NFING
  - axes as 1-D datasets: L, VGS, VDS, VSB
  - variables: ID VT IGD IGS GM GMB GDS CGG CGS CSG CGD CDG CGB CDD CSS
               FT GM_ID GAIN VDSAT   (+ noise STH SFL)

Usage:
  # Full extraction
  python extract_new.py --fine --outdir <path> --workers 5

  # Incremental: only L range
  python extract_new.py --fine --outdir <path> --L-range 0 5 --workers 5

  # Append to existing .h5
  python extract_new.py --fine --outdir <path> --L-range 9 13 --workers 5
"""
import time, sys, os, numpy as np
from multiprocessing import Pool


def save_mat(filename, data_dict):
    import h5py
    with h5py.File(filename + '.h5', 'w') as f:
        for k, v in data_dict.items():
            if isinstance(v, str):
                f.create_dataset(k, data=v, dtype=h5py.special_dtype(vlen=str))
            elif np.isscalar(v):
                f.create_dataset(k, data=v)
            else:
                arr = np.asarray(v)
                if arr.ndim >= 4:
                    chunks = arr.shape[:3] + (1,)
                    f.create_dataset(k, data=arr, dtype=np.float32, chunks=chunks,
                                     compression='gzip', compression_opts=6,
                                     shuffle=True)
                else:
                    f.create_dataset(k, data=arr)
    print('  -> %s.h5' % filename)


def extract_one(raw_dir, c):
    from psf_reader import _read_all_traces, _parse_noise_psf

    dc_traces = _read_all_traces(raw_dir)

    dc_n, dc_p = {}, {}
    for (sig, units, coeffs) in c['n']:
        if sig in dc_traces:
            vals = dc_traces[sig].copy()
            for m, vname in enumerate(c['outvars']):
                if coeffs[m] != 0:
                    if vname not in dc_n:
                        dc_n[vname] = np.zeros_like(vals)
                    dc_n[vname] += vals * coeffs[m]

    for (sig, units, coeffs) in c['p']:
        if sig in dc_traces:
            vals = dc_traces[sig].copy()
            for m, vname in enumerate(c['outvars']):
                if coeffs[m] != 0:
                    if vname not in dc_p:
                        dc_p[vname] = np.zeros_like(vals)
                    dc_p[vname] += vals * coeffs[m]

    noise_n, noise_p = {}, {}
    noise_files = sorted([f for f in os.listdir(raw_dir) if f.endswith(".noise")])
    if noise_files:
        noise_data = [_parse_noise_psf(os.path.join(raw_dir, f)) for f in noise_files]
        nfiles = len(noise_data)

        for k, (sig, _) in enumerate(c['n_noise']):
            trace_name = sig.split(":")[0]
            field_name = sig.split(":")[1] if ":" in sig else ""
            nv = []
            for nd in noise_data:
                if field_name and trace_name in nd and isinstance(nd[trace_name], dict):
                    field_vals = nd[trace_name].get(field_name, [])
                    if isinstance(field_vals, list):
                        nv.append(field_vals)
                    else:
                        nv.append([field_vals])
                elif trace_name in nd and isinstance(nd[trace_name], list):
                    nv.append(nd[trace_name])
                elif trace_name in nd and isinstance(nd[trace_name], (int, float)):
                    nv.append([nd[trace_name]])
            if len(nv) == nfiles:
                arr = np.array(nv, dtype=np.float64)
                noise_n[c['outvars_noise'][k]] = arr

        for k, (sig, _) in enumerate(c['p_noise']):
            trace_name = sig.split(":")[0]
            field_name = sig.split(":")[1] if ":" in sig else ""
            nv = []
            for nd in noise_data:
                if field_name and trace_name in nd and isinstance(nd[trace_name], dict):
                    field_vals = nd[trace_name].get(field_name, [])
                    if isinstance(field_vals, list):
                        nv.append(field_vals)
                    else:
                        nv.append([field_vals])
                elif trace_name in nd and isinstance(nd[trace_name], list):
                    nv.append(nd[trace_name])
                elif trace_name in nd and isinstance(nd[trace_name], (int, float)):
                    nv.append([nd[trace_name]])
            if len(nv) == nfiles:
                arr = np.array(nv, dtype=np.float64)
                noise_p[c['outvars_noise'][k]] = arr

    return dc_n, dc_p, noise_n, noise_p


def make_mat_data(c, L_arr, device):
    nL = len(L_arr)
    nVGS = len(c['VGS'])
    nVDS = len(c['VDS'])
    nVSB = len(c['VSB'])
    d = {
        'INFO': c['modelinfo'], 'CORNER': c['corner'], 'DEVICE': device,
        'TEMP': c['temp'], 'NFING': c['NFING'], 'W': c['WIDTH'],
        'L': np.asarray(L_arr),
        'VGS': np.asarray(c['VGS']),
        'VDS': np.asarray(c['VDS']),
        'VSB': np.asarray(c['VSB']),
    }
    for v in c['outvars']:
        d[v] = np.zeros((nVSB, nVDS, nVGS, nL), dtype=np.float32)
    for v in c['outvars_noise']:
        d[v] = np.full((nVSB, nVDS, nVGS, nL), np.nan, dtype=np.float32)
    return d


def extract_corner(corner, fine, outdir, l_range, voltage, srcdir=None):
    from config_tsmcN40 import get_config

    c = get_config(corner, coarse=not fine)
    nVGS = len(c['VGS'])
    nVDS = len(c['VDS'])
    nVSB = len(c['VSB'])
    rundir_base = os.path.join(srcdir if srcdir else outdir, c['rundir_base'])

    L_n = c['LENGTH']
    L_p = c['LENGTH']
    all_L = list(L_n)
    idx_n = {v: i for i, v in enumerate(all_L)}
    idx_p = {v: i for i, v in enumerate(all_L)}

    l_start, l_end = l_range
    l_start = max(0, l_start)
    l_end = min(len(all_L), l_end)

    print()
    print('===== Corner: %s =====' % corner)
    print('L range: %d..%d (%.3f..%.3f um)' %
          (l_start, l_end - 1, all_L[l_start], all_L[l_end - 1]))

    fn_n = os.path.join(outdir, c['savefilen'])
    fn_p = os.path.join(outdir, c['savefilep'])
    import h5py

    mat_n_exists = os.path.exists(fn_n + '.h5')
    mat_p_exists = os.path.exists(fn_p + '.h5')
    append_mode = mat_n_exists and mat_p_exists

    if not append_mode:
        print('  Creating: %s.h5' % fn_n)
        print('  Creating: %s.h5' % fn_p)
        os.makedirs(outdir, exist_ok=True)
        save_mat(fn_n, make_mat_data(c, L_n, 'nch'))
        save_mat(fn_p, make_mat_data(c, L_p, 'pch'))
    else:
        print('  Appending to existing .h5 files')

    t0 = time.time()
    processed = 0

    nout = len(c['outvars'])
    nout_noise = len(c['outvars_noise'])

    with h5py.File(fn_n + '.h5', 'r+') as f_n, \
         h5py.File(fn_p + '.h5', 'r+') as f_p:

        for ii in range(l_start, l_end):
            lval = all_L[ii]

            buf_n = {v: np.zeros((nVSB, nVDS, nVGS), dtype=np.float32) for v in c['outvars']}
            buf_p = {v: np.zeros((nVSB, nVDS, nVGS), dtype=np.float32) for v in c['outvars']}
            buf_n_noise = {v: np.full((nVSB, nVDS, nVGS), np.nan, dtype=np.float32) for v in c['outvars_noise']}
            buf_p_noise = {v: np.full((nVSB, nVDS, nVGS), np.nan, dtype=np.float32) for v in c['outvars_noise']}

            for j in range(nVSB):
                raw_dir = '%s/L%03d_%.3fum_VSB%03d_%+.2fV.raw' % (
                    rundir_base, ii, lval, j, c['VSB'][j])

                t1 = time.time()
                try:
                    dc_n, dc_p, noise_n, noise_p = extract_one(raw_dir, c)
                except Exception as e:
                    print(' [%s] L=%.3fum VSB=%+.2fV ... SKIPPED (%s)' %
                          (c['corner'], lval, c['VSB'][j], str(e)[:50]))
                    continue
                elapsed = time.time() - t1

                if lval in idx_n:
                    for vname, vals in dc_n.items():
                        buf_n[vname][j, :, :] = vals
                    for vname, vals in noise_n.items():
                        buf_n_noise[vname][j, :, :] = vals

                if lval in idx_p:
                    for vname, vals in dc_p.items():
                        buf_p[vname][j, :, :] = vals
                    for vname, vals in noise_p.items():
                        buf_p_noise[vname][j, :, :] = vals

                processed += 1
                print(' [%s] L=%.3fum VSB=%+.2fV ... OK (%.1fs)' %
                      (c['corner'], lval, c['VSB'][j], elapsed))

            # Batch write all VSB at once per L
            if lval in idx_n:
                in_n = idx_n[lval]
                for vname in c['outvars']:
                    f_n[vname][:, :, :, in_n] = buf_n[vname]
                for vname in c['outvars_noise']:
                    f_n[vname][:, :, :, in_n] = buf_n_noise[vname]

            if lval in idx_p:
                ip = idx_p[lval]
                for vname in c['outvars']:
                    f_p[vname][:, :, :, ip] = buf_p[vname]
                for vname in c['outvars_noise']:
                    f_p[vname][:, :, :, ip] = buf_p_noise[vname]

    t_total = time.time() - t0
    print('  [%s] %d combos in %.0fs' % (corner, processed, t_total))

    # Verify
    ok = _verify(fn_n + '.h5', l_start, l_end, all_L, idx_n, nVSB)
    ok &= _verify(fn_p + '.h5', l_start, l_end, all_L, idx_p, nVSB)
    if ok:
        print('  Verify OK')
    else:
        print('  Verify FAILED')
    return ok


def _verify(fname, l_start, l_end, all_L, idx_map, nVSB):
    import h5py
    try:
        with h5py.File(fname, 'r') as f:
            id_data = np.array(f['ID'])
            for ii in range(l_start, l_end):
                lval = all_L[ii]
                if lval not in idx_map:
                    continue
                idx = idx_map[lval]
                a = id_data[0, :, :, idx].sum()
                b = id_data[nVSB - 1, :, :, idx].sum()
                if abs(a) < 1e-20 and abs(b) < 1e-20:
                    print('  WARN: L=%.3fum has near-zero ID' % lval)
                    return False
                if np.allclose(a, b):
                    print('  WARN: L=%.3fum same VSB (cache bug?)' % lval)
                    return False
        return True
    except Exception as e:
        print('  Verify error: %s' % e)
        return False


def _worker(args):
    return extract_corner(*args)


def main():
    args = sys.argv[1:]
    voltage = 'n40'
    fine = '--fine' in args
    if fine:
        args.remove('--fine')
    outdir = ''
    if '--outdir' in args:
        idx = args.index('--outdir')
        args.pop(idx)
        outdir = args.pop(idx) if idx < len(args) else ''
    srcdir = ''
    if '--srcdir' in args:
        idx = args.index('--srcdir')
        args.pop(idx)
        srcdir = args.pop(idx) if idx < len(args) else ''
    if '--voltage' in args:
        idx = args.index('--voltage')
        args.pop(idx)
        voltage = args.pop(idx) if idx < len(args) else '18'
    l_range = None
    if '--L-range' in args:
        idx = args.index('--L-range')
        args.pop(idx)
        if idx < len(args):
            s = int(args.pop(idx))
            e = int(args.pop(idx)) if idx < len(args) else s + 1
            l_range = (s, e)
    workers = 1
    if '--workers' in args:
        idx = args.index('--workers')
        args.pop(idx)
        workers = int(args.pop(idx)) if idx < len(args) else 1

    corners = args if args else ['tt', 'ff', 'ss', 'fs', 'sf']
    if not l_range:
        from config_tsmcN40 import get_config
        c0 = get_config(corners[0], coarse=not fine)
        nL = len(c0['LENGTH'])
        l_range = (0, nL)

    print('Fine: %s  Outdir: %s  Srcdir: %s  Workers: %d  L-range: %d-%d' %
          (fine, outdir if outdir else '.', srcdir if srcdir else '(same)', workers,
           l_range[0], l_range[1] - 1))

    t0 = time.time()
    tasks = [(c, fine, outdir, l_range, voltage, srcdir if srcdir else None) for c in corners]

    if workers <= 1:
        results = [_worker(t) for t in tasks]
    else:
        with Pool(workers) as pool:
            results = pool.map(_worker, tasks)

    failed = [corners[i] for i, r in enumerate(results) if not r]
    if failed:
        print('FAILED corners:', failed)
        sys.exit(1)
    print('All done in %.0fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
