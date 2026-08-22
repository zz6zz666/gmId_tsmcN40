#!/usr/bin/env python3
"""1.1V core device Spectre simulation runner (no data extraction)."""
import time, subprocess, sys, os
from multiprocessing import Pool
from paths import PATHS
from config_tsmcN40 import get_config, write_netlist


def run_sim_corner(args):
    corner, fine_flag, outdir, l_range = args
    import tempfile
    os.makedirs(PATHS['TMPDIR'], exist_ok=True)
    workdir = tempfile.mkdtemp(prefix='spec_work_', dir=PATHS['TMPDIR'])
    cwd0 = os.getcwd()
    os.chdir(workdir)
    try:
        return _run_sim_corner(corner, fine_flag, outdir, l_range)
    finally:
        os.chdir(cwd0)
        import shutil
        shutil.rmtree(workdir, ignore_errors=True)


def _run_sim_corner(corner, fine_flag, outdir, l_range):
    from config_tsmcN40 import get_config, write_netlist
    c = get_config(corner, coarse=not fine_flag)

    nL = len(c['LENGTH'])
    nVSB = len(c['VSB'])
    l_start, l_end = l_range
    l_start = max(0, l_start)
    l_end = min(nL, l_end)

    c['rundir_base'] = os.path.join(outdir, c['rundir_base'])
    os.makedirs(c['rundir_base'], exist_ok=True)

    n_jobs = (l_end - l_start) * nVSB
    print()
    print("===== Corner: " + corner + " =====")
    print("  L range %d..%d (%.3f..%.3f um) x %d VSB = %d jobs" %
          (l_start, l_end - 1, c['LENGTH'][l_start], c['LENGTH'][l_end - 1],
           nVSB, n_jobs))

    t_start = time.time()
    skipped = []
    for i in range(l_start, l_end):
        for j in range(nVSB):
            raw_dir = "%s/L%03d_%.3fum_VSB%03d_%+.2fV.raw" % (
                c['rundir_base'], i, c['LENGTH'][i], j, c['VSB'][j])

            write_netlist(c, raw_dir)

            with open(c['paramfile'], 'w') as f:
                f.write('parameters length = %.10g\n' % c['LENGTH'][i])
                f.write('parameters sb = %.10g\n' % c['VSB'][j])

            t0 = time.time()
            print(' [%s] L=%.3fum VSB=%+.2fV ... ' %
                  (c['corner'], c['LENGTH'][i], c['VSB'][j]),
                  end='', flush=True)

            result = subprocess.run(c['simcmd'], shell=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    timeout=300)
            if result.returncode != 0:
                print('SKIPPED (crash)')
                skipped.append('L=%.3fum VSB=%+.2fV' % (c['LENGTH'][i], c['VSB'][j]))
                continue
            print('OK (%.1fs)' % (time.time() - t0))
    if skipped:
        os.makedirs(PATHS['LOGDIR'], exist_ok=True)
        with open(os.path.join(PATHS['LOGDIR'], 'sim_skipped.log'), 'a') as fsk:
            for s in skipped:
                fsk.write('%s  %s\n' % (c['corner'], s))
        print("  [%s] SKIPPED %d VSB points: %s" % (c['corner'], len(skipped), skipped))
    print("  [%s] done in %.0fs" % (c['corner'], time.time() - t_start))
    return True


def main():
    args = sys.argv[1:]
    fine = '--fine' in args
    if fine:
        args.remove('--fine')
    outdir = ''
    if '--outdir' in args:
        idx = args.index('--outdir')
        args.pop(idx)
        outdir = args.pop(idx) if idx < len(args) else ''
    l_range = None
    if '--L-range' in args:
        idx = args.index('--L-range')
        args.pop(idx)
        if idx < len(args):
            s = int(args.pop(idx))
            e = int(args.pop(idx)) if idx < len(args) else s + 1
            l_range = (s, e)
    corners = args if args else ['tt', 'ff', 'ss', 'fs', 'sf']

    if l_range is None:
        from config_tsmcN40 import get_config
        c = get_config(corners[0], coarse=not fine)
        l_range = (0, len(c['LENGTH']))

    if outdir:
        print('Output dir:', outdir)
    t0 = time.time()
    with Pool(len(corners)) as pool:
        results = pool.map(run_sim_corner, [(c, fine, outdir, l_range) for c in corners])
    failed = [corners[i] for i, r in enumerate(results) if not r]
    if failed:
        print('FAILED corners:', failed)
    print('Total: %.0fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
