"""
Configuration for TSMC N40 (CRN40LP) 1.1V core device sweep.
Output: one .h5 lookup table per device + corner (nch_tt.h5, pch_tt.h5, ...).

NOTE: Per the PDK doc (T-N40-CM-SP-001), do NOT additionally apply the 0.9
"shrink" (setscale scalefactor=0.9) on top of the base model file — the base
model card already carries its own scale factor. Stacking the shrink0d9 wrapper
would double-scale (0.9 x 0.9 = 0.81). So the base file is included directly
with the stat / global / {corner} sections.
"""
import numpy as np
from paths import require


def get_config(corner="tt", coarse=True):
    c = {}

    c['corner'] = corner
    base = require('PDK_MODEL_FILE')
    c['modelfile'] = (
        f'include "{base}" section=stat\n'
        f'include "{base}" section=global\n'
        f'include "{base}" section={corner}\n'
    )
    c['modelinfo'] = 'tsmcN40, BSIM4, 1.1V'
    c['temp'] = 300
    c['modeln'] = 'nch'
    c['modelp'] = 'pch'
    c['simcmd'] = (require('SPECTRE_BIN') + ' -64 '
                   f'techsweep_tsmcN40_{corner}.scs '
                   f'+log techsweep_tsmcN40_{corner}.out')
    c['rundir_base'] = f'raw_tsmcN40_{corner}'
    c['paramfile'] = f'techsweep_params_tsmcN40_{corner}.scs'
    c['sweep'] = 'sweepvds_sweepvgs-sweep'
    c['sweep_noise'] = 'sweepvds_noise_sweepvgs_noise-sweep'

    if coarse:
        c['VGS_step'] = 0.02
        c['VDS_step'] = 0.02
        c['VSB_step'] = 0.15
        c['VGS_max'] = 1.1
        c['VDS_max'] = 1.1
        c['VSB_max'] = 0.8
        c['LENGTH'] = np.array([0.04, 0.05, 0.06, 0.08, 0.1, 0.13, 0.18,
                                0.25, 0.35, 0.5, 0.7, 1.0, 1.5, 2.5, 5.0])
        c['savefilen'] = f'nch_{corner}'
        c['savefilep'] = f'pch_{corner}'
    else:
        c['VGS_step'] = 0.02
        c['VDS_step'] = 0.02
        c['VSB_step'] = 0.1
        c['VGS_max'] = 1.1
        c['VDS_max'] = 1.1
        c['VSB_max'] = 0.8
        c['LENGTH'] = np.concatenate([
            np.arange(4, 15) / 100,          # 0.04 - 0.14 um, step 0.01 (11 pts)
            np.arange(16, 51, 2) / 100,      # 0.16 - 0.50 um, step 0.02 (18 pts)
            np.arange(55, 101, 5) / 100,     # 0.55 - 1.00 um, step 0.05 (10 pts)
            np.array([1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0]),  # 9 pts
        ])
        c['savefilen'] = f'nch_{corner}'
        c['savefilep'] = f'pch_{corner}'

    c['VGS'] = np.arange(0, c['VGS_max'] + c['VGS_step'] / 2, c['VGS_step'])
    c['VDS'] = np.arange(0, c['VDS_max'] + c['VDS_step'] / 2, c['VDS_step'])
    c['VSB'] = np.arange(0, c['VSB_max'] + c['VSB_step'] / 2, c['VSB_step'])
    c['WIDTH'] = 5   # um per finger, used in the netlist (w=)
    c['NFING'] = 2   # finger count (m=)

    # 19 output variables. Stored into .h5 (float32). Names are self-explanatory:
    # FT = fT (transit frequency), GM_ID = gm/id, GAIN = gm/gds (intrinsic gain).
    c['outvars'] = [
        'ID','VT','IGD','IGS','GM','GMB','GDS',
        'CGG','CGS','CSG','CGD','CDG','CGB','CDD','CSS',
        'FT','GM_ID','GAIN','VDSAT',
    ]
    # Coeffs: ID,VT,IGD,IGS,GM,GMB,GDS,CGG,CGS,CSG,CGD,CDG,CGB,CDD,CSS,FT,GM_ID,GAIN,VDSAT
    c['n'] = [
        ('mn:ids','A',        [1, 0,0,0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0]),
        ('mn:vth','V',        [0, 1,0,0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0]),
        ('mn:igd','A',        [0, 0,1,0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0]),
        ('mn:igs','A',        [0, 0,0,1,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0]),
        ('mn:gm', 'S',        [0, 0,0,0,1,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0]),
        ('mn:gmbs','S',       [0, 0,0,0,0,1,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0]),
        ('mn:gds','S',        [0, 0,0,0,0,0,1, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0]),
        ('mn:cgg','F',        [0, 0,0,0,0,0,0, 1, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0]),
        ('mn:cgs','F',        [0, 0,0,0,0,0,0, 0,-1, 0, 0, 0, 0, 0, 0, 0,0,0,0]),
        ('mn:cgd','F',        [0, 0,0,0,0,0,0, 0, 0, 0,-1, 0, 0, 0, 0, 0,0,0,0]),
        ('mn:cgb','F',        [0, 0,0,0,0,0,0, 0, 0, 0, 0, 0,-1, 0, 0, 0,0,0,0]),
        ('mn:cdd','F',        [0, 0,0,0,0,0,0, 0, 0, 0, 0, 0, 0, 1, 0, 0,0,0,0]),
        ('mn:cdg','F',        [0, 0,0,0,0,0,0, 0, 0, 0, 0,-1, 0, 0, 0, 0,0,0,0]),
        ('mn:css','F',        [0, 0,0,0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 1, 0,0,0,0]),
        ('mn:csg','F',        [0, 0,0,0,0,0,0, 0, 0,-1, 0, 0, 0, 0, 0, 0,0,0,0]),
        ('mn:cjd','F',        [0, 0,0,0,0,0,0, 0, 0, 0, 0, 0, 0, 1, 0, 0,0,0,0]),
        ('mn:cjs','F',        [0, 0,0,0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 1, 0,0,0,0]),
        ('mn:fug','Hz',       [0, 0,0,0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 1,0,0,0]),
        ('mn:gmoverid','V',   [0, 0,0,0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,1,0,0]),
        ('mn:self_gain','rall',[0,0,0,0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,1,0]),
        ('mn:vdsat','V',      [0, 0,0,0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,1]),
    ]

    c['p'] = [
        ('mp:ids','A',        [-1,0, 0, 0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0]),
        ('mp:vth','V',        [0,-1, 0, 0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0]),
        ('mp:igd','A',        [0, 0,-1, 0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0]),
        ('mp:igs','A',        [0, 0, 0,-1,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0]),
        ('mp:gm', 'S',        [0, 0, 0, 0,1,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0]),
        ('mp:gmbs','S',       [0, 0, 0, 0,0,1,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0]),
        ('mp:gds','S',        [0, 0, 0, 0,0,0,1, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0]),
        ('mp:cgg','F',        [0, 0, 0, 0,0,0,0, 1, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0]),
        ('mp:cgs','F',        [0, 0, 0, 0,0,0,0, 0,-1, 0, 0, 0, 0, 0, 0, 0,0,0,0]),
        ('mp:cgd','F',        [0, 0, 0, 0,0,0,0, 0, 0, 0,-1, 0, 0, 0, 0, 0,0,0,0]),
        ('mp:cgb','F',        [0, 0, 0, 0,0,0,0, 0, 0, 0, 0, 0,-1, 0, 0, 0,0,0,0]),
        ('mp:cdd','F',        [0, 0, 0, 0,0,0,0, 0, 0, 0, 0, 0, 0, 1, 0, 0,0,0,0]),
        ('mp:cdg','F',        [0, 0, 0, 0,0,0,0, 0, 0, 0, 0,-1, 0, 0, 0, 0,0,0,0]),
        ('mp:css','F',        [0, 0, 0, 0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 1, 0,0,0,0]),
        ('mp:csg','F',        [0, 0, 0, 0,0,0,0, 0, 0,-1, 0, 0, 0, 0, 0, 0,0,0,0]),
        ('mp:cjd','F',        [0, 0, 0, 0,0,0,0, 0, 0, 0, 0, 0, 0, 1, 0, 0,0,0,0]),
        ('mp:cjs','F',        [0, 0, 0, 0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 1, 0,0,0,0]),
        ('mp:fug','Hz',       [0, 0, 0, 0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 1,0,0,0]),
        ('mp:gmoverid','V',   [0, 0, 0, 0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,1,0,0]),
        ('mp:self_gain','rall',[0,0, 0, 0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,1,0]),
        ('mp:vdsat','V',      [0, 0, 0, 0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,1]),
    ]

    c['outvars_noise'] = ['STH', 'SFL']
    c['n_noise'] = [('mn:id',''), ('mn:fn','')]
    c['p_noise'] = [('mp:id',''), ('mp:fn','')]

    c['_netlist_tmpl'] = (
        f'//techsweep_tsmcN40_{corner}.scs' '\n'
        '%s'
        f'include "{c["paramfile"]}"' '\n'
        'save mn:ids mn:vth mn:igd mn:igs mn:gm mn:gmbs mn:gds '
        'mn:cgg mn:cgs mn:cgd mn:cgb mn:cdd mn:cdg mn:css mn:csg mn:cjd mn:cjs '
        'mn:fug mn:gmoverid mn:self_gain mn:vdsat\n'
        'save mp:ids mp:vth mp:igd mp:igs mp:gm mp:gmbs mp:gds '
        'mp:cgg mp:cgs mp:cgd mp:cgb mp:cdd mp:cdg mp:css mp:csg mp:cjd mp:cjs '
        'mp:fug mp:gmoverid mp:self_gain mp:vdsat\n'
        'parameters gs=0 ds=0 \n'
        'vnoi     (vx  0)         vsource dc=0  \n'
        'vdsn     (vdn vx)        vsource dc=ds  \n'
        'vgsn     (vgn 0)         vsource dc=gs  \n'
        'vbsn     (vbn 0)         vsource dc=-sb \n'
        'vdsp     (vdp vx)        vsource dc=-ds \n'
        'vgsp     (vgp 0)         vsource dc=-gs \n'
        'vbsp     (vbp 0)         vsource dc=sb  \n'
        '\n'
        'mn       (vdn vgn 0 vbn) %s  l=length*1e-6 w=%de-6 m=%d \n'
        'mp       (vdp vgp 0 vbp) %s  l=length*1e-6 w=%de-6 m=%d \n'
        '\n'
        'simOptions options gmin=1e-13 reltol=1e-4 vabstol=1e-6 iabstol=1e-10 '
        'temp=%d tnom=27 rawfmt=psfbin rawfile="%s" \n'
        'sweepvds sweep param=ds start=0 stop=%.10g step=%.10g { \n'
        '   sweepvgs dc param=gs start=0 stop=%.10g step=%.10g \n'
        '}\n'
        'sweepvds_noise sweep param=ds start=0 stop=%.10g step=%.10g { \n'
        '   sweepvgs_noise noise freq=1 oprobe=vnoi param=gs start=0 stop=%.10g step=%.10g \n'
        '}\n'
    )

    return c


def write_netlist(c, raw_dir):
    netlist = c['_netlist_tmpl'] % (
        c['modelfile'],
        c['modeln'], c['WIDTH'], c['NFING'],
        c['modelp'], c['WIDTH'], c['NFING'],
        c['temp'] - 273, raw_dir,
        c['VDS_max'], c['VDS_step'],
        c['VGS_max'], c['VGS_step'],
        c['VDS_max'], c['VDS_step'],
        c['VGS_max'], c['VGS_step'],
    )
    netlist_file = f'techsweep_tsmcN40_{c["corner"]}.scs'
    with open(netlist_file, 'w') as fid:
        fid.write(netlist)
    return netlist_file
