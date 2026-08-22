(*
  tsmcN40_lookup.wl
  =================
  Mathematica loader + interpolator for TSMC N40 gm/ID lookup tables (.h5).

  Files are produced by extract_new.py, e.g.  nch_tt.h5, pch_tt.h5, ...
  Layout (HDF5):
      metadata datasets : CORNER, DEVICE, INFO, TEMP, W, NFING
      axis datasets     : L, VGS, VDS, VSB
      4-D variable data : ID VT IGD IGS GM GMB GDS CGG CGS CSG CGD CDG CGB
                          CDD CSS FT GM_ID GAIN VDSAT  (+ noise STH SFL)
      array order       : (VSB, VDS, VGS, L)

  Usage
  -----
    (* load a single device+corner *)
    data = LoadTsmcN40["D:\\tsmcN40_lookup\\nch_tt.h5"];

    (* 4-D interpolant over (VSB, VDS, VGS, L), order-3 spline *)
    f = N40Interpolant[data, "GM_ID"];
    f[0.0, 0.7, 0.6, 0.1]           (* gm/id at that bias point *)

    (* quick grid values are also directly available *)
    data["GM_ID"];                  (* (VSB,VDS,VGS,L) numeric array *)
    data["GM"]/data["ID"];          (* same as "GM_ID" *)
    data["GM"]/data["CGG"];         (* gm/Cgg ratio *)

    (* fix VSB/VDS/L and plot gm/id vs VGS *)
    cur = SliceVGS[data, "GM_ID", 0.0, 0.7, 0.1];
    ListLinePlot[cur, AxesLabel -> {"VGS (V)", "gm/id (S/A)"}]
*)

ClearAll[LoadTsmcN40, N40Interpolant, SliceVGS];

(* Import every dataset of an .h5 file into an Association keyed by name.
   Quiet suppresses benign LibraryFunction::fpexc warnings raised when an
   array contains NaN (the noise datasets STH/SFL are NaN where the device
   is off). *)
LoadTsmcN40[file_String] := Module[{names},
  names = Import[file, "Datasets"];
  Association[StringTrim[#, "/"] -> Quiet[Import[file, {"Datasets", #}]] & /@ names]
];

(* Build a 4-D InterpolatingFunction over (VSB, VDS, VGS, L). *)
N40Interpolant[data_Association, varKey_String] := Module[
  {vsb, vds, vgs, l, arr},
  vsb = data["VSB"]; vds = data["VDS"];
  vgs = data["VGS"]; l  = data["L"];
  arr = data[varKey];
  ListInterpolation[arr, {vsb, vds, vgs, l}, InterpolationOrder -> 3]
];

(* 1-D sweep of a variable vs VGS at a fixed (VSB, VDS, L). *)
SliceVGS[data_Association, varKey_String, vsb_, vds_, l_] :=
  Module[{f}, f = N40Interpolant[data, varKey];
    Table[{g, f[vsb, vds, g, l]}, {g, data["VGS"]}]
  ];
