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
      array order       : (L, VGS, VDS, VSB)   (* official Murmann layout *)

  Usage
  -----
    (* load a single device+corner *)
    data = LoadTsmcN40["D:\\tsmcN40_lookup\\nch_tt.h5"];

    (* 4-D interpolant over (L, VGS, VDS, VSB), order-3 spline *)
    f = N40Interpolant[data, "GM_ID"];
    f[0.04, 0.6, 0.7, 0.0]           (* gm/id at that bias point *)

    (* quick grid values are also directly available *)
    data["GM_ID"];                  (* (L,VGS,VDS,VSB) numeric array *)
    data["GM"]/data["ID"];          (* same as "GM_ID" *)
    data["GM"]/data["CGG"];         (* gm/Cgg ratio *)

    (* fix L/VDS/VSB and plot gm/id vs VGS *)
    cur = SliceVGS[data, "GM_ID", 0.04, 0.7, 0.0];
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

(* Build a 4-D InterpolatingFunction over (L, VGS, VDS, VSB). *)
N40Interpolant[data_Association, varKey_String] := Module[
  {l, vgs, vds, vsb, arr},
  l   = data["L"]; vgs = data["VGS"];
  vds = data["VDS"]; vsb = data["VSB"];
  arr = data[varKey];
  ListInterpolation[arr, {l, vgs, vds, vsb}, InterpolationOrder -> 3]
];

(* 1-D sweep of a variable vs VGS at a fixed (L, VDS, VSB). *)
SliceVGS[data_Association, varKey_String, l_, vds_, vsb_] :=
  Module[{f}, f = N40Interpolant[data, varKey];
    Table[{g, f[l, g, vds, vsb]}, {g, data["VGS"]}]
  ];
