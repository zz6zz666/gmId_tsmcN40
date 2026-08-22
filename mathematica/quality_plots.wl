(* ::Package:: *)

(* =====================================================================
   quality_plots.wl  \[LongDash]  N40 \:67e5\:627e\:8868\:6570\:636e\:8d28\:91cf\:5206\:6790\:4f5c\:56fe
   \:4f7f\:7528\:524d\:5148\:8fd0\:884c 0 \:53f7 Setup cell\:ff1b\:5404\:8282\:53ef\:5355\:72ec\:590d\:5236\:5230 Notebook \:4e2d\:8fd0\:884c\:3002
   \:6570\:636e\:4f4d\:7f6e: D:\tsmcN40_lookup\*.h5  (10 \:4e2a: nch/pch x tt/ff/ss/fs/sf)
   ===================================================================== *)

(* =====================================================================
   0. Setup \[LongDash] \:8f7d\:5165\:6570\:636e
   ===================================================================== *)
dataDir = "D:\\tsmcN40_lookup";
Get[FileNameJoin[{dataDir, "mathematica", "tsmcN40_lookup.wl"}]];

nch[corner_] := LoadTsmcN40[FileNameJoin[{dataDir, "nch_" <> corner <> ".h5"}]];
pch[corner_] := LoadTsmcN40[FileNameJoin[{dataDir, "pch_" <> corner <> ".h5"}]];

tt = nch["tt"];
corners = {"tt", "ff", "ss", "fs", "sf"};

(* --- \:6570\:636e\:7f51\:683c\:4fe1\:606f --- *)
vars4 = Select[Keys[tt], (ArrayQ[tt[#]] && Length[Dimensions[tt[#]]] == 4) &];
Print["===== \:6570\:636e\:7f51\:683c\:4fe1\:606f ====="];
Print["  \:5668\:4ef6: ", tt["DEVICE"], " ", tt["CORNER"], "  |  W=", tt["W"],
  " um/\:6307 x ", tt["NFING"], "  |  TEMP=", tt["TEMP"], " K"];
Print["  VSB[", Length[tt["VSB"]], "] = ", tt["VSB"]];
Print["  VDS[", Length[tt["VDS"]], "] = ", Min[tt["VDS"]], " .. ", Max[tt["VDS"]],
  " V, \:6b65\:957f ", tt["VDS"][[2]] - tt["VDS"][[1]]];
Print["  VGS[", Length[tt["VGS"]], "] = ", Min[tt["VGS"]], " .. ", Max[tt["VGS"]],
  " V, \:6b65\:957f ", tt["VGS"][[2]] - tt["VGS"][[1]]];
Print["  L[", Length[tt["L"]], "] = ", tt["L"]];
Print["  4-D \:53d8\:91cf\:4e2a\:6570 = ", Length[vars4],
  ", \:7ef4\:5ea6 (VSB,VDS,VGS,L) = ", Dimensions[tt[vars4[[1]]]]];
Print["  ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===="];

Ls = {0.04, 0.1, 0.18, 0.35, 1.0};
cols = ColorData["Rainbow"][#] & /@ Subdivide[0, 0.9, Length[Ls] - 1];
vgs = tt["VGS"];

(* \:56fa\:5b9a (VSB,VDS,L) \:6cbf VGS \:626b\:51fa\:5404 L \:7684\:4e00\:7ec4\:66f2\:7ebf *)
fam[data_, var_, vsb_, vds_, Ls_] :=
  With[{fi = N40Interpolant[data, var]},
    Table[Map[fi[vsb, vds, #, L] &, vgs], {L, Ls}]];

legend[Ls_] := LineLegend[cols, Map["L=" <> ToString[#] <> " um" &, Ls],
  LegendLabel -> "L"];

(* =====================================================================
   1. gm/ID vs VGS  \[LongDash] \:5f31\:53cd\:578b\:5cf0\:503c\:5e94 \[TildeTilde] 25~35\:ff0c\:5f3a\:53cd\:578b\:5355\:8c03\:4e0b\:964d\:3001\:65e0\:6bdb\:523a
   ===================================================================== *)
ListLinePlot[fam[tt, "GM_ID", 0, 0.55, Ls], DataRange -> {0, 1.1},
  PlotStyle -> cols, PlotLegends -> legend[Ls],
  AxesLabel -> {"VGS (V)", "gm/ID (S/A)"},
  PlotLabel -> "nch tt: gm/ID vs VGS  @ VDS=0.55V VSB=0",
  GridLines -> Automatic]

(* =====================================================================
   2. ID vs VGS\:ff08\:5bf9\:6570\:ff09 \[LongDash] \:4e9a\:9608\:533a\:5e94\:5448\:76f4\:7ebf\:3001\:65e0\:6298\:70b9/\:53cc\:5cf0\:ff1bVGS=0 \:4e3a\:5173\:65ad\:7535\:6d41
   ===================================================================== *)
ListLogPlot[fam[tt, "ID", 0, 0.55, Ls], DataRange -> {0, 1.1},
  PlotStyle -> cols, PlotLegends -> legend[Ls],
  AxesLabel -> {"VGS (V)", "ID (A)"},
  PlotLabel -> "nch tt: ID vs VGS (log)  @ VDS=0.55V VSB=0",
  GridLines -> Automatic]

(* =====================================================================
   3. ID vs VDS\:ff08\:8f93\:51fa\:7279\:6027\:65cf\:ff09 \[LongDash] \:9971\:548c\:533a\:5e94\:5e73\:5766\:3001\:65e0\:626d\:7ed3(kink)\:3001\:65e0\:5f02\:5e38\:53cd\:5f39
   ===================================================================== *)
vgsO = {0.4, 0.6, 0.8, 1.0};
With[{fi = N40Interpolant[tt, "ID"], vd = tt["VDS"]},
  ListLinePlot[
    Table[{vds, fi[0, vds, vg, 0.1]}, {vg, vgsO}, {vds, vd}],
    AxesLabel -> {"VDS (V)", "ID (A)"},
    PlotLabel -> "nch tt: ID vs VDS  L=0.1um VSB=0",
    PlotLegends -> LineLegend[Automatic, Map["VGS=" <> ToString[#] <> "V" &, vgsO]],
    GridLines -> Automatic]]

(* =====================================================================
   4. \:672c\:5f81\:589e\:76ca GAIN = gm/gds vs gm/ID
   ===================================================================== *)
ListLogLinearPlot[
  Transpose[#] & /@ Transpose[{fam[tt, "GM_ID", 0, 0.55, Ls],
     fam[tt, "GAIN", 0, 0.55, Ls]}],
  PlotStyle -> cols, PlotLegends -> legend[Ls],
  AxesLabel -> {"gm/ID (S/A)", "gm/gds (V/V)"},
  PlotLabel -> "nch tt: Intrinsic Gain  @ VDS=0.55V VSB=0",
  GridLines -> Automatic]

(* =====================================================================
   5. fT vs gm/ID
   ===================================================================== *)
ListLogLogPlot[
  Transpose[#] & /@ Transpose[{fam[tt, "GM_ID", 0, 0.55, Ls],
     fam[tt, "FT", 0, 0.55, Ls]}],
  PlotStyle -> cols, PlotLegends -> legend[Ls],
  AxesLabel -> {"gm/ID (S/A)", "fT (Hz)"},
  PlotLabel -> "nch tt: fT  @ VDS=0.55V VSB=0",
  GridLines -> Automatic]

(* =====================================================================
   6. \:7535\:5bb9 Cgg / Cgs / Cgd vs VGS  \[LongDash] \:5e94>0 \:4e14\:5e73\:6ed1\:3001\:91cf\:7ea7\:5408\:7406
   ===================================================================== *)
With[{fg = N40Interpolant[tt, "CGG"], fs = N40Interpolant[tt, "CGS"],
   fd = N40Interpolant[tt, "CGD"]},
  ListLinePlot[{
    Map[fg[0, 0.55, #, 0.1] &, vgs],
    Map[fs[0, 0.55, #, 0.1] &, vgs],
    Map[fd[0, 0.55, #, 0.1] &, vgs]},
   DataRange -> {0, 1.1}, PlotRange -> All,
   PlotStyle -> {Blue, Orange, Red},
   PlotLegends -> LineLegend[{"Cgg", "Cgs", "Cgd"}],
   AxesLabel -> {"VGS (V)", "Cap (F)"},
   PlotLabel -> "nch tt: capacitances  L=0.1um VDS=0.55V VSB=0",
   GridLines -> Automatic]]

(* =====================================================================
   7. VDSAT vs VGS  \[LongDash] \:5f3a\:53cd\:578b VDSAT \[TildeTilde] 2/(gm/ID)
   ===================================================================== *)
ListLinePlot[fam[tt, "VDSAT", 0, 0.55, Ls], DataRange -> {0, 1.1},
  PlotStyle -> cols, PlotLegends -> legend[Ls],
  AxesLabel -> {"VGS (V)", "VDSAT (V)"},
  PlotLabel -> "nch tt: VDSAT  @ VDS=0.55V VSB=0",
  GridLines -> Automatic]

(* =====================================================================
   8. \:4f53\:6548\:5e94 \[LongDash] ID vs VGS, VSB=0..0.8 \:ff08VT \:5e94\:968f |VSB| \:589e\:5927\:3001\:4e9a\:9608\:659c\:7387\:7565\:53d8\:ff09
   ===================================================================== *)
With[{fi = N40Interpolant[tt, "ID"], vsbs = {0.0, 0.2, 0.4, 0.6, 0.8}},
  ListLogPlot[
    Table[Map[fi[vsb, 0.55, #, 0.1] &, vgs], {vsb, vsbs}],
    DataRange -> {0, 1.1}, PlotStyle -> Automatic,
    PlotLegends -> LineLegend[Automatic, Map["VSB=" <> ToString[#] <> "V" &, vsbs]],
    AxesLabel -> {"VGS (V)", "ID (A)"},
    PlotLabel -> "nch tt: body effect  L=0.1um VDS=0.55V",
    GridLines -> Automatic]]

(* =====================================================================
   9. \:5de5\:827a\:89d2\:5bf9\:6bd4 \[LongDash] \:56fa\:5b9a\:504f\:7f6e\:4e0b gm/ID \:4e0e ID \:968f corner \:7684\:6563\:5e03
   (\:6ce8\:610f: \:6784\:9020\:5f85\:503c\:5f62\:5f0f\:5f97\:7528 With \:5148\:7ed1\:5b9a\:ff0c\:5426\:5219\:4f1a\:5bf9\:6bcf\:4e2a VGS \:70b9\:91cd\:590d\:52a0\:8f7d/\:91cd\:5efa\:63d2\:503c)
   ===================================================================== *)
ListLinePlot[
  Table[
    With[{fi = N40Interpolant[nch[c], "GM_ID"]},
      Map[fi[0, 0.55, #, 0.1] &, vgs]],
    {c, corners}],
  DataRange -> {0, 1.1},
  PlotStyle -> ColorData["DarkRainbow"][#] & /@ Subdivide[0, 1, Length[corners] - 1],
  PlotLegends -> LineLegend[Automatic, corners],
  AxesLabel -> {"VGS (V)", "gm/ID (S/A)"},
  PlotLabel -> "nch: gm/ID corners  L=0.1um VDS=0.55V VSB=0",
  GridLines -> Automatic]

(* =====================================================================
   10. gm/ID \:4e8c\:7ef4\:4e91\:56fe (VGS, VDS) \[LongDash] \:68c0\:67e5\:7f51\:683c\:8986\:76d6/\:8fde\:7eed\:6027
   ===================================================================== *)
With[{f = N40Interpolant[tt, "GM_ID"]},
  ListContourPlot[
    Flatten[Table[{vds, vgs, f[0, vds, vgs, 0.1]},
      {vds, 0, 1.1, 0.02}, {vgs, 0, 1.1, 0.02}], 1],
   FrameLabel -> {"VGS (V)", "VDS (V)"},
   PlotLabel -> "nch tt: gm/ID contour  L=0.1um VSB=0",
   ColorFunction -> "Rainbow", PlotLegends -> Automatic]]

(* =====================================================================
   11. \:566a\:58f0 STH / SFL vs VGS \[LongDash] \:5173\:65ad\:533a\:4e3a NaN(\:7a7a)\:3001\:5f00\:542f\:533a\:5e94\:6709\:5408\:7406\:91cf\:7ea7
   ===================================================================== *)
With[{fs = N40Interpolant[tt, "STH"], fl = N40Interpolant[tt, "SFL"]},
  ListLinePlot[{
    Map[fs[0, 0.55, #, 0.1] &, vgs],
    Map[fl[0, 0.55, #, 0.1] &, vgs]},
   DataRange -> {0, 1.1}, PlotRange -> All,
   PlotLegends -> LineLegend[{"STH (id, A^2/Hz)", "SFL (fn, A^2/Hz)"}],
   AxesLabel -> {"VGS (V)", "Noise PSD (A^2/Hz)"},
   PlotLabel -> "nch tt: noise  L=0.1um VDS=0.55V VSB=0",
   GridLines -> Automatic]]

(* =====================================================================
   12. \:6570\:636e\:8d28\:91cf\:62a5\:544a
   ===================================================================== *)
nanFrac[arr_] := N[Total[Flatten@Boole[Map[! NumericQ[#] &, arr, {-1}]]] /
    Max[1, Length[Flatten[arr]]]];
Print["=== \:6570\:636e\:8d28\:91cf\:62a5\:544a: nch_tt.h5 ==="];
Print["  NaN/\:975e\:6570\:503c\:5360\:6bd4 (\:5e94\:4e3a 0, \:566a\:58f0\:9664\:5916):"];
Table[With[{a = tt[v]}, Print["    ", v, ": ", nanFrac[a]]],
  {v, {"ID", "GM", "FT", "GAIN", "CGG", "STH", "SFL"}}];

idarr = tt["ID"];
neg = Flatten@Table[Select[Differences[idarr[[iv, id, All, il]]], # < 0 &],
    {iv, 9}, {id, 56}, {il, 48}];
Print["  ID \:6cbf VGS \:5355\:8c03\:6027: \:8d1f\:589e\:91cf\:4e2a\:6570 = ", Length[neg],
  "  \:6700\:5927\:8d1f\:589e\:91cf = ", If[Length[neg] > 0, Min[neg], 0.]];

off = Max[Flatten@Table[idarr[[iv, id, 1, il]], {iv, 9}, {id, 56}, {il, 48}]];
Print["  \:6700\:5927\:5173\:65ad\:7535\:6d41 ID@VGS=0 (\:5e94\:4e3a\:6781\:5c0f\:4e9a\:9608\:6f0f\:7535) = ", off, " A"];

Print["=== \:5b8c\:6210 ==="]



