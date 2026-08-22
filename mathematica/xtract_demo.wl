(*
  xtract_demo.wl
  Mathematica port of lookup_funs/xtract_demo.py — EKV extraction demo.
  Reproduces:
    (a) drain current density  : lookup vs reconstructed basic EKV
    (b) transconductance efficiency : lookup vs reconstructed
    (c) percentage difference between the two

  Run:  wolframscript -script xtract_demo.wl
*)

dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{dir, "tsmcN40_lookup.wl"}]];
Get[FileNameJoin[{dir, "ekv_extract.wl"}]];

file = "D:\\tsmcN40_lookup\\nch_tt.h5";
nch = LoadTsmcN40[file];

(* N40 example point *)
L = 0.1;      (* um *)
VDS = 0.55;   (* V, VDD/2 *)
VSB = 0.0;    (* V *)
rho = 0.6;
TEMP = 300.0;

y = XTRACT[nch, L, VDS, VSB, rho, TEMP];

Print["========================================"];
Print["XTRACT output (L = 0.10 um, VDS = 0.55 V)"];
Print["========================================"];
Print["  VDS           = ", y[[1]], " V"];
Print["  n             = ", y[[2]]];
Print["  VT            = ", y[[3]], " V"];
Print["  JS            = ", y[[4]]*1.*^6, " uA/um"];
Print["  dn/dVDS       = ", y[[5]]];
Print["  dVT/dVDS      = ", y[[6]]];
Print["  dlogJS/dVDS   = ", y[[7]]];
Print["  d2n/dVDS2     = ", y[[8]]];
Print["  d2VT/dVDS2    = ", y[[9]]];
Print["  d2logJS/dVDS2 = ", y[[10]]];
Print["========================================"];

(* Reconstruct the basic EKV model *)
k = 1.380649*^-23; qe = 1.602176634*^-19;
UT = k*TEMP/qe;
n = y[[2]]; VT = y[[3]]; JS = y[[4]];

q = 10.^Subdivide[-3, 1, 199];
i = q^2 + q;
VP = UT*(2*(q - 1) + Log[q]);
VGSekv = n*VP + VT;
JDekv = i*JS;
gmIDekv = 1/(n*UT*(1 + q));

(* Original lookup-table data *)
vgs = nch["VGS"];
fid = N40Interpolant[nch, "ID"];
fgmid = N40Interpolant[nch, "GM_ID"];
JDLut = Map[fid[VSB, VDS, #, L] &, vgs]/nch["W"];
gmIDLut = Map[fgmid[VSB, VDS, #, L] &, vgs];

(* Reference point on reconstructed curves *)
ref = rho*Max[gmIDLut];
fref = Interpolation[Transpose[{Reverse[gmIDekv], Reverse[VGSekv]}],
  InterpolationOrder -> 1];
VGSr = fref[ref];
fJDekv = Interpolation[Transpose[{VGSekv, JDekv}], InterpolationOrder -> 1];
JDr = fJDekv[VGSr];

(* (c) percentage difference  (np.interp clamps outside the EKV VGS range) *)
clampEval[f_, xs_List, lo_, hi_] := Map[f[Clip[#, {lo, hi}]] &, xs];
JDekvOnLut = clampEval[fJDekv, vgs, Min[VGSekv], Max[VGSekv]];
Dpct = Table[
  If[JDLut[[k2]] > 0, (JDekvOnLut[[k2]] - JDLut[[k2]])/JDLut[[k2]]*100.,
   Indeterminate],
  {k2, Length[vgs]}];

p1 = ListLogPlot[Transpose[{vgs, Abs[JDLut]*1.*^6}], Joined -> True,
  PlotStyle -> {Blue, Thick}, PlotRange -> All, AxesLabel -> {"VGS (V)", "JD (uA/um)"},
  PlotLabel -> "(a) Drain Current Density"];
p1 = Show[p1, ListLogPlot[Transpose[{VGSekv, JDekv*1.*^6}], Joined -> True,
  PlotStyle -> {Red, Dashed}],
  ListPlot[{{VGSr, JDr*1.*^6}}, PlotStyle -> {Black, PointSize[0.02]}]];
p1 = Show[p1, GridLines -> Automatic];

p2 = ListPlot[Transpose[{vgs, gmIDLut}], Joined -> True,
  PlotStyle -> {Blue, Thick}, PlotRange -> All, AxesLabel -> {"VGS (V)", "gm/id (S/A)"},
  PlotLabel -> "(b) Transconductance Efficiency"];
p2 = Show[p2, ListPlot[Transpose[{VGSekv, gmIDekv}], Joined -> True,
  PlotStyle -> {Red, Dashed}],
  ListPlot[{{VGSr, ref}}, PlotStyle -> {Black, PointSize[0.02]}]];

p3 = ListPlot[Select[Transpose[{vgs, Dpct}], NumberQ[#[[2]]] &],
  Joined -> True, PlotStyle -> Green, PlotRange -> {-20, 20},
  AxesLabel -> {"VGS (V)", "D (%)"}, PlotLabel -> "(c) Percentage Difference"];
p3 = Show[p3, Plot[0, {x, 0, 1.1}, PlotStyle -> {Black, Thin}]];

out = "D:\\tsmcN40_lookup\\xtract_demo.png";
Export[out, GraphicsGrid[{{p1, p2, p3}}, ImageSize -> 1500]];
Print["plot saved to ", out];
Print["DONE"];
