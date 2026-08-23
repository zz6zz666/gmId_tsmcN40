(*
  xtract_demo.wl
  Mathematica port of lookup_funs/xtract_demo.py — EKV extraction demo.
  Reproduces:
    (a) drain current density  : lookup vs reconstructed basic EKV
    (b) transconductance efficiency : lookup vs reconstructed
    (c) percentage difference between the two

  Run:  wolframscript -script xtract_demo.wl
*)

dataDir = "D:\\tsmcN40_lookup";
scriptDir = FileNameJoin[{dataDir, "mathematica"}];
Get[FileNameJoin[{scriptDir, "tsmcN40_lookup.wl"}]];
Get[FileNameJoin[{scriptDir, "ekv_extract.wl"}]];

ClearAll[invq, EKVParams, EkVIdVds, EkVGds, EKVAInt];
invq[x_?NumericQ] := Module[{f},
  f[q_?NumericQ] := 2*(q - 1) + Log[q] - x;
  q /. FindRoot[f[q], {q, 0.5}]];
invq[x_List] := invq /@ x;
EKVParams[data_Association, l_?NumericQ, vds_, vsb_?NumericQ,
  rho_ : 0.6, temp_ : 300.0] := XTRACT[data, l, vds, vsb, rho, temp];
EkVIdVds[data_Association, l_?NumericQ, vds_List, vgs_?NumericQ,
  vsb_ : 0, rho_ : 0.6, temp_ : 300.0] := Module[
  {ut = kB*temp/qe, y, n, vt, js, qs, id},
  y = EKVParams[data, l, vds, vsb, rho, temp];
  n = y[[All, 2]]; vt = y[[All, 3]]; js = y[[All, 4]];
  qs = invq[((vgs - vt)/n)/ut];
  id = data["W"]*js*(qs^2 + qs);
  {id, qs, y}];
EkVGds[data_Association, l_?NumericQ, vds_List, vgs_?NumericQ,
  vsb_ : 0, rho_ : 0.6, temp_ : 300.0] := Module[
  {ut = kB*temp/qe, y, n, vt, js, qs, id, gm, x},
  y = EKVParams[data, l, vds, vsb, rho, temp];
  n = y[[All, 2]]; vt = y[[All, 3]]; js = y[[All, 4]];
  qs = invq[((vgs - vt)/n)/ut];
  id = data["W"]*js*(qs^2 + qs);
  gm = data["W"]*js/(n*ut)*qs;
  x = ((vgs - vt)/n)/ut;
  -gm*(y[[All, 6]] + x*y[[All, 5]]) + id*y[[All, 7]]];
EKVAInt[data_Association, l_?NumericQ, vds_List, vgs_?NumericQ,
  vsb_ : 0, rho_ : 0.6, temp_ : 300.0] := Module[
  {ut = kB*temp/qe, y, n, vt, qs, gmid},
  y = EKVParams[data, l, vds, vsb, rho, temp];
  n = y[[All, 2]]; vt = y[[All, 3]];
  qs = invq[((vgs - vt)/n)/ut];
  gmid = 1/(n*ut*(1 + qs));
  x = ((vgs - y[[All, 3]])/y[[All, 2]])/ut;
  1/(y[[All, 7]]/gmid - y[[All, 6]] - x*y[[All, 5]])];

file = FileNameJoin[{dataDir, "nch_tt.h5"}];
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
JDLut = Map[fid[L, #, VDS, VSB] &, vgs]/nch["W"];
gmIDLut = Map[fgmid[L, #, VDS, VSB] &, vgs];

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
  ListLogPlot[{{VGSr, JDr*1.*^6}}, PlotStyle -> {Black, PointSize[0.02]}]];
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

Print[p1];
Print[p2];
Print[p3];

vdsVec = Select[nch["VDS"], # >= 0.3 &];
vgDemo = 0.4;
{idEkv, qsv, params} = EkVIdVds[nch, L, vdsVec, vgDemo, VSB, rho, TEMP];
idLut = Map[fid[L, vgDemo, #, VSB] &, vdsVec];
p4 = ListLinePlot[{Transpose[{vdsVec, idLut*10^6}],
     Transpose[{vdsVec, idEkv*10^6}]},
   PlotStyle -> {Blue, {Red, Dashed}}, PlotLegends -> {"Lookup", "EKV"},
   AxesLabel -> {"VDS (V)", "ID (uA)"}, PlotLabel -> "ID vs VDS",
   GridLines -> Automatic];

gdsEkv = EkVGds[nch, L, vdsVec, vgDemo, VSB, rho, TEMP];
gdsLut = Map[N40Interpolant[nch, "GDS"][L, vgDemo, #, VSB] &, vdsVec];
p5 = ListLinePlot[{Transpose[{vdsVec, gdsLut*10^6}],
     Transpose[{vdsVec, gdsEkv*10^6}]},
   PlotStyle -> {Blue, {Red, Dashed}}, PlotLegends -> {"Lookup", "EKV"},
   AxesLabel -> {"VDS (V)", "gds (uS)"}, PlotLabel -> "gds vs VDS",
   GridLines -> Automatic];

vdsDer = Select[nch["VDS"], 0.5 <= # <= 0.7 &];
der = XTRACT[nch, L, vdsDer, VSB, rho, TEMP];
iDer = First@Ordering[Abs[vdsDer - 0.6], 1];
 nDer = der[[iDer, 2]]; vtDer = der[[iDer, 3]];
 dnDer = der[[iDer, 5]]; svt = der[[iDer, 6]]; sis = der[[iDer, 7]];
gmidLut = Map[fgmid[L, #, 0.6, VSB] &, vgs];
gdsIdLut = Map[N40Interpolant[nch, "GDS"][L, #, 0.6, VSB] &, vgs]/
   Map[fid[L, #, 0.6, VSB] &, vgs];
qGain = Map[Max[1/(nDer*UT*#) - 1, $MachineEpsilon] &, gmidLut];
xGain = 2*(qGain - 1) + Log[qGain];
gdsIdEkv = -gmidLut*(svt + xGain*dnDer) + sis;
p6 = Show[
    ListPlot[Transpose[{gmidLut, gdsIdLut}], Joined -> True,
     PlotStyle -> Blue],
    ListPlot[Transpose[{gmidLut, gdsIdEkv}], Joined -> True,
     PlotStyle -> {Red, Dashed}],
   Frame -> True, FrameLabel -> {"gm/ID (S/A)", "gds/ID (1/V)"},
   PlotLabel -> "gds/ID vs gm/ID", GridLines -> Automatic];

gainLut = Map[N40Interpolant[nch, "GAIN"][L, #, 0.6, VSB] &, vgs];
gainEkv = 1/(sis/gmidLut - svt - xGain*dnDer);
p7 = Show[
   ListLogPlot[Transpose[{gmidLut, gainLut}], Joined -> True,
    PlotStyle -> Blue],
   ListLogPlot[Transpose[{gmidLut, gainEkv}], Joined -> True,
    PlotStyle -> {Red, Dashed}],
   Frame -> True, FrameLabel -> {"gm/ID (S/A)", "Intrinsic gain"},
   PlotLabel -> "Intrinsic gain", GridLines -> Automatic];

Print[p4];
Print[p5];
Print[p6];
Print[p7];
Print["DONE"];
