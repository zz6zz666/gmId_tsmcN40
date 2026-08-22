(*
  ekv_extract.wl
  Mathematica port of lookup_funs/ekv_extract.py — EKV parameter extraction
  (Murmann, "Systematic Design of Analog CMOS Circuits", Appendix A.1).

  Requires tsmcN40_lookup.wl (LoadTsmcN40 / N40Interpolant) to be loaded.

  XTRACT[data, L, VDS, VSB, rho, TEMP]   extract from a lookup table
  XTRACT2[VGS, ID, rho, TEMP]            extract from a raw ID(VGS) curve
*)

ClearAll[XTRACT, XTRACT2, NPGradient];

kB = 1.380649*^-23;   (* Boltzmann constant  [J/K] *)
qe = 1.602176634*^-19; (* elementary charge   [C] *)

(* Numeric gradient matching numpy.gradient (non-uniform spacing).
   Interior: central difference; edges: forward/backward difference. *)
NPGradient[y_List, x_List] := Module[{n = Length[x], h},
  h = Differences[x];
  Table[
    Which[
      i == 1, (y[[2]] - y[[1]])/h[[1]],
      i == n, (y[[n]] - y[[n - 1]])/h[[n - 1]],
      True, (y[[i + 1]] - y[[i - 1]])/(x[[i + 1]] - x[[i - 1]])],
    {i, n}]
];

(* Core: given sorted, finite (vgs, gm_id, jd), return {n, VT, JS}. *)
ekvCore[vgs_List, gmId_List, jd_List, rho_, UT_] := Module[
  {idxMax, M, n, gmIdRef, gmSide, vgSide, u, gmInc, vgInc, vgsO, jdO,
   qo, vpO, vt, io, js},
  idxMax = Ordering[gmId, -1][[1]];
  M = gmId[[idxMax]];
  If[M <= 0 || !NumberQ[M], Indeterminate,
    n = 1/(M*UT);
    gmIdRef = rho*M;
    (* right side of the peak, reversed so gm_id is increasing *)
    gmSide = Reverse[gmId[[idxMax ;;]]];
    vgSide = Reverse[vgs[[idxMax ;;]]];
    u = DeleteDuplicatesBy[Transpose[{gmSide, vgSide}], #[[1]] &];
    gmInc = u[[All, 1]]; vgInc = u[[All, 2]];
    If[Length[gmInc] < 2 || gmIdRef < Min[gmInc] || gmIdRef > Max[gmInc],
      Indeterminate,
      vgsO = Interpolation[Transpose[{gmInc, vgInc}], InterpolationOrder -> 1][gmIdRef];
      If[!NumberQ[vgsO], Indeterminate,
        jdO = Interpolation[Transpose[{vgs, jd}], InterpolationOrder -> 1][vgsO];
        If[!NumberQ[jdO] || jdO <= 0, Indeterminate,
          qo = 1/rho - 1;
          vpO = UT*(2*(qo - 1) + Log[qo]);
          vt = vgsO - n*vpO;
          io = qo^2 + qo;
          js = jdO/io;
          {n, vt, js}
        ]
      ]
    ]
  ]
];

(*
  XTRACT[data, L, VDS, VSB, rho, TEMP]
  Scalar VDS -> {VDS, n, VT, JS, dn/dVDS, dVT/dVDS, dlogJS/dVDS,
                 d2n/dVDS2, d2VT/dVDS2, d2logJS/dVDS2}
  Vector VDS -> matrix with one row per VDS value.
*)
XTRACT[data_Association, L_?NumericQ, VDS_, VSB_?NumericQ,
       rho_ : 0.6, TEMP_ : 300.0] := Module[
  {UT, vdsList, fid, fgmid, rows, nVec, vtVec, jsVec, vdsVec, derivs},
  UT = kB*TEMP/qe;
  fid = N40Interpolant[data, "ID"];
  fgmid = N40Interpolant[data, "GM_ID"];
  vdsList = Flatten[{VDS}];
  rows = Table[
    Module[{vgs, jd, gmid, sel, core},
      vgs = data["VGS"];
      jd = Map[fid[VSB, vds, #, L] &, vgs]/data["W"];
      gmid = Map[fgmid[VSB, vds, #, L] &, vgs];
      sel = Select[Transpose[{vgs, jd, gmid}],
        NumberQ[#[[2]]] && NumberQ[#[[3]]] &];
      core = If[sel == {}, Indeterminate,
        ekvCore[sel[[All, 1]], sel[[All, 3]], sel[[All, 2]], rho, UT]];
      If[core === Indeterminate,
        {vds, Indeterminate, Indeterminate, Indeterminate},
        {vds, core[[1]], core[[2]], core[[3]]}
      ]
    ],
    {vds, vdsList}];

  If[Length[rows] == 1, Join[rows[[1]], ConstantArray[0., 6]],
    Module[{ln, lvt, ljs},
      vdsVec = rows[[All, 1]];
      nVec = rows[[All, 2]];
      vtVec = rows[[All, 3]];
      jsVec = rows[[All, 4]];
      If[MemberQ[rows[[All, 2 ;;]], Indeterminate, Infinity],
        Join[rows, ConstantArray[Indeterminate, {Length[rows], 6}], 2],
        ln = NPGradient[nVec, vdsVec];
        lvt = NPGradient[vtVec, vdsVec];
        ljs = NPGradient[Log[jsVec], vdsVec];
        derivs = Join[
          Transpose[{ln, lvt, ljs}],
          Transpose[{NPGradient[ln, vdsVec], NPGradient[lvt, vdsVec],
                     NPGradient[ljs, vdsVec]}], 2];
        Join[rows, derivs, 2]
      ]
    ]
  ]
];

(*
  XTRACT2[VGS, ID, rho, TEMP]
  ID may be a 1-D curve or a 2-D matrix (each column one curve).
  Returns {n, VT, IS} for a single curve, or a matrix (rows per curve).
*)
XTRACT2[VGS_List, ID_, rho_ : 0.6, TEMP_ : 300.0] := Module[
  {UT, idMat, ncols, out},
  UT = kB*TEMP/qe;
  idMat = If[MatrixQ[ID],
    If[Length[ID] == Length[VGS], ID, Transpose[ID]],
    {ID}];
  ncols = Length[idMat];
  out = Table[
    Module[{vgs, idv, valid, gmId, core},
      vgs = VGS;
      idv = idMat[[All, c]];
      valid = Select[Transpose[{vgs, idv}],
        NumberQ[#[[2]]] && #[[2]] > 0 &];
      If[Length[valid] < 3, {Indeterminate, Indeterminate, Indeterminate},
        vgs = valid[[All, 1]];
        idv = valid[[All, 2]];
        gmId = NPGradient[Log[idv], vgs];
        core = ekvCore[vgs, gmId, idv, rho, UT];
        If[core === Indeterminate, {Indeterminate, Indeterminate, Indeterminate}, core]
      ]
    ],
    {c, ncols}];
  If[ncols == 1, out[[1]], out]
];
