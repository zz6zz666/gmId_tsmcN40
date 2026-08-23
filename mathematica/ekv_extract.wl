(*)
  ekv_extract.wl
  Mathematica port of the two core EKV extraction functions.

  Requires tsmcN40_lookup.wl to be loaded.
*)

ClearAll[XTRACT, XTRACT2, NPGradient, ekvCore];

kB = 1.380649*^-23;
qe = 1.602176634*^-19;

NPGradient[y_List, x_List] := Module[{n = Length[x], h},
  h = Differences[x];
  Table[
    Which[
      i == 1, (y[[2]] - y[[1]])/h[[1]],
      i == n, (y[[n]] - y[[n - 1]])/h[[n - 1]],
      True, (y[[i + 1]] - y[[i - 1]])/(x[[i + 1]] - x[[i - 1]])],
    {i, n}]
];

ekvCore[vgs_List, gmId_List, jd_List, rho_, UT_] := Module[
  {idxMax, M, n, gmIdRef, gmSide, vgSide, u, gmInc, vgInc, vgsO, jdO,
   qo, vpO, vt, io, js},
  idxMax = Ordering[gmId, -1][[1]];
  M = gmId[[idxMax]];
  If[M <= 0 || !NumberQ[M], Indeterminate,
    n = 1/(M*UT);
    gmIdRef = rho*M;
    gmSide = Reverse[gmId[[idxMax ;;]]];
    vgSide = Reverse[vgs[[idxMax ;;]]];
    u = DeleteDuplicatesBy[Transpose[{gmSide, vgSide}], #[[1]] &];
    gmInc = u[[All, 1]]; vgInc = u[[All, 2]];
    If[Length[gmInc] < 2 || gmIdRef < Min[gmInc] || gmIdRef > Max[gmInc],
      Indeterminate,
       (* Cubic interpolation keeps JS(VDS) smooth across VGS grid nodes. *)
       vgsO = Interpolation[Transpose[{gmInc, vgInc}], InterpolationOrder -> 3][gmIdRef];
      If[!NumberQ[vgsO], Indeterminate,
         jdO = Interpolation[Transpose[{vgs, jd}], InterpolationOrder -> 3][vgsO];
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

XTRACT[data_Association, L_?NumericQ, VDS_, VSB_?NumericQ,
  rho_ : 0.6, TEMP_ : 300.0] := Module[
  {UT, vdsList, vdsGrid, fid, fgmid, rows, validRows, vdsVec, nVec, vtVec,
   jsVec, vdsVec1, vdsVec2, d1n, d1vt, d1js, d2n, d2vt, d2js, interp,
   eval, out},
  UT = kB*TEMP/qe;
  fid = N40Interpolant[data, "ID"];
  fgmid = N40Interpolant[data, "GM_ID"];
  vdsList = Flatten[{VDS}];
  vdsGrid = Rest[data["VDS"]];
  rows = Table[
    Module[{vgs, jd, gmid, sel, core},
      vgs = data["VGS"];
      jd = Map[fid[L, #, vds, VSB] &, vgs]/data["W"];
      gmid = Map[fgmid[L, #, vds, VSB] &, vgs];
      sel = Select[Transpose[{vgs, jd, gmid}],
        NumberQ[#[[2]]] && NumberQ[#[[3]]] &];
      core = If[sel == {}, Indeterminate,
        ekvCore[sel[[All, 1]], sel[[All, 3]], sel[[All, 2]], rho, UT]];
      If[core === Indeterminate,
        {vds, Indeterminate, Indeterminate, Indeterminate},
        {vds, core[[1]], core[[2]], core[[3]]}]
    ],
    {vds, vdsGrid}];
  validRows = Select[rows,
    And @@ (NumberQ /@ #[[2 ;; 4]]) && #[[4]] > 0 &];
  If[Length[validRows] < 4,
    out = ConstantArray[Indeterminate, {Length[vdsList], 10}];
    out[[All, 1]] = vdsList;
    If[Length[vdsList] == 1, First[out], out],
    vdsVec = validRows[[All, 1]];
    nVec = validRows[[All, 2]];
    vtVec = validRows[[All, 3]];
    jsVec = validRows[[All, 4]];
    vdsVec1 = MovingAverage[vdsVec, 2];
    vdsVec2 = MovingAverage[vdsVec1, 2];
    d1n = Differences[nVec]/Differences[vdsVec];
    d1vt = Differences[vtVec]/Differences[vdsVec];
    d1js = Differences[Log[jsVec]]/Differences[vdsVec];
    d2n = Differences[d1n]/Differences[vdsVec1];
    d2vt = Differences[d1vt]/Differences[vdsVec1];
    d2js = Differences[d1js]/Differences[vdsVec1];
    interp[x_, y_] := Interpolation[Transpose[{x, y}], InterpolationOrder -> 3];
    eval[x_, y_, z_] := interp[x, y][Clip[z, {First[x], Last[x]}]];
    out = Join[
      Transpose[{eval[vdsVec, nVec, #] & /@ vdsList,
        eval[vdsVec, vtVec, #] & /@ vdsList,
        eval[vdsVec, jsVec, #] & /@ vdsList}],
      Transpose[{eval[vdsVec1, d1n, #] & /@ vdsList,
        eval[vdsVec1, d1vt, #] & /@ vdsList,
        eval[vdsVec1, d1js, #] & /@ vdsList,
        eval[vdsVec2, d2n, #] & /@ vdsList,
        eval[vdsVec2, d2vt, #] & /@ vdsList,
        eval[vdsVec2, d2js, #] & /@ vdsList}], 2];
    out = Join[Transpose[{vdsList}], out, 2];
    If[Length[vdsList] == 1, First[out], out]
  ]
];

XTRACT2[VGS_List, ID_, rho_ : 0.6, TEMP_ : 300.0] := Module[
  {UT, idMat, ncols, out},
  UT = kB*TEMP/qe;
  idMat = Which[
    VectorQ[ID], List /@ ID,
    MatrixQ[ID] && Length[ID] == Length[VGS], ID,
    MatrixQ[ID] && Length[First[ID]] == Length[VGS], Transpose[ID],
    True, Return[$Failed]
  ];
  ncols = Dimensions[idMat][[2]];
  out = Table[
    Module[{vgs, idv, valid, gmId, core},
      vgs = VGS;
      idv = idMat[[All, c]];
      valid = Select[Transpose[{vgs, idv}],
        NumberQ[#[[2]]] && #[[2]] > 0 &];
      If[Length[valid] < 3,
        {Indeterminate, Indeterminate, Indeterminate},
        vgs = valid[[All, 1]];
        idv = valid[[All, 2]];
        gmId = NPGradient[Log[idv], vgs];
        core = ekvCore[vgs, gmId, idv, rho, UT];
        If[core === Indeterminate,
          {Indeterminate, Indeterminate, Indeterminate}, core]
      ]
    ],
    {c, ncols}];
  If[ncols == 1, out[[1]], out]
];
