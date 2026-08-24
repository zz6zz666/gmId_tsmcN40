(*
  tsmcN40_lookup.wl
  =================
  Mathematica API for TSMC N40 gm/ID lookup tables (.h5).

  Files are produced by extract_new.py, e.g.  nch_tt.h5, pch_tt.h5, ...
  Layout (HDF5):
      metadata datasets : CORNER, DEVICE, INFO, TEMP, W, NFING
      axis datasets     : L, VGS, VDS, VSB
      4-D variable data : ID VT IGD IGS GM GMB GDS CGG CGS CSG CGD CDG CGB
                          CDD CSS VDSAT  (+ noise STH SFL)
      derived on the fly (textbook naming): GM_ID = GM/ID,
                          GM_CGG = GM/CGG (=> fT = GM_CGG/(2 Pi)),
                          GM_GDS = GM/GDS (=> intrinsic gain), ID_W = ID/W,
                          ... (see n40Array)
      array order       : (L, VGS, VDS, VSB)   (* official Murmann layout *)

  Usage
  -----
    (* load a single device+corner *)
    data = LoadTsmcN40["D:\\tsmcN40_lookup\\nch_tt.h5"];

    (* Appendix/Python-compatible public API *)
    lookup[data, "GM_ID", "VGS", 0.6, "VDS", 0.7, "L", 0.04]
    lookup[data, "GM_CGG", "GM_ID", {5, 10, 15}, "VDS", 0.7]
    lookupVGS[data, "GM_ID", 15, "VDS", 0.7, "L", 0.04]

    (* quick grid values are also directly available *)
    data["VDSAT"];                  (* stored 4-D array *)
    data["GM"]/data["ID"];          (* = GM_ID, computed on the fly *)
    data["GM"]/data["CGG"];         (* gm/Cgg ratio *)

    (* fix L/VDS/VSB and plot gm/id vs VGS *)
    cur = SliceVGS[data, "GM_ID", 0.04, 0.7, 0.0];
    ListLinePlot[cur, AxesLabel -> {"VGS (V)", "gm/id (S/A)"}]
*)

ClearAll[LoadTsmcN40, lookup, lookupVGS, N40Interpolant, SliceVGS,
  n40ParseArgs, n40Array, n40GridValues, n40GridPointValues, n40Squeeze,
  n40DirectLookup, n40DirectLookupScalarVector, n40DirectLookupPointVectors,
  n40DirectLookupScalarPair, n40DirectLookupPairBatch, n40PackedRealArrayQ, n40Compiled4D,
  n40Compiled4DPair, n40Compiled4DPairBatch, n40CompiledInputs,
  n40CompiledValues, n40CompiledPairValues, n40CompiledPairBatchValues,
  n40CompiledPchipBatch, n40CurveLookupBatch,
  n40Bracket, n40LinearValue, n40PchipEndpoint, n40PchipSlopes,
  n40PchipValue, n40CurveLookup, n40LookupVGSOne];

(* Import every dataset of an .h5 file into an Association keyed by name.
   Quiet suppresses benign LibraryFunction::fpexc warnings raised when an
   array contains NaN (the noise datasets STH/SFL are NaN where the device
   is off). *)
LoadTsmcN40[file_String] := Module[{names},
  names = Import[file, "Datasets"];
  Association[StringTrim[#, "/"] -> Quiet[Import[file, {"Datasets", #}]] & /@ names]
];

(* Selective loading avoids importing every large 4-D array when a workflow
   needs only a few variables, for example a GM_ID corner comparison. *)
LoadTsmcN40[file_String, keys_List] := Association[
  (ToUpperCase[#] -> Quiet[Import[file, {"Datasets", "/" <> ToUpperCase[#]}]]) & /@ keys
];

(* Build a lightweight callable over (L, VGS, VDS, VSB). Multidimensional
   lookup is linear, matching the Appendix implementation. *)
N40Interpolant[data_Association, varKey_String] := Module[
  {axes, arr},
  axes = data /@ {"L", "VGS", "VDS", "VSB"};
  arr = n40Array[data, varKey];
  Function[{l, vgs, vds, vsb}, n40LinearValue[axes, arr, {l, vgs, vds, vsb}]]
];

(* Convert the Matlab/Python-style alternating name/value arguments to an
   Association. Parameter names are case-insensitive. *)
n40ParseArgs[args_List] := Module[{pairs},
  If[OddQ[Length[args]], Return[$Failed]];
  pairs = Partition[args, 2];
  If[!And @@ (StringQ[First[#]] & /@ pairs), Return[$Failed]];
  Association[(ToUpperCase[First[#]] -> Last[#]) & /@ pairs]
];

(* Resolve a stored variable or a ratio such as GM_CGG or GM_ID.
   Textbook naming: GM_ID = GM/ID, GM_CGG = GM/CGG (=> fT = GM_CGG/(2 Pi)),
   GM_GDS = GM/GDS (=> intrinsic gain), ID_W = ID/W. VDSAT is stored. *)
n40Array[data_Association, name_String] := Module[
  {key = ToUpperCase[name], parts, numerator, denominator},
  If[KeyExistsQ[data, key], Return[data[key]]];
  parts = StringSplit[key, "_"];
  If[Length[parts] != 2, Return[$Failed]];
  numerator = Lookup[data, parts[[1]], Missing["NotFound"]];
  denominator = If[parts[[2]] == "W", Lookup[data, "W", Missing["NotFound"]],
    Lookup[data, parts[[2]], Missing["NotFound"]]];
  If[MissingQ[numerator] || MissingQ[denominator], $Failed,
    Quiet[numerator/denominator]]
];

(* Resolve only selected grid entries. This avoids constructing a full 4-D
   ratio array when a query needs one VGS sweep at fixed L/VDS/VSB. *)
n40GridValues[data_Association, name_String, il_, ig_, id_, is_] := Module[
  {key = ToUpperCase[name], stored, parts, numerator, denominator},
  stored = Lookup[data, key, Missing["NotFound"]];
  If[ArrayQ[stored, 4], Return[stored[[il, ig, id, is]]]];
  parts = StringSplit[key, "_"];
  If[Length[parts] != 2, Return[$Failed]];
  numerator = Lookup[data, parts[[1]], Missing["NotFound"]];
  denominator = Lookup[data, parts[[2]], Missing["NotFound"]];
  If[parts[[2]] == "W",
    denominator = First[Flatten[{Lookup[data, "W", Missing["NotFound"]]}]]];
  If[!ArrayQ[numerator, 4] || MissingQ[denominator], Return[$Failed]];
  Quiet[If[ArrayQ[denominator, 4],
    numerator[[il, ig, id, is]]/denominator[[il, ig, id, is]],
    numerator[[il, ig, id, is]]/denominator]]
];

(* Paired form for a path through 4-D bias space. Each row of indices selects
   one point, unlike Part with several list indices, which forms a product. *)
n40GridPointValues[data_Association, name_String, indices_List] := Module[
  {key = ToUpperCase[name], stored, parts, numerator, denominator},
  stored = Lookup[data, key, Missing["NotFound"]];
  If[ArrayQ[stored, 4], Return[Extract[stored, indices]]];
  parts = StringSplit[key, "_"];
  If[Length[parts] != 2, Return[$Failed]];
  numerator = Lookup[data, parts[[1]], Missing["NotFound"]];
  denominator = Lookup[data, parts[[2]], Missing["NotFound"]];
  If[parts[[2]] == "W",
    denominator = First[Flatten[{Lookup[data, "W", Missing["NotFound"]]}]]];
  If[!ArrayQ[numerator, 4] || MissingQ[denominator], Return[$Failed]];
  Quiet[If[ArrayQ[denominator, 4],
    Extract[numerator, indices]/Extract[denominator, indices],
    Extract[numerator, indices]/denominator]]
];

(* Constant-time eligibility check for the compiled path. PackedArrayQ avoids
   scanning a full 4-D table; one element establishes that the packed type is
   machine real rather than integer. *)
n40PackedRealArrayQ[array_, depth_Integer] :=
  Developer`PackedArrayQ[array] && Length[Dimensions[array]] == depth &&
    MachineNumberQ[Extract[array, ConstantArray[1, depth]]];

(* Numeric multilinear interpolation kernel. It targets the Wolfram VM, so no
   external C compiler is required. Column 2 marks a zero ratio denominator. *)
n40Compiled4D = Compile[{
    {lAxis, _Real, 1}, {gAxis, _Real, 1}, {dAxis, _Real, 1},
    {sAxis, _Real, 1}, {numerator, _Real, 4}, {denominator, _Real, 4},
    {points, _Real, 2}, {divide, _Integer}},
  Module[{n = Length[points], output, p, il, ig, id, is, a, b, c, d,
    l, g, vd, vs, wl, wg, wd, ws, value, valid, den, corner, weight},
    output = Table[0.0, {Length[points]}, {2}];
    For[p = 1, p <= n, p++,
      l = points[[p, 1]]; g = points[[p, 2]];
      vd = points[[p, 3]]; vs = points[[p, 4]];
      il = 1; While[il < Length[lAxis] - 1 && lAxis[[il + 1]] <= l, il++];
      ig = 1; While[ig < Length[gAxis] - 1 && gAxis[[ig + 1]] <= g, ig++];
      id = 1; While[id < Length[dAxis] - 1 && dAxis[[id + 1]] <= vd, id++];
      is = 1; While[is < Length[sAxis] - 1 && sAxis[[is + 1]] <= vs, is++];
      wl = (l - lAxis[[il]])/(lAxis[[il + 1]] - lAxis[[il]]);
      wg = (g - gAxis[[ig]])/(gAxis[[ig + 1]] - gAxis[[ig]]);
      wd = (vd - dAxis[[id]])/(dAxis[[id + 1]] - dAxis[[id]]);
      ws = (vs - sAxis[[is]])/(sAxis[[is + 1]] - sAxis[[is]]);
      value = 0.0; valid = 1.0;
      For[a = 0, a <= 1, a++, For[b = 0, b <= 1, b++,
        For[c = 0, c <= 1, c++, For[d = 0, d <= 1, d++,
          den = denominator[[il + a, ig + b, id + c, is + d]];
          If[divide == 1,
            If[den == 0.0, valid = 0.0; corner = 0.0,
              corner = numerator[[il + a, ig + b, id + c, is + d]]/den],
            corner = numerator[[il + a, ig + b, id + c, is + d]]];
          weight = If[a == 0, 1.0 - wl, wl] If[b == 0, 1.0 - wg, wg]
            If[c == 0, 1.0 - wd, wd] If[d == 0, 1.0 - ws, ws];
          value += weight corner;
        ]]]];
      output[[p, 1]] = value; output[[p, 2]] = valid;
    ];
    output],
  RuntimeOptions -> "Speed"];

(* Fused kernel for mode 3. Bracket search and 16-corner weights are shared
   while the input and output ratios are accumulated independently. *)
n40Compiled4DPair = Compile[{
    {lAxis, _Real, 1}, {gAxis, _Real, 1}, {dAxis, _Real, 1},
    {sAxis, _Real, 1}, {numerator1, _Real, 4}, {denominator1, _Real, 4},
    {divide1, _Integer}, {numerator2, _Real, 4}, {denominator2, _Real, 4},
    {divide2, _Integer}, {points, _Real, 2}},
  Module[{n = Length[points], output, p, il, ig, id, is, a, b, c, d,
    l, g, vd, vs, wl, wg, wd, ws, value1, value2, valid1, valid2,
    den1, den2, corner1, corner2, weight},
    output = Table[0.0, {Length[points]}, {4}];
    For[p = 1, p <= n, p++,
      l = points[[p, 1]]; g = points[[p, 2]];
      vd = points[[p, 3]]; vs = points[[p, 4]];
      il = 1; While[il < Length[lAxis] - 1 && lAxis[[il + 1]] <= l, il++];
      ig = 1; While[ig < Length[gAxis] - 1 && gAxis[[ig + 1]] <= g, ig++];
      id = 1; While[id < Length[dAxis] - 1 && dAxis[[id + 1]] <= vd, id++];
      is = 1; While[is < Length[sAxis] - 1 && sAxis[[is + 1]] <= vs, is++];
      wl = (l - lAxis[[il]])/(lAxis[[il + 1]] - lAxis[[il]]);
      wg = (g - gAxis[[ig]])/(gAxis[[ig + 1]] - gAxis[[ig]]);
      wd = (vd - dAxis[[id]])/(dAxis[[id + 1]] - dAxis[[id]]);
      ws = (vs - sAxis[[is]])/(sAxis[[is + 1]] - sAxis[[is]]);
      value1 = 0.0; value2 = 0.0; valid1 = 1.0; valid2 = 1.0;
      For[a = 0, a <= 1, a++, For[b = 0, b <= 1, b++,
        For[c = 0, c <= 1, c++, For[d = 0, d <= 1, d++,
          den1 = denominator1[[il + a, ig + b, id + c, is + d]];
          den2 = denominator2[[il + a, ig + b, id + c, is + d]];
          If[divide1 == 1,
            If[den1 == 0.0, valid1 = 0.0; corner1 = 0.0,
              corner1 = numerator1[[il + a, ig + b, id + c, is + d]]/den1],
            corner1 = numerator1[[il + a, ig + b, id + c, is + d]]];
          If[divide2 == 1,
            If[den2 == 0.0, valid2 = 0.0; corner2 = 0.0,
              corner2 = numerator2[[il + a, ig + b, id + c, is + d]]/den2],
            corner2 = numerator2[[il + a, ig + b, id + c, is + d]]];
          weight = If[a == 0, 1.0 - wl, wl] If[b == 0, 1.0 - wg, wg]
            If[c == 0, 1.0 - wd, wd] If[d == 0, 1.0 - ws, ws];
          value1 += weight corner1; value2 += weight corner2;
        ]]]];
      output[[p, 1]] = value1; output[[p, 2]] = valid1;
      output[[p, 3]] = value2; output[[p, 4]] = valid2;
    ];
    output],
  RuntimeOptions -> "Speed"];

(* Batch form for Cartesian mode-3 sweeps. L/VDS/VSB brackets are computed
   once per bias combination and reused for every VGS sample in that row. *)
n40Compiled4DPairBatch = Compile[{
    {lAxis, _Real, 1}, {gAxis, _Real, 1}, {dAxis, _Real, 1},
    {sAxis, _Real, 1}, {numerator1, _Real, 4}, {denominator1, _Real, 4},
    {divide1, _Integer}, {numerator2, _Real, 4}, {denominator2, _Real, 4},
    {divide2, _Integer}, {combinations, _Real, 2}, {vgs, _Real, 1}},
  Module[{nk = Length[combinations], ng = Length[vgs], output, p, q, row,
    il, ig, id, is, a, b, c, d, l, g, vd, vs, wl, wg, wd, ws, value1,
    value2, valid1, valid2, den1, den2, corner1, corner2, weight},
    output = Table[0.0, {nk ng}, {4}];
    For[p = 1, p <= nk, p++,
      l = combinations[[p, 1]]; vd = combinations[[p, 2]];
      vs = combinations[[p, 3]];
      il = 1; While[il < Length[lAxis] - 1 && lAxis[[il + 1]] <= l, il++];
      id = 1; While[id < Length[dAxis] - 1 && dAxis[[id + 1]] <= vd, id++];
      is = 1; While[is < Length[sAxis] - 1 && sAxis[[is + 1]] <= vs, is++];
      wl = (l - lAxis[[il]])/(lAxis[[il + 1]] - lAxis[[il]]);
      wd = (vd - dAxis[[id]])/(dAxis[[id + 1]] - dAxis[[id]]);
      ws = (vs - sAxis[[is]])/(sAxis[[is + 1]] - sAxis[[is]]);
      For[q = 1, q <= ng, q++,
        g = vgs[[q]];
        ig = 1; While[ig < Length[gAxis] - 1 && gAxis[[ig + 1]] <= g, ig++];
        wg = (g - gAxis[[ig]])/(gAxis[[ig + 1]] - gAxis[[ig]]);
        value1 = 0.0; value2 = 0.0; valid1 = 1.0; valid2 = 1.0;
        For[a = 0, a <= 1, a++, For[b = 0, b <= 1, b++,
          For[c = 0, c <= 1, c++, For[d = 0, d <= 1, d++,
            den1 = denominator1[[il + a, ig + b, id + c, is + d]];
            den2 = denominator2[[il + a, ig + b, id + c, is + d]];
            If[divide1 == 1,
              If[den1 == 0.0, valid1 = 0.0; corner1 = 0.0,
                corner1 = numerator1[[il + a, ig + b, id + c, is + d]]/den1],
              corner1 = numerator1[[il + a, ig + b, id + c, is + d]]];
            If[divide2 == 1,
              If[den2 == 0.0, valid2 = 0.0; corner2 = 0.0,
                corner2 = numerator2[[il + a, ig + b, id + c, is + d]]/den2],
              corner2 = numerator2[[il + a, ig + b, id + c, is + d]]];
            weight = If[a == 0, 1.0 - wl, wl] If[b == 0, 1.0 - wg, wg]
              If[c == 0, 1.0 - wd, wd] If[d == 0, 1.0 - ws, ws];
            value1 += weight corner1; value2 += weight corner2;
          ]]]];
        row = (p - 1) ng + q;
        output[[row, 1]] = value1; output[[row, 2]] = valid1;
        output[[row, 3]] = value2; output[[row, 4]] = valid2;
      ];
    ];
    output],
  RuntimeOptions -> "Speed"];

(* Resolve a stored variable or A_B ratio for the compiled kernel. Unsupported
   or unpacked arrays return $Failed and transparently fall back to the general
   Wolfram Language implementation. Packed machine NaNs propagate normally. *)
n40CompiledInputs[data_Association, name_String] := Module[
  {key = ToUpperCase[name], stored, parts, numerator, denominator, divide,
   scale = 1., w},
  stored = Lookup[data, key, Missing["NotFound"]];
  If[n40PackedRealArrayQ[stored, 4],
    numerator = stored; denominator = stored; divide = 0,
    parts = StringSplit[key, "_"];
    If[Length[parts] != 2, Return[$Failed]];
    numerator = Lookup[data, parts[[1]], Missing["NotFound"]];
    If[!n40PackedRealArrayQ[numerator, 4], Return[$Failed]];
    If[parts[[2]] == "W",
      w = Flatten[{Lookup[data, "W", Missing["NotFound"]]}];
      If[Length[w] != 1 || !NumericQ[First[w]] || First[w] == 0,
        Return[$Failed]];
      denominator = numerator; divide = 0; scale = 1./N[First[w]],
      denominator = Lookup[data, parts[[2]], Missing["NotFound"]];
      If[!n40PackedRealArrayQ[denominator, 4], Return[$Failed]];
      divide = 1]];
  {numerator, denominator, divide, scale}
];

n40CompiledValues[data_Association, name_String, points_List] := Module[
  {inputs, axes, packedPoints, raw},
  If[points == {}, Return[{}]];
  inputs = n40CompiledInputs[data, name];
  If[inputs === $Failed, Return[$Failed]];
  axes = data /@ {"L", "VGS", "VDS", "VSB"};
  If[!And @@ (n40PackedRealArrayQ[#, 1] & /@ axes), Return[$Failed]];
  packedPoints = Developer`ToPackedArray[N[points]];
  If[!n40PackedRealArrayQ[packedPoints, 2], Return[$Failed]];
  raw = n40Compiled4D[Sequence @@ axes, inputs[[1]], inputs[[2]], packedPoints,
    inputs[[3]]];
  MapThread[If[#2 == 1., inputs[[4]] #1, Indeterminate] &,
    {raw[[All, 1]], raw[[All, 2]]}]
];

n40CompiledPairValues[data_Association, name1_String, name2_String,
    points_List] := Module[{inputs1, inputs2, axes, packedPoints, raw},
  If[points == {}, Return[{{}, {}}]];
  inputs1 = n40CompiledInputs[data, name1];
  inputs2 = n40CompiledInputs[data, name2];
  If[inputs1 === $Failed || inputs2 === $Failed, Return[$Failed]];
  axes = data /@ {"L", "VGS", "VDS", "VSB"};
  If[!And @@ (n40PackedRealArrayQ[#, 1] & /@ axes), Return[$Failed]];
  packedPoints = Developer`ToPackedArray[N[points]];
  If[!n40PackedRealArrayQ[packedPoints, 2], Return[$Failed]];
  raw = n40Compiled4DPair[Sequence @@ axes,
    inputs1[[1]], inputs1[[2]], inputs1[[3]],
    inputs2[[1]], inputs2[[2]], inputs2[[3]], packedPoints];
  {MapThread[If[#2 == 1., inputs1[[4]] #1, Indeterminate] &,
     {raw[[All, 1]], raw[[All, 2]]}],
   MapThread[If[#2 == 1., inputs2[[4]] #1, Indeterminate] &,
     {raw[[All, 3]], raw[[All, 4]]}]}
];

n40CompiledPairBatchValues[data_Association, name1_String, name2_String,
    combinations_List, vgs_List] := Module[
  {inputs1, inputs2, axes, packedCombinations, packedVgs, raw},
  If[combinations == {} || vgs == {}, Return[{{}, {}}]];
  inputs1 = n40CompiledInputs[data, name1];
  inputs2 = n40CompiledInputs[data, name2];
  If[inputs1 === $Failed || inputs2 === $Failed, Return[$Failed]];
  axes = data /@ {"L", "VGS", "VDS", "VSB"};
  If[!And @@ (n40PackedRealArrayQ[#, 1] & /@ axes), Return[$Failed]];
  packedCombinations = Developer`ToPackedArray[N[combinations]];
  packedVgs = Developer`ToPackedArray[N[vgs]];
  If[!n40PackedRealArrayQ[packedCombinations, 2] ||
      !n40PackedRealArrayQ[packedVgs, 1], Return[$Failed]];
  raw = n40Compiled4DPairBatch[Sequence @@ axes,
    inputs1[[1]], inputs1[[2]], inputs1[[3]],
    inputs2[[1]], inputs2[[2]], inputs2[[3]], packedCombinations, packedVgs];
  {MapThread[If[#2 == 1., inputs1[[4]] #1, Indeterminate] &,
     {raw[[All, 1]], raw[[All, 2]]}],
   MapThread[If[#2 == 1., inputs2[[4]] #1, Indeterminate] &,
     {raw[[All, 3]], raw[[All, 4]]}]}
];

(* Sort and invert each mode-3 curve inside the VM. Output column 2 is a
   validity flag so out-of-range targets can be restored as Indeterminate. *)
n40CompiledPchipBatch = Compile[{
    {xRows, _Real, 2}, {yRows, _Real, 2}, {targets, _Real, 1},
    {branchMode, _Integer}},
  Module[{nr = Length[xRows], ng = Length[xRows[[1]]], nt = Length[targets],
    output, xs, ys, h, delta, slopes, p, j, k, i, row, start, finish,
    maxIndex, n, unique, xv, yv, target, h1, h2, d1, d2, slope, w1, w2,
    interval, u, u2, u3, value},
    output = Table[0.0, {nr nt}, {2}];
    xs = Table[0.0, {ng}]; ys = Table[0.0, {ng}];
    h = Table[0.0, {ng - 1}]; delta = Table[0.0, {ng - 1}];
    slopes = Table[0.0, {ng}];
    For[p = 1, p <= nr, p++,
      maxIndex = 1;
      For[j = 2, j <= ng, j++,
        If[xRows[[p, j]] > xRows[[p, maxIndex]], maxIndex = j]];
      start = If[branchMode == 1, maxIndex, 1];
      finish = If[branchMode == 2, maxIndex, ng];
      n = 0;
      For[j = start, j <= finish, j++,
        xv = xRows[[p, j]]; yv = yRows[[p, j]];
        If[xv == xv && yv == yv,
          n++;
          k = n;
          While[k > 1 && xs[[k - 1]] > xv,
            xs[[k]] = xs[[k - 1]]; ys[[k]] = ys[[k - 1]]; k--];
          xs[[k]] = xv; ys[[k]] = yv;
        ];
      ];
      unique = 0;
      For[j = 1, j <= n, j++,
        If[unique == 0 || xs[[j]] != xs[[unique]],
          unique++; xs[[unique]] = xs[[j]]; ys[[unique]] = ys[[j]]];
      ];
      If[unique >= 2,
        For[j = 1, j < unique, j++,
          h[[j]] = xs[[j + 1]] - xs[[j]];
          delta[[j]] = (ys[[j + 1]] - ys[[j]])/h[[j]];
        ];
        If[unique == 2,
          slopes[[1]] = delta[[1]]; slopes[[2]] = delta[[1]],
          h1 = h[[1]]; h2 = h[[2]]; d1 = delta[[1]]; d2 = delta[[2]];
          slope = ((2.0 h1 + h2) d1 - h1 d2)/(h1 + h2);
          If[Sign[slope] != Sign[d1], slope = 0.0,
            If[Sign[d1] != Sign[d2] && Abs[slope] > 3.0 Abs[d1],
              slope = 3.0 d1]];
          slopes[[1]] = slope;
          For[j = 2, j < unique, j++,
            If[delta[[j - 1]] delta[[j]] > 0.0,
              w1 = 2.0 h[[j]] + h[[j - 1]];
              w2 = h[[j]] + 2.0 h[[j - 1]];
              slopes[[j]] = (w1 + w2)/(w1/delta[[j - 1]] + w2/delta[[j]]),
              slopes[[j]] = 0.0];
          ];
          h1 = h[[unique - 1]]; h2 = h[[unique - 2]];
          d1 = delta[[unique - 1]]; d2 = delta[[unique - 2]];
          slope = ((2.0 h1 + h2) d1 - h1 d2)/(h1 + h2);
          If[Sign[slope] != Sign[d1], slope = 0.0,
            If[Sign[d1] != Sign[d2] && Abs[slope] > 3.0 Abs[d1],
              slope = 3.0 d1]];
          slopes[[unique]] = slope;
        ];
      ];
      For[k = 1, k <= nt, k++,
        row = (p - 1) nt + k; target = targets[[k]];
        If[unique >= 2 && target >= xs[[1]] && target <= xs[[unique]],
          interval = 1;
          While[interval < unique - 1 && xs[[interval + 1]] <= target,
            interval++];
          h1 = xs[[interval + 1]] - xs[[interval]];
          u = (target - xs[[interval]])/h1; u2 = u u; u3 = u2 u;
          value = (2.0 u3 - 3.0 u2 + 1.0) ys[[interval]] +
            (u3 - 2.0 u2 + u) h1 slopes[[interval]] +
            (-2.0 u3 + 3.0 u2) ys[[interval + 1]] +
            (u3 - u2) h1 slopes[[interval + 1]];
          output[[row, 1]] = value; output[[row, 2]] = 1.0,
          output[[row, 1]] = 0.0; output[[row, 2]] = 0.0];
      ];
    ];
    output],
  RuntimeOptions -> "Speed"];

n40CurveLookupBatch[xRows_List, yRows_List, targets_List, variable_String,
    method_String] := Module[
  {packedX, packedY, packedTargets, branchMode, raw, rowValid, validIndices,
   fallbackIndices, output},
  If[ToLowerCase[method] != "pchip", Return[$Failed]];
  packedTargets = Developer`ToPackedArray[N[targets]];
  If[!n40PackedRealArrayQ[packedTargets, 1], Return[$Failed]];
  rowValid = MapThread[Function[{x, y},
      And @@ MapThread[NumericQ[#1] && NumericQ[#2] &, {x, y}]],
    {xRows, yRows}];
  validIndices = Flatten[Position[rowValid, True]];
  fallbackIndices = Flatten[Position[rowValid, False]];
  output = ConstantArray[{}, Length[xRows]];
  branchMode = Which[variable == "GM_ID", 1,
    MemberQ[{"GM_CGG", "GM_CGS"}, variable], 2, True, 0];
  If[validIndices =!= {},
    packedX = Developer`ToPackedArray[N[xRows[[validIndices]]]];
    packedY = Developer`ToPackedArray[N[yRows[[validIndices]]]];
    If[!n40PackedRealArrayQ[packedX, 2] || !n40PackedRealArrayQ[packedY, 2],
      Return[$Failed]];
    raw = n40CompiledPchipBatch[packedX, packedY, packedTargets, branchMode];
    output[[validIndices]] = Partition[
      MapThread[If[#2 == 1., #1, Indeterminate] &,
        {raw[[All, 1]], raw[[All, 2]]}], Length[targets]]];
  Do[output[[i]] = n40CurveLookup[xRows[[i]], yRows[[i]], targets, variable,
      method], {i, fallbackIndices}];
  output
];

n40Squeeze[value_] := Module[{dims},
  If[!ArrayQ[value], Return[value]];
  dims = Select[Dimensions[value], # != 1 &];
  Which[dims == {}, First[Flatten[value]], dims == Dimensions[value], value,
    True, ArrayReshape[value, dims]]
];

n40Bracket[axis_List, value_?NumericQ] := Module[{n, i},
  n = Length[axis];
  If[n < 2 || value < First[axis] || value > Last[axis], Return[$Failed]];
  i = Clip[Total[Boole[Thread[axis <= value]]], {1, n - 1}];
  {i, (value - axis[[i]])/(axis[[i + 1]] - axis[[i]])}
];

n40LinearValue[axes_List, arr_, point_List] := Module[
  {brackets, index, weight},
  brackets = MapThread[n40Bracket, {axes, point}];
  If[MemberQ[brackets, $Failed], Return[Indeterminate]];
  index = brackets[[All, 1]];
  weight = brackets[[All, 2]];
  Sum[
    arr[[index[[1]] + a, index[[2]] + b, index[[3]] + c, index[[4]] + d]]*
      (If[a == 0, 1 - weight[[1]], weight[[1]]])*
      (If[b == 0, 1 - weight[[2]], weight[[2]]])*
      (If[c == 0, 1 - weight[[3]], weight[[3]]])*
      (If[d == 0, 1 - weight[[4]], weight[[4]]]),
    {a, 0, 1}, {b, 0, 1}, {c, 0, 1}, {d, 0, 1}]
];

(* Vectorized path for one VGS sweep. Only the eight L/VDS/VSB neighbor lines
   are read, and ratio arithmetic is performed only on those grid entries. *)
n40DirectLookupScalarVector[data_Association, outvar_String, l_, vgs_List,
    vds_, vsb_] := Module[
  {axes, bl, bd, bs, il, id, is, wl, wd, ws, valid, safeVgs, ig, wg,
   result, compiled, a, c, d, lower, upper, weight},
  axes = data /@ {"L", "VGS", "VDS", "VSB"};
  bl = n40Bracket[axes[[1]], l];
  bd = n40Bracket[axes[[3]], vds];
  bs = n40Bracket[axes[[4]], vsb];
  If[MemberQ[{bl, bd, bs}, $Failed],
    Return[ConstantArray[Indeterminate, Length[vgs]]]];
  {il, wl} = bl; {id, wd} = bd; {is, ws} = bs;
  valid = NumericQ[#] && First[axes[[2]]] <= # <= Last[axes[[2]]] & /@ vgs;
  safeVgs = MapThread[If[#2, #1, First[axes[[2]]]] &, {vgs, valid}];
  compiled = n40CompiledValues[data, outvar,
    Transpose[{ConstantArray[l, Length[vgs]], safeVgs,
      ConstantArray[vds, Length[vgs]], ConstantArray[vsb, Length[vgs]]}]];
  If[compiled =!= $Failed,
    Return[MapThread[If[#2, #1, Indeterminate] &, {compiled, valid}]]];
  ig = Clip[Total[Boole[Thread[axes[[2]] <= #]]],
      {1, Length[axes[[2]]] - 1}] & /@ safeVgs;
  wg = (safeVgs - axes[[2]][[ig]])/
    (axes[[2]][[ig + 1]] - axes[[2]][[ig]]);
  result = ConstantArray[0., Length[vgs]];
  Do[
    lower = n40GridValues[data, outvar, il + a, ig, id + c, is + d];
    upper = n40GridValues[data, outvar, il + a, ig + 1, id + c, is + d];
    If[lower === $Failed || upper === $Failed, Return[$Failed]];
    weight = If[a == 0, 1 - wl, wl] If[c == 0, 1 - wd, wd]
      If[d == 0, 1 - ws, ws];
    result += weight ((1 - wg) lower + wg upper),
    {a, 0, 1}, {c, 0, 1}, {d, 0, 1}];
  MapThread[If[#2, #1, Indeterminate] &, {result, valid}]
];

(* Prepare one scalar-bias VGS sweep and evaluate two variables together. *)
n40DirectLookupScalarPair[data_Association, name1_String, name2_String, l_,
    vgs_List, vds_, vsb_] := Module[
  {axes, bl, bd, bs, valid, safeVgs, pair},
  axes = data /@ {"L", "VGS", "VDS", "VSB"};
  bl = n40Bracket[axes[[1]], l];
  bd = n40Bracket[axes[[3]], vds];
  bs = n40Bracket[axes[[4]], vsb];
  If[MemberQ[{bl, bd, bs}, $Failed],
    Return[ConstantArray[Indeterminate, {2, Length[vgs]}]]];
  valid = NumericQ[#] && First[axes[[2]]] <= # <= Last[axes[[2]]] & /@ vgs;
  safeVgs = MapThread[If[#2, #1, First[axes[[2]]]] &, {vgs, valid}];
  pair = n40CompiledPairValues[data, name1, name2,
    Transpose[{ConstantArray[l, Length[vgs]], safeVgs,
      ConstantArray[vds, Length[vgs]], ConstantArray[vsb, Length[vgs]]}]];
  If[pair === $Failed, Return[$Failed]];
  Map[MapThread[If[#2, #1, Indeterminate] &, {#, valid}] &, pair]
];

(* Evaluate every bias combination and VGS point in one compiled call. The
   returned pair contains two matrices with one VGS sweep per row. *)
n40DirectLookupPairBatch[data_Association, name1_String, name2_String,
    combinations_List, vgs_List] := Module[
  {axes, combinationValid, vgsValid, safeCombinations, safeVgs, pair, rows,
   valid},
  If[combinations == {} || vgs == {}, Return[{{}, {}}]];
  axes = data /@ {"L", "VGS", "VDS", "VSB"};
  combinationValid = Function[c,
      Length[c] == 3 && And @@ MapThread[
        NumericQ[#2] && First[#1] <= #2 <= Last[#1] &,
        {axes[[{1, 3, 4}]], c}]] /@ combinations;
  vgsValid = NumericQ[#] && First[axes[[2]]] <= # <= Last[axes[[2]]] & /@ vgs;
  safeCombinations = MapThread[
    If[#2, N[#1], First /@ axes[[{1, 3, 4}]]] &,
    {combinations, combinationValid}];
  safeVgs = MapThread[If[#2, N[#1], First[axes[[2]]]] &, {vgs, vgsValid}];
  pair = n40CompiledPairBatchValues[data, name1, name2, safeCombinations, safeVgs];
  If[pair === $Failed, Return[$Failed]];
  valid = Flatten[Outer[And, combinationValid, vgsValid]];
  pair = Map[MapThread[If[#2, #1, Indeterminate] &, {#, valid}] &, pair];
  rows = Length[vgs];
  Partition[#, rows] & /@ pair
];

(* Vectorized interpolation of paired coordinates along an arbitrary path.
   Used by lookupVGS mode 2, where VGS/VDS/VSB change together. *)
n40DirectLookupPointVectors[data_Association, outvar_String, l_, vgs_List,
    vds_List, vsb_List] := Module[
  {axes, n, coordinates, validAxes, safe, indices, weights, valid, result,
   compiled, offsets, pointIndices, values, cornerWeight, a, b, c, d},
  n = Length[vgs];
  If[Length[vds] != n || Length[vsb] != n, Return[$Failed]];
  axes = data /@ {"L", "VGS", "VDS", "VSB"};
  coordinates = {ConstantArray[l, n], vgs, vds, vsb};
  validAxes = MapThread[
    Function[{axis, coordinateValues},
      NumericQ[#] && First[axis] <= # <= Last[axis] & /@ coordinateValues],
    {axes, coordinates}];
  safe = MapThread[
    Function[{axis, coordinateValues, flags},
      MapThread[If[#2, #1, First[axis]] &, {coordinateValues, flags}]],
    {axes, coordinates, validAxes}];
  valid = And @@@ Transpose[validAxes];
  compiled = n40CompiledValues[data, outvar, Transpose[safe]];
  If[compiled =!= $Failed,
    Return[MapThread[If[#2, #1, Indeterminate] &, {compiled, valid}]]];
  indices = MapThread[
    Function[{axis, coordinateValues},
      Clip[Total[Boole[Thread[axis <= #]]], {1, Length[axis] - 1}] & /@
        coordinateValues],
    {axes, safe}];
  weights = MapThread[
    Function[{axis, coordinateValues, indexValues},
      (coordinateValues - axis[[indexValues]])/
        (axis[[indexValues + 1]] - axis[[indexValues]])],
    {axes, safe, indices}];
  result = ConstantArray[0., n];
  Do[
    offsets = {a, b, c, d};
    pointIndices = Transpose[MapThread[Plus, {indices, offsets}]];
    values = n40GridPointValues[data, outvar, pointIndices];
    If[values === $Failed, Return[$Failed]];
    cornerWeight = Times @@ MapThread[
      If[#2 == 0, 1 - #1, #1] &, {weights, offsets}];
    result += cornerWeight values,
    {a, 0, 1}, {b, 0, 1}, {c, 0, 1}, {d, 0, 1}];
  MapThread[If[#2, #1, Indeterminate] &, {result, valid}]
];

n40DirectLookup[data_Association, outvar_String, params_Association] := Module[
  {axes, query, arr, dims, values},
  axes = {data["L"], data["VGS"], data["VDS"], data["VSB"]};
  query = Flatten[{Lookup[params, #]}] & /@ {"L", "VGS", "VDS", "VSB"};
  dims = Length /@ query;
  If[Length[query[[1]]] == 1 && Length[query[[3]]] == 1 &&
      Length[query[[4]]] == 1,
    values = n40DirectLookupScalarVector[data, outvar, First[query[[1]]],
      query[[2]], First[query[[3]]], First[query[[4]]]];
    If[values === $Failed, Return[$Failed]];
    Return[n40Squeeze[ArrayReshape[values, dims]]]];
  arr = n40Array[data, outvar];
  If[arr === $Failed, Return[$Failed]];
  values = n40LinearValue[axes, arr, #] & /@ Tuples[query];
  n40Squeeze[ArrayReshape[values, dims]]
];

n40PchipEndpoint[h1_, h2_, delta1_, delta2_] := Module[{slope},
  slope = ((2*h1 + h2)*delta1 - h1*delta2)/(h1 + h2);
  Which[
    Sign[slope] != Sign[delta1], 0.,
    Sign[delta1] != Sign[delta2] && Abs[slope] > 3*Abs[delta1], 3*delta1,
    True, slope]
];

n40PchipSlopes[x_List, y_List] := Module[
  {n = Length[x], h, delta, slopes, k, w1, w2},
  h = Differences[x];
  delta = Differences[y]/h;
  If[n == 2, Return[{First[delta], First[delta]}]];
  slopes = ConstantArray[0., n];
  slopes[[1]] = n40PchipEndpoint[h[[1]], h[[2]], delta[[1]], delta[[2]]];
  slopes[[-1]] = n40PchipEndpoint[h[[-1]], h[[-2]], delta[[-1]], delta[[-2]]];
  For[k = 2, k <= n - 1, k++,
    If[delta[[k - 1]]*delta[[k]] > 0,
      w1 = 2*h[[k]] + h[[k - 1]];
      w2 = h[[k]] + 2*h[[k - 1]];
      slopes[[k]] = (w1 + w2)/(w1/delta[[k - 1]] + w2/delta[[k]])]
  ];
  slopes
];

n40PchipValue[pairs_List, targets_List] := Module[
  {x, y, slopes, derivativeData, interpolant, xmin, xmax},
  x = pairs[[All, 1]]; y = pairs[[All, 2]];
  slopes = n40PchipSlopes[x, y];
  derivativeData = MapThread[{{#1}, #2, #3} &, {x, y, slopes}];
  interpolant = Interpolation[derivativeData, Method -> "Hermite"];
  xmin = First[x]; xmax = Last[x];
  If[NumericQ[#] && xmin <= # <= xmax, Quiet[interpolant[#]], Indeterminate] & /@
    targets
];

n40PchipValue[pairs_List, target_?NumericQ] :=
  First[n40PchipValue[pairs, {target}]];

(* Final one-dimensional interpolation used by cross-lookup and lookupVGS. *)
n40CurveLookup[x_List, y_List, target_, variable_String, method_String] := Module[
  {pairs, idx, normalizedMethod, order, targets, scalarTarget, result,
   interpolant},
  scalarTarget = !ListQ[target];
  targets = Flatten[{target}];
  pairs = Select[Transpose[{x, y}], NumericQ[#[[1]]] && NumericQ[#[[2]]] &];
  If[Length[pairs] < 2, Return[If[scalarTarget, Indeterminate,
    ConstantArray[Indeterminate, Length[targets]]]]];
  idx = First@Ordering[pairs[[All, 1]], -1];
  pairs = Which[
    variable == "GM_ID", pairs[[idx ;;]],
    MemberQ[{"GM_CGG", "GM_CGS"}, variable], pairs[[;; idx]],
    True, pairs];
  pairs = DeleteDuplicatesBy[SortBy[pairs, First], First];
  If[Length[pairs] < 2, Return[If[scalarTarget, Indeterminate,
    ConstantArray[Indeterminate, Length[targets]]]]];
  normalizedMethod = ToLowerCase[method];
  result = If[normalizedMethod == "pchip",
    n40PchipValue[pairs, targets],
    order = If[normalizedMethod == "linear", 1, Min[3, Length[pairs] - 1]];
    interpolant = Interpolation[pairs, InterpolationOrder -> order];
    If[!NumericQ[#] || # < First[pairs][[1]] || # > Last[pairs][[1]],
      Indeterminate, Quiet[interpolant[#]]] & /@ targets];
  If[scalarTarget, First[result], result]
];

lookup[data_Association, outvar_String, args___] := Module[
  {params, defaults, method, warning, nonCoordinates, inputVar, targets,
   l, vds, vsb, combinations, rows, result, pair, vgs},
  params = n40ParseArgs[{args}];
  If[params === $Failed, Return[$Failed]];
  defaults = <|"L" -> Min[data["L"]], "VGS" -> data["VGS"],
    "VDS" -> Max[data["VDS"]]/2, "VSB" -> 0.,
    "METHOD" -> "pchip", "WARNING" -> "off"|>;
  params = Join[defaults, params];
  method = ToString[params["METHOD"]];
  warning = ToLowerCase[ToString[params["WARNING"]]];
  nonCoordinates = Complement[Keys[params], Keys[defaults]];

  If[nonCoordinates == {}, Return[n40DirectLookup[data, outvar, params]]];
  If[Length[nonCoordinates] != 1, Return[$Failed]];

  inputVar = First[nonCoordinates];
  If[!StringContainsQ[outvar, "_"] || !StringContainsQ[inputVar, "_"],
    Return[$Failed]];
  targets = Flatten[{params[inputVar]}];
  l = Flatten[{params["L"]}];
  vds = Flatten[{params["VDS"]}];
  vsb = Flatten[{params["VSB"]}];
  combinations = Tuples[{l, vds, vsb}];
  vgs = Flatten[{params["VGS"]}];
  pair = n40DirectLookupPairBatch[data, inputVar, outvar, combinations, vgs];
  rows = If[pair === $Failed,
    Table[Module[{x, y},
        x = Flatten[{n40DirectLookup[data, inputVar,
            <|"L" -> c[[1]], "VGS" -> params["VGS"], "VDS" -> c[[2]],
              "VSB" -> c[[3]]|>]}];
        y = Flatten[{n40DirectLookup[data, outvar,
            <|"L" -> c[[1]], "VGS" -> params["VGS"], "VDS" -> c[[2]],
              "VSB" -> c[[3]]|>]}];
        n40CurveLookup[x, y, targets, inputVar, method]], {c, combinations}],
    result = n40CurveLookupBatch[pair[[1]], pair[[2]], targets, inputVar, method];
    If[result === $Failed,
      MapThread[n40CurveLookup[#1, #2, targets, inputVar, method] &,
        {pair[[1]], pair[[2]]}], result]];
  result = n40Squeeze[ArrayReshape[Flatten[rows], {Length[l], Length[vds],
      Length[vsb], Length[targets]}]];
  If[warning == "on" && !FreeQ[result, Indeterminate],
    Print["lookup warning: ", inputVar, " input out of range (Indeterminate returned)."]];
  result
];

n40LookupVGSOne[data_Association, targetVar_String, targets_List, l_, vds_, vsb_, method_] :=
  Module[{x, vgs = data["VGS"]},
    x = Flatten[{n40DirectLookup[data, targetVar,
        <|"L" -> l, "VGS" -> vgs, "VDS" -> vds, "VSB" -> vsb|>]}];
    n40CurveLookup[x, vgs, targets, targetVar, method]
  ];

lookupVGS[data_Association, args___] := Module[
  {params, method, warning, targetVar, targets, l, vds, vsb, vdb, vgb, rows,
   sourceBias, vgsSearch, vdsSearch, x, valid, combinations, vectorCount,
   result, step},
  params = n40ParseArgs[{args}];
  If[params === $Failed, Return[$Failed]];
  method = ToString[Lookup[params, "METHOD", "pchip"]];
  warning = ToLowerCase[ToString[Lookup[params, "WARNING", "off"]]];
  targetVar = Which[KeyExistsQ[params, "ID_W"], "ID_W",
    KeyExistsQ[params, "GM_ID"], "GM_ID", True, Return[$Failed]];
  targets = Flatten[{params[targetVar]}];
  l = Flatten[{Lookup[params, "L", Min[data["L"]]]}];

  If[KeyExistsQ[params, "VDB"] || KeyExistsQ[params, "VGB"],
    If[!KeyExistsQ[params, "VDB"] || !KeyExistsQ[params, "VGB"], Return[$Failed]];
    If[Length[l] != 1, Return[$Failed]];
    vdb = params["VDB"]; vgb = params["VGB"];
    If[!NumericQ[vdb] || !NumericQ[vgb], Return[$Failed]];
    step = Abs[data["VGS"][[2]] - data["VGS"][[1]]];
    sourceBias = Range[Max[data["VSB"]], Min[data["VSB"]], -step];
    vgsSearch = vgb - sourceBias;
    vdsSearch = vdb - sourceBias;
    x = n40DirectLookupPointVectors[data, targetVar, First[l], vgsSearch,
      vdsSearch, sourceBias];
    valid = Select[Transpose[{x, vgsSearch}], NumericQ[First[#]] &];
    If[Length[valid] < 2, Return[ConstantArray[Indeterminate, Length[targets]]]];
    result = n40Squeeze[n40CurveLookup[valid[[All, 1]], valid[[All, 2]], targets,
        targetVar, method]];
    If[warning == "on" && targetVar == "GM_ID" && !FreeQ[result, Indeterminate],
      Print["lookupVGS: GM_ID input larger than maximum!"]];
    Return[result]
  ];

  vds = Flatten[{Lookup[params, "VDS", Max[data["VDS"]]/2]}];
  vsb = Flatten[{Lookup[params, "VSB", 0.]}];
  vectorCount = Count[{Length[targets], Length[l], Length[vds], Length[vsb]},
    n_ /; n > 1];
  If[vectorCount > 1, Return[$Failed]];
  combinations = Tuples[{l, vds, vsb}];
  rows = n40LookupVGSOne[data, targetVar, targets, #[[1]], #[[2]], #[[3]],
      method] & /@ combinations;
  result = n40Squeeze[rows];
  If[warning == "on" && targetVar == "GM_ID" && !FreeQ[result, Indeterminate],
    Print["lookupVGS: GM_ID input larger than maximum!"]];
  result
];

(* 1-D sweep of a variable vs VGS at a fixed (L, VDS, VSB). *)
SliceVGS[data_Association, varKey_String, l_, vds_, vsb_] :=
  Module[{f}, f = N40Interpolant[data, varKey];
    Table[{g, f[l, g, vds, vsb]}, {g, data["VGS"]}]
  ];
