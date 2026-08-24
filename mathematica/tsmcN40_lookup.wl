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
  n40PackedRealArrayQ, n40Compiled4D, n40CompiledValues,
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

(* Resolve a stored variable or A_B ratio for the compiled kernel. Unsupported
   or unpacked arrays return $Failed and transparently fall back to the general
   Wolfram Language implementation. Packed machine NaNs propagate normally. *)
n40CompiledValues[data_Association, name_String, points_List] := Module[
  {key = ToUpperCase[name], stored, parts, numerator, denominator, divide,
   scale = 1., w, axes, packedPoints, raw},
  If[points == {}, Return[{}]];
  stored = Lookup[data, key, Missing["NotFound"]];
  If[n40PackedRealArrayQ[stored, 4],
    numerator = stored; denominator = stored; divide = 0,
    parts = StringSplit[key, "_"];
    If[Length[parts] != 2, Return[$Failed]];
    numerator = Lookup[data, parts[[1]], Missing["NotFound"]];
    If[!n40PackedRealArrayQ[numerator, 4], Return[$Failed]];
    If[parts[[2]] == "W",
      w = Flatten[{Lookup[data, "W", Missing["NotFound"]]}];
      If[Length[w] != 1 || !MachineNumberQ[First[w]] || First[w] == 0.,
        Return[$Failed]];
      denominator = numerator; divide = 0; scale = 1./First[w],
      denominator = Lookup[data, parts[[2]], Missing["NotFound"]];
      If[!n40PackedRealArrayQ[denominator, 4], Return[$Failed]];
      divide = 1]];
  axes = data /@ {"L", "VGS", "VDS", "VSB"};
  If[!And @@ (n40PackedRealArrayQ[#, 1] & /@ axes), Return[$Failed]];
  packedPoints = Developer`ToPackedArray[N[points]];
  If[!n40PackedRealArrayQ[packedPoints, 2], Return[$Failed]];
  raw = n40Compiled4D[Sequence @@ axes, numerator, denominator, packedPoints,
    divide];
  MapThread[If[#2 == 1., scale #1, Indeterminate] &,
    {raw[[All, 1]], raw[[All, 2]]}]
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
   l, vds, vsb, combinations, rows, result},
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
  rows = Table[
    Module[{x, y},
      x = Flatten[{n40DirectLookup[data, inputVar,
          <|"L" -> c[[1]], "VGS" -> params["VGS"], "VDS" -> c[[2]],
            "VSB" -> c[[3]]|>]}];
      y = Flatten[{n40DirectLookup[data, outvar,
          <|"L" -> c[[1]], "VGS" -> params["VGS"], "VDS" -> c[[2]],
            "VSB" -> c[[3]]|>]}];
      n40CurveLookup[x, y, targets, inputVar, method]
    ], {c, combinations}];
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
