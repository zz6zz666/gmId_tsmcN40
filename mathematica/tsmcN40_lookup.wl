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
  n40ParseArgs, n40Array, n40Squeeze, n40DirectLookup,
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

n40DirectLookup[data_Association, outvar_String, params_Association] := Module[
  {axes, query, arr, dims, values},
  axes = {data["L"], data["VGS"], data["VDS"], data["VSB"]};
  query = Flatten[{Lookup[params, #]}] & /@ {"L", "VGS", "VDS", "VSB"};
  arr = n40Array[data, outvar];
  If[arr === $Failed, Return[$Failed]];
  dims = Length /@ query;
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

n40PchipValue[pairs_List, target_?NumericQ] := Module[
  {x = pairs[[All, 1]], y = pairs[[All, 2]], slopes, bracket, i, t, h},
  bracket = n40Bracket[x, target];
  If[bracket === $Failed, Return[Indeterminate]];
  i = bracket[[1]]; t = bracket[[2]]; h = x[[i + 1]] - x[[i]];
  slopes = n40PchipSlopes[x, y];
  (2*t^3 - 3*t^2 + 1)*y[[i]] + (t^3 - 2*t^2 + t)*h*slopes[[i]] +
    (-2*t^3 + 3*t^2)*y[[i + 1]] + (t^3 - t^2)*h*slopes[[i + 1]]
];

(* Final one-dimensional interpolation used by cross-lookup and lookupVGS. *)
n40CurveLookup[x_List, y_List, target_, variable_String, method_String] := Module[
  {pairs, idx, normalizedMethod, order},
  pairs = Select[Transpose[{x, y}], NumericQ[#[[1]]] && NumericQ[#[[2]]] &];
  If[Length[pairs] < 2, Return[Indeterminate]];
  idx = First@Ordering[pairs[[All, 1]], -1];
  pairs = Which[
    variable == "GM_ID", pairs[[idx ;;]],
    MemberQ[{"GM_CGG", "GM_CGS"}, variable], pairs[[;; idx]],
    True, pairs];
  pairs = DeleteDuplicatesBy[SortBy[pairs, First], First];
  If[Length[pairs] < 2 || target < First[pairs][[1]] || target > Last[pairs][[1]],
    Return[Indeterminate]];
  normalizedMethod = ToLowerCase[method];
  If[normalizedMethod == "pchip", Return[n40PchipValue[pairs, target]]];
  order = If[normalizedMethod == "linear", 1, Min[3, Length[pairs] - 1]];
  Quiet[Interpolation[pairs, InterpolationOrder -> order][target]]
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
      n40CurveLookup[x, y, #, inputVar, method] & /@ targets
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
    n40CurveLookup[x, vgs, #, targetVar, method] & /@ targets
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
    x = MapThread[
      n40DirectLookup[data, targetVar,
        <|"L" -> First[l], "VGS" -> #1, "VDS" -> #2, "VSB" -> #3|>] &,
      {vgsSearch, vdsSearch, sourceBias}];
    valid = Select[Transpose[{x, vgsSearch}], NumericQ[First[#]] &];
    If[Length[valid] < 2, Return[ConstantArray[Indeterminate, Length[targets]]]];
    result = n40Squeeze[n40CurveLookup[valid[[All, 1]], valid[[All, 2]], #,
        targetVar, method] & /@ targets];
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
