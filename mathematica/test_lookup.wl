(*
  test_lookup.wl — verify tsmcN40_lookup.wl against a real N40 .h5 table.
  Run:  wolframscript -script test_lookup.wl
*)
Get[FileNameJoin[{DirectoryName[$InputFileName], "tsmcN40_lookup.wl"}]];

file = "D:\\tsmcN40_lookup\\nch_tt.h5";

data = LoadTsmcN40[file];
Print["datasets: ", Keys[data]];
Print["CORNER = ", data["CORNER"], "  DEVICE = ", data["DEVICE"]];
Print["axes lengths: VSB=", Length[data["VSB"]], " VDS=", Length[data["VDS"]],
  " VGS=", Length[data["VGS"]], " L=", Length[data["L"]]];
Print["ID dims = ", Dimensions[data["ID"]]];

(* 4-D interpolant at a known bias point *)
f = N40Interpolant[data, "GM_ID"];
Print["GM_ID @ (VSB=0, VDS=0.7, VGS=0.6, L=0.04) = ", f[0.0, 0.7, 0.6, 0.04]];

fid = N40Interpolant[data, "ID"];
Print["ID   @ (VSB=0, VDS=0.7, VGS=0.6, L=0.04) = ", fid[0.0, 0.7, 0.6, 0.04]];

(* 1-D VGS slice at fixed (VSB, VDS, L) *)
cur = SliceVGS[data, "GM_ID", 0.0, 0.7, 0.04];
Print["SliceVGS: n=", Length[cur], "  first=", First[cur], "  last=", Last[cur]];

Print["DONE"];
