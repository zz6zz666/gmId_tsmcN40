# TSMC N40 查找表使用指南

## 1. 查找表文件（.h5）

由 `extract_new.py` 生成，HDF5 格式，一个文件 = 器件 + 工艺角：

```text
D:\tsmcN40_lookup\
├── nch_tt.h5      (nch, typical corner)
├── nch_ff.h5      (nch, fast corner)
├── nch_fs.h5
├── nch_sf.h5
├── nch_ss.h5
├── pch_tt.h5
├── pch_ff.h5
├── pch_fs.h5
├── pch_sf.h5
└── pch_ss.h5
```

正常设计使用 `tt` corner。所有查表值都以**参考宽度归一化**：元数据 W = 10 µm。

## 2. 工艺参数概览（fine 模式）

| 参数           | 数值                                             |
| -------------- | ------------------------------------------------ |
| **电源电压**   | **1.1 V**（核心器件）                            |
| **默认 VDS**   | **0.55 V**（VDD/2）                              |
| **默认 VSB**   | 0 V                                              |
| VGS 范围       | 0 ~ 1.1 V，步长 0.02 V（56 点）                 |
| VDS 范围       | 0 ~ 1.1 V，步长 0.02 V（56 点）                 |
| VSB 范围       | 0 ~ 0.8 V，步长 0.1 V（9 点）        |
| L 范围         | 0.04 ~ 5 µm（48 点），详细值见下方               |
| 仿真宽度 W     | 10 µm                                        |
| 单变量网格     | 48 × 56 × 56 × 9（(L, VGS, VDS, VSB)）= 1,354,752 点 |

**L 详细值（48 点）**：

- 0.04 ~ 0.14（步长 0.01）：`0.04 ... 0.14`
- 0.16 ~ 0.50（步长 0.02）：`0.16 ... 0.50`
- 0.55 ~ 1.00（步长 0.05）：`0.55 ... 1.00`
- 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0

## 3. 数据结构

- 4-D 数据数组 `(L, VGS, VDS, VSB)`（与官方 Murmann 查找表布局一致），float32 + gzip 压缩。
- 变量：`ID VT IGD IGS GM GMB GDS CGG CGS CSG CGD CDG CGB CDD CSS FT GM_ID GAIN VDSAT`，外加噪声 `STH SFL`。
- 元数据数据集：`CORNER DEVICE INFO TEMP W NFING`。
- 坐标轴数据集：`L VGS VDS VSB`。
- 其余比例量（如 `GM_CGG = GM/CGG`、`ID_W = ID/W`）由原始量现算。

## 4. Python API

```python
from lookup_table import loadmat, lookup, lookupVGS

nch = loadmat(r"D:\tsmcN40_lookup\nch_tt.h5")
pch = loadmat(r"D:\tsmcN40_lookup\pch_tt.h5")

# 基础查表
gm_id = lookup(nch, 'GM_ID', 'VGS', 0.6, 'VDS', 0.7, 'VSB', 0, 'L', 0.1)
gain  = lookup(nch, 'GAIN',  'VGS', 0.6, 'VDS', 0.7, 'VSB', 0, 'L', 0.1)
ft    = lookup(nch, 'FT',    'VGS', 0.6, 'VDS', 0.7, 'VSB', 0, 'L', 0.1)

# 反查 VGS（给定 gm/id 或电流密度）
VGS = lookupVGS(nch, 'GM_ID', 15, 'VDS', 0.7, 'VSB', 0, 'L', 0.1)

# 交叉查表（gm/Cgg vs gm/id）
wt = lookup(nch, 'GM_CGG', 'GM_ID', [5, 10, 15], 'VDS', 0.7, 'VSB', 0, 'L', 0.1)
```

> pch 文件同样使用**正电压域**（|VSG|、|VSD|）查询，与 nch 一致。

## 5. Mathematica

HDF5 可直接用 Mathematica `Import` 读取。见项目 `mathematica/tsmcN40_lookup.wl`：

```mathematica
<< "mathematica/tsmcN40_lookup.wl";
data = LoadTsmcN40["D:\\tsmcN40_lookup\\nch_tt.h5"];
f = N40Interpolant[data, "GM_ID"];
f[0.1, 0.6, 0.7, 0.0]              (* 4-D 插值查表 (L,VGS,VDS,VSB) *)
cur = SliceVGS[data, "GM_ID", 0.1, 0.7, 0.0];   (* 固定 (L,VDS,VSB) 扫 VGS *)
```

也可直接取数组后用 `ListInterpolation` 自建插值。
