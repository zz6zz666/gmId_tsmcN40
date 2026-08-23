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

### 附录缺省值

| 函数 | 缺省值 |
| --- | --- |
| `lookup` | `L=min(data.L)`，`VGS=data.VGS`，`VDS=max(data.VDS)/2`，`VSB=0` |
| `lookup` 交叉查表 | 最终一维插值 `METHOD='pchip'`，提示 `WARNING='on'`；多维插值固定为线性 |
| `lookupVGS` | `L=min(data.L)`，`VDS=max(data.VDS)/2`，`VSB=0`，`METHOD='pchip'` |
| `XTRACT` | `rho=0.6`，`TEMP=300 K` |
| `XTRACT2` | `rho=0.6`，温度按 `300 K`；输入 `ID` 不做宽度归一化，返回 `IS` |

`lookupVGS` 按附录限制为至多一个向量输入。交叉查表显式传入 `VGS` 时，该向量用作
交点搜索范围；未传入时使用完整的 `data.VGS`。

`XTRACT2` 额外保留 `TEMP` 和 `W` 扩展参数。两者均不影响缺省调用：`TEMP` 缺省为
300 K，`W=None` 表示直接使用总电流并返回 `IS`；只有显式传入 `W`（µm）时才使用
`ID/W` 并返回电流密度 `JS`。

```python
from ekv_extract import XTRACT, XTRACT2

y = XTRACT(nch, 0.1, 0.55, 0)       # rho=0.6, TEMP=300 K
p = XTRACT2(vgs, ids)                 # [n, VT, IS], rho=0.6, 300 K
p_density = XTRACT2(vgs, ids, W=10)   # [n, VT, JS]
```

## 5. Mathematica

Mathematica API 与附录及本项目 Python API 保持相同的函数划分和参数顺序。见项目
`mathematica/tsmcN40_lookup.wl`：

```mathematica
<< "mathematica/tsmcN40_lookup.wl";
data = LoadTsmcN40["D:\\tsmcN40_lookup\\nch_tt.h5"];

gmId = lookup[data, "GM_ID", "VGS", 0.6, "VDS", 0.7, "VSB", 0, "L", 0.1];
vgs = lookupVGS[data, "GM_ID", 15, "VDS", 0.7, "VSB", 0, "L", 0.1];
wt = lookup[data, "GM_CGG", "GM_ID", {5, 10, 15}, "VDS", 0.7,
  "VSB", 0, "L", 0.1];
```

参数使用与 Matlab/Python 相同的交替名称和值形式，名称不区分大小写，未提供的
`L/VGS/VDS/VSB` 使用附录规定的默认值。`N40Interpolant` 和 `SliceVGS` 仍保留为底层接口。

```mathematica
<< "mathematica/ekv_extract.wl";
y = XTRACT[data, 0.1, 0.55, 0];
p = XTRACT2[vgs, ids];
pDensity = XTRACT2[vgs, ids, 0.6, 300., 10.];
```

Mathematica 的 `XTRACT2` 宽度扩展缺省为 `Automatic`，语义与 Python 的 `W=None`
相同：不做宽度归一化。
