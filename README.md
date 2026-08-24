# gmId_tsmcN40

TSMC N40 (CRN40LP) PDK 的 1.1V 核心 MOS 器件（nch / pch，BSIM4）gm/ID 查找表生成工具集：Spectre 扫描仿真 + 数据提取，输出 HDF5（.h5）查找表。

---

## 项目结构

```
.
├── config_tsmcN40.py     # 1.1V 核心器件扫描参数配置（nch / pch）
├── run_sim.py            # Spectre 并行仿真
├── extract_new.py        # 数据提取，生成 .h5 查找表
├── psf_reader.py         # PSF 格式仿真结果读取器
├── pipeline.sh           # 逐 L 流水线（仿真 → 提取 → 清理）
├── get_nl.py             # 输出 fine 模式 L 点数（供 pipeline.sh 使用）
├── lookup_funs/          # Python 查找表查询 API 与使用文档
├── mathematica/          # Mathematica 加载 / 插值脚本
│   └── tsmcN40_lookup.wl
└── README.md
```

## 参考器件与扫描网格（config_tsmcN40.py, fine 模式）

| 参数       | 值                                                          |
| ---------- | ----------------------------------------------------------- |
| 参考宽度 W | 10 µm                                                    |
| VGS        | 0 ~ 1.1 V，步长 0.02 → 56 点                                |
| VDS        | 0 ~ 1.1 V，步长 0.02 → 56 点                                |
| VSB        | 0 ~ 0.8 V，步长 0.1 → 9 点                                |
| L          | 0.04 ~ 5 µm → 48 点（0.04~0.14 步 0.01；0.16~0.50 步 0.02；0.55~1.0 步 0.05；1.2,1.4,1.6,1.8,2,2.5,3,4,5） |
| 单变量网格 | 48 × 56 × 56 × 9（(L, VGS, VDS, VSB)）= 1,354,752 点        |

每文件 21 个变量（19 DC + 2 噪声），共 10 文件（nch/pch × tt/ff/ss/fs/sf）。

## 存储格式（.h5，HDF5）

- 每文件 = 器件+工艺角，如 `nch_tt.h5`、`pch_tt.h5`。
- 4-D 数据数组 `(L, VGS, VDS, VSB)`，
  **float32 + gzip(level 6) + shuffle**，chunk = 一个 L 的整片 `(1, VGS, VDS, VSB)`。
- 元数据存为标量/字符串数据集：`CORNER DEVICE INFO TEMP W NFING`。
- 坐标轴存为 1-D 数据集：`L VGS VDS VSB`。
- 变量名：`ID VT IGD IGS GM GMB GDS CGG CGS CSG CGD CDG CGB CDD CSS VDSAT`（+ 噪声 `STH SFL`）。
  `VDSAT`=BSIM 饱和电压，唯一存储的复合量；其余复合量均为比例，现算：
  `GM_ID`=gm/id、`GM_CGG`=gm/Cgg（fT=`GM_CGG/(2π)`）、`GM_GDS`=gm/gds（固有增益）、`ID_W`=ID/W。

每文件约 180 MB（float32 + gzip6），10 文件合计约 1.8 GB。

## 关键配置说明（config_tsmcN40.py）

- **模型文件**：直接使用基础模型文件（非 shrink0d9 包装），路径由
  `machine.env` 中的 `PDK_MODEL_FILE` 指定（默认 `$HOME/PDKs/tsmc40nm/models/spectre/crn40lp_2d5_v2d0_2.scs`）。
  分别以 `section=stat` / `section=global` / `section={corner}` 包含三次。
- **不做额外 0.9 缩放**：PDK 文档（T-N40-CM-SP-001）说明基础模型卡已自带缩放因子，
  若再叠加 `setscale scalefactor=0.9` 会导致双重缩放（0.9×0.9=0.81）。
- **工艺角**：`tt ss ff fs sf`（对应模型 section）。
- **Spectre**：`machine.env` 中的 `SPECTRE_BIN`（默认 `/opt/eda/cadence/SPECTRE181/tools.lnx86/bin/spectre`）。

## 运行环境（VM）

| 组件                 | 环境要求                                     |
| -------------------- | -------------------------------------------- |
| Spectre 扫描仿真     | Cadence SPECTRE 18.1（`/opt/eda/cadence/SPECTRE181`） |
| 数据提取 & .h5 生成  | Python ≥ 3.10 + numpy, h5py, **psf-parser**（二进制 PSF 解析） |
| 查找表 API           | Python 3.x + numpy, scipy, h5py              |
| Pipeline             | Bash + Python                                     |

机器路径（PDK 模型、Spectre、venv Python、数据目录）集中在项目根目录的 `machine.env`
（模板 `machine.env.example`，已 git 忽略），值可用 `$HOME` 或引用前面的键（如 `$SIMDIR/raw`）。
Python 经 `paths.py` 读取，`pipeline.sh` 经 `source machine.env` 读取，脚本不硬编码目录或用户名；
环境变量同名优先。

原始仿真数据与查找表放在 `$HOME/simulation`：`RAWDIR`（原始 PSF，提取后即删）、
`MATDIR`（最终 .h5）、`LOGDIR`（日志）。提取需 Python ≥ 3.10 的 venv（含 psf-parser 解析
二进制 PSF；`PYTHON_BIN` 默认 `$HOME/venvs/gmId313/bin/python`），否则 `psf_reader.py`
回退到 ASCII 解析。

## 快速开始

### 1. 配置 machine.env

```bash
cp machine.env.example machine.env   # 然后填入本机路径
```

### 2. Pipeline 一键运行（推荐）

```bash
bash pipeline.sh
```

逐 L 增量处理，每级栅长 L 独立走完 仿真→提取→清理，5 个工艺角并行。
原始仿真数据放 `$RAWDIR`（默认 `$HOME/simulation/raw`），提取完即删。
.h5 输出到 `$MATDIR`（默认 `$HOME/simulation/out`）。启动时会自动清理上次
中断残留的 raw 数据、日志和过期 spectre 工作目录。

### 3. 分步运行

```bash
# Step 1 — 扫描仿真（指定 corner 和输出目录）
$PYTHON_BIN run_sim.py tt --fine --outdir "$RAWDIR"

# Step 2 — 数据提取
$PYTHON_BIN extract_new.py --fine --srcdir "$RAWDIR" --outdir "$MATDIR" --workers 5

# Step 3 — 使用查找表 API（示例见下文"Python 使用"）
from lookup_funs.lookup_table import loadmat, lookup, lookupVGS
nch = loadmat(r"D:\tsmcN40_lookup\nch_tt.h5")
```

## Python 使用

`lookup_funs/lookup_table.py` 提供查询 API（numpy + scipy + h5py），完整说明见
`lookup_funs/LOOKUP_README.md`：

```python
from lookup_funs.lookup_table import loadmat, lookup, lookupVGS

nch = loadmat(r"D:\tsmcN40_lookup\nch_tt.h5")
pch = loadmat(r"D:\tsmcN40_lookup\pch_tt.h5")

gm_id = lookup(nch, 'GM_ID', 'VGS', 0.6, 'VDS', 0.7, 'VSB', 0, 'L', 0.1)     # 查表
VGS   = lookupVGS(nch, 'GM_ID', 15, 'VDS', 0.7, 'VSB', 0, 'L', 0.1)          # 反查 VGS
wt    = lookup(nch, 'GM_CGG', 'GM_ID', [5, 10, 15], 'VDS', 0.7, 'VSB', 0, 'L', 0.1)  # 交叉查表

nch.VDSAT                                # 直接取存储的 4-D 数组 (L,VGS,VDS,VSB)
nch.GM / nch.CGG                         # 任意比例量现算（= GM_CGG，fT=GM_CGG/(2π)）
```

- `loadmat` 返回 `LookupTable`，存储变量名为属性（`ID GM VDSAT ...`），4-D 数组
  `(L, VGS, VDS, VSB)`（与官方 .mat 布局一致；也可直接加载官方 `.mat` 文件）。
- `lookup(data, outvar, ...)`：`outvar` 可为存储变量或比例量（`GM_ID`、`ID_W`、`GM_CGG`、`GM_GDS` 等，
  由原始量现算），支持 4-D 插值、交叉查表；`lookupVGS(...)` 按 gm/id 或电流密度反查 VGS。
- pch 文件用正电压域（|VSG|、|VSD|）查询，与 nch 一致。

## Mathematica 使用

HDF5 可直接用 Mathematica `Import` 读取。见 `mathematica/tsmcN40_lookup.wl`：

```mathematica
<< "mathematica/tsmcN40_lookup.wl"
data = LoadTsmcN40["D:\\tsmcN40_lookup\\nch_tt.h5"];
gmId = lookup[data, "GM_ID", "VGS", 0.6, "VDS", 0.7, "VSB", 0, "L", 0.1];
vgs = lookupVGS[data, "GM_ID", 15, "VDS", 0.7, "VSB", 0, "L", 0.1];
wt = lookup[data, "GM_CGG", "GM_ID", {5, 10, 15}, "VDS", 0.7, "L", 0.1];

(* MapThread 参数按位置配对；原生 VDS 列表仍是笛卡尔维度 *)
gain = lookup[data, "GM_GDS", "GM_ID", #1, "L", #2,
    "VDS", {0.5, 0.7, 0.9}, "VSB", 0] & ~MapThread~
  {{8, 10, 12}, {0.04, 0.06, 0.08}};
Dimensions[gain]                         (* {3, 3} = 配对维度 × VDS维度 *)
```

- `lookup`、`lookupVGS` 的函数划分、交替名称/值参数及默认值与附录和 Python API 对齐。
- 加载查表文件后，函数体直接为 `lookup` 或 `lookupVGS` 的 `/@` 与 `MapThread`
  会自动合并为分块批量查询。`MapThread` 中的列表按位置配对，函数体内直接传给查表函数的
  普通列表仍形成笛卡尔积；其他 Mathematica 函数的 `Map`/`MapThread` 行为不变。
- 可用 `DisableTsmcN40MapOptimization[]` 禁用上述系统规则，并用
  `EnableTsmcN40MapOptimization[]` 再次启用。复杂 Slot 表达式或不支持的查表模式自动回退到逐项求值。
- `XTRACT`、`XTRACT2` 位于 `mathematica/ekv_extract.wl`，保持附录中的位置参数接口。
- 附录缺省值为：`lookup` 使用最小 `L`、完整 `VGS`、`VDD/2` 和 `VSB=0`；
  `lookupVGS` 另使用 PCHIP；`XTRACT/XTRACT2` 使用 `rho=0.6` 和 300 K。
- `XTRACT2` 缺省直接使用总电流并返回 `IS`。显式提供扩展宽度参数时才返回 `JS=IS/W`。
- `N40Interpolant`、`SliceVGS` 保留为需要直接操作插值函数时使用的底层接口。
- 验证：`wolframscript -script mathematica/xtract_demo.wl`（读取真实 .h5，
  打印与 Python `lookup` 一致的插值结果）。

## 备注

- N40 核心器件最小 L 为 0.04 µm（40nm），L 网格从 0.04 µm 起步。
- 如需扩展 2.5V IO（nch_25 / pch_25）或 LVT / HVT，在 `config_tsmcN40.py` 中
  增加相应 device 分组即可。
