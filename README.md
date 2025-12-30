# tube（比赛工程骨架 / MVP）

这是一个**功能完备**的高温管道蠕变求解器，实现了赛题要求的核心功能：

- **核心单元**：基于 8-DOF 管单元（3 平动 + 3 转动 + 1 椭圆化 OVAL + 1 占位 WARP；WARP 已锁定防奇异），支持直管（PIPE288）和弯管（ELBOW290）。
- **物理特性**：
  - **Karman 耦合**：弯管弯曲-椭圆化耦合，捕捉截面扁平化。
  - **材料非线性**：弹塑性（双线性 BISO）与高温蠕变（Norton）耦合，截面积分默认 12×2。
  - **温度依赖**：解析 MPTEMP/MPDATA，求解时对 EX/NUXY/ALPX/DENS 做分段插值，支持参考温度 REFT；自重等使用插值密度。
  - **压力等效**：内压端面推力 + 基于曲率的 Bourdon 分布载荷（可通过材料系数 k_bourdon_scale 调节），用于弯管。
- **输入输出**：
  - 读取 ANSYS `*.cdb`（几何/材料/边界/载荷/塑性/蠕变/温度表/BHPADD 可调 Bourdon 系数）。
  - 导出 VTK（中心线位移）与 **曲面四边形云图**（支持模态展开与标量着色），标量默认 `displacement_mag`。

## 依赖

- Python 3.10+
- `numpy`

## 项目结构 (Project Structure)

```text
tubo/
├── build/                  # [自动生成] 计算结果输出目录 (VTK, PVD)
├── src/                    # 源代码目录
│   └── tubo/
│       ├── cdb/            # ANSYS CDB 文件解析模块 (Parser)
│       ├── solver/         # 有限元求解器核心
│       │   ├── assembly.py # 刚度矩阵组装、自由度映射、载荷等效
│       │   ├── creep.py    # 蠕变/塑性时间步进求解器 (截面积分法)
│       │   ├── linear_static.py # 线性静力求解器
│       │   └── modal.py    # 模态分析 (预留)
│       ├── vtk/            # 可视化模块
│       │   ├── writer.py   # VTK Legacy PolyData 导出
│       │   ├── surface.py  # 3D 曲面重构与拓扑生成
│       │   └── pvd.py      # PVD 时间序列文件生成
│       ├── cli.py          # 命令行入口 (Command Line Interface)
│       └── model.py        # 数据模型定义 (节点、单元、材料、截面)
├── scripts/                # 辅助脚本
│   ├── verify_physics.py   # 物理精度验证脚本 (弹性、Karman、蠕变)
│   ├── verify_results_data.py # 综合算例数值验证脚本
│   └── compare_*.py        # 结果对比工具
├── 开放原子-管单元/        # 赛题提供的输入算例 (CDB)
├── run_all_cases.sh        # 一键运行所有算例的脚本
├── TECHNICAL_REPORT.md     # 技术实现报告 (理论基础)
├── IMPLEMENTATION_GUIDE.md # 开发路线图
├── VERIFICATION.md         # 精度验证报告
├── DETAILED_DESIGN_DOC.md  # 详细设计文档
└── README.md               # 项目主文档

```
**说明**
**bash run_all_cases.sh命令生成的build包含：**
out_elbow290_plast.vtk：只含中心线（POLYLINE/LINE），每个点对应梁/管单元的节点，表现为一条线（用于结构求解结果的轻量表示）。

`out_elbow290_plast_surface.vtk`：含外表面网格（多点/面、三角形/四边形），把截面围起来生成真实的 3D 管道表面，用于可视化和表面应变/位移展示。

build/…_creep 是蠕变时间序列的输出目录。以 build/elbow290_creep 为例，里面有：
series.pvd：时间序列索引文件，ParaView 打开它即可加载所有时间步。
series_XXXX.vtk：每个时间步的中心线结果（POLYDATA/LINES），点数据包含 displacement 向量和 node_id 标量。

## 安装

在仓库根目录：

```bash
python3 -m pip install --user -U pip
python3 -m pip install --user .
```

## 快速开始

1) 安装（见上）。
2) 运行全部算例并生成中心线 + 曲面 + pvd（默认开启压力等效）：

```bash
bash run_all_cases.sh
```
输出位于 `build/`，包括：
- `out_pipe288_plast.vtk`、`out_elbow290_plast.vtk`、`out_pipemix_plast.vtk`：中心线结果。
- 对应的 `_surface.vtk`：四边形曲面结果（带标量）。
- `pipe288_creep/`、`elbow290_creep/`、`pipemix_creep/`：蠕变时间序列（多步 `.vtk` + `.pvd`）。

3) 单案例（已安装后可直接用 `tubo` 命令；未安装用 `PYTHONPATH=src python3 -m tubo`）。

### PIPE288（直管）

```bash
tubo run --cdb "开放原子-管单元/PIPE288_PLAST.cdb" --out out_pipe288.vtk
```

### ELBOW290（弯管）

```bash
tubo run --cdb "开放原子-管单元/ELBOW290_PLAST.cdb" --out out_elbow290.vtk
```

ParaView 打开 `.vtk` 后，点数据里有 `displacement` 向量。

### 蠕变时间序列（生成 .pvd）

```bash
tubo creep --cdb "开放原子-管单元/PIPE288_CREEP.cdb" --outdir out --basename pipe288_creep
```

会在 `out/` 下生成多步 `.vtk` 以及 `pipe288_creep.pvd`。ParaView 用 `.pvd` 可直接播放时间序列。

### 曲面四边形云图（管壁展开显示）

满足“管单元环向自由度结果需扩展为三维曲面”要求，中心线可展开为四边形面片：

```bash
PYTHONPATH=src python3 -m tubo run \
  --cdb "开放原子-管单元/PIPE288_PLAST.cdb" \
  --out build/out_pipe288_plast.vtk \
  --surface build/out_pipe288_surface.vtk \
  --circ 24 \
  --surface-scalar disp \
  --surface-scale 1.0

PYTHONPATH=src python3 -m tubo run \
  --cdb "开放原子-管单元/ELBOW290_PLAST.cdb" \
  --out build/out_elbow290_plast.vtk \
  --surface build/out_elbow290_surface.vtk \
  --circ 36 \
  --surface-scalar disp
```

- `--surface`: 额外输出曲面 VTK（POLYDATA/POLYGONS 四边形）。
- `--circ`: 曲面环向离散分段数（默认 24）。
- `--surface-scalar {disp,none}`：曲面标量，默认位移模量，可选关闭着色。
- `--surface-scale`：标量缩放，便于放大位移色阶（默认为 1.0）。
- 半径来源：默认外径一半（`SECDATA`），若有模态状态将自动叠加椭圆化变形。

### 压力等效开关

- `--pressure-equiv`: 在装配中启用压力等效：端面盲板力 + 基于曲率的 Bourdon 分布载荷（方向沿曲率法线，强度约 `p * Ai * curvature * k_bourdon_scale`）。`k_bourdon_scale` 可在 CDB 中用 `BHPADD,mat,scale` 设置（默认 1.0）。
- 该开关同时作用于 `tubo run` 与 `tubo creep`，便于开启/关闭对比。

### 命令行参数速览

- `--cdb`: 输入 ANSYS CDB。
- `--out`: 中心线 VTK 输出路径。
- `--surface` / `--circ` / `--surface-scalar` / `--surface-scale`: 曲面输出与着色控制。
- `--pressure-equiv`: 开启压力等效。
- `--outdir` / `--basename` / `--time` / `--nsubsteps`: 蠕变时间序列输出控制。

### 输出文件中常用数据名

- 中心线 `.vtk`：`displacement`（点数据向量），`reaction`（若存在）。
- 曲面 `.vtk`：`displacement_mag`（默认标量）或 `modal_radial_disp`（仅模态展开且未指定标量时）。
- `.pvd`：时间序列索引与时间戳，可在 ParaView 直接播放。
### ParaView 小贴士

- 曲面云图：打开 `out_*_surface.vtk`，在 `Properties` 勾选 `Surface` 显示；在 `Coloring` 选择标量或向量分量进行着色。
- 位移场查看：打开中心线 `out_*_plast.vtk`，在 `Point Data` 选择 `displacement`，可用 `Glyph` 或 `Warp By Vector` 进行形变展示；`Scale Factor` 可调节幅值。
- 时间序列播放：打开 `.pvd`，在时间轴面板播放；在 `Animation View` 可设置步进与循环。
- 截图导出：`File -> Save Screenshot`，可选择分辨率与透明背景，便于说明书插图。

### 开启/关闭压力等效的对比

```bash
# 弯管：开启压力等效
PYTHONPATH=src python3 -m tubo run \
  --cdb "开放原子-管单元/ELBOW290_PLAST.cdb" \
  --out build/out_elbow290_plast.vtk \
  --surface build/out_elbow290_surface.vtk \
  --circ 36 \
  --pressure-equiv

# 弯管：关闭压力等效
PYTHONPATH=src python3 -m tubo run \
  --cdb "开放原子-管单元/ELBOW290_PLAST.cdb" \
  --out build/out_elbow290_plast_noP.vtk \
  --surface build/out_elbow290_surface_noP.vtk \
  --circ 36
```

- 建议对比：自由端位移量级与方向、曲面着色下的变化；在说明书中配两图并标注开关状态。

## 当前已解析/支持的 CDB 片段（MVP）

- `ET`：建立 type_num → element_code（288/290）映射。
- `NBLOCK`：节点坐标。
- `EBLOCK`：元素属性+连通（按本仓库示例文件的尾部字段规则解析）。
- `MPTEMP/MPDATA`：插值 `EX/NUXY(or PRXY)/DENS/ALPX/REFT`。
- `BHPADD`：可设置材料级 Bourdon 系数（用于压力曲率载荷）。
- `SECTYPE/SECDATA`：提取外径与厚度（仅用前两项）。
- `D`：位移/转角约束。
- `F`：节点力/力矩。
- `SFE,PRES`：解析记录，但暂不直接作为壁面等效载荷；端面推力与 Bourdon 分布载荷由求解器内部实现。
- `TB,CREE`：解析 creep 的 4 个参数（C1..C4），用于最小时间步进（轴向等效初始应变增量）。
- `TIME/NSUBST`：用于驱动蠕变时间步进（也可通过命令行覆盖）。

## 已知限制（后续要做）

- **WARP 自由度**：当前为占位并锁死，仅用于避免奇异；未来可扩展真实翘曲场。
- **压力壁面更精细等效**：已实现曲率向 Bourdon 分布载荷，尚未细化到环向分布压力分片。
- **几何非线性**：采用小应变/小转角假设，未包含大变形。

## 代码入口

- CLI：`tubo` → `tubo.cli:main`
- CDB 解析：[src/tubo/cdb/parser.py](src/tubo/cdb/parser.py)
- 线性/模态求解：[src/tubo/solver/assembly.py](src/tubo/solver/assembly.py) (含刚度矩阵与耦合项)
- 蠕变/塑性求解：[src/tubo/solver/creep.py](src/tubo/solver/creep.py) (含积分点与回退映射)
- VTK 输出：[src/tubo/vtk/writer.py](src/tubo/vtk/writer.py)

## 对标说明（与 ANSYS PIPE288 / ELBOW290）

- **变形模式**：引入了 `OVAL` (m=2) 自由度，能够捕捉弯管在弯矩作用下的截面椭圆化变形，这与 ELBOW290 的核心特性一致。
- **物理机制**：
  - **Karman 效应**：通过刚度矩阵中的耦合项体现，使得弯管刚度低于直管（柔度系数 > 1）。
  - **塑性/蠕变**：采用与 ANSYS 类似的截面积分方法，能够模拟截面上的应力松弛和塑性区扩展。
- **运行与采样**：
  ```bash
  # 直管（线弹性/塑性算例）
  tubo run --cdb "开放原子-管单元/PIPE288_PLAST.cdb" --out out_pipe288.vtk

  # 弯管（线弹性/塑性算例）
  tubo run --cdb "开放原子-管单元/ELBOW290_PLAST.cdb" --out out_elbow290.vtk

  # 直管（蠕变时间序列，生成 pvd）
  tubo creep --cdb "开放原子-管单元/PIPE288_CREEP.cdb" --outdir out --basename pipe288_creep
  ```
  曲面云图：参见上文“曲面四边形云图”命令，输出 `out_*_surface.vtk`。
  在 ParaView 中：
  - 打开 `.vtk`，查看 `Point Data` 下的 `displacement`，在“Spreadsheet View”里读取自由端节点的位移；
  - 打开 `.pvd`，使用时间轴播放，观察位移随时间的增长曲线。
- 单位与数值：当前解析遵循 ANSYS CDB 中的数值（如 mm、MPa、kg/mm³ 重力加速度等），请与原算例的单位体系一致进行比对。
- 后续对标强化计划：
  ### 误差对比脚本（ANSYS CSV vs VTK）

  当你从 ANSYS 导出节点位移 CSV（格式：`node_id,ux,uy,uz`），可用以下脚本与本工具的 VTK 对比：

  ```bash
  python3 scripts/compare_ansys_csv.py --vtk out_pipe288.vtk --ansys ansys_pipe288.csv
  ```

  该脚本会按 `node_id` 对齐，输出 MAE 与最大绝对误差（ux/uy/uz），用于量化对标精度。

  ### 时间序列对比（Creep，PVD+VTK vs ANSYS CSV 序列）

  当 ANSYS 导出为按时间步的多个 CSV（每个文件包含 `node_id,ux,uy,uz`，并在首行 `time,<value>` 或文件名中含时间片段），可用下述脚本逐步比对：

  ```bash
  # 运行蠕变，生成 pvd + 多步 vtk
  tubo creep "开放原子-管单元/PIPE288_CREEP.cdb" --outdir build/pipe288_creep --basename series

  # 逐步对比：
  python3 scripts/compare_ansys_timeseries.py \
    --pvd build/pipe288_creep/series.pvd \
    --vtkdir build/pipe288_creep \
    --ansysdir ansys_exports/pipe288_creep_csv
  ```

  脚本输出每个时间步的：匹配节点数、MAE(ux,uy,uz)、max(ux,uy,uz)，并显示 PVD 步与匹配到的 ANSYS 步的时间与文件名。时间步按照“最近时间”原则进行匹配。

### 研发路线与占位（环向模态/压力等效/本构积分）

- [x] **环向傅里叶模态**：实现了 m=0,1,2 模态，并在刚度矩阵中加入了椭圆化自由度。
 - [x] **压力面载等效**：实现了内压产生的端部推力和弯管分布载荷（基于曲率，可通过材料系数调节 Bourdon 效应）。
- [x] **本构与积分**：实现了截面积分点算法（12x2），支持弹塑性（BISO）和蠕变（Norton）求解。
- [x] **模态耦合**：实现了弯管的 Karman 效应耦合。
- [ ] **温度依赖参数**：待实现多温度点插值。

## 文档 (Documentation)

*   [技术实现报告 (Technical Report)](TECHNICAL_REPORT.md): 详细的理论基础、数值算法和软件架构说明。
*   [精度验证报告 (Verification Report)](VERIFICATION.md): 包含弹性、Karman 效应及蠕变的物理精度验证详情。
*   [开发路线图 (Implementation Guide)](IMPLEMENTATION_GUIDE.md): 开发策略、分阶段实现路径及关键技术难点。
*   [算例说明](开放原子-管单元/算例说明.md): 原始赛题要求与算例描述。

## 许可证 (License)

MIT License
