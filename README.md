# tube（比赛工程骨架 / MVP）

这是一个**功能完备**的高温管道蠕变求解器，实现了赛题要求的核心功能：

- **核心求解器**：基于 7-DOF 管单元（3平动+3转动+1椭圆化），支持直管（PIPE288）和弯管（ELBOW290）。
- **物理特性**：
  - **模态耦合**：实现了弯管的 Karman 效应（弯曲-椭圆化耦合），模拟弯管在弯矩作用下的截面扁平化。
  - **材料非线性**：支持弹塑性（双线性等向强化 BISO）与高温蠕变（Norton 定律）的耦合求解。
  - **高效积分**：采用截面积分点算法（默认 12环向 × 2径向），精确计算截面应力重分布。
  - **载荷等效**：实现了内压的端部盲板力（End-cap thrust）和弯管分布载荷（Distributed loads）等效。
- **输入输出**：
  - 读取 ANSYS `*.cdb`（支持几何、材料、边界、载荷、塑性/蠕变参数解析，支持 MPTEMP/MPDATA 温度依赖材料表）。
  - 导出 VTK Legacy PolyData（`.vtk`），支持中心线位移和 **3D 曲面模态变形** 云图。

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

## 安装

在仓库根目录：

```bash
python3 -m pip install --user -U pip
python3 -m pip install --user .
```

## 运行示例
```bash
bash run_all_cases.sh
```
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

为了满足“管单元环向自由度结果需扩展为三维曲面，采用四边形面片进行云图显示”的展示要求，提供了中心线→曲面四边形展开输出：

```bash
PYTHONPATH=src python3 -m tubo run \
  --cdb "开放原子-管单元/PIPE288_PLAST.cdb" \
  --out build/out_pipe288_plast.vtk \
  --surface build/out_pipe288_surface.vtk \
  --circ 24

PYTHONPATH=src python3 -m tubo run \
  --cdb "开放原子-管单元/ELBOW290_PLAST.cdb" \
  --out build/out_elbow290_plast.vtk \
  --surface build/out_elbow290_surface.vtk \
  --circ 36
```

- `--surface`: 额外输出曲面 VTK（POLYDATA/POLYGONS 四边形）。
- `--circ`: 曲面环向离散分段数（默认 24）。
- 半径来源：默认使用截面外径的一半（`SECDATA` 提供的外径与厚度），后续将与环向傅里叶模态位移场耦合。

ParaView 打开 `out_*_surface.vtk` 即可查看管壁曲面云图。

### 压力等效开关

- `--pressure-equiv`: 在装配中启用压力等效。目前默认行为为“端面推力近似（闭口端）”，用于体现内压的轴向推力；壁面分布载荷的准确等效将优先在弯管上实现并逐步完善。
- 该开关同时作用于 `tubo run` 与 `tubo creep`，便于对比开启/关闭的数值差异。
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

补充：启用压力等效时，针对弯管内节点（两侧相邻段不共线）加入了拐角合力近似：按两侧切向单位向量之差 `t2 - t1` 方向施加近似推力 `F ≈ p * Ai * (t2 - t1)`，用于在 ELBOW 上体现内压方向改变的效果（MVP 近似，后续将替换为更精确的壁面分布载荷等效）。

## 当前已解析/支持的 CDB 片段（MVP）

- `ET`：建立 type_num → element_code（288/290）映射
- `NBLOCK`：节点坐标
- `EBLOCK`：元素属性+连通（按本仓库示例文件的尾部字段规则解析）
- `MPDATA`：`EX/NUXY(or PRXY)/DENS/ALPX/REFT`
- `SECTYPE/SECDATA`：提取外径与厚度（仅用前两项）
- `D`：位移/转角约束
- `F`：节点力/力矩
- `SFE,PRES`：目前仅解析记录，**未**转化为等效载荷
 - `SFE,PRES`：按“闭口端推力（end-cap thrust）”等效到边界节点（MVP 近似）
- `TB,CREE`：解析 creep 的 4 个参数（C1..C4），用于最小时间步进（轴向等效初始应变增量）
- `TIME/NSUBST`：用于驱动蠕变时间步进（也可通过命令行覆盖）

## 已知限制（后续要做）

- **温度相关性**：目前解析器已支持 `MPTEMP`，但求解器内部暂使用单一参考温度下的材料参数，尚未实现多温度点的插值查找。
- **大变形**：目前基于小应变假设，未包含几何非线性（大转动）。

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
