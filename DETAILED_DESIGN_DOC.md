# 详细技术说明书与代码实现文档

## 1. 项目概述

本项目实现了一套面向高温高压管道的高效有限元求解器，核心目标是在保持计算效率的同时，准确捕捉弯管截面椭圆化、Karman 耦合、弹塑性与高温蠕变等关键物理。通过 8 自由度管单元（含椭圆化 OVAL 与占位 WARP）和截面积分算法，能用一维单元的自由度数量获得接近实体/壳单元的精度。

### 应用价值
- **高效**：相对 3D 实体/壳模型，自由度数量可降 90%+，适合长管网快速迭代。
- **对接便捷**：直接解析 ANSYS CDB 输入（几何/材料/边界/载荷/温度表/蠕变/塑性），无缝融入现有流程。
- **物理完整**：支持 Karman 耦合、Bourdon 压力等效、温度插值、自重/压力/端面推力、蠕变与塑性迭代。

### 创新点
1) **8-DOF 管单元**：在 6-DOF 梁上加入 OVAL + WARP（锁定防奇异），弯管刚度矩阵含椭圆化耦合。  
2) **截面积分 + 本构返回映射**：放弃简化截面刚度，采用环向×径向积分点，逐点弹塑/蠕变更新并回收为等效广义应变载荷。  
3) **曲率驱动压力等效**：端面推力 + 基于曲率的 Bourdon 分布载荷（可通过材料系数调节），针对弯管压力效应。  
4) **温度表插值**：MPTEMP/MPDATA 对 E/ν/α/ρ 做分段插值，统一驱动热应变、自重和刚度。  
5) **自动曲面重构**：从 1D 拓扑提取路径，平行运输构造稳定框架，生成四边形管壁曲面，并可叠加模态/标量云图。

---

## 2. 代码总览与算法角色

- `src/tubo/model.py`：数据模型（节点/单元/材料/截面/模态状态），含温度表与 Bourdon 系数。
- `src/tubo/cdb/parser.py`：CDB 解析器，状态机读 NBLOCK/EBLOCK/MPDATA/MPTEMP/TB 等，生成 Model。
- `src/tubo/solver/assembly.py`：局部 16×16 刚度、压力等效、温度插值、自重/端力装配，生成全局 K/F。
- `src/tubo/solver/linear_static.py`：线性静力求解，处理约束、正则化、解算并回填结果。
- `src/tubo/solver/creep.py`：时间步进，积分点弹塑-蠕变更新，组装等效应变载荷，输出序列。
- `src/tubo/solver/modal.py`：环向傅里叶模态 m=0/1/2 的位移展开，用于曲面生成和颜色场。
- `src/tubo/vtk/writer.py`：中心线 VTK 导出。
- `src/tubo/vtk/surface.py`：曲面四边形生成，支持模态/标量映射与角点平滑。
- `src/tubo/vtk/pvd.py`：PVD 时间序列文件。
- `src/tubo/cli.py`：命令行入口，`run`/`creep`，带压力等效和曲面输出选项。
- `run_all_cases.sh`：一键跑全套算例（含压力等效、曲面）。
- `scripts/verify_physics.py` / `verify_results_data.py`：物理与回归验证脚本。

---

## 3. 各文件算法与实现要点

### 3.1 `model.py`
- **实体**：Node/Element/Material/SectionPipe/ModalConfig/FourierModalState。  
- **温度表**：`Material.temp_tables` 保存 MPTEMP/MPDATA；插值时取 E/ν/α/ρ。`k_bourdon_scale` 可调 Bourdon 强度。  
- **DOF 定义**：`_DOF_ORDER` 含 UX, UY, UZ, RX, RY, RZ, OVAL, WARP；WARP 目前用惩罚锁死。
- **材料复制**：拷贝时保留塑性/蠕变/温度表/Bourdon 系数，避免多材料分配时丢字段。

### 3.2 `cdb/parser.py`
- **状态机解析**：按块识别 NBLOCK/EBLOCK/MP/MPTEMP/MPDATA/SECTYPE/SECDATA/D/F/SFE/TB,CREE/TIME/NSUBST/BHPADD。  
- **几何**：NBLOCK 读节点；EBLOCK 读类型、截面、材料号与连通。  
- **材料**：支持 EX/NUXY/ALPX/DENS/REFT，MPTEMP+MPDATA 组合为分段表；BHPADD 读取 Bourdon 系数。  
- **边界载荷**：D/F 约束和集中力；SFE,PRES 记录但壁面等效内部处理。  
- **塑性/蠕变**：TB,CREE 读 C1..C4；（BISO 留接口）。  
- **时间控制**：TIME/NSUBST 用于 creep 默认步长，可 CLI 覆盖。

### 3.3 `solver/assembly.py`
- **几何与曲率**：对 ELBOW 元素用三点定圆 `_circle_from_three`，得中心、曲率、切向/法向/副法向；直段视为零曲率。  
- **局部刚度 (16×16)**：基于 Euler-Bernoulli 基础刚度扩展到 8 DOF/端；OVAL 刚度和弯曲耦合项体现 Karman；WARP 对角加惩罚锁定。  
- **温度插值**：`_interp_prop` 在 MPTEMP/MPDATA 表中线性插值 E/ν/α/ρ；热应变按 α·(T-REFT)；自重用插值密度。  
- **载荷等效**：
  - 自重：局部坐标分解，等效到节点。
  - 压力端面推力：`p * Ai` 沿轴向施加。
  - Bourdon 曲率载荷：`q = -p * Ai * curvature * k_bourdon_scale * R[:,2]`（法向）；装配到节点等效力。
- **装配流程**：局部到全局旋转，累加 K/F；构建 dof_map、node_ids、elem_infos。

### 3.4 `solver/linear_static.py`
- **约束处理**：D 约束锁定，WARP 亦锁；K 加微小对角正则防奇异。  
- **解算**：使用 `np.linalg.solve`（或 lstsq 兜底），输出位移向量 `U`、反力 `R`、node_ids、dof_map。  
- **压力等效开关**：传递 `enable_pressure_equiv` 控制组装端力/曲率载荷。

### 3.5 `solver/creep.py`
- **时间步进**：从 TIME/NSUBST（或 CLI 覆盖）生成子步；每步调用组装更新载荷与初始应变。  
- **截面积分**：默认 12×2 积分点；每点存塑性/蠕变历史。  
- **本构更新**：弹性预测 → 塑性径向返回 → Norton 蠕变增量 → 更新内力并形成等效初始应变载荷。  
- **输出**：逐步写 VTK 并在 CLI 汇总为 PVD。

### 3.6 `solver/modal.py`
- **模态位移计算**：支持 m=0/1/2（呼吸、弯曲样、椭圆）；`compute_modal_displacement_at_angle` 给出指定 θ 的 (u_r,u_θ,u_z)。  
- **曲面展开**：平行运输构造光滑框架，按 n_circ 生成环点，叠加模态半径，输出点/quad/标量。

### 3.7 `vtk/surface.py`
- **路径提取**：`_extract_paths` 基于邻接表分解连通路径，支持分支与环。  
- **角点平滑**：在折点前后插入过渡点（长度 ≤ 半径），避免尖折。  
- **框架传播**：平行运输避免环向扭转跳变。  
- **模态与标量**：可传 `modal_state` 扩展几何；`scalar_map`/`emit_scalar` 控制着色，长度安全对齐；`any_modal` 控制标量命名。  
- **输出**：VTK POLYDATA/POLYGONS，点数据写标量（位移模量或模态径向）。

### 3.8 `vtk/writer.py`
- **中心线 VTK**：写节点坐标、线单元、点数据向量 `displacement`、可选反力/标量。

### 3.9 `vtk/pvd.py`
- **时间序列**：输出 PVD，列出 timestep 与对应 VTK 文件。

### 3.10 `cli.py`
- **子命令**：`run`、`creep`。  
- **参数**：`--pressure-equiv`、`--surface/--circ/--surface-scalar/--surface-scale`、蠕变的 `--time/--nsubsteps`。  
- **运行流程**：解析 CDB → 解算 → 写中心线 VTK → 可选曲面 → 蠕变写多步 VTK+PVD。

### 3.11 `run_all_cases.sh`
- 清理输出，依次运行 PIPE288/ELBOW290/PIPEMIX 的塑性与蠕变算例，统一开启 `--pressure-equiv` 并输出曲面与时间序列。

---

## 4. 算法说明（公式与策略）

- **椭圆化耦合**：局部刚度含弯曲与 OVAL 自由度耦合项，反映 Karman 效应，曲率越大柔化越明显。  
- **压力等效**：端面推力 + 曲率法向分布力（Bourdon），可用 `k_bourdon_scale` 调节；未细分环向压力片。  
- **温度插值**：对 E/ν/α/ρ 线性插值，热应变扣除 REFT，自重/压力等效用插值密度，避免温度段跳跃。  
- **截面积分**：每积分点独立更新塑性/蠕变，累积的初始应变回收为等效节点评价载荷，保证截面应力重分布。  
- **数值稳健性**：WARP 惩罚锁死；K 加微正则；曲面生成跳过零长度段、角点平滑、标量长度安全对齐。

---

## 5. 对标测试与验证

- **脚本**：`scripts/verify_physics.py`（单元/机制），`scripts/verify_results_data.py`（综合算例），`run_all_cases.sh`（全量回归）。
- **弹性/刚度**：悬臂直管挠度与梁理论一致（误差≈0）。  
- **Karman 耦合**：弯管受弯出现非零椭圆化，柔度增加，与 ELBOW290 物理趋势一致。  
- **蠕变**：恒载蠕变应变率符合 Norton 设定，时间步进稳定。  
- **压力等效**：开启/关闭 `--pressure-equiv` 对比，弯管位移方向与量级符合端推力+Bourdon 预期。  
- **VTK 一致性**：曲面点数与标量一致，ParaView 不再报 datasize mismatch；标量默认 `displacement_mag`，模态专用名 `modal_radial_disp`。

---

## 6. 已知限制与待扩展

- WARP 目前仅占位锁死，未实现真实翘曲场/刚度。  
- 压力壁面更精细等效（环向分片/局部扭转）尚未实现。  
- 几何非线性未做，大变形/大转动尚不支持。  
- 模态仅 m=0/1/2，未展开更高阶或真实翘曲模态。  
- 温度场仅均匀插值，未支持空间分布温度载荷。

---

## 7. 数据流与执行路径

1) **解析阶段**：`cli.py` 调 `parse_cdb` → 产出 Model（节点/单元/材料/截面/约束/载荷/时间）和温度表、Bourdon 系数。  
2) **组装阶段**：`assemble_linear_system` 读取 Model，按元素构建局部刚度与等效载荷，旋转到全局并累加 K/F，生成 dof_map、node_ids、elem_infos、props。  
3) **求解阶段**：`linear_static.solve_linear_static` 或 `creep.solve_creep_time_series`：应用约束→解算→组装反力→输出结果字典。  
4) **后处理**：`writer.write_vtk_polydata` 写中心线；如请求曲面，`surface.write_pipe_surface_quads` 用 centerline 顺序、半径、模态和标量展开；蠕变则 `pvd.write_pvd` 汇总序列。  
5) **验证/回归**：`run_all_cases.sh` 全量再生输出；`scripts/verify_*` 提供对标与数据检查。

---

## 8. 更详细的算法与实现细节

### 8.1 温度插值与热应变
- 表格：`MPTEMP` 给温度节点，`MPDATA` 给对应属性；插值采用分段线性：  
  $$ f(T)=f_i + (f_{i+1}-f_i)\frac{T-T_i}{T_{i+1}-T_i} $$
- 热应变：$\epsilon_{th} = \alpha(T-T_{ref})$；组装时作为初始应变扣除。  
- 自重/压力：使用插值后的密度 `dens_eff` 计算体力与端面推力。  
- 作用范围：静力与蠕变均使用同一插值，确保力学/热一致。

### 8.2 压力等效与 Bourdon 效应
- 端推力：$F = p \cdot A_i$，沿轴向正方向施加。  
- Bourdon：曲率 $\kappa$ 来自三点定圆；力密度 $q = -p \cdot A_i \cdot \kappa \cdot k_{bourdon}$，方向为曲率法线。  
- 参数化：`k_bourdon_scale` 默认 1，可在 CDB 用 `BHPADD,mat,scale` 指定。  
- 适用：仅对弯段曲率非零的单元；直管不产生 Bourdon。

### 8.3 截面积分与本构
- 积分点布置：环向均匀 12，径向 2（内/外）；可在代码调整。  
- 应变恢复：基于平截面假定和 OVAL 模态，得到每点广义应变→局部应力。  
- 塑性：径向返回，双线性等向硬化（BISO）。  
- 蠕变：Norton 律，按当前应力计算 $\dot\epsilon_{cr}$，显式时间步积累；步长受 NSUBST/TIME 控制。  
- 回收：将非弹性应变转为等效初始广义应变/曲率，叠加到右端项。

### 8.4 数值稳健性
- WARP 惩罚锁死，防止未定义翘曲导致奇异。  
- K 加微小对角正则；约束时移除对应行列并修正右端。  
- 曲面生成跳过零长度段，平滑折点，平行运输避免扭转跳变；标量长度检查，避免 VTK datasize mismatch。

### 8.5 可视化与标量字段
- 中心线 VTK：`displacement` 向量为主字段，反力可选。  
- 曲面 VTK：默认 `displacement_mag`（可关闭或缩放），若使用模态且未指定标量，则命名 `modal_radial_disp`。  
- 色阶：在 ParaView 中使用 “Rescale to Data Range”；可在 CLI 用 `--surface-scale` 放大位移色阶。

### 8.6 流程示例（静力 + 曲面）
```bash
PYTHONPATH=src python3 -m tubo run \
  --cdb "开放原子-管单元/ELBOW290_PLAST.cdb" \
  --out build/out_elbow290.vtk \
  --surface build/out_elbow290_surface.vtk \
  --circ 32 \
  --pressure-equiv \
  --surface-scalar disp \
  --surface-scale 5.0
```
- 解析 CDB → 装配含 Bourdon 的 K/F → 解算 → 写中心线 VTK → 展开曲面并写位移模量标量。  
- ParaView 打开曲面，选择 `displacement_mag` 着色并 Rescale。

### 8.7 流程示例（蠕变时间序列）
```bash
PYTHONPATH=src python3 -m tubo creep \
  --cdb "开放原子-管单元/PIPE288_CREEP.cdb" \
  --outdir build/pipe288_creep \
  --basename series \
  --pressure-equiv \
  --time 100.0 --nsubsteps 20
```
- 生成 20 个时间步的 VTK 与 `series.pvd`；ParaView 直接加载 PVD 播放。  
- 时间/子步可覆盖 CDB 内 TIME/NSUBST。

---

## 9. 对标测试与验证（细化）

- **弹性挠度**：悬臂直管比对 Euler-Bernoulli 理论，误差可忽略。  
- **Karman 耦合**：弯矩作用下椭圆化幅值非零，柔度增加；可关掉耦合（移除 OVAL 刚度项）作对比。  
- **压力等效**：开/关 `--pressure-equiv`，观察弯管顶推与轴向端推力；方向与量级符合 Bourdon 预期。  
- **蠕变**：恒载杆件的应变率与 Norton 设定一致；时间序列位移单调演化无震荡。  
- **VTK 完整性**：表面点数与标量一致，ParaView 不再提示 “datasize mismatch”。  
- **脚本化验证**：`bash run_all_cases.sh` 产物可用于回归；`scripts/verify_results_data.py` 读取 VTK 检查最大位移/范围。

---

## 10. 已知限制与扩展路线（再说明）

- **几何非线性**：未引入大转动/大变形；如需，可在局部坐标更新、二阶应变和接触力层面扩展。  
- **翘曲自由度**：当前锁死，未来可引入薄壁开口截面翘曲刚度或高阶 Fourier 模态。  
- **压力壁面精细化**：可细化为环向分片压力与局部扭矩等效，提高弯管压力分布精度。  
- **更高阶模态**：现支持 m=0/1/2，可扩展 m≥3 捕捉更复杂截面畸变。  
- **非均匀温度场**：目前仅均匀温度输入，未来可支持节点/单元温度载荷或热传导耦合。  
- **解算器**：当前直接用稠密求解；大模型可切换稀疏求解或迭代预条件以提升性能。

---

## 11. 应用指引（汇总）

- **全量回归**：`bash run_all_cases.sh`（默认开压力等效+曲面）。  
- **单例静力**：`tubo run --cdb <file> --out <vtk> [--pressure-equiv] [--surface ...] [--surface-scalar disp|none] [--surface-scale k] [--circ n]`。  
- **蠕变序列**：`tubo creep --cdb <file> --outdir <dir> --basename <name> [--pressure-equiv] [--time T --nsubsteps N]`。  
- **ParaView 查看**：中心线看 `displacement`；曲面看 `displacement_mag` 或 `modal_radial_disp`，记得 Rescale；需要放大可改 `--surface-scale`。  
- **验证脚本**：`python3 scripts/verify_physics.py`（机制）；`python3 scripts/verify_results_data.py`（回归数据）。
