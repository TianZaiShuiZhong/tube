# 求解器实现思路与开发路线图 (Implementation Strategy & Roadmap)

本文档面向开发者，给出从解析到求解再到可视化的实现思路、阶段路线和关键算法要点，并梳理尚未完成的扩展方向。

---

## 1. 总体策略

### 1.1 核心理念：降维 + 物理保真
- 用 1D 管单元替代 3D 实体/壳，加入椭圆化（OVAL）与占位翘曲（WARP）自由度，保持关键物理（Karman、Bourdon、蠕变/塑性）。
- 截面积分点取代宏观截面刚度，逐点进行本构更新，提升材料非线性精度。

### 1.2 解耦架构
- 解析层：CDB → Model（几何/材料/温度表/载荷/时间）。
- 组装层：局部 16×16 刚度、温度插值、压力/自重等效 → 全局 K/F。
- 求解层：线性静力或时间步进（蠕变），约束/正则化后解方程。
- 可视化层：中心线 VTK + 曲面四边形（可叠加模态/标量）。

---

## 2. 阶段路线图（已完成与未完成功能）

### Phase 1 基础通路
1) CDB 解析（NBLOCK/EBLOCK/MP/SECTYPE/SECDATA/D/F 等）。  
2) 6-DOF 梁刚度、局部-全局旋转，中心线 VTK 导出。

### Phase 2 管单元增强
1) 扩展为 8 DOF（OVAL + WARP 占位锁死）。  
2) Karman 耦合：弯曲-椭圆化刚度耦合；弯管曲率提取（定圆）。  
3) 压力等效：端面推力 + 曲率法向 Bourdon 载荷（`k_bourdon_scale` 可调）。  
4) 角点平滑 + 曲面重构，四边形管壁输出。

### Phase 3 材料非线性
1) 截面积分点（默认 12×2），应变恢复、内力回收。  
2) 塑性：BISO 径向返回；蠕变：Norton 律显式步进。  
3) 蠕变时间序列输出（VTK+PVD）。

### Phase 4 温度与标量/模态可视化
1) MPTEMP/MPDATA 分段插值 E/ν/α/ρ，热应变扣 REFT，自重/压力用插值密度。  
2) 曲面标量控制：`--surface-scalar {disp,none}`、`--surface-scale k`；模态叠加时安全对齐标量。

### Phase 5 验证与脚本
1) `run_all_cases.sh` 一键回归（全算例 + 压力等效 + 曲面）。  
2) `scripts/verify_physics.py` 机制验证；`verify_results_data.py` 回归检查。

### 待办/扩展
- 真实 WARP 刚度/模态；压力壁面更精细等效（环向分片/扭矩）；几何非线性（大转动）；高阶傅里叶模态；空间温度场。

---

## 3. 模块与代码导读（按执行顺序）

### 3.1 解析层 `src/tubo/cdb/parser.py`
- 状态机读取：NBLOCK/EBLOCK/MP/MPTEMP/MPDATA/SECTYPE/SECDATA/D/F/SFE/TB,CREE/TIME/NSUBST/BHPADD。
- 结果：Model（节点/单元/材料/截面/约束/集中力/时间控制），材料含温度表、Bourdon 系数。

### 3.2 数据模型 `src/tubo/model.py`
- 实体：Node/Element/Material/SectionPipe/ModalConfig/FourierModalState。  
- DOF 顺序：UX, UY, UZ, RX, RY, RZ, OVAL, WARP（WARP 惩罚锁死）。  
- 材料复制保留塑性/蠕变/温度表/Bourdon。

### 3.3 组装层 `src/tubo/solver/assembly.py`
- 几何：三点定圆求曲率与 Frenet 框架；直段曲率零。  
- 刚度：16×16 局部矩阵，含弯曲-OVAL 耦合；WARP 对角惩罚。  
- 温度：线性插值 E/ν/α/ρ；热应变 α·(T-REFT)；自重/压力用插值密度。  
- 载荷：自重；端面推力；Bourdon 力密度 `-p*Ai*curvature*k_scale` 沿曲率法线。  
- 装配：局部→全局旋转，累加 K/F，生成 dof_map、node_ids、elem_infos。

### 3.4 静力求解 `src/tubo/solver/linear_static.py`
- 处理约束（含 WARP 锁定），K 加微正则；求解 U、反力 R；压力等效开关透传给组装。

### 3.5 蠕变求解 `src/tubo/solver/creep.py`
- 时间步：TIME/NSUBST 或 CLI 覆盖。  
- 截面积分：12×2 默认；每点存塑性/蠕变历史。  
- 本构：弹性预测 → BISO 返回 → Norton 蠕变增量；回收等效初始应变载荷。  
- 输出：多步 VTK + PVD。

### 3.6 模态与曲面 `src/tubo/solver/modal.py` / `src/tubo/vtk/surface.py`
- 模态：m=0/1/2 位移展开，平行运输构造稳定框架。  
- 曲面：提取路径、角点平滑、生成四边形；可叠加模态几何与标量，长度对齐防 VTK 警告；标量名依场景选择。

### 3.7 VTK 与 CLI
- `vtk/writer.py`：中心线 VTK（`displacement`）。  
- `vtk/pvd.py`：时间序列。  
- `cli.py`：`run`/`creep`，参数含 `--pressure-equiv`、`--surface*`、蠕变时间控制。  
- `run_all_cases.sh`：全算例回归（塑性+蠕变，带曲面与压力等效）。

---

## 4. 关键算法要点

- **椭圆化耦合**：局部刚度含弯曲-OVAL 耦合项，模拟 Karman 柔化。  
- **压力等效**：端推力 + 曲率法向 Bourdon（可调系数），直管无 Bourdon。  
- **温度插值**：MPTEMP/MPDATA 线性插值 E/ν/α/ρ；热应变扣 REFT；自重/压力用插值密度。  
- **截面积分**：环向×径向积分点，逐点本构更新并回收等效载荷。  
- **稳健性**：WARP 锁死，K 正则；曲面跳过零长度段、角点平滑、标量长度对齐。

---

## 5. 对标与验证

- `run_all_cases.sh`：回归输出中心线/曲面/蠕变序列。  
- 弹性挠度：悬臂直管 vs 梁理论。  
- Karman：弯管弯曲出现非零 OVAL 幅值，柔度增加。  
- 压力等效：开/关 `--pressure-equiv` 对比弯管响应。  
- 蠕变：恒载应变率符合 Norton；时间序列单调演化。  
- VTK 完整性：曲面点/标量一致，ParaView 无 datasize mismatch。

---

## 6. 待扩展路线（建议优先级）

1) **几何非线性**：Co-rotational/更新坐标，支持大转动/大变形。  
2) **真实 WARP/高阶模态**：引入翘曲刚度或 m≥3 傅里叶模态。  
3) **压力壁面精细等效**：环向分片压力与局部扭矩。  
4) **非均匀温度场/热传导耦合**：节点/单元温度载荷或瞬态温度。  
5) **稀疏/迭代解算器**：大规模模型性能优化。  
6) **更多验证脚本**：对标 ANSYS 导出 CSV，批量误差统计。

---

## 7. 快速使用与开发提示

- 全量：`bash run_all_cases.sh`。  
- 静力：`tubo run --cdb ... --out ... [--pressure-equiv] [--surface ...] [--surface-scalar disp|none] [--surface-scale k] [--circ n]`。  
- 蠕变：`tubo creep --cdb ... --outdir ... --basename ... [--pressure-equiv] [--time T --nsubsteps N]`。  
- ParaView：中心线看 `displacement`；曲面看 `displacement_mag` 或 `modal_radial_disp`，必要时 Rescale；色阶不足用 `--surface-scale` 放大。  
- 开发调试：优先跑小算例（PIPE288/ELBOW290），配合 `verify_physics.py` 检查机制。  
- 提交前：运行 `run_all_cases.sh` 确认无异常警告（尤其 VTK datasize）。
