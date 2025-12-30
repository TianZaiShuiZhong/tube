# 求解器精度验证报告 (Solver Verification Report)

本文档说明求解器的验证策略、案例设置、对标方法与当前结果，确保弹性、Karman 耦合、压力等效、温度插值以及蠕变/塑性机制均按预期工作。

---

## 1. 验证策略

- **理论对标**：使用有解析/半解析解的经典模型（悬臂梁挠度、恒载蠕变率、热膨胀量级）。  
- **机制验证**：针对关键特性（Karman 耦合、Bourdon 压力等效、温度插值、截面积分本构）构造最小案例。  
- **系统级回归**：对赛题完整 CDB 模型跑全流程（含压力等效与曲面），比对位移/色阶/VTK 完整性。  
- **脚本化检查**：`scripts/verify_physics.py`（机制单测），`scripts/verify_results_data.py`（回归数据检查），`run_all_cases.sh`（全量再生）。

---

## 2. 自动化脚本与使用

- 机制验证：`PYTHONPATH=src python3 scripts/verify_physics.py`  
  - 覆盖悬臂挠度、Karman 椭圆化响应、恒载蠕变率等。  
- 回归检查：`PYTHONPATH=src python3 scripts/verify_results_data.py`  
  - 读取 `build/` 中生成的 VTK，提取最大位移/范围，用于对标和异常检测。  
- 全量再生：`bash run_all_cases.sh`  
  - 生成中心线 VTK、曲面 VTK、蠕变 PVD/VTK，默认开启压力等效。

---

## 3. 测试案例详解

### 案例 1：悬臂管弹性弯曲（线性刚度）
- **目标**：验证 6-DOF 基础刚度矩阵与局部-全局变换。  
- **设置**：L=1000 mm，D=100 mm，t=5 mm，一端固支，端横向力 F=1000 N。  
- **理论值**：$\delta = \frac{F L^3}{3 E I}$，I = π/64 (D^4-(D-2t)^4)。  
- **结果**：数值挠度与理论一致（误差≈0），线性刚度正确。

### 案例 2：弯管 Karman 耦合
- **目标**：验证弯曲-椭圆化耦合（Karman 柔化）。  
- **设置**：90° 弯管，曲率半径 R=500 mm，施加面内弯矩。  
- **标准**：出现非零椭圆化幅值 a，且随弯矩线性放大；关掉耦合则 a→0。  
- **结果**：a 非零且线性随弯矩变化，证明耦合项生效。

### 案例 3：压力等效（端推力 + Bourdon）
- **目标**：验证压力等效开关与曲率载荷。  
- **设置**：弯管受内压 p，启用 `--pressure-equiv`。  
- **标准**：端面产生轴向推力，弯管沿曲率法线出现 Bourdon 推力；关掉开关则消失。  
- **结果**：开启后位移方向/量级符合端推力+Bourdon 预期；直管无 Bourdon。

### 案例 4：温度插值与热应变
- **目标**：验证 MPTEMP/MPDATA 插值与热应变扣除。  
- **设置**：多温度点给出 E/ν/α/ρ，设 REFT。  
- **标准**：热应变 α (T-REFT) 与刚度/自重密度一致；跨温度段平滑无突跳。  
- **结果**：插值平滑，热应变方向正确，自重随密度变化。

### 案例 5：蠕变本构积分
- **目标**：验证 Norton 蠕变率实现与时间步进稳定性。  
- **设置**：单轴恒应力杆件或赛题蠕变算例，TIME/NSUBST 或 CLI 覆盖。  
- **标准**：应变率与设定 C1/C2 一致，时间序列位移单调演化无震荡。  
- **结果**：单元测试与时间序列均符合预期。

### 系统级算例（赛题 CDB）
- `PIPE288_PLAST/ELBOW290_PLAST/PIPEMIX_PLAST`：塑性/热膨胀/压力等效，检查最大位移、色阶合理性。  
- `PIPE288_CREEP/ELBOW290_CREEP/PIPEMIX_CREEP`：蠕变时间序列，检查位移随时间单调演化。

---

## 4. 主要修正记录（与验证相关）

1) 压力等效：弯管改用曲率法向 Bourdon 分布载荷（保留端推力），可调 `k_bourdon_scale`。  
2) 温度：MPTEMP/MPDATA 线性插值；热应变扣 REFT；自重/压力等效用插值密度。  
3) VTK 完整性：曲面点数与标量对齐，消除 ParaView “datasize mismatch” 警告。  
4) WARP 占位锁定：防止未定义翘曲导致奇异。  
5) 曲面标量：支持 `--surface-scalar/--surface-scale`，可关闭或放大色阶。

---

## 5. 当前验证结论

- **弹性与刚度**：悬臂挠度对标理论误差可忽略。  
- **Karman/Bourdon**：弯管出现椭圆化与压力引起的法向推力，开关可对比。  
- **温度**：插值平滑，热应变与载荷一致。  
- **蠕变**：Norton 率实现正确，时间序列稳定。  
- **VTK 输出**：中心线与曲面文件一致可用，ParaView 无尺寸警告。

---

## 6. 使用建议（验证相关）

- 回归：`bash run_all_cases.sh`（默认开压力等效+曲面）。  
- 机制单测：`python3 scripts/verify_physics.py`。  
- 数据检查：`python3 scripts/verify_results_data.py`（查看最大位移/范围）。  
- ParaView：载入曲面后 “Rescale to Data Range”；色阶过暗可用 `--surface-scale` 放大。  
- 对比压力等效：同一算例开/关 `--pressure-equiv`，观察弯管位移变化。  
- 若需 ANSYS 对标，可在 VTK 输出上用自备 CSV 对齐节点位移进行误差统计。
