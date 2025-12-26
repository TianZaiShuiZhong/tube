# 求解器实现思路与开发路线图 (Implementation Strategy & Roadmap)

本文档详细阐述了“基于管单元的高温管道蠕变高效求解器”的开发策略、分阶段实现路径以及关键技术难点的攻克思路。

## 1. 开发策略 (Development Strategy)

### 1.1 核心理念：降维求解 (Dimensionality Reduction)
*   **问题**：高温管道蠕变通常需要使用 3D 实体单元或壳单元进行精细化建模，计算量巨大，难以满足快速评估需求。
*   **策略**：采用 **广义梁理论 (Generalized Beam Theory)**。将 3D 变形特征（如截面椭圆化、翘曲）凝聚为 1D 梁节点上的附加自由度（DOF）。
*   **优势**：自由度数量减少 2-3 个数量级，同时保留了关键的结构非线性特征。

### 1.2 解耦设计 (Decoupled Architecture)
*   **几何非线性**（大变形/Karman效应）通过 **刚度矩阵修正** 处理。
*   **材料非线性**（塑性/蠕变）通过 **截面积分点 (Integration Points)** 处理。
*   两者在全局迭代循环中解耦，便于独立开发和调试。

## 2. 分阶段实现路线 (Phased Implementation Roadmap)

### 第一阶段：基础设施与线性求解 (Phase 1: Infrastructure & Linear Solver)
**目标**：打通数据流，实现基本的 6自由度梁单元求解。
1.  **CDB 解析器开发**：
    *   实现 `NBLOCK`, `EBLOCK` 解析，构建网格拓扑。
    *   解析 `MP`, `SECTYPE` 获取材料与截面属性。
2.  **6-DOF 梁单元刚度矩阵**：
    *   基于 Euler-Bernoulli 梁理论，构建 $12 \times 12$ 刚度矩阵。
    *   实现坐标变换（局部坐标系 -> 全局坐标系）。
3.  **VTK 可视化管道**：
    *   初步实现将计算结果（中心线位移）导出为 `.vtk` 文件。

### 第二阶段：几何非线性与管单元增强 (Phase 2: Geometric Nonlinearity & Pipe Enhancement)
**目标**：捕捉弯管特有的力学行为（Karman 效应）。
1.  **引入第 7 自由度 (Ovalization)**：
    *   扩展节点自由度至 7 个 ($u_x, u_y, u_z, \theta_x, \theta_y, \theta_z, a$)。
    *   推导环向刚度 $k_{oval}$。
2.  **Karman 耦合实现**：
    *   推导弯曲-椭圆化耦合项 $K_{c}$。
    *   组装 $14 \times 14$ 单元刚度矩阵。
3.  **压力载荷修正**：
    *   实现端部盲板力（轴向拉伸）。
    *   实现弯管分布力（Bourdon 效应，抵抗弯曲）。

### 第三阶段：材料非线性架构 (Phase 3: Material Nonlinearity Architecture)
**目标**：处理塑性和蠕变，引入“积分点”概念。
1.  **截面离散化**：
    *   开发 `IntegrationPoint` 类。
    *   实现截面网格生成（环向 $\times$ 径向）。
2.  **应变恢复与应力积分**：
    *   实现 `节点位移 -> 截面广义应变 -> 积分点应变` 的映射。
    *   实现 `积分点应力 -> 截面内力 -> 节点力` 的逆映射。

### 第四阶段：本构模型与时间步进 (Phase 4: Constitutive Models & Time Stepping)
**目标**：实现完整的物理场耦合求解。
1.  **塑性回路 (Plasticity Loop)**：
    *   实现 J2 流动理论。
    *   实现 **径向返回映射 (Radial Return Mapping)** 算法处理屈服面。
    *   解析 ANSYS `TB,BISO` 参数。
2.  **蠕变回路 (Creep Loop)**：
    *   实现 Norton 蠕变律 ($\dot{\epsilon} = C_1 \sigma^{C_2}$)。
    *   实现时间步进算法（显式/半隐式 Euler）。
3.  **全局牛顿迭代**：
    *   处理非线性残差力，确保每个时间步的平衡。

### 第五阶段：高级可视化与验证 (Phase 5: Advanced Visualization & Validation)
**目标**：直观展示截面变形，验证精度。
1.  **曲面重构**：
    *   利用计算出的椭圆化幅值 $a$，在后处理中动态修正管壁网格坐标。
    *   生成包含真实变形的 3D 表面 VTK。
2.  **对标验证**：
    *   与 ANSYS 商业软件结果（ELBOW290/PIPE288）进行对比。

## 3. 关键技术难点与解决方案 (Key Challenges & Solutions)

| 难点 | 解决方案 |
| :--- | :--- |
| **弯管刚度过大** | 引入 Karman 耦合项，模拟截面扁平化导致的刚度折减。 |
| **塑性迭代不收敛** | 采用径向返回映射算法（Radial Return Mapping），保证应力状态始终位于屈服面上。 |
| **蠕变时间步长敏感** | 采用自适应时间步长或较小的固定步长，结合 Predictor-Corrector 策略。 |
| **可视化抽象** | 开发专门的 Surface 重构模块，将 1D 单元数据映射回 3D 几何。 |

## 4. 未来优化方向 (Future Optimization)

1.  **温度相关性**：完善 `MPTEMP` 插值逻辑，支持变温场下的材料属性更新。
2.  **大转动理论**：引入 Co-rotational 坐标系，支持几何大变形（不仅是小应变）。
3.  **多模态扩展**：引入 $m=3, 4$ 等更高阶傅里叶模态，模拟更复杂的截面翘曲。
