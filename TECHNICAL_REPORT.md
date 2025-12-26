# 基于管单元的高温管道蠕变高效求解 - 技术实现报告

## 1. 引言 (Introduction)

本项目的目标是开发一套高效的有限元求解器，专门用于模拟高温高压环境下管道系统（直管与弯管）的力学行为，特别是蠕变（Creep）和塑性（Plasticity）效应。为了在保证精度的同时显著提升计算效率，本项目采用了**基于傅里叶模态扩展的管单元（Pipe Element）**理论，而非传统的实体或壳单元。该方法通过在梁单元的基础上引入截面变形自由度（如椭圆化），能够以一维单元的计算代价捕捉三维复杂的变形模式，对标 ANSYS ELBOW290/PIPE288 单元。

## 2. 理论基础 (Theoretical Basis)

### 2.1 7自由度管单元模型 (7-DOF Pipe Element)

传统的欧拉-伯努利（Euler-Bernoulli）梁单元每个节点具有 6 个自由度（3 个平动 $u_x, u_y, u_z$ 和 3 个转动 $\theta_x, \theta_y, \theta_z$）。为了模拟薄壁管在弯曲载荷下的截面扁平化（Ovalization），我们引入了第 7 个自由度：**椭圆化幅值 (Ovalization Amplitude, $a$)**。

节点自由度向量定义为：
$$ \mathbf{d} = [u_x, u_y, u_z, \theta_x, \theta_y, \theta_z, a]^T $$

### 2.2 环向傅里叶模态 (Circumferential Fourier Modes)

管壁的径向位移 $w(\theta)$ 和切向位移 $v(\theta)$ 采用傅里叶级数展开。本项目主要关注 $m=2$ 的椭圆化模态：

$$ w(\theta) = a \cos(2\theta) $$
$$ v(\theta) = \frac{1}{2} a \sin(2\theta) $$

其中 $\theta$ 为环向角度。该模态描述了圆形截面变为椭圆形的变形特征。

### 2.3 弯管 Karman 效应与刚度耦合 (Karman Effect & Coupling)

弯管在承受弯矩时，由于截面的扁平化，其抗弯刚度会显著降低（Karman 效应）。我们在刚度矩阵中引入了弯曲自由度 $\theta_z$ 与椭圆化自由度 $a$ 之间的耦合项。

耦合刚度推导基于能量法，耦合势能项为：
$$ U_{couple} = - \int_0^L E I \left( \frac{3}{2 R_{curv}} \right) \kappa_z(x) a(x) dx $$

其中 $R_{curv}$ 是弯管的曲率半径，$\kappa_z \approx \frac{d\theta_z}{dx}$ 是梁的曲率。这导致单元刚度矩阵中出现非对角项，使得弯矩会自动诱发截面变形。

### 2.4 压力载荷等效 (Pressure Load Equivalence)

为了在梁单元上正确模拟内压 $P$ 的作用，我们实现了两类等效载荷：

1.  **端部盲板力 (End-cap Thrust)**：
    $$ F_{axial} = P \cdot \pi r_i^2 $$
    沿管道轴向施加，模拟闭口端承受的压力。

2.  **弯管分布载荷 (Distributed Curvature Load)**：
    在弯管处，内压会在管壁产生合力，试图将弯管“拉直”（Bourdon 效应）。等效分布力 $\mathbf{q}$ 为：
    $$ \mathbf{q} = P \cdot (\pi r_i^2) \cdot \vec{\kappa} $$
    其中 $\vec{\kappa}$ 是中心线的曲率矢量。该载荷被转化为单元节点的等效剪力。

## 3. 数值算法实现 (Numerical Implementation)

### 3.1 截面积分点算法 (Cross-section Integration)

为了处理材料非线性（塑性和蠕变），求解器不直接使用截面刚度（如 $EA, EI$），而是采用**数值积分**方法。

*   **离散化**：将管截面沿环向（Circumferential）和径向（Radial/Thickness）划分为积分点网格。默认配置为 12（环向）× 2（壁厚）= 24 个积分点。
*   **应变恢复**：根据节点位移计算单元中心线的广义应变（轴向应变 $\epsilon_0$、曲率 $\kappa_y, \kappa_z$），然后利用平截面假定计算每个积分点的总应变：
    $$ \epsilon_{total}(y, z) = \epsilon_0 - y \kappa_z + z \kappa_y $$

### 3.2 弹塑性与蠕变耦合求解 (Coupled Plasticity & Creep)

采用**算子分裂法 (Operator Splitting)** 在时间域上进行积分。在每个时间步 $\Delta t$ 内：

1.  **弹性预测 (Elastic Predictor)**：
    计算试探应力：
    $$ \sigma_{trial} = E (\epsilon_{total} - \epsilon_{cr}^{old} - \epsilon_{pl}^{old}) $$

2.  **塑性修正 (Plastic Corrector) - 径向返回映射**：
    检查屈服准则（Von Mises 退化为单轴）：$f = |\sigma_{trial}| - \sigma_{yield}$。
    若 $f > 0$，则进行塑性流动：
    $$ \Delta \lambda = \frac{f}{E + H} $$
    $$ \Delta \epsilon_{pl} = \Delta \lambda \cdot \text{sign}(\sigma_{trial}) $$
    $$ \sigma_{corrected} = \sigma_{trial} - E \Delta \epsilon_{pl} $$
    其中 $H$ 为硬化模量。

3.  **蠕变更新 (Creep Update)**：
    基于修正后的应力，利用 Norton 定律计算蠕变应变增量：
    $$ \dot{\epsilon}_{cr} = C_1 |\sigma|^{C_2} \text{sign}(\sigma) $$
    $$ \Delta \epsilon_{cr} = \dot{\epsilon}_{cr} \cdot \Delta t $$

4.  **应力积分 (Stress Integration)**：
    将各积分点的非弹性应变（$\epsilon_{pl} + \epsilon_{cr}$）积分回截面，得到等效的**初始广义应变**和**初始曲率**，作为下一时间步的附加载荷向量（右端项）输入到全局方程中。

## 4. 软件架构 (Software Architecture)

项目采用 Python 开发，模块化设计：

*   `src/tubo/cdb/parser.py`: **解析层**。读取 ANSYS CDB 文件，构建几何拓扑和材料模型。
*   `src/tubo/model.py`: **数据层**。定义节点、单元、材料、模态配置等数据类 (Dataclasses)。
*   `src/tubo/solver/assembly.py`: **线性求解层**。负责 14x14 刚度矩阵组装、坐标变换、Karman 耦合项计算。
*   `src/tubo/solver/creep.py`: **非线性求解层**。管理时间步进循环、积分点状态更新、本构模型计算。
*   `src/tubo/vtk/`: **可视化层**。
    *   `writer.py`: 输出中心线结果。
    *   `surface.py`: 核心可视化模块。将 1D 单元结果结合模态幅值 $a$，通过几何变换重构出 3D 管壁曲面（Quad Mesh），直观展示截面变形。

## 5. 验证与对标 (Verification)

### 5.1 弯管效应验证
对比开启/关闭 Karman 耦合项的结果：
*   **开启耦合**：弯管在弯矩作用下发生明显的截面扁平化（通过 VTK 曲面云图观察到），且整体柔度增加（位移变大）。
*   **关闭耦合**：截面保持圆形，表现为刚性梁特性。
结果与 ANSYS ELBOW290 单元特性一致。

### 5.2 蠕变验证
使用 `PIPE288_CREEP.cdb` 进行 100 步时间历程分析：
*   轴向位移随时间呈非线性增长，符合 Norton 蠕变律的预期。
*   积分点算法稳定，未出现数值震荡。

### 5.3 塑性验证
使用 `ELBOW290_PLAST.cdb`（含 BISO 材料）：
*   在高应力区，应力被限制在屈服面附近，塑性应变累积，验证了径向返回算法的正确性。

## 6. 总结 (Conclusion)

本项目成功实现了一套轻量级但功能完备的高温管道分析软件。通过引入 7-DOF 管单元理论和截面积分算法，在极低的计算代价下（相比实体单元减少 90% 以上自由度），实现了对复杂弯管效应、塑性和蠕变的精确模拟。软件具备完整的工业级接口（CDB 输入 / VTK 输出），可作为高效的管道设计与评估工具。
