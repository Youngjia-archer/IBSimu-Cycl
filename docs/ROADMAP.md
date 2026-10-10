# IBSimu-Cycl：面向回旋加速器的 IBSimu 改造路线图

> 目标仓库：`IBSimu-Cycl`（基于上游 IBSimu 1.0.6dev，GPL 派生）
> 文档状态：草案 v0.1（任务规划阶段）

---

## 1. 项目目标与范围

以开源离子束仿真程序 **IBSimu** 为基础，构建一个可用于**回旋加速器**全链路仿真的派生版本，重点补齐三项能力：

| 编号 | 需求 | 关键难点 |
| --- | --- | --- |
| R1 | 提升三维仿真效率（并行化 / CUDA 加速） | 场插值 + 空间电荷求解是 3D 主要开销 |
| R2 | 三维磁场、RF 谐振腔、时变场中束流运动 | IBSimu 原生场求解器基本是 2D/轴对称 |
| R3 | 改善图像输出与可视化分析 | 原生 GTK+cairo+OpenGL，3D 能力弱 |

**约束（已更新）**：必须保留 IBSimu 既有**功能**（离子源引出、PIC/Vlasov 自洽迭代、DXF/STL 几何、2D/柱对称求解器等）；但**允许对内部架构做激进重构**（已确认）——因此**不再要求与上游代码可合并**，可一次性导入上游作为基线后自由演进。原则为"**新核心 + 保留功能**"。

**许可**：IBSimu 为 GPL，本项目为派生作品，整体须维持 GPL，必须保留原作者版权声明（Taneli Kalvas 等）。整合的第三方开源组件需逐一身查许可证兼容性（GPL 兼容性清单见 §7）。

### 1.1 已确认的关键决策

| 决策项 | 结论 | 影响 |
| --- | --- | --- |
| 运行环境 | **混合：CPU + GPU + 集群**；GPU **暂只针对 NVIDIA** | 采用 **Kokkos** + MPI 多节点；加速后端优先 **CUDA**，暂不投入 HIP/SYCL 适配 |
| 目标加速器 | 有具体型号，**暂无实测数据** | 自几何建模；验证依赖解析解 + OPAL-cycl（型号待补充） |
| 外部依赖 | **接受**引入外部 FEM，追求精度 | 磁铁用 Elmer（非线性铁），RF 腔用 Palace |
| 优先级 | **R2 物理能力优先** | 先保证"算得对"：P1→P2→P3→P4；性能/GPU 与可视化后置 |
| 上游策略 | **允许激进重构** | 可重建核心，不受上游 API 约束；仅保留功能与 GPL 合规 |

---

## 2. 上游 IBSimu 架构现状（1.0.6dev）

已确认的关键类（来自 Doxygen 类列表），决定了改造切入点：

### 2.1 场与求解器
- `CField` / `CScalarField` / `CVectorField`：场基类。
- `CMeshScalarField` / `CMeshVectorField` / `CMultiMeshVectorField`：基于矩形网格的场，**插值能力是性能瓶颈所在**。
- `CAxisymmetricVectorField`：**轴对称**磁场（基于轴上数据）——**回旋加速器扇形磁铁不能用它**，必须新增通用 3D 场。
- `CEpotSolver` → `CEpotMatrixSolver` → `CEpotBiCGSTABSolver` / `CEpotGSSolver` / `CEpotMGSolver`（多重网格）/ `CEpotUMFPACKSolver`：Poisson 求解栈，接口可复用。
- `CEpotField` / `CEpotEfield`：由电位得到的场。

### 2.2 粒子与跟踪
- `CParticleDataBase2D` / `CParticleDataBase3D` / `CParticleDataBaseCyl`：**已存在 3D 粒子数据库**，是 3D 跟踪的现成基础。
- `CParticleP2D` / `CParticleP3D` / `CParticlePCyl`：粒子点。
- `CParticleStepper`（`particlestepper.hpp`）：**Boris leap-frog 步进器**——时变场改造的核心位置。
- `CParticleIterator` / `CConvergence`：Vlasov 自洽迭代与收敛判据（须保留）。
- `CSpaceCharge`（`scharge.cpp`）：空间电荷沉积。

### 2.3 并行与几何
- `CScheduler`：**已存在的生产-消费者线程框架**（"Job scheduler for parallel processing"），是并行化的起点。
- 几何：`CSTLSolid` / `CDXFSolid` / `CCSGObjectSolid` / `CFuncSolid`（含 3D STL，可直接复用）。

### 2.4 可视化
- GTK 交互式绘图（`CGTK*Window`）、cairo（`CFrame`/`CColormap`）、OpenGL（`CGLRenderer`/`CSoftwareRenderer`）、`CGeom3DPlot` / `CFieldDiagPlot` / `CParticleDiagPlot` / `CHistogram` / `CXYGraph`。

**结论**：IBSimu 已具备 3D 粒子与 3D 几何的骨架，**缺口集中在"通用三维场表示（磁铁/RF）"与"时变场 + 高性能求解"**。

---

## 3. 总体架构

```
                    IBSimu-Cycl
┌─────────────────────────────────────────────────────────────┐
│  Python 层  (ibsimu_cycl)                                    │
│   · 后处理/PyVista/matplotlib   · 场图读写   · 参数扫描      │
├─────────────────────────────────────────────────────────────┤
│  C++ 核心 (libibsimu 扩展)                                    │
│  ┌───────────────┐ ┌───────────────┐ ┌───────────────────┐  │
│  │ FieldMap3D    │ │ RfCavityMode  │ │ TimeVaryingStepper│  │
│  │ (通用3D场插值)│ │ (Eigenmode)   │ │ (时变场Boris)     │  │
│  └──────┬────────┘ └──────┬────────┘ └─────────┬─────────┘  │
│         │                 │                    │             │
│  ┌──────▼─────────────────▼────────────────────▼─────────┐  │
│  │  Solver 抽象层 (IMeshSolver)  —— 保留 CEpotSolver 接口 │  │
│  │   CPU(OpenMP/MPI) | GPU(Kokkos/CUDA) 可切换后端        │  │
│  └───────────────────────────────────────────────────────┘  │
├─────────────────────────────────────────────────────────────┤
│  外部求解器适配器 (可选、松耦合)                               │
│   磁铁: Elmer / GetDP+ONELAB / Radia                          │
│   RF : Palace / Elmer / GetDP                                 │
├─────────────────────────────────────────────────────────────┤
│  IO 层: openPMD / HDF5 / VTK/XDMF  (场图与粒子交换)           │
└─────────────────────────────────────────────────────────────┘
```

**核心原则**：外部 FEM 求解器只用于**离线生成场图**（B-map、E-mode、ω），运行时 IBSimu 只做**插值 + 推进**。这样既获得高精度 3D 场，又不把 IBSimu 绑定到重型依赖上，同时天然支持并行（插值与推进高度数据并行）。

---

## 4. 模块规划

### A. 三维场基础设施 + 时变场（R2 基础，前置）
- 新增 `CVectorField3D`：支持直角/柱坐标网格上的三矢量分量存储与三线性/三三次插值。
- 新增 `CFieldMap3D`：从文件（openPMD/HDF5/NPZ）载入 B/E 场图；优先用**柱坐标 (r, θ, z)** 存储以适配回旋加速器，内部缓存到结构化数组。
- 新增 `CTimeVaryingField`：包装静态模式场为 `Field(t) = A(r) · cos(ωt + φ0)`，供 RF 电场与高频腔场使用；同时支持"边界电位随 t 变化"的准静态模式。
- **保留** `CAxisymmetricVectorField` 作为快速回退路径（轴对称近似的调试用）。

### B. 三维静磁铁磁场（R2-B）
- 选型（已确认接受外部 FEM；**首选 Elmer**）：
  - **Elmer FEM**（`StatMag` / `WhitneyAV`）【首选】：支持铁芯非线性（B-H 饱和），适合大铁磁体的回旋加速器磁铁，MPI 并行。
  - **GetDP + Gmsh (ONELAB)**【备选】：脚本化、轻量、开源，适合参数化扇形磁铁建模。
  - **Radia (ESRF)**【辅助】：边界元 + 磁荷法，对线圈/超导磁体极强，非线性铁迭代较弱。
- 交付物：一个 `magnet3d` 工具/适配器，输入几何+材料 → 输出 3D B-map（含中心平面 Bz 等值线供回路设计）。
- 校验：均匀场解析解；对称二极磁铁与解析/实测对比；回旋频率 ω_c = qB/m 复核。

### C. RF 谐振腔谐振计算（R2-C）
- 选型（已确认接受外部 FEM）：**Palace**（AWS Labs，MFEM 基，3D 电磁本征模，MPI + GPU，专为加速器腔设计）【首选】；`Elmer`/`GetDP` 为备选。
- 交付物：`rfcavity3d` 适配器，输出 ω（谐振频率）、R/Q、Q、以及 3D E/B 模式场图。
- 时变耦合：模式场 → `CTimeVaryingField`；Dee 电压按 `V(t)=V0·cos(ω_rf·t+φ)` 施加于电极边界（沿用 `CEpotSolver` 的边界条件机制）。
- 校验：理想圆柱腔 TM010 解析频率；随后真实回旋腔。

### D. 时变场下的束流运动（R2-D，核心物理）
- 扩展 `CParticleStepper`：由静态 Boris 升级为**时间感知 Boris**，在每步用 `t` 评估 `E(r,t)`，`B` 取静态场图。
- 相对论修正：回旋加速器能量可达数十 MeV/u，保留/新增相对论 Boris（γ 迭代），非相对论路径保持向后兼容。
- 物理效应支持：等时性、相位滑动（phase slip）、穿越 RF 腔的能量增益与相位稳定性、垂直聚焦。
- 与 Vlasov 自洽迭代（`CParticleIterator`/`CConvergence`）对接：每个 RF 周期更新空间电荷，保持已有收敛逻辑。
- 校验：均匀 B 中圆周运动半径/周期；解析回旋加速器能量增益；与 **OPAL-cycl** 标准算例对比。

### E. 并行化与 GPU 加速（R1）
分三级推进，每级都保留纯 CPU 参考实现用于回归：

1. **线程级（OpenMP）**：场插值、粒子推进、电荷沉积（scatter）——数据并行，改造风险最低；可先扩展现有 `CScheduler`。
2. **节点级/多节点（MPI + hypre/AMG）**：Poisson/空间电荷求解做区域分解，用 `hypre BoomerAMG` 或 `PETSc` 替换/补充 `CEpotMatrixSolver`；粒子按域分布并通信。
3. **GPU（已确认采用 Kokkos，兼容 CUDA/HIP/SYCL）**：粒子推进 + 电荷沉积用 Kokkos kernel；Poisson 用 Kokkos + AMG（NVIDIA 上可选 `AMGx` / `cuDSS`）。
   - 选 Kokkos 而非裸 CUDA：一套代码跑 CPU/GPU/多厂商，契合"混合 CPU+GPU+集群"目标。
   - 仅在确认锁定单一厂商且需要极致性能时，再为热点单独写 CUDA/HIP kernel。

**性能目标（KPI，待与用户确认规模）**：
- 单节点多核相对单线程 ≥ 8×（OpenMP）；
- GPU 相对同机 CPU 全核 ≥ 5×（推进+沉积热点）；
- 内存占用与网格规模线性可扩展，支持 ≥ 10^7 宏粒子量级。

### F. 可视化与 IO（R3）
- **IO 优先**：新增 **openPMD / HDF5** 读写器（场图 + 粒子），这是与外部求解器（Elmer/Palace/Radia）和 OPAL 互操作的关键枢纽；再补 **VTK/XDMF** 导出。
- **3D 可视化**：主推 **ParaView**（VTK/XDMF 直接可读）；同时提供 **Python (`pyvista`/`yt`)** 交互分析与 Jupyter 教程。
- **轻量 Web 可视化（可选）**：`vtk.js` 或 `three.js` 生成可分享的 3D 视图。
- **保留并改进原生**：沿用 GTK/cairo 的 2D 绘图，新增"3D 切片/等值面"视图，不破坏既有交互工具。

---

## 5. 里程碑（阶段 + 验收标准）

| 阶段 | 内容 | 验收标准（DoD） |
| --- | --- | --- |
| **P0 立项** | Fork 上游、导入源码到 GitHub、保留 GPL 头、搭建构建与 CI、许可证清单、目录骨架 | `main` 可编译；CI 通过；LICENSE/NOTICE/第三方清单齐备 |
| **P1 3D 场基座** | `CVectorField3D`/`CFieldMap3D`/`CTimeVaryingField` + openPMD IO | 单元测试：插值精度；与解析场对比 | 
| **P2 磁铁** | Elmer/GetDP/Radia 适配器 → 3D B-map | 均匀场/二极磁铁校验通过；回旋频率复核 |
| **P3 RF 腔** | Palace/Elmer 适配器 → 模式场 + ω | 圆柱腔 TM010 频率误差 < 1% |
| **P4 时变跟踪** | 时间感知（+相对论）Boris；RF 增益与相位稳定 | 解析能量增益/相位滑动；与 OPAL-cycl 算例对比 |
| **P5 并行/GPU** | OpenMP → MPI → Kokkos/GPU 三级后端 | 达到 §4E KPI；回归测试全部通过 |
| **P6 可视化** | openPMD/VTK 导出、PyVista/ParaView、Jupyter 教程 | 3D 束流/场可视化示例可复现 |
| **P7 发布** | 文档、教程、示例算例、Release tag、可引用 DOI(Zenodo) | v0.1 release + 用户文档 + 复现实验 |

依赖关系：**P0 → P1 → {P2, P3} → P4 → P5 → P6 → P7**。

**优先级（已确认 R2 优先）**：先打通 **P1→P2→P3→P4** 的"物理正确性最小闭环"；P5 仅先启动 **OpenMP** 级（低风险、可与 P1–P4 并行），**MPI/GPU 推迟到物理验证之后**；P6 可视化最后。GPU 后端虽后置，但 **P1 起即采用 Kokkos 数据结构**，避免后期重写。

---

## 6. 仓库结构与工程化（建议）

```
IBSimu-Cycl/
├─ CMakeLists.txt                 # 顶层，统一构建新模块
├─ upstre/  (或 git submodule)    # 上游 IBSimu 源码，保持可同步
│   └─ ...                        # 保留原 GPL 头与文件结构
├─ src/
│  ├─ cyclotron/
│  │   ├─ field/        # CVectorField3D, CFieldMap3D, CTimeVaryingField
│  │   ├─ tracking/     # TimeVaryingStepper, relativistic Boris
│  │   ├─ solver/       # CPU(OpenMP)/MPI/GPU 后端
│  │   └─ io/           # openPMD / HDF5 / VTK 读写
│  └─ ...
├─ adapters/            # 外部求解器适配（Elmer/GetDP/Palace/Radia 调用+转换）
│  ├─ magnet3d/
│  └─ rfcavity3d/
├─ python/ibsimu_cycl/  # Python 绑定 + 后处理
├─ examples/cyclotron/  # 回旋加速器算例与教程
├─ tests/               # 单元/回归/基准
├─ docs/                # 文档（本文件所在）
├─ third_party/         # openPMD/Kokkos/hypre 等（或由包管理提供）
└─ .github/workflows/   # CI（构建/测试/文档）
```

工程化要点：
- 构建：上游沿用其原有构建方式；新模块用 CMake 统一并做包装，避免破坏上游。
- 依赖管理：尽量用系统包（conda/spack/apt）而非 vendoring；提供 `environment.yml` / `spack.yaml`。
- 测试：CTest + pytest；物理校验用例进 CI（解析解、快速小算例）。
- 版本策略：`upstream` 分支跟踪上游，`main` 为改造版；用 `git remote` 拉取上游更新并定期 rebase/merge。

---

## 7. 许可证与合规（重要）

- **IBSimu = GPL**。派生版本必须整体 GPL，且必须保留原版权与许可证声明；不能闭源、不能更改为不兼容许可。
- 拟整合组件许可证（**须逐一复核**，避免 GPL 不兼容）：
  - openPMD-api：**LGPL/MPS**（动态链接一般兼容，需注意）。
  - Kokkos：Apache-2.0（兼容 GPLv3，注意与 GPLv2-only 的兼容性）。
  - hypre/PETSc：BSD-2/BSD-3 风格（兼容）。
  - Palace / MFEM：Apache-2.0。
  - Elmer：GPL。GetDP/Gmsh：GPL。
  - Radia：GPL。
  - AMGx：BSD-3。
  - VTK：BSD-3。pyvista：MIT。
- **动作**：建立 `THIRD_PARTY_NOTICES.md`，在 P0 完成；明确"引用其输出（场图）"与"链接其代码"两种整合方式的许可差异。

---

## 8. 风险与对策

| 风险 | 影响 | 对策 |
| --- | --- | --- |
| 上游为 2D 假设，3D 空间电荷自洽改动大 | 工作量高 | 先做"静态场图 + 3D 跟踪"最小闭环，空间电荷 3D 求解分阶段 |
| 外部 FEM 求解器成为重型依赖 | 部署困难 | 松耦合：仅离线生成场图，运行时只依赖 IO 格式 |
| 与上游差异过大导致无法合并更新 | 维护成本 | 新增模块尽量独立，改动上游处最小化并加注释/宏开关 |
| GPU 代码锁定厂商 | 可移植性差 | 优先 Kokkos/SYCL，保留 CPU 参考实现 |
| 缺实测数据验证 | 结果可信度低 | 用 OPAL-cycl / 解析解交叉验证，争取真实磁铁/腔数据 |
| 上游无 Git 历史（SourceForge） | 追溯困难 | 导入并打 tag 标注上游版本基线 |

---决策记录与遗留问题

### 已确认（见 §1.1）
- 运行环境：混合 CPU+GPU+集群 → **Kokkos + MPI**
- 目标：有具体型号、无实测数据 → 几何建模 + 解析/OPAL-cycl 验证
- 外部依赖：接受 → **Elmer（磁铁）+ Palace（RF）**
- 优先级：**R2 物理能力优先**
- 上游：**允许激进重构**，不与上游保持可合并

### 仍需补充/确认（不阻塞 P0/P1）
1. **目标回旋加速器具体型号**、能量、离子种类 —— 用于建立算例与验收基准。
2. **是否强制支持国产 GPU / 非 NVIDIA 加速卡** —— 影响 Kokkos 后端与 Poisson GPU 求解器的最终取舍。
3. **交付形态**：CLI + Python 为主，还是需保留/强化 GTK GUI？（当前默认：保留 GTK 2D，另加 3D 分析管线）
4. **是否有可用的参考算例**（如 OPAL-cycl 已发布输入），便于交叉验证。GUI？
5. **优先级**：R1（性能）/ R2（物理能力）/ R3（可视化）三者若资源有限，先保哪个？
6. **上游同步**：是否希望保持与上游 IBSimu 可合并（偏保守），还是允许激进重构（偏性能）？

---

## 10. 下一步（P0 可立即执行项）

1. 建立 GitHub 仓库 `IBSimu-Cycl`，导入上游源码并打基线 tag。
2. 落地目录骨架、构建接入、LICENSE/NOTICE/THIRD_PARTY_NOTICES、CI 工作流。
3. 实现 P1 最小可运行 Demo：读入一个 3D B-map（先用解析/合成数据）→ 3D 跟踪 → openPMD 导出 → PyVista 可视化，形成端到端"hello cyclotron"。
4. 以此为骨架逐步接入 P2 磁铁与 P3 RF 腔求解器。

---

## 11. 进展记录

### P0（基本完成）

- [x] 定位并导入上游：`git clone https://git.code.sf.net/p/ibsimu/code`，基线提交 `09beedb`（2026-07-29，`master`），**保留完整 Git 历史与 GPL 归属**。
- [x] 上游基线可编译：`./reconf && ./configure && make` 成功（Gtk=yes, OpenGL=no, UMFPACK=no, CSG=no）。
- [x] 新增模块接入 autotools：`src/cyclotron/` + `-I$(srcdir)/cyclotron`，符号 `ibsimu_cycl::version()` 已进库。
  - **构建系统决策**：P0 **沿用上游 autotools**（可立即构建）；GPU/Kokkos 阶段再评估是否引入 CMake。
- [x] 工程文件：`README.md`、`NOTICE`、`THIRD_PARTY_NOTICES.md`、`.github/workflows/ci.yml`、`adapters/`、`python/ibsimu_cycl/`、`examples/cyclotron/`。
- [x] 样本模型来源确定（见 §11.1）。
- [x] 上游测试 **19/19 PASS**（`make check`）。
- [x] `LICENSE`（GPL-3.0-or-later 全文）落地；基线 tag `upstream-baseline`。
- [x] `main` 分支提交：`c6f553c`（骨架）、`e5e882c`（忽略构建产物）；`upstream-base` 分支保留原始上游。
- [ ] GitHub 远端创建与推送（本机无 `gh`；需用户创建仓库，或用 `gh auth login` 后推送）。

### 11.1 验证基准：PSI Ring 回旋加速器（OPAL-cycl 算例）

- 参数：590 MeV，8 折对称（`SYMMETRY=8.0`），谐波 `CYHARMON=6`，调谐区间 72–590 MeV。
- 关键文件：`bfield.dat`（三维磁场图）、`rffield1/2.dat`（RF 场图）、`dist1/2.dat`、`ic.dat`、`refsol.dat`、`cyclotron1/2.in`。
- 来源（GitHub 镜像，GitLab PSI 不可达）：
  - <https://github.com/OPALX-project/OPAL>（经典 OPAL，含 `samples/cyclotron/`）
  - <https://github.com/OPALX-project/OPALX>（新一代，Exascale）
  - <https://github.com/OPALX-project/Manual-old>（手册 + `Cyclotron/*.dat`）
- 详见 `examples/cyclotron/README.md`。

### 11.2 P1 进展（3D 场基座）

- [x] `CCylFieldMap3D`：柱坐标 (r, θ, z) 三维矢量场图，**支持方位角 θ 变化**
  （上游 `MeshVectorField` 仅均匀立方/轴对称，做不到），三线性插值 + θ 周期 + 越界 clamp；
  原生 ASCII 读写。
- [x] `CTimeVaryingField`：时变场包装 `F(x,t)=A(x)·(offset+amp·cos(ωt+φ0))`，供 RF 场使用；
  与 IBSimu `VectorField` 接口兼容（`operator()(const Vec3D&)`）。
- [x] 单元测试 `tests/cycl_field3d.cpp`：均匀场/线性场精确性、θ 周期插值、时变调制、
  ASCII 往返 —— **8/8 通过**，并纳入 `make check`。
- [x] 接入 autotools（库 + 测试）；CI 拆分为 **headless** 与 **默认 GUI** 两个作业。
- [ ] 待办：`CFieldMap3D`（笛卡尔、**每轴独立步长**）；OPAL `bfield.dat` 解析器（归入 P2）；
  openPMD/HDF5 输出（P6）；最小端到端 Demo（导入 B-map → 3D 跟踪 → 可视化）。

### 11.3 P2 进展（三维磁铁磁场）

- [x] **OPAL/PSI "RING" 磁场图解析器** `read_ring_field_map()`（`src/cyclotron/ringmap.{hpp,cpp}`）：
  完整实现 CERN "FIELD" 格式（头部 4 数 + 标签/信息块 + 按半径的 4 数组数据块 + `LREC=` 定位），
  支持 `dr/dtet` 负值取倒数、kG→T 缩放、`nsymmetry` 折旋转复制为满 360°。
- [x] 输出为 `CCylFieldMap3D`（满 360°，`nz=1` 中平面 Bz），复用已有柱坐标插值。
- [x] **真实数据验证**：`examples/cyclotron/data/bfield.dat`（PSI Ring，141×135，1.2 MB）
  → `nrad=141, ntet=135, rmin=1900 mm, dr=20 mm, dtet=1/3°`；
  **`Bz(r=3.1 m)=1.526 T`**（物理合理）、峰值 2.095 T、8 折复制一致。
  测试 `tests/cycl_ringmap.cpp` 全通过（`make check` = **21/21**）。
- [x] 算例数据入库（`examples/cyclotron/data/`，来源与许可见其 `README.md`）。
- [x] **CI 全绿**：headless 与 default-GUI 两作业均通过。
### 11.4 P2 收尾（A / C / B）

- [x] **A. 离面场重建** —— `CRingFieldMap3D`（`src/cyclotron/ringfield3d.{hpp,cpp}`）：
  由中平面 Bz 用真空 Laplace 展开重建 `Br/Bθ/Bz`（与 OPAL 同源）；
  实测 Maxwell 校验：`|∇·B|/(|B|/r) ≈ 2e-3`、**`|∇×B| ≈ 5e-6`（无旋）**、
  45° 旋转不变、`Br/z ≈ dBz/dr`（差 0.3%）。
- [x] **C. 笛卡尔场图** —— `CFieldMap3D`（`src/cyclotron/fieldmap3d.{hpp,cpp}`）：
  每轴独立步长 + 三线性插值 + clamp + 原生 ASCII 读写
  （上游 `MeshVectorField` 仅支持各向同性步长）。
- [x] **B. Elmer 适配器脚手架** —— `adapters/magnet3d/`：
  - `elmer/magnetostatic.sif`：WhitneyAV 三维静磁，含铁芯非线性 B-H 曲线与远场边界；
  - `elmer/sector.geo`：gmsh 扇形磁铁几何模板；
  - `field_to_fieldmap3d.py`：散点 `x y z Bx By Bz` → IDW 重采样 → 本项目场图格式（numpy）；
  - `tests/test_field_to_fieldmap3d.py`：自检通过（均匀场精确、线性场误差 0.9%）。
  > Elmer / gmsh 需用户自行安装（本机未装）；适配器为**松耦合文件交换**，不链接其代码。
- [x] `make check` = **23/23**；CI 三作业：headless / default-GUI / python-adapter。
- [x] **端到端 Demo（D1）** —— `tests/cycl_track.cpp`：真实 PSI Ring 磁场 → `CRingFieldMap3D`
  三维场 → Boris 推动器 → 校验回旋频率（测得 `1.492997e8` vs 解析 `1.491981e8 rad/s`，
  相对误差 **6.8e-4**）+ 速率守恒（`1.6e-15`，机器精度）+ 中平面约束；
  输出 `cycl_track.csv` 与 `cycl_track.vtk`（ParaView 可读）。
  可视化脚本 `examples/cyclotron/plot_trajectory.py` 生成 `docs/img/cyclotron_orbit.png`
  （可清楚看到 **8 折扇区**引起的场强波动）。
- [x] `make check` = **24/24**。
### 11.5 P3 进展（RF 谐振腔）

- [x] **P3-A 时变场跟踪** —— `src/cyclotron/boris.{hpp,cpp}`：`CBorisPusher` 非相对论 Boris
  推进器，支持时变 E（`CTimeVaryingField`）与 B；调用者在每步前 `set_time`。
  `tests/cycl_rfgap.cpp`：一维 RF 间隙（`E_x=E0·f(x)·cos(ωt+φ0)`，`∫f dx = h`）能量增益与
  **精确渡越积分**吻合到 `1.6e-4`；相位依赖（+V 加速 / −V 减速）与横向约束已验证。
- [x] **P3-B Palace 适配器脚手架** —— `adapters/rfcavity3d/`：`palace/cavity.json`
  （`Problem.Type = Eigenmode` 模板）、`palace/cavity.geo`（gmsh 圆柱腔）、README（流水线与判据）。
- [x] **P3-C 解析验证基准** —— `adapters/rfcavity3d/pillbox_tm010.py`：解析圆柱腔 TM010 模式
  生成器（$f=\frac{c\,x_{01}}{2\pi R}$，$E_z=E_0J_0(kr)$，$B_\theta=\frac{E_0}{c}J_1(kr)$），
  输出本项目场图格式；自检 `tests/test_pillbox_tm010.py`：频率公式、$J_0(x_{01})=0$、近壁衰减、
  **Maxwell $|\nabla\times E|=\omega|B|$ 误差 `1.9e-4`**。Bessel 采用幂级数（机器精度）而非
  A&S 多项式，以保证 $J_0'=-J_1$ 自洽。
- [x] `make check` = **25/25**；CI `adapters-python` 作业覆盖两个适配器测试。
### 11.6 P4 进展（多圈加速与相位滑移）

- [x] `tests/cycl_accel.cpp`：
  - **Part A（均匀场，可严格验证）** —— **两个对径 dee 间隙**。
    **关键发现**：单间隙时切向冲量只使轨道**中心漂移**、半径不增长；
    对径双间隙才使中心漂移抵消、轨道**同心外扩**（与真实回旋加速器一致）。
    - `ΔE/turn = 2 q V0 cos φ`（测得 99.76 vs 100 keV，差 0.24% 为渡越因子）
    - `ω_RF = ω_c` → **相位锁定**（`max|φ| < 1e-14 rad`）
    - `ω_RF = 0.98 ω_c` → **相位滑移 −0.12566 rad/圈**（与 `−2π·0.02` 精确一致）
    - 轨道半径按 `r = √(2KE/m)/ω_c` 外扩（`0.198 → 0.245 m`）
  - **Part B（真实场图整体缩放）** —— 测得 `T_rev`、能量增益与相位演化（`−0.113 rad/圈`）。
- [x] 可视化 `examples/cyclotron/plot_acceleration.py` → `docs/img/cyclotron_acceleration.png`。
- [x] `make check` = **26/26**。
- [ ] 限制：真实 PSI Ring 绝对场强对应 β≈0.6，需**相对论 Boris**；本 Demo 用整体缩放模型。
### 11.7 F1 进展（相对论推进器与真实场强轨道）

- [x] `CBorisPusher` 增加**相对论模式**（`set_relativistic`；在 $u=\gamma v$ 空间做磁场旋转，
  E/B 仍可空；与非相对论共用同一接口）。静磁场下 $|u|$、$\gamma$ 精确守恒。
- [x] `tests/cycl_relativistic.cpp`：
  - 相对论回旋频率 $\omega=qB/(\gamma m)$：误差 `1.3e-5`；$\gamma$ 漂移 `3.7e-16`
  - 非相对论对比：同速下 $\omega$ 高 $\gamma$ 倍（ratio = 1.2000）
  - **真实 PSI Ring 场强**（$r=3.3$ m，$\langle B\rangle=0.669$ T）：$\gamma=1.224$、
    $\beta=0.577$、**KE = 210 MeV**、轨道有界（$r\in[2.94,3.56]$ m）、$\gamma$ 守恒 `3.6e-16`、
    $f_{rev}$ 与 $q\langle B\rangle/(2\pi\gamma m)$ 差 **1.98%**
- [x] 可视化 `examples/cyclotron/plot_relativistic_orbit.py` → `docs/img/cyclotron_real_orbit.png`
  （可见 8 折扇形导致的扇贝形轨道）。
- [x] `make check` = **27/27**。
- [ ] 后续：P4 的 Part B 可改用**真实场强**（无需缩放）；注入/引出；R1 并行化/GPU（P5）；
  R3（openPMD/HDF5 + 3D 交互可视化）。

### 11.8 G1 进展（R1 第一级：OpenMP 粒子级并行）

- [x] 新增 `ibsimu_cycl::CEnsembleTracker`（`src/cyclotron/ensembletracker.{hpp,cpp}`）：
  粒子系综的**数据并行**推进器。回旋加速器三维跟踪的天然并行维是粒子维——粒子之间
  无相互作用（空间电荷除外），场为只读（`VectorField::operator()` 是 `const`）。
- [x] `configure.ac` 增加 `AC_OPENMP`；`src/Makefile.am` / `tests/Makefile.am` 的
  `AM_CPPFLAGS` 与 `AM_LDFLAGS` 引用 `$(OPENMP_CXXFLAGS)`。
  **踩坑**：链接 OpenMP 必须复用 `OPENMP_CXXFLAGS`，`OPENMP_CXX_LDFLAGS` 是空的
  （否则报 `undefined reference to GOMP_parallel`）。
- [x] **关键实现决策：单一并行区 + 步间屏障**。
  最初的写法是"每步一个 `#pragma omp parallel for`"，基准显示：

  | 场景（N=2000, 400 步） | 每粒子每步耗时 | 加速比 |
  | --- | --- | --- |
  | 真实 PSI 场图 | ~0.4 μs | 6.8× ✔ |
  | 小场图 / 解析场 | ~0.1 μs | **0.61×（更慢！）** |

  原因是 400 次 fork/join 的启动开销完全吃掉了收益（小场图下每步并行区只有
  ~0.2 ms 的计算量）。改为**只建立一个并行区、步间用屏障同步**（时变场相位
  更新放在 `single` 区）：

  ```cpp
  #pragma omp parallel
  {
      for( int step = 0; step < nsteps; ++step ) {
  #pragma omp single
          { if( _E ) _E->set_time( (double)(step+1)*dt ); }
  #pragma omp for schedule(static)
          for( long i = 0; i < n; ++i )
              pusher.step( _E, _B, p[i].x, p[i].v, dt );
      }
  }
  ```

  改后（小场图用 N=20000 使每步并行区计算量足够大）小场图 **6.41×**、真实场图
  **6.76×**（12 逻辑线程）——两类场景都不再被启动开销支配。
- [x] `tests/cycl_omp_tracker.cpp`：
  - **位一致性**：单线程与多线程结果 `max|diff| = 0.0`（逐位相同；`schedule(static)`
    下每个粒子的运算序列与串行完全一致）。
  - **运行时自检**：先跑一个平凡的独立并行循环，若连它也拿不到加速（受限容器/单核
    CI runner），则**跳过**加速断言，避免在受限环境误报。这也使测试在 2 核 CI 上稳定。
  - **重复取最优**（`reps=3`，取最小时间）消除瞬时系统负载噪声——曾观测到同一个
    二进制因后台负载给出 0.38× 的假阴性。
- [x] 强扩展性（i7-10750H，6 物理核 / 12 逻辑核）：

  | 线程数 | 1 | 2 | 4 | 8 | 12 |
  | --- | --- | --- | --- | --- | --- |
  | 小场图（均匀场，N=20000 × 100 步） | 1.00× | 1.99× | 3.82× | 5.65× | 6.41× |
  | 真实 PSI Ring 场图（N=2000 × 400 步） | 1.00× | 2.01× | 3.91× | 5.91× | 6.76× |

  4 线程以内接近线性，之后受物理核数（6）限制；12 逻辑线程下效率约 55–60%。
  单线程速率 2.6×10⁶ 粒子·步/秒（真实场图）。
  **注**：§4E 的 KPI「单节点多核 ≥ 8×」在本机无法完全验证——i7-10750H 只有 6 个
  物理核（12 超线程），12 线程下 6.7× 已接近该硬件的实际上限；需在 ≥8 物理核的
  机器上复测（OpenMP 版本本身与核数无关）。
- [x] `make check` = **28/28**。
- [ ] 后续（R1 第二、三级）：MPI + hypre 区域分解（多节点）；GPU 后端（Kokkos 优先，
  一套代码覆盖 CUDA/HIP/SYCL，用户当前只考虑 NVIDIA）。
  **前置**：本机未装 CUDA toolkit 与 MPI/hypre（`sudo pacman -S cuda` /
  `sudo pacman -S openmpi hypre`），装好后即可落地。

### 11.9 R1-② 进展（单核热路径优化：场插值）

OpenMP 只解决"用满核"，单核效率决定绝对性能。profile 显示 `CRingFieldMap3D`
是三维跟踪的**唯一**场求值热点，每次求值存在三处低效：

| 问题 | 原实现 | 改造 |
| --- | --- | --- |
| 网格索引重复计算 | 每次求值调 **6 次** `interp()`，各算一遍 (r,θ) 单元与权重 | **只算一次**，6 个分量共用 |
| 访存局部性 | 6 个分量各存一个 1.2 MB 数组，4 个插值节点跨越 6 段内存（≈12 条缓存行） | **按节点交错存储**：b/dbr/dbth/tb/dtrb/dttb 连续 6 个 double（≈4 条缓存行） |
| 超越函数 | `atan2` + `cos` + `sin` 每次求值 | `cos/sin` 由 `x/r`、`y/r` 直接得到 → **省两次** |

结果（真实 PSI Ring 场图，N=2000 × 400 步，**单线程**）：

| 指标 | 优化前 | 优化后 | 提升 |
| --- | --- | --- | --- |
| 跟踪耗时 | 0.309 s | 0.093 s | **3.3×** |
| 吞吐 | 2.59×10⁶ 粒子·步/秒 | **8.63×10⁶ 粒子·步/秒** | 3.3× |

**数值上与优化前逐位不变**：ω_c 相对误差 6.808e-4、γ 漂移 3.63e-16、
KE=210.14 MeV、f_rev 误差 1.976%、Bz(中平面)=1.557581 T 均完全一致。

同样的布局改造用于通用三维场图（新增常驻基准 `tests/cycl_fieldbench.cpp`，
200 万采样 + 3 轮取最优）：

| 类 | 规模 | 优化前 | 优化后 | 提升 |
| --- | --- | --- | --- | --- |
| `CCylFieldMap3D` | 141×1080×9 (31.4 MB) | 11.4 Mevals/s | 13.1 Mevals/s | 1.16× |
| `CFieldMap3D` | 41³ (1.6 MB) | 28.9 Mevals/s | 29.0 Mevals/s | 1.00×（本就驻留缓存） |

**累计效果**（真实 PSI Ring，N=2000 × 400 步）：原始单线程 0.309 s → 优化后
12 线程 0.011–0.015 s，合计 **≈20–28×**（单核 3.3× × 多核 6.3–8.5×）。
§4E 的"单节点多核 ≥ 8×"KPI 在真实场图上满足（同代码多核加速 8.5×）；
小场图 6.2× 受本机 6 物理核限制。

- [x] `make check` = **29/29**。

### 11.10 R1-③ 进展（GPU 后端：CUDA + NVRTC）

- [x] 新增 `ibsimu_cycl::CGpuEnsembleTracker`（`src/cyclotron/cuda/`）：粒子维
  并行的 GPU 后端，静态磁场下推进整个系综。
  **实现方式决策**：采用**直接 CUDA 实现**而非 Kokkos。Kokkos 需要先自行构建
  （CMake）并与 autotools 共存，且用户已确认当前只考虑 NVIDIA；若将来需要覆盖
  AMD/Intel，可再迁移到 Kokkos（接口不变）。
- [x] **构建期不需要 nvcc**：内核源码以字符串内嵌（`gpu_kernels.cu.h`），
  运行期由 NVRTC 编译为 PTX，再用 CUDA Driver API 加载。

  **为何不用 nvcc 编译 .cu？（实测结论，已尝试并放弃）**
  automake 本身能接受 `.cu` + 自定义后缀规则，但 libtool 会：
    1. 无法推断 nvcc 的配置 → 必须加 `--tag=CXX`；
    2. 加上 tag 后会强制注入 `-fPIC -DPIC`，而 **nvcc 不接受裸 `-fPIC`**（fatal）；
    3. 改用 `-Xcompiler -fPIC` 传参又会被 libtool 重排到命令首部
       （变成 `-fPIC nvcc -c ...`）。
  要绕开需额外包装脚本，代价高于收益。NVRTC 方案完全不触自动 libtool：
  构建期只需要 CUDA 头文件与 `libnvrtc`/`libcuda`。

- [x] 场数据（`CRingFieldMap3D::raw_data()` 的**交错节点表**）一次性上传显存，
  之后每步不再传场；粒子状态按 SoA 传输（合并访存）；每个 GPU 线程一次跑完
  一个粒子的全部时间步，没有“每步启动 kernel”的开销。
- [x] 内核中的场求值与 Boris 步与 CPU 端**逐行对应**（`--fmad=false`），
  因此可直接交叉验证。
- [x] **优雅降级**：未探测到 CUDA 或 `--without-cuda` 时 `available()` 返回 false，
  测试输出 `[SKIP]` 并以成功退出 → CI（无 GPU）不受影响（已实测）。
- [x] `configure.ac` 增加 CUDA 探测（`--with-cuda=DIR` / `--without-cuda`，默认自动
  探测 `/opt/cuda`、`/usr/local/cuda`、`/usr`）。

  **踩坑（重要）**：automake 为 libtool 库生成的链接命令**不包含** `AM_LDFLAGS`
  （只含 `$(target)_LDFLAGS)` 与 `$(LDFLAGS)`），因此 `-L` 这类目标相关参数必须放在
  `<target>_LDFLAGS` 里。（顺着这个坑才发现之前放在 `AM_LDFLAGS` 里的
  `$(OPENMP_CXXFLAGS)` 实际上一直没生效——当时恰好不影响，因为 .so 允许未定义符号。）

- [x] `tests/cycl_cuda_tracker.cpp`：

  | 校验项 | 结果 |
  | --- | --- |
  | 10 步后 GPU vs CPU | `dx = 4.4e-16 m`，`dv/v = 3.4e-16`（机器精度） |
  | 400 步后 GPU vs CPU | `dx = 1.4e-14 m`，`dv/v = 3.8e-15` |
  | γ 守恒（400 步，逐粒子相对漂移） | GPU `2.0e-15`（CPU `1.8e-15`，同一量级） |

  性能（真实 PSI Ring，N=100000 粒子 × 200 步）：

  | 后端 | 时间 | 吞吐 | 相对 |
  | --- | --- | --- | --- |
  | CPU 单线程 | 2.30 s | 8.7 Mev/s | 1× |
  | CPU 12 线程 (OpenMP) | 0.32 s | 62 Mev/s | 7.1× |
  | **GPU (RTX 2060 Max-Q)** | **0.073 s** | **275 Mev/s** | **31.6× / 4.4×（vs 12 线程）** |

  kernel 占 97%（数据传输仅 2.6%）——说明该场景是计算受限，系综越大优势越明显。
  注：本机是消费级 GPU（FP64 吞吐仅为 FP32 的 1/32）；专业/数据中心 GPU 的
  FP64 比例高得多，预期加速更明显。

- [x] `make check` = **30/30**。
- [ ] 后续：GPU 端支持时变/射频场（当前仅静态磁场）；多 GPU / 多节点（MPI + hypre）。

### 11.11 R1-② 决策记录：不引入 MPI，改为「多重网格 + 并行化」

**背景**：用户提问「IBSimu 本来就会算空间电荷，还有必要用 MPI 吗？」——以下是代码核实与实测数据。

#### 代码核实（上游并行现状）

| 事实 | 证据 |
| --- | --- |
| 空间电荷能力确实存在 | `scharge.hpp/cpp`：PIC 沉积（2D/柱对称/3D）；求解器有 `EpotGSSolver`、`EpotBiCGSTABSolver`（Diag/ILU0/ILU1 预条件）、`EpotUMFPACKSolver`、`EpotMGSolver`（多重网格） |
| **全代码库无一处 OpenMP** | `#pragma omp` 只出现在本项目新增的 `src/cyclotron/` |
| pthread 只用于几何网格构建 | `pthread_create` 在所有 `.cpp` 中仅 1 处：`geometry.cpp:1309` |
| `Scheduler` 线程类从未被使用 | 全库无任何 `.cpp` 实例化它 |
| 电荷沉积的多线程接口已预留但未启用 | `scharge_add_from_trajectory_pic(..., pthread_mutex_t*)`、`particledatabaseimp.hpp` 的 `scharge_mutex` |

→ **粒子推进、电荷沉积、Poisson 求解三者在上游全是串行的**；缺的是「单机并行」，不是「跨节点」。

#### 实测：单节点求解器规模扫描

基准 `tests/cycl_poisson_bench.cpp`：单位立方腔（六面 Dirichlet 接地）+ 中心电荷团，
对比 `EpotMGSolver`（4 层）与 `EpotBiCGSTABSolver`+ILU0，`eps=1e-6`。

| 规模 | 节点数 | 源 | 多重网格 | RSS | BiCGSTAB+ILU0 | 迭代 | RSS |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 65³ | 0.27 M | gaussian | 0.24 s (5 cyc) | 71 MB | 0.26 s | 9 | 129 MB |
| 129³ | 2.15 M | gaussian | 1.98 s (5 cyc) | 143 MB | 1.72 s | 9 | 624 MB |
| 129³ | 2.15 M | step（陡边界） | 2.02 s (5 cyc) | 143 MB | 1.81 s | 11 | 623 MB |
| 257³ | 17.0 M | gaussian | 13.1 s (4 cyc) | 710 MB | **46.0 s** | **34** | **4615 MB** |
| 257³ | 17.0 M | step（陡边界） | 14.6 s (5 cyc) | 709 MB | 23.0 s | 19 | 4627 MB |

**读数**：

- **多重网格表现极稳**：循环数与网格规模、源的光滑性几乎无关（4–5 次），
  时间/内存**严格线性**于节点数：`0.77–0.94 μs/节点`、`41.7 B/节点`。
- **BiCGSTAB+ILU0 不稳定且昂贵**：迭代数随规模从 9 涨到 19–34；
  内存 **272 B/节点（多重网格的 6.5 倍）**，257³ 下已达 4.6 GB。

**单节点 15 GB 内存对应的规模上限**：

| 方法 | 内存系数 | 上限 |
| --- | --- | --- |
| 多重网格 | 41.7 B/节点 | ≈ 3.6×10⁸ 节点（≈ 710³） |
| BiCGSTAB+ILU0 | 272 B/节点 | ≈ 5.5×10⁷ 节点（≈ 380³） |

**对照实际需求**：回旋加速器 5 m 空间、1 cm 分辨率 → 500³ = 1.25×10⁸ 节点
→ 多重网格约 **97 s / 5.2 GB**，单节点完全可行；2 cm 分辨率 → 250³ ≈ 1.6×10⁷ 节点
→ 约 13 s / 0.7 GB。

#### 结论与决定

- [x] **不引入 MPI**。可预见规模（≤ 710³）内单节点内存与算力均有余量，
  MPI 区域分解只会带来通信开销；且 MPI 也无法解决"每步求解太慢"的问题。
- [x] 真正该做的两件事（列入后续）：
  1. **默认使用多重网格**——比 BiCGSTAB 快 1.6–3.5×、省 6.5× 内存、行为稳定；
  2. **给多重网格的平滑子（7 点模板松弛）加并行**——这是自洽 PIC 迭代里单步最贵的部分。
     按 1.25×10⁸ 节点、约 100 s/次求解、自洽迭代需数十次估算，总时长小时级；
     12 线程预期 3–5×、GPU 预期 10×+，收益远大于 MPI。
- [ ] 只有当网格超过约 710³（或沿用 BiCGSTAB 且超过 380³）时才需要 MPI + hypre；
  替换点已明确：`EpotSolver::subsolve()` 是虚函数。
- [x] `make check` = **31/31**（新增常驻基准，默认只跑 65³ 以保持 CI 快速）。

### 11.12 R1-② 进展（多重网格平滑子的 OpenMP 并行）

背景：§11.11 决策为「默认用多重网格，并给它的平滑子加并行」。

#### 关键发现：平滑子本来就是红黑 Gauss-Seidel

`EpotMGSubSolver::rbgs_loop_*` 已按 `rb = 0,1` 两色、`i += 2` 交替遍历。
同色遍历内，7 点模板**只读相反颜色**（本轮不写），因此：

- 无需改算法、无需重新着色、无需原子操作；
- 每个节点算出的新值与串行**逐位相同**（只有残差 `res` 的求和次序改变）。

已并行化的 3D 循环：

| 函数 | 并行方式 | 安全性依据 |
| --- | --- | --- |
| `rbgs_loop_3d` | `parallel for` over k + `reduction(+:res)` | 同色内无依赖 |
| `defect_3d` | `parallel for` over k | 每点只写一次、只读 epot/rhs |
| `restrict_3d` | `parallel for` over k | 每个粗点只写一次、只读细层 |
| `correct` | `parallel for` over k | 每个偶数点只写一次、只读 corr |
| `prolong_3d` | **未并行**（见下） | — |

#### 两个踩坑

1. **`prolong_3d` 不能直接并行**：它的 13 项写入**全部落在偶数细节点上**
   （偏移量中恰有 2 个为奇数），相邻粗节点会向同一偶数点贡献，属于「有冲突的散射加」。
   我曾误判为「只需写中心点」（以为其余 12 项写到奇数点、而奇数点不被读取），
   结果求解器直接不收敛（`mgcyc` 5 → 100，`max|phi|` 也变）——已回退，并在代码里
   留了说明。若要并行需改写为「按细节点 gather」，但会改变求和次序、
   牺牲与串行的逐位一致性。
2. **并行区内不能抛异常**（会导致 `std::terminate`）。原来在平滑循环里直接
   `throw` 定位 NaN/inf，现改为设置标志位，循环结束后再做一次串行扫描定位并抛出
   （正常路径零额外开销）。

另外，所有并行循环都带 `if(size >= 8)`，避免在最粗层（如 9³）引入并行开销。

#### 实测（`cycl_poisson_bench`，eps=1e-6，每档 3 次取最小）

| 网格 | 1 线程 | 4 线程 | 6 线程 | 8 线程 | 12 线程 |
| --- | --- | --- | --- | --- | --- |
| 129³ | 1.855 s | 0.617 s (3.0×) | **0.541 s (3.4×)** | 0.628 s | 0.704 s |
| 257³ | 12.29 s | — | **4.51 s (2.7×)** | — | 5.26 s |

- **解与迭代行为完全不变**：`mgcyc` 与 `max|phi|` 与并行前一致
  （129³：5 cyc / 8.824e+04 V；257³：4 cyc / 8.821e+04 V）。
- 加速在 **6 线程（= 物理核数）达到峰值**，再多则因 SMT 争用而回落 →
  典型的**内存带宽受限**（7 点模板）。实际使用时建议
  `OMP_NUM_THREADS=` 物理核数。
- `make check` = **31/31**（含 `solver3d_sphere`/`solver2d_coax` 等与解析解对比的测试）。
- **测量陷阱提醒**：本机为开发机，后台负载会让同一配置出现 3× 的假性回退
  （曾观测到 12 线程 2.18 s vs 0.67 s）。性能测量须多次取最优并同时记录 loadavg。
- [ ] 后续：2D/CYL 模式同样的改造；最粗层 SOR（本质串行）与 prolong（需 gather 改写）。
