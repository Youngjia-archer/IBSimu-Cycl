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
- [ ] 后续（P3）：Palace RF 腔本征模适配器；`CTimeVaryingField` 接入推进器实现**时变场跟踪**；
  柱坐标场图直转工具；更多可视化（R3）。
