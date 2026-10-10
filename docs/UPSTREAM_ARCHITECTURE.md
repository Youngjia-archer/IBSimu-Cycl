# IBSimu 上游架构总结（1.0.6dev）

> 目的：让任何人在不重新通读 249 个源文件的情况下，快速掌握 IBSimu 的结构、
> 数据流与关键约定。本文档由 IBSimu-Cycl 项目在集成过程中逐项核实编写，
> 与上游代码一一对应（基线提交 `09beedb`，2026-07-29）。
>
> 相关文档：`ROADMAP.md`（本项目的规划与进展）、`WORK_LOG.md`（工作日志）。

---

## 1. 概述

| 项 | 内容 |
| --- | --- |
| 定位 | **离子源 / 低能束光学模拟器**（Ion Beam Simulator），非通用加速器程序 |
| 语言 | C++（无 OpenMP；并行基于 pthread） |
| 构建 | GNU autotools（`configure.ac` + `Makefile.am` + `reconf`），无 CMake |
| 布局 | **扁平**：249 个源文件全在 `src/`；`doc/` 为 Doxygen；`tests/` 为 20+ 测试程序 |
| 许可 | GPL-2.0-or-later（派生项目可合法升级为 GPL-3.0-or-later） |
| 依赖 | 必需 pkg-config/zlib/libpng/fontconfig/freetype/cairo/gsl；可选 GTK3、gtkglext、UMFPACK |
| 入口 | 全局单例 `IBSimu ibsimu;`（`ibsimu.cpp`），用户代码以库的方式调用 |

---

## 2. 模块地图

### 2.1 基础设施

| 文件 | 职责 |
| --- | --- |
| `vec3d/vec4d/mat3d/mvector/transformation` | 向量、矩阵、坐标变换 |
| `matrix` + `ccolmatrix/crowmatrix/coordmatrix` | 稀疏矩阵三种存储（列主/行主/坐标），可互转 |
| `precond` + `diag_precond/ilu0_precond/ilu1_precond/empty_precond` | 预条件子 |
| `error` | 异常体系（`Error`/`ErrorNoMem`/`ErrorUnimplemented`…），`ERROR_LOCATION` 宏 |
| `ibsimu` | 全局配置单例：消息等级、随机数类型、**线程数**、输出流 |
| `timer/statusprint/comptime/compmath/constants/types/sort/cfifo` | 杂项 |
| `callback` | 回调函数对象（`CallbackFunctorB_V` 等），用于等离子体区域/固体判定 |

### 2.2 几何与网格

| 文件 | 职责 |
| --- | --- |
| `geometry` | 核心：`Geometry(MODE_?, Int3D size, Vec3D origo, double h)`，边界条件、固体标记、网格构建（**唯一使用 pthread 的上游代码**） |
| `solid/func_solid/csgobject_solid/dxf_solid/stl_solid/vtriangle` | 实体几何：函数式、CSG、DXF、STL |
| `mesh/meshscalarfield/meshvectorfield/interpolation/coordmapper` | 网格场容器与插值；`MeshScalarField::clear()` 即 `memset(0)` |
| `axisymmetricvectorfield/multimeshvectorfield` | 特殊矢量场包装 |

**关键约定**（写扩展代码时必须遵守）：

- 模式：`MODE_1D / MODE_2D / MODE_CYL / MODE_3D`。**2D/CYL 时 `size(2) == 1`**——
  这会让"按 k 并行"之类的优化静默失效（本项目踩过此坑）。
- `Geometry::mesh(i,j,k)` 返回位标记：`SMESH_NODE_ID_MASK` 取节点类型
  （`..._PURE_VACUUM / _NEAR_SOLID / _NEUMANN / Dirichlet`），`SMESH_NODE_FIXED` 表示固定节点，
  `SMESH_NEAR_SOLID_INDEX_MASK` 配合 `nearsolid_ptr()` 取"到固体的分数距离"表。
- CYL 模式下 **`j` 是径向**（`1/j` 项、`EPOT_SOLVER_BYMIN` 即轴 r=0）；
  轴上的模板是特例（只读 `(i,j+1)`）。

### 2.3 场与求解器

| 文件 | 职责 |
| --- | --- |
| `field/scalarfield/vectorfield` | 抽象基类；`VectorField::operator()(const Vec3D&)` 返回直角坐标分量 |
| `epot_field/epot_efield` | 电位场、电位→电场 |
| `solver/epot_solver` | 求解器基类：`solve(MeshScalarField &epot, const ScalarField &scharge)`、纯虚 `subsolve()`、`preprocess/postprocess`、等离子体模型 |
| `epot_gssolver` | 红黑 Gauss-Seidel / SOR（简单但慢） |
| `epot_matrixsolver` → `epot_bicgstabsolver` / `epot_umfpacksolver` | 矩阵法：BiCGSTAB（可配预条件）/ UMFPACK 直接解 |
| `epot_mgsolver` + `epot_mgsubsolver` | **几何多重网格**（本项目的主力）。`EpotMGSolver` 管层级与网格传输；`EpotMGSubSolver` 管单层平滑（红黑 GS）与缺陷计算 |

**等离子体模型**（离子源场景的核心）：`PLASMA_PEXP`（正离子）、`PLASMA_NSIMP`（负离子）、
`PLASMA_SHIELD`；通过 `set_pexp_plasma/set_nsimp_plasma/set_shield_plasma` 配置，
在节点上解非线性 Newton 方程（`solve_*_potential`）。

**实测选型结论**（本项目，见 `ROADMAP.md` §11.11）：多重网格远优于 BiCGSTAB+ILU0——
循环数与网格规模无关（4–5 次 vs 9→34 次），内存 41.7 B/节点 vs 272 B/节点（6.5 倍）。

### 2.4 粒子与自洽迭代

| 文件 | 职责 |
| --- | --- |
| `particles/particlestepper` | 粒子数据结构与**固定步长 Boris 步进**（`ParticleStepper`，仅 MODE_3D） |
| `particleiterator` | **GSL 自适应 ODE 积分**的轨迹追踪；`_scharge_mutex` 形参用于多线程沉积 |
| `particledatabase` / `particledatabaseimp` | 3D/2D/Cyl 数据库；**多线程追踪的实际实现**（`Scheduler` + 每线程一个 `ParticleIterator`） |
| `scheduler` | 线程池：`run(iterators)`、`get_next_problem()`、`wait_finish()`、错误收集 |
| `scharge` | 电荷沉积：PIC 与 linear 两种；2D/CYL/3D 三套重载；`scharge_finalize_*` 做边界修正 |
| `trajectory/trajectorydiagnostics/particlestatistics` | 轨迹存储与诊断（发射度等） |
| `convergence` | **自洽迭代的收敛判据**：跟踪 epot/scharge/emittance 的逐代变化 |

**自洽迭代的典型形态**（**循环写在用户代码里**，库不提供驱动）：

```cpp
EpotBiCGSTABSolver solver( geom );          // 或 EpotMGSolver（更快）
for( size_t i = 0; i < N; i++ ) {
    solver.solve( epot, scharge );                        // ① 每代求解空间电荷
    efield.recalculate();
    pdb.clear();
    pdb.add_cylindrical_beam_with_energy( ... );           // 放入粒子
    pdb.iterate_trajectories( scharge, efield, bfield );   // ② 每代追踪（含沉积）
    // convergence.evaluate_iteration(); → 判断是否退出
}
```

#### 2.4.1 两条粒子推进路径（步长选择）

`ParticleDataBase` **同时提供两条推进路径**。它们不是互斥的替代实现，而是面向不同用法的两个入口：

| 入口 | 步长 | 终止条件 | 实现 |
| --- | --- | --- | --- |
| `iterate_trajectories( scharge, efield, bfield )` | **GSL 自适应**（`epsabs`/`epsrel`） | 边界 / 表面 / `max_time` / `max_steps` | `ParticleIterator<PP>`（`particleiterator.hpp`），每线程一个实例 |
| `step_particles( scharge, efield, bfield, dt )` | **用户给定固定 `dt`**（Boris） | 无——调用一次推进一个 `dt` | `ParticleStepper<PP>`（`particlestepper.hpp`） |

要点：

- `step_particles()` 每次调用**内部先 `scharge.clear()`**，因此它只沉积「这一步」的电荷；
  调用者应对每个时间步调用一次。这与需要固定时间栅格（RF 同步）的加速器仿真天然契合。
- `step_particles()` 只实现 **`MODE_3D`**；2D/CYL 会抛 `ErrorUnimplemented`。
- 首次调用时用 `ParticleStepper::initialize()` 对速度做半步反踢（leapfrog 启动），
  因此第一次 `step_particles()` 之后粒子恰好前进一个 `dt`，且 `t` 严格按 `dt` 累加。
- 粒子状态数组布局为 `[t, x, vx, y, vy, z, vz]`（注意 x/y/z 与速度分量**交错**）。

**相对论语义**：`set_relativistic(true)` 在两个入口下语义一致，都得到正确的
相对论回旋频率 $\omega_c = qB/(\gamma m)$：

- 自适应路径：`ParticleP3D::get_derivatives()` 使用**质量矩阵**形式（`src/particles.cpp`）。
- 固定步长路径：`ParticleStepper` 在 $u=\gamma v$ 空间做 Boris 旋转，旋转向量
  $t = \dfrac{qB}{2\gamma m}\Delta t$ **必须含 $1/\gamma$**。漏掉它会得到 $\omega=qB/m$，
  即偏大 $\gamma$ 倍（$\gamma=1.22$ 时约 22%；两圈后相位偏差 158°，位置完全错）。

两条路径的一致性由 `tests/cycl_stepper_modes.cpp` 在均匀 Bz 场中对解析解验证。

### 2.5 IO 与可视化

| 文件 | 职责 |
| --- | --- |
| `readascii/file/hbio` | ASCII / 二进制 IO |
| `mydxf*`（15 个文件） | DXF 读写（CAD 几何导入） |
| `stlfile` | STL 读写 |
| `plotter/renderer/softwarerenderer/glrenderer` | 绘图后端抽象（cairo 软件渲染 / OpenGL） |
| `graph/xygraph/graph3d/eqpotgraph/fieldgraph/fielddiagplot/meshgraph/histogram/colormap/palette/legend/ruler/frame/label` | 二维/三维图元 |
| `geomplot/geomplotter/fielddiagplotter/particlediagplotter` | 高层绘图组合 |
| `gtk*`（14 个文件） | GTK3 交互界面（可选编译） |

---

## 3. 并行现状（本项目逐项核实）

> 初次核实曾遗漏头文件并得出"上游没有任何并行"的错误结论，已更正
> （见 `ROADMAP.md` §11.11 的更正说明）。以下为**已核实**的事实。

| 环节 | 状态 | 证据 |
| --- | --- | --- |
| 几何网格构建 | ✅ pthread 并行 | `geometry.cpp:1309` `build_mesh_parallel_entry` |
| **轨迹追踪** | ✅ **已支持多线程**（默认关闭） | `IBSimu::set_thread_count(N)`；`particledatabaseimp.hpp` 为 N 个 `ParticleIterator` 各建一个 pthread，用 `Scheduler` 线程池派发粒子 |
| 电荷沉积 | ⚠️ 单一全局互斥锁 | `scharge_add_from_trajectory_pic(..., pthread_mutex_t*)`：临界区仅 4 次 double 累加，每条轨迹片段都要 lock/unlock |
| 空间电荷求解 | ❌ 串行（本项目已并行化） | `EpotMGSubSolver`/`EpotMGSolver` 中的平滑与网格传输 |
| OpenMP | ❌ 上游无 | `#pragma omp` 仅出现在本项目新增的 `src/cyclotron/` |

**实测（本项目）**：自洽循环中 trace 占 **92%**、solve 占 8%（1 线程）；
开启多线程后 trace 6 核扩展 **5.14×**、端到端 **4.4×**。
沉积锁**不是**瓶颈（锁开销被 GSL 自适应积分的场求值摊薄）。

---

## 4. 本项目的扩展与改动（相对上游）

| 类别 | 内容 |
| --- | --- |
| **新增模块** `src/cyclotron/` | `CCylFieldMap3D`（柱坐标三维场图，支持方位角变化）、`CRingFieldMap3D`（由中平面 Bz 重建三维场）、`CFieldMap3D`（笛卡尔）、`CTimeVaryingField`（时变场）、`CBorisPusher`（相对论 Boris）、`find_closed_orbit()`（单圈映射不动点＝闭合轨道求解）、`CEnsembleTracker`（OpenMP 系综）、`CGpuEnsembleTracker`（CUDA） |
| **并行化改动** | `epot_mgsubsolver.cpp` / `epot_mgsolver.cpp`：红黑 GS 平滑、缺陷、限制、插值、修正的 OpenMP 化（3D/2D/CYL）；`ibsimu.cpp`：默认线程数改为硬件并发数（可用 `IBSIMU_THREADS` 覆盖） |
| **新增基准** | `tests/cycl_fieldbench.cpp`（场求值）、`tests/cycl_poisson_bench.cpp`（求解器规模扫描）、`tests/cycl_pic_bench.cpp`（自洽循环时间构成）、`tests/cycl_omp_tracker.cpp`（系综并行）、`tests/cycl_cuda_tracker.cpp`（GPU） |
| **适配器** `adapters/` | `magnet3d`（Elmer + 场图转换）、`rfcavity3d`（Palace + 解析圆柱腔） |

---

## 5. 阅读建议（最小路径）

想快速理解 IBSimu，按此顺序读即可：

1. `ibsimu.hpp` — 全局单例与配置
2. `geometry.hpp` + `meshscalarfield.hpp` — 网格与标记约定
3. `epot_solver.hpp` — 求解器接口与等离子体模型
4. `epot_mgsolver.cpp` 的 `mg_recurse()` — 多重网格循环（**本项目并行化的主战场**）
5. `scharge.hpp` + `particleiterator.hpp` — 电荷沉积与轨迹追踪
6. `particledatabaseimp.hpp` — 多线程追踪的实际装配
7. `convergence.hpp` — 自洽迭代的收敛判据
