# IBSimu-Cycl 工作日志

> **用途**：跨会话交接。记录已完成的工作、关键决策及其理由、踩过的坑与当前状态。
> 详细技术叙述见 `ROADMAP.md` 的 §11.x 进展记录；上游结构见 `UPSTREAM_ARCHITECTURE.md`。
>
> 最近更新：2026-10-10（提交 `8f0b39c`）

---

## 1. 项目定位

基于 IBSimu（GPL-2.0-or-later，基线 `09beedb`）的派生项目，针对**回旋加速器**场景做
针对性改造，整体以 **GPL-3.0-or-later** 发布：
<https://github.com/Youngjia-archer/IBSimu-Cycl>

三大目标（用户原始需求）：
1. 提升三维仿真效率（并行化 / CUDA）
2. 补足回旋加速器特有的物理：三维磁场、射频腔、时变场下的束流运动
3. 改善可视化与图形输出

---

## 2. 时间线（阶段 → 提交 → 结果）

| 阶段 | 提交 | 内容 | 验证 |
| --- | --- | --- | --- |
| **P0** 立项 | `c6f553c` `e5e882c` `a913f71` | 导入上游、autotools、GPL-3.0、CI、发布 GitHub | `make check` 19/19 |
| **P1** 三维场基座 | `353bcc9` | `CCylFieldMap3D`（支持方位角变化）+ `CTimeVaryingField` | 8 项检查 |
| **P2** 磁场 | `43e248b` | OPAL/PSI RING 场图解析 → 真实 PSI Ring（Bz(3.1m)=1.526 T） | 13 项检查 |
| **P2-A** 离面重建 | `07b22ef` | `CRingFieldMap3D`：由中平面 Bz 用 Laplace 展开重建 Br/Bθ | Maxwell \|∇·B\|≈2e-3 |
| **P2-C** 笛卡尔 | `fad83ab` | `CFieldMap3D`（每轴独立步长） | 4 项检查 |
| **P2-B** 适配器 | `564d122` | Elmer 静磁 + 场图转换器 | 自检通过 |
| **D1** 端到端 | `616e95f` | 真实 B 场中三维跟踪 | ω_c 误差 6.8e-4 |
| **P3-A** Boris | `66aa1ef` | `CBorisPusher` + RF 间隙渡越积分 | 吻合 1.6e-4 |
| **P3-B/C** 射频腔 | `f590485` | Palace 脚手架 + 解析 Pillbox TM010 | Maxwell 1.9e-4 |
| **P4/E1** 多圈加速 | `0074b5c` | 双对径 dee 间隙：ΔE=2qV₀cosφ、相位滑移 | 滑移 −0.12566 rad/圈 |
| **F1** 相对论 | `249cf5f` | 相对论 Boris + 真实场强轨道 | γ=1.224，KE=210 MeV，f_rev 1.98% |
| **G1** OpenMP | `79bf675` `ff41642` `938017c` | 粒子系综 OpenMP 并行 | 小场 6.4×、真实场 6.8× |
| **R1-②** 单核优化 | `a7c54c7` | 场插值交错存储 + 单次索引计算 | **3.3×** 单线程 |
| **R1-②** MPI 决策 | `6a0bf1b` | 求解器规模基准 → **决定不引入 MPI** | 见 §3.2 |
| **R1-②** MG 并行 | `77bb989` `0bb1a5a` | 多重网格平滑/缺陷/限制/修正 OpenMP 化（3D+2D+CYL） | 2.4–3.7× |
| **R1-③** GPU | `9cf2d5d` | CUDA 后端（NVRTC + Driver API） | **31.6×** vs 单线程 |
| **R1-③** PIC 修正 | `d57c719` `8f0b39c` | 更正上游并行现状的误判；PIC 基准 + **默认开启多线程** | 端到端 **4.4×** |

`make check` 当前 **32/32**；CI 三作业（headless / default-GUI / adapter tooling）全绿。

---

## 3. 关键工程决策与理由

### 3.1 许可证升到 GPL-3.0-or-later
利用上游"or later"条款，便于将来链接 Apache-2.0 组件（Kokkos / Palace）。

### 3.2 不引入 MPI（数据驱动决策）
实测 `cycl_poisson_bench`（多重网格 vs BiCGSTAB+ILU0）：

| 规模 | 多重网格 | RSS | BiCGSTAB+ILU0 | 迭代 | RSS |
| --- | --- | --- | --- | --- | --- |
| 65³ | 0.24 s (5 cyc) | 71 MB | 0.26 s | 9 | 129 MB |
| 257³ | 13.1 s (4 cyc) | 710 MB | 46.0 s | 34 | 4615 MB |

单节点 15 GB 内存上限：多重网格 ≈ 710³；而 5 m 回旋加速器按 1 cm 分辨率仅 500³
（≈97 s / 5.2 GB）→ **单机足够**。且 MPI 解决不了"单次求解 100 s × 自洽迭代数十次"
这一**单机算力**问题。结论：改用多重网格 + 并行化，MPI 留到 >710³ 再说
（替换点 `EpotSolver::subsolve()` 已是虚函数）。

### 3.3 GPU 后端用 NVRTC 而非 nvcc（构建系统约束）
automake 接受 `.cu` + 自定义规则，但 libtool 会：
1. 无法推断 nvcc 配置 → 必须 `--tag=CXX`；
2. 加 tag 后强制注入 `-fPIC -DPIC`，而 **nvcc 拒绝裸 `-fPIC`**（fatal）；
3. 改用 `-Xcompiler -fPIC` 又会被 libtool 重排到命令首部。

→ 改为**内核源码内嵌字符串 + 运行期 NVRTC 编译为 PTX + Driver API 加载**：
构建期只需 CUDA 头文件与 `libnvrtc`/`libcuda`，完全不触自动 libtool，
且无 GPU 时可优雅降级（`available()` 返回 false，测试 `[SKIP]`）。

### 3.4 并行粒度：单一并行区 + 步间屏障
粒子系综每步开一个 `parallel for` 会因 fork/join 开销反噬：
小场图下每步并行区仅 ~0.2 ms，实测仅 **0.61×（更慢）**。改为只建立一个并行区、
步间用屏障同步后 → **6.4×**。

### 3.5 场插值：按节点交错存储
6 个导数各存一个 1.2 MB 数组 → 一次插值跨约 12 条缓存行；交错后约 4 条。
配合"单元索引只算一次"，真实 PSI Ring 场图上单线程 **3.3×**（数值逐位不变）。

### 3.6 自洽 PIC 迭代：先测量，再动手（重要）
原以为"电荷沉积用一把全局锁保护 4 次累加"必然是瓶颈并准备重写。实测后：

| 线程 | solve | trace | 总 | trace 扩展 |
| --- | --- | --- | --- | --- |
| 1 | 0.855 s (8%) | 10.307 s (92%) | 11.163 s | 1.00× |
| 6 | 0.562 s | 2.006 s | 2.568 s | **5.14×** |

**两个结论推翻了先前假设**：
1. **trace 占 92%**，solve 只占 8% → 求解器并行化只影响约十分之一的运行时间；
2. **沉积锁不是瓶颈**——每条轨迹片段之间隔着 GSL 自适应积分的几十次场求值，
   锁开销被充分摊薄（6 核 5.14×，86% 效率）。

真正的浪费是：`ibsimu.cpp` 把 `_threadcount` 硬编码为 1，且**全项目无人调用
`set_thread_count()`** → 整个循环白白慢 4.4×。已改为默认硬件并发数
（`IBSIMU_THREADS=N` 可覆盖；=1 时严格复现）。

---

## 4. 踩坑清单

### 4.1 构建系统
- **automake 为 libtool 库生成的链接命令不含 `AM_LDFLAGS`**（只含 `$(target)_LDFLAGS`
  与 `$(LDFLAGS)`）→ `-L` 等目标相关参数必须放 `<target>_LDFLAGS`。
  （此前放在 `AM_LDFLAGS` 的 OpenMP 链接标志其实一直空转，只是恰好不影响。）
- 链接 OpenMP 必须复用 `$(OPENMP_CXXFLAGS)`；`OPENMP_CXX_LDFLAGS` 是空的。
- `tests/Makefile.am` 的 `check_PROGRAMS` 续行写成 `\ \` 会导致 configure 失败。

### 4.2 代码约定
- **2D/CYL 模式下 `Geometry::size(2) == 1`** → "按 k 并行"之类的优化会**静默失效**。
  曾因此在 `correct()` 里写了一个永不生效的 `parallel for`（3D 测得好掩盖了问题）。
- `MeshScalarField::clear()` 是 `memset(0)`；`correct()` 只读写**偶数点**
  （`(i+j+k)` 偶），奇数点由第一次 RBGS 覆盖。
- **`prolong_3d` 的 13 项写入全部落在偶数细节点**（偏移量恰有 2 个奇数分量）→
  是完整三线性插值、属于有冲突的散射加，**不能直接并行**。曾误判为"只写中心点"，
  导致求解器不收敛（`mgcyc` 5→100），已回退。要并行需改写为按细节点 gather，
  但会牺牲与串行的逐位一致性。
- **OpenMP 并行区内不能抛异常**（会 `std::terminate`）→ 改为标志位 + 循环后串行定位。
- **CYL 模式下 `j` 是径向**；轴上（`EPOT_BYMIN`）模板是特例（只读 `(i,j+1)`）。
- **`MeshVectorField( geom, fout, ... )` 里 `fout[i]=true` 表示该分量「存在」**，
  `false` 表示恒为零。写 `{false,false,false}` 得到的是**零场**（曾因此让粒子走直线，
  两条推进路径一起“错”，差点误判为步进器 bug）。
- **`Vec3D::norm2()` 返回的是 2-范数本身，不是范数的平方**（平方和另有 `ssqr()`）。
  把 `norm2()/SPEED_C2` 当 $\beta^2$ 用会让 $\gamma$ 恒等于 1，
  从而**掩盖真实误差**（相对论算例里表现为“误差存在但 γ 正常”的矛盾现象）。
- **相对论 Boris 的旋转向量必须含 $1/\gamma$**：$t=qB\Delta t/(2\gamma m)$，
  不是教科书记忆里的 $qB\Delta t/(2m)$。

### 4.3 检索与工具
- **检索并行/线程相关代码必须包含 `*.hpp`**：IBSimu 把大量实现放在头文件里，
  只 grep `.cpp` 会得出错误结论（我曾据此误判"上游无任何并行"）。
- **fish shell**：不能用 `VAR=value`、`$$`、bash 的 `for..do/done`；
  `$"` 在双引号内也会被当变量展开（正则要用单引号）。
- `replace_string_in_file` 的 `oldString` 必须**逐字节精确**（含制表符与行尾空格）；
  建议先用 `sed -n 'X,Yp' | cat -A` 核对空白字符。大批量替换后务必读回文件确认。
- CI 的 YAML：`name:` 中未加引号的 `&`（锚点）与 `: `（映射）会让整个 workflow 失效。

### 4.4 测量方法学
- **开发机后台负载会造成数倍假性回退**（曾观测 12 线程 2.18 s vs 0.67 s）。性能测量
  必须**多次取最优**并同时记录 loadavg。
- 基准必须足够长（如 200 万次求值）才能分辨 1.2× 级别的差异。
- 用 `git stash push -- <指定文件>` 做单变量前后对比非常有效。
- **先有基准再优化**：`mgcyc`/`max|phi|` 这类现成指标当场就能暴露语义破坏。
- **端到端“导出真实轨道”这种带场景的验证，比单元测试更容易暴露物理缺陷**：
  R3 的可视化导出测试本来是验证 IO，却抳出了 `CBorisPusher` 的 $O(\omega\Delta t)$
  启动瞬态（半径振荡 4.88e-4）。**测得的振荡幅度与理论值符合到 4 位有效数字**
  后才敢确认成因，而不是拍脑袋归因于“数值噪声”。

---

## 5. 当前状态

**代码结构**

```
src/cyclotron/          本项目新增（命名空间 ibsimu_cycl）
  cylfmap3d.*           柱坐标三维场图（支持方位角变化）
  ringmap.*             OPAL/PSI RING 场图格式解析
  ringfield3d.*         由中平面 Bz 重建三维场（交错存储）
  fieldmap3d.*          笛卡尔三维场图
  timevaryingfield.*    时变场包装
  boris.*               相对论 Boris 推进器
  ensembletracker.*     OpenMP 粒子系综
  cuda/                 CGpuEnsembleTracker（NVRTC + Driver API）
src/io/                 本项目新增（命名空间 ibsimu_cycl）
  vtkwriter.*           VTK XML 导出（.vti 场图 / .vtp 轨迹，零外部依赖）
python/ibsimu_cycl/     Python 侧：vtk_io 读取器（交叉验证 + 后处理）
adapters/               外部求解器适配器（Elmer / Palace），松耦合文件交换
examples/cyclotron/     真实算例数据、绘图与三维查看脚本
tests/cycl_*.cpp        本项目新增的测试与基准（共 15 个）
docs/                   ROADMAP / DEVELOPMENT / UPSTREAM_ARCHITECTURE / WORK_LOG
```

**性能成果汇总**（真实 PSI Ring 场图）

| 环节 | 加速 | 说明 |
| --- | --- | --- |
| 粒子系综 · 多核 | 6.3–8.5× | 单一并行区 + 步间屏障 |
| 粒子系综 · GPU | 31.6× vs 单线程 | NVRTC；与 CPU 一致到机器精度 |
| 场求值 · 单核 | 3.3× | 交错存储 + 单次索引 |
| 空间电荷求解 · 多核 | 2.4–3.7× | 红黑 GS 同色并行（逐位一致） |
| 自洽 PIC 循环 · 多核 | 4.4× | 默认开启多线程 |

**明确的未完成项**

- 2D/CYL 的 `prolong_*` 与最粗层 SOR 仍未并行（前者需 gather 改写，后者本质串行且网格很小）。
- GPU 后端只支持**静态磁场**（时变/射频场仍在 CPU）。
- 电荷沉积的累加顺序随线程数变化 → 结果有 ~1e-15 差异，严格复现需 `IBSIMU_THREADS=1`。
- 空间电荷的自洽迭代尚未与 GPU/Boris 路径打通（属**物理模型选择**问题）。

---

## 6. 步长选择：两条推进路径的实现与验证

**背景**：回旋加速器需要「与 RF 周期同步的等时间步」，而 IBSimu 默认的
`iterate_trajectories()` 是 GSL 自适应步长。目标是把步长选择做成可选参数，
**且两种方案都真正可用**（而不是只留个接口）。

### 6.1 查证：接口已经存在

上游**已经**在 `ParticleDataBase` 上提供两条路径，
所以“步长选择”在 API 层面已是可选参数。缺的不是接口，而是 **(a) 验证、(b) 相对论支持、(c) 文档**：

| 入口 | 步长 | 用法 |
| --- | --- | --- |
| `iterate_trajectories(scharge, efield, bfield)` | GSL 自适应 | 一次把轨道追到底 |
| `step_particles(scharge, efield, bfield, dt)` | Boris 固定 `dt` | 自己控制时间栅格，每次一个 `dt` |

**验证状态**：改动前仓库里 **18 处调用全部是 `iterate_trajectories`**，
`step_particles` **零测试、零文档**。

### 6.2 新增验证（`tests/cycl_stepper_modes.cpp`）

以均匀 Bz 场中的单粒子回旋运动（有解析解）为基准，两条路径 × 两种束流条件各追踪 2 圈：

| 路径 | 束流 | 轨道误差 | 速率漂移 | $\gamma$ 守恒 |
| --- | --- | --- | --- | --- |
| 自适应 | 非相对论 | 4.4e-5 | 4.4e-5 | 1.000000 |
| 固定步长 | 非相对论 | 4.1e-5 | 4.9e-6 | 1.000000 |
| 自适应 | 相对论 $\gamma$=1.22 | 6.4e-5 | 2.0e-5 | 1.220012 |
| 固定步长 | 相对论 $\gamma$=1.22 | 6.2e-5 | 3.3e-6 | 1.220002 |

（固定步长取 `dt = T/1000`，每圈 1000 步。）

### 6.3 修复的实质缺陷

`ParticleStepper`（固定步长）原先**完全没有相对论分支**：
`set_relativistic(true)` 对它无效，回旋频率会偏大 $\gamma$ 倍。已补齐：

- `src/particlestepper.hpp`：新增 `_relativistic` 成员、`set_relativistic()`、
  `v_to_u()`/`u_to_v()` 以及 $u=\gamma v$ 空间的 Boris 旋转；
- `src/particledatabaseimp.hpp`：`step_particles()` 构造 `ParticleStepper` 时传入 `_relativistic`；
- **非相对论分支保持与原实现逐位一致**（显式分叉，不共用表达式）。

### 6.4 一个“想当然”错误的代价

我最初按记忆把相对论旋转向量写成 $t = qB\Delta t/(2m)$（无 $\gamma$），
实测轨道误差 **1.61**（位置完全错位）。正确形式是 $t=qB\Delta t/(2\gamma m)$。

关键的修复线索是：本项目 `CBorisPusher`（已独立验证过 $\omega_c=qB/(\gamma m)$）
的实现里明确写了 `ghalf` 除法。**两个独立实现互相印证**比重新推导更快也更可靠——
以后遇到此类问题，优先去查仓库里已有的、被验证过的同类实现。

### 6.5 结论

- 两条路径均可用，误差量级相同（~5e-5，受各自设置支配）。
- **选择建议**：需要固定时间栅格 / 与 RF 同步 → `step_particles(dt)`；
  需要自适应精度、非 3D 网格、或一次性追到底 → `iterate_trajectories()`。
- 已知限制：`step_particles` 仅支持 **MODE_3D**（2D/CYL 抛 `ErrorUnimplemented`）；
  且它目前只接受静态 `bfield`，与 RF 时变场的对接见 ROADMAP 时变场条目。
- 另一处不对称：`ParticleStepper` **不更新 `ParticleStatistics`**
  （边界碰撞/电流统计只由 `iterate_trajectories` 填充），
  依赖 `get_statistics()` 的诊断代码在固定步长路径下拿不到数据。

---

## 7. 三维可视化：零依赖 VTK 导出（R3-①/②）

**背景**：用户原始三大需求（效率 / 回旋加速器物理 / 可视化）中，可视化是**唯一尚未推进**
的一项。R3 的规划路线是先 IO、后 3D 交互。

### 7.1 关键取舍：先不做 HDF5

查证发现 `configure.ac` 与 `src/Makefile.am` 里**没有任何 HDF5 接线**
（之前文档里“HDF5 已就绪”指的是系统装了库，不是项目链接了它）。
而 **VTK XML 的 `.vti`/`.vtp` 是纯文本格式**，不需链接任何库，
ParaView 与 PyVista 都原生支持 → 可视化可以在**零新依赖**下先落地，
HDF5/openPMD 留给大规模数据的二进制后端（R3-③）。

### 7.2 交付物

| 文件 | 作用 |
| --- | --- |
| `src/io/vtkwriter.{hpp,cpp}` | `.vti`（规则网格标量+矢量）、`.vtp`（轨迹折线，每点带 `t`） |
| `tests/cycl_vtk_export.cpp` | 导出 → 极简解析器逐值回读；真实螺旋轨道端到端 |
| `python/ibsimu_cycl/vtk_io.py` | **独立实现**的读取器（仅 numpy） |
| `python/tests/test_vtk_io.py` | 跨语言交叉验证（CI 两个构建作业均运行） |
| `examples/cyclotron/view_3d.py` | PyVista 查看/出图；无 pyvista 时退化为文本摘要 |

**为何要两套实现**：C++ 写出器与 Python 读取器互相独立，用后者校验前者才叫验证；
Python 侧还刻意**不复用 C++ 的常数**，而用文件内部的自洽关系（`Bmag == |B|`、
`z/t` 恒定、半径恒定、场的对称性）做判据——任一侧出错都会被发现。

### 7.3 意外收获：端到端验证掉出一个物理缺陷

螺旋轨道（`CFieldMap3D` 均匀场 + `CBorisPusher`，2 圈 800 步）的回旋半径并不恒定：

| | 半径振荡振幅 |
| --- | --- |
| 未初始化 | 4.8839e-4 |
| 理论 $(\omega\Delta t/2)\cdot r$ | 4.8840e-4 |
| **已初始化** | **3.8359e-6（抑制 127×）** |

标准 leapfrog/Boris 的速度定义在半步 $v^{n-1/2}$ 上，用 $t=0$ 的速度直接起步会留下
$O(\omega\Delta t)$ 初始瞬态（磁场**相位**推进仍正确，只有半径受影响）。
上游 `ParticleStepper` 有 `initialize()` 做半步反踢，`CBorisPusher` 没有。

→ 新增 `CBorisPusher::initialize()`（非相对论/相对论两分支）。
测试同时断言“未初始化振幅 = 理论值”，**把瞬态本身变成可验证的物理量**，
而不是放宽容差掩盖它。

### 7.4 结论与后续

- [x] `make check` = **34/34**（新增 `cycl_vtk_export`）。
- [x] CI 两个构建作业增加独立 Python 交叉验证步骤。
- [x] **接入真实算例**：`tests/cycl_track.cpp`（真实 PSI Ring 磁场图）现额外导出
  `cycl_track.vtp`（2001 点，带 `t` 标量可按时间着色）、
  `cycl_track_field.vti`（41×41×5 轨道附近三维 B）与
  `cycl_track_map.vti`（95×95 **整机中场图**，中平面 360°，r∈[1.9,4.7] m）。
  原 legacy `.vtk` 保留（但**不含时间标量**，是本次补充的动机之一）。
- [x] **成图**：新增 `examples/cyclotron/plot_field_map.py` →
  `docs/img/cyclotron_field_map.png`（整机中场 |B| + 五圈轨道、轨道区域放大、
  |B| 随方位角（45° 周期的 8 折结构）、|B| 随半径（扇区中心 vs 边界））。
  图中标签统一用英文——与既有 `plot_*.py` 一致，避免 matplotlib 中文字体依赖。

  **端到端自检**：Python 从 `.vti` 读回的 $r=3.3$ m、$\theta=90^\circ$ 处
  $|B| = 1.557581$ T，与 `cycl_track` 的解析判据**逐位一致**——
  即“C++ 解析真实场图 → 导出 VTK → Python 读回”全链条无失真。
- [ ] R3-③：HDF5/openPMD 二进制后端（大网格/大粒子数）；`.pvti` + 二进制追加段。

---

## 8. 真实场轨道为什么不闭合（诊断与修正）

**触发**：用户指出 `docs/img/cyclotron_real_orbit.png` 的轨道“重合且不稳定”，
并追问为何 `cyclotron_orbit.png` 只算了一个 0.5 m 范围的局部回旋。**两个问题都成立。**

### 8.1 局部 0.5 m 轨道：故意的，但要标清楚

那是 0.1c 质子在 $B=1.5576$ T 下的**拉莫尔回旋**：
$r_L=mv/(qB)=0.201$ m → 直径 0.4 m。选这么小是**故意的**——局部磁场近似均匀，
才能干净地测出 $\omega_c=qB/m$（这就是该算例的目的）。但它**不是**回旋加速器轨道，
原 README 的措辞容易让人误读，已在文档中改正。

### 8.2 真实场轨道：确实不闭合（已定量确认）

| 观测量 | 实测 | 判读 |
| --- | --- | --- |
| 径向摆动 | $r\in[2.94,3.56]$ m（±9.5%） | 比物理扇贝调制大一个量级 |
| 逐圈重复性 | 相邻圈同方位角半径差 **0.24 m（RMS）** | 每次回旋走全新曲线 → “轨道重合” |
| 45° 折叠后 $r$ 的范围 | 单个扇区内跨度即达 0.5 m | 不具 8 折周期性 ⇒ 不是闭合轨道 |

原测试的判据太松（“均值在 15% 内”、“$f_{rev}$ 在 5% 内”）掩盖了这个问题。

### 8.3 根因：初值不是闭合轨道

朴素初值 = 在 $r=3.3$ m 处**纯切向**发射、动量取 $p=q\langle B_z\rangle r$。
这不满足闭合条件。数值求出的不动点在**同一方位角**上带
$v_r=-1.60\times10^{7}$ m/s（偏离切向约 5°）。

物理解释：PSI Ring 是**强分离扇区** ring cyclotron——磁极只占每个 45° 扇区的
约 19°（极场 1.557 T × 0.43 = 平均场 0.669 T ✓ 自洽，与场图实测的方位分布一致），
其余是近零场的漂移段。在这种场里不动点的径向速度分量本就非零，
按切向发射会激励出大幅 betatron 振荡。

### 8.4 修正与验证（`tests/cycl_closed_orbit.cpp`，新增）

| 指标 | 朴素发射 | 闭合轨道 |
| --- | --- | --- |
| 逐圈重复性 | 0.24 m | **1.17e-2 m（改善 21×）** |
| 径向摆动 | 0.62 m | 0.13 m（±2.2%） |
| 整圈闭合性 $\Delta r$ | — | 3.4e-3 m |
| $\nu_r$ | — | 1.263（非整数 ⇒ 稳定） |

- [x] `tests/cycl_closed_orbit.cpp` 已纳入 `make check`（R2-D）；
- [x] `make check` = **35/35**。

### 8.5 后续修正（用户复查）

用户复查后指出三件事，均已处理：

1. **删除“局部拉莫尔回旋”**：0.1c 质子在该场中回旋半径仅 0.2 m，那是局部回旋、
   不是回旋加速器轨道。已删除 `tests/cycl_track.cpp` 的运行 B、`plot_trajectory.py`
   与图 `cyclotron_orbit.png`；`cycl_track` 现在只做 ω_c 校验 + 导出整机中场图。
2. **闭合轨道图不再做“未匹配”对比**：`plot_closed_orbit.py` 只画闭合轨道，
   并明确标注**这是单能量稳态轨道、不含 RF 加速**——真实回旋加速器里束流每圈
   获得能量、半径按 $r=\sqrt{2KE/m}/\omega_c$ 外扩，是一条向外盘旋的螺旋。
3. **加速度图的真实场数据已修正**（原 Part B「已知限制」，现已解决）：

   根因：初值是「在 $r=r_{ref}$ 处切向发射」，而在强分离扇区机器里这**不是**
   闭合轨道，轨迹被大幅 betatron 摆动主导（与 §8.3 同源）。修法：先在**同一
   能量**下求出闭合轨道（单圈映射的不动点）再起步，为此新增库模块
   `src/cyclotron/closedorbit.{hpp,cpp}`（`find_closed_orbit()`，求法同 §8.6）。

   | 量 | 修正前（切向发射） | 修正后（闭合轨道起步） |
   | --- | --- | --- |
   | 半径演化 | 3.50 → **3.46 m（反向）** | 3.295 → **3.541 m（外扩 ✓）** |
   | 能量增益 | 3.72 keV/圈（9%） | **24.22 keV/圈（61%）** |
   | 相位滑移 | — | **−0.246 rad/圈** |

   修正后半径随能量**增大**，与 $r\propto\sqrt{KE}$ 方向一致 ✓。`cycl_accel.cpp`
   的 Part B 现在有两条**真正的验收判据**：半径必须增大、增益必须高于设计值的一半。

   剩余的 −0.246 rad/圈 相位滑移是**真实的物理**：整体缩放后的场在当前能量下
   $\langle B_z\rangle(r)$ 上升得不够快，**不满足等时性**
   $\langle B_z\rangle(r)\propto\gamma(r)$，因此第 10 圈后相位滑过加速相位、
   能量增益趋于饱和（图中平台段）。这正是下一个待办项。

- [x] `make check` = **35/35**（删图不影响测试）。
- [x] 真实场中的 RF 加速（闭合轨道起步）—— 见上表。
- [ ] 待办：等时性场（$\langle B_z\rangle(r)\propto\gamma$ 使 $\omega_{RF}$ 保持常数，
  消除 −0.246 rad/圈 滑移）；R3-④ HDF5/openPMD 大规模数据后端。

### 8.6 踩坑记录（数值方法，重要）

- **数值 Jacobian 的步长必须匹配响应量纲**：$\partial r'/\partial v_r\sim1/\omega_c\sim2\times10^{-8}$ s，
  取 $\Delta v_r=0.1$ m/s 时响应被浮点噪声淹没，Jacobian 在 $v_r$ 方向退化为**奇异**，
  Newton 彻底失效。必须取 $\Delta v_r\sim10^{3}$ m/s。
- **线搜索的范数要与 Newton 方向一致**：用加权 1-范数（$|F_r|+10^{-3}|F_v|$）时
  Newton 方向已不是下降方向；改成按自然尺度归一的 2-范数（$V_0=\omega_c R_0$）后才正常。
- **速度是半步量**：Boris 的 $v_{n\pm1/2}$ 若按整步线性插值，会引入
  $\sim a\,dt/2=1.4\times10^{5}$ m/s 的系统偏差——那正是 Newton 残差地板的来源。
  按 $(f+0.5)$ 插值后地板下降。
- **利用对称性缩短积分**：求整圈不动点改为求 45° 单扇区不动点，噪声与条件数都显著变好。
- **先扫描再 Newton**：残差 0.2 m 时已超出线性区，先粗扫描定位再迭代。
- 最终地板（$\sim10^{-3}$ m）由**场图的双线性插值间断**决定，不是算法问题。

---

## 9. 建议的下一步（按性价比）

1. **R3 可视化（已启动，见 §7）**：`.vti`/`.vtp` 导出 + ParaView/PyVista 通路已完成。
   下一步是 R3-③：openPMD/HDF5 二进制后端（支持 >1e8 点）、把导出接进真实算例脚本
   （PSI Ring 场图 + 多圈轨道）、以及 Jupyter 教程。
2. **物理侧扩展**：注入/引出（螺旋偏转板）、Palace 真实腔模式场接入、
   相对论 RF 渡越。
3. **GPU 端支持时变/射频场**：让多圈加速也能跑在 GPU 上。
4. 若将来出现超大网格（>710³）：再考虑 MPI + hypre 区域分解。

---

## 10. 常用命令

```bash
./reconf && ./configure && make -j12     # 构建（含 CUDA 自动探测）
make check                               # 全部测试（当前 35/35）
./tests/cycl_poisson_bench 129 257       # 求解器规模扫描
./tests/cycl_poisson_bench --2d 1025     # 2D 模式
./tests/cycl_poisson_bench --cyl 513 257 # CYL 模式
./tests/cycl_pic_bench 4 1 6 12          # 自洽循环时间构成
IBSIMU_THREADS=1 make check              # 强制串行（严格复现）
./configure --without-cuda               # 关闭 GPU 后端
python3 python/tests/test_vtk_io.py      # VTK 导出器的 Python 交叉验证
python3 examples/cyclotron/view_3d.py a.vti b.vtp   # 三维查看（需 pyvista）
```
