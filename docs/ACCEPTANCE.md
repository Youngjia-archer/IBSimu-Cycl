# 验收矩阵：需求 → 判据 → 测试 → 证据

> 一键运行：`make check`（全量，含本测试）或 `./tests/cycl_acceptance`（仅验收矩阵，< 2 s）。
> 本文件中的实测数字来自 2026-10-10 开发机（i7-10750H 6P12T + RTX 2060 Max-Q），
> 由 `tests/cycl_acceptance` 与各专项测试的当前输出摘录。

---

## 1. 需求与验收判据（来源：`docs/ROADMAP.md` §1/§5）

| 编号 | 需求（用户原始目标） | 验收判据（DoD） |
| --- | --- | --- |
| P0 | 立项与工程化 | `main` 可编译；LICENSE/NOTICE/第三方清单/CI 齐备 |
| R1 | 提升三维仿真效率（并行化 / CUDA） | 多核并行有效且结果与串行一致；GPU 通路可用（可选后端） |
| R2 | 三维磁场、RF 谐振腔、时变场中的束流运动 | 真实场图解析与离面重建（真空 Maxwell）；相对论回旋；RF 渡越增益；真实场闭合轨道 |
| R3 | 改善可视化与 IO | VTK/HDF5 写出并可独立回读；示例与成图资产齐备 |

## 2. 顶层门禁：`tests/cycl_acceptance`

把上述判据压缩成 **42 项检查**（< 2 s）输出“需求 → 结论”矩阵；任何一项回归
都会让 `make check` / CI 失败。CI 的两个构建作业中另有同名显式步骤
（Acceptance gate），使门禁在 CI 界面上可见。

| 分区 | 覆盖内容 | 本次实测证据 |
| --- | --- | --- |
| §0 P0 | 11 个工程文件 + LICENSE 全文抽查；可选后端探测 | 12/12 通过（OpenMP 12 线程 / HDF5 2.2.0 / CUDA 已编译） |
| §1 R1 | OpenMP 可用；默认线程数 = 硬件并发；并行/串行逐位一致；加速比自检；GPU vs CPU | `max\|diff\| = 0`；均匀场 12 线程 **6.05×**；GPU `max\|dx\|=1.3e-15 m` |
| §2 R2 | 场图解析；真空 Maxwell；$\omega_c=qB/(\gamma m)$；RF 渡越；闭合轨道 | Bz(3.1 m)=1.526363 T；`\|div B\|`=8.7e-4（`\|B\|/r`=0.472）；$\omega_c$ 偏差 1.3e-5；RF 渡越 1.6e-4；闭合轨道 r=3.2890 m |
| §3 R3 | VTK `.vti`/`.vtp` 写读；HDF5 写读；10 个可视/示例资产 | 全部逐位/逐值一致；资产 10/10 |

### 2.1 运行时自检与降级路径（CI 安全设计）

- **CUDA / HDF5 为可选后端**：未编译或运行期不可用时输出 `[SKIP]` 并记通过，
  无 GPU 的 CI runner（`sudo apt` 不装 CUDA）不受影响；
- **加速比断言带环境自检**：仅当 `nthreads >= 4` 且独立平凡并行循环实测加速
  `>= 2×` 时才断言（2 核 CI runner 自动跳过数值断言，只验正确性）；
- **测量抗抖动**：大窗口（serial ≈ 0.2 s）+ 3 轮取最优 + 未达标自动重测一轮。
  （踩坑：0.01 s 级窗口在长套件/后台负载下曾给出 0.32× 的假性回退——
  屏障同步会放大瞬时负载的影响。）

## 3. 与专项测试的分工

`cycl_acceptance` 是**快而浅**的顶层门禁（回答“需求是否还立着”）；
深度验证由各专项测试承担（回答“数值对不对”）。两者都在 `make check`
（当前 **37/37**）内运行：

| 需求领域 | 专项测试（深度验证） | 关键证据 |
| --- | --- | --- |
| R1 并行 | `cycl_omp_tracker` / `cycl_fieldbench` / `cycl_poisson_bench` / `cycl_pic_bench` | 真实场 6.3–8.5×；场求值单核 3.3×；多重网格 2.4–3.7×；PIC 循环 4.4× |
| R1 GPU | `cycl_cuda_tracker` | 31.6× vs 单线程；与 CPU 一致到机器精度 |
| R2 磁场 | `cycl_ringmap` / `cycl_ringfield3d` / `cycl_fieldmap3d` / `cycl_track` / `cycl_closed_orbit` | 141×135 实图；Maxwell 残差 2e-3；$\omega_c$ 6.8e-4；闭合残差 1.2e-2 m |
| R2 RF/加速 | `cycl_rfgap` / `cycl_accel` / `cycl_relativistic` / `cycl_stepper_modes` | 渡越 1.6e-4；双隙 ΔE=2qV₀cosφ；相位滑移 −0.0058 rad/圈；γ 守恒 |
| R3 IO | `cycl_vtk_export` / `cycl_hdf5_export` + `python/tests/*` | 独立解析器交叉验证；HDF5 2.7× 更小 / 2× 更快 |

## 4. 运行方式

```bash
make check                                   # 全量套件（含验收测试）
./tests/cycl_acceptance                      # 仅验收矩阵（<2 s）
./tests/cycl_acceptance /path/to/bfield.dat  # 指定场图数据文件
```

两条运行目录约定均受支持：`make check`（CWD = `tests/`）与仓库根手工运行。
