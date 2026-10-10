# 回旋加速器验证算例

本目录用于存放**回旋加速器**算例与验证基准。P0 阶段先完成算例来源登记，
P2/P3/P4 阶段再导入实际输入文件与场图。

## 基准模型：PSI Ring 回旋加速器（来自 OPAL-cycl）

我们采用 **OPAL-cycl** 的标准回旋加速器算例作为物理验证基准（PSI Ring，590 MeV，
8 扇形区，谐波数 h=6，调谐范围 72–590 MeV）。

### 算例文件清单

| 文件 | 说明 | 对应我们的模块 |
| --- | --- | --- |
| `cyclotron1.in` | 调谐（tune）计算输入 | P4 跟踪 |
| `cyclotron2.in` | 加速轨道计算输入 | P4 跟踪 |
| **`bfield.dat`** | **三维磁场图** | **P2 磁铁** |
| `rffield1.dat` / `rffield2.dat` | RF 场图 | P3 RF 腔 |
| `dist1.dat` / `dist2.dat` | 初始粒子分布 | P4 |
| `ic.dat` / `refsol.dat` | tune 计算初值与参考解 | 验证 |
| `cyclotron1.gpl` / `cyclotron2.gpl` | gnuplot 绘图脚本 | P6 可视化 |
| `plotTunes.py` | 调谐曲线分析（Python3） | P6 可视化 |

### 绘图脚本（本项目新增，均在仓库根目录下运行）

| 脚本 | 输入（由 `make check` 或单独运行产生） | 输出 |
| --- | --- | --- |
| `plot_field_map.py` | `tests/cycl_track_map.vti` | `docs/img/cyclotron_field_map.png`（整机中场图，**不含任何轨迹**） |
| `plot_closed_orbit.py` | `tests/cycl_closed_orbit.vtp`（+ 可选 `cycl_track_map.vti` 作底图） | `docs/img/cyclotron_real_orbit.png`（**闭合轨道，单能量**） |
| `plot_acceleration.py` | `tests/cycl_accel_locked.csv`、`cycl_accel_real.csv` | `docs/img/cyclotron_acceleration.png`（Part B 为**已知限制**） |
| `view_3d.py` | 任意 `.vti` / `.vtp` | 三维交互视图或 PNG（需 `pyvista`） |

典型流程::

```bash
make -C tests cycl_track && ./tests/cycl_track                 # 校验 ω_c + 整机中场图
python3 examples/cyclotron/plot_field_map.py                    # 场图
make -C tests cycl_closed_orbit && ./tests/cycl_closed_orbit    # 闭合轨道
python3 examples/cyclotron/plot_closed_orbit.py \
    --orbit tests/cycl_closed_orbit.vtp --map tests/cycl_track_map.vti \
    -o docs/img/cyclotron_real_orbit.png
python3 examples/cyclotron/view_3d.py tests/cycl_track_map.vti  # 三维交互（需 pyvista）
```

> ⚠️ 已删除的脚本：`plot_trajectory.py` 与图 `cyclotron_orbit.png`。它们展示的是一条
> 半径仅 0.2 m 的**局部拉莫尔回旋**（0.1c 质子），不是回旋加速器轨道，属误导。

> 注：`tests/cycl_track` 默认从 `tests/` 目录运行（它按 `../examples/cyclotron/data/` 找场图）。

### 关键输入参数（摘自 `cyclotron2.in` 等）

```
RingCycl: CYCLOTRON, TYPE=RING, CYHARMON=6,
          PHIINIT=0.0, PRINIT=pr0, RINIT=r0, SYMMETRY=8.0,
          RFFREQ=f1, FMAPFN="bfield.dat", FMLOWE=72, FMHIGHE=590;
```

- `SYMMETRY=8.0` → 8 折对称扇形磁铁（**非轴对称**，IBSimu 原生磁场无法处理）
- `FMAPFN="bfield.dat"` → 外部磁场图
- `FMLOWE/FMHIGHE` → 72–590 MeV

## 算例来源（可访问的 GitHub 镜像）

OPAL 官方源码托管在 `gitlab.psi.ch`（部分网络环境不可达）。可用 GitHub 上的镜像：

- <https://github.com/OPALX-project/OPAL> — 经典 OPAL（C++），含 `samples/cyclotron/`
- <https://github.com/OPALX-project/OPALX> — 新一代 OPAL（Exascale）
- <https://github.com/OPALX-project/Manual-old> — 手册与算例文件（`examples/cyclotron.md`、`Cyclotron/*.dat`）
- 在线文档：<https://opalx-project.github.io/Examples/cyclotron>

## 数据（已导入 `data/`）

`data/` 下已导入来自 OPAL 的 PSI Ring 算例数据（来源与许可见 `data/README.md`）：

| 文件 | 内容 |
| --- | --- |
| `bfield.dat` | CERN FIELD 格式三维（中平面）磁场图，141 × 135 网格，1/8 扇区，1.2 MB |
| `rffield1.dat` / `rffield2.dat` | 一维径向 RF 场剖面 |
| `cyclotron1.in` / `cyclotron2.in` | OPAL 输入文件（参考） |

复现来源：

```bash
git clone --depth 1 https://github.com/OPALX-project/Manual-old /tmp/opal-manual
ls /tmp/opal-manual/examples/Cyclotron/
```

> OPAL 与 IBSimu 均为 GPL 系许可证；我们只使用其**算例数据**作为验证基准，不使用其代码。
> 来源、格式说明与许可见 `data/README.md` 与 `THIRD_PARTY_NOTICES.md`。
