# 第三方组件与许可证清单

本项目整体以 **GPL-3.0-or-later** 发布。下表列出**已使用**或**计划整合**的第三方组件
及其许可证。凡计划整合的组件，在正式链接/分发前必须完成许可证兼容性复核（标记为"待核实"）。

> 关键区别：**调用外部程序处理其输出**（如用 Elmer 生成磁场图 `.dat`）通常不构成"链接"，
> 不产生许可证组合问题；只有**链接库代码**（如链接 Kokkos、openPMD）才需要严格的许可证兼容。
> 本项目优先采用"外部程序 + 文件交换"的松耦合方式。

## 已使用

| 组件 | 版本/来源 | 许可证 | 整合方式 | 兼容性 |
| --- | --- | --- | --- | --- |
| IBSimu | 上游 `09beedb` | GPL-2.0-or-later | 源码派生 | ✅ 本项目升级为 GPL-3.0-or-later |
| zlib | 系统 | Zlib | 链接 | ✅ 宽松许可，兼容 |
| libpng | 系统 | libpng-2.0 | 链接 | ✅ |
| fontconfig | 系统 | MIT-like | 链接 | ✅ |
| freetype2 | 系统 | FTL / GPL-2.0-or-later | 链接 | ✅ |
| GTK+3（可选） | 系统 | LGPL-2.1-or-later | 链接（可选） | ✅ |
| cairo（可选） | 系统 | LGPL-2.1 / MPL-1.1 | 链接（可选） | ✅ |

## 计划整合（离线生成场图，外部程序）

| 组件 | 用途 | 许可证 | 整合方式 | 兼容性 |
| --- | --- | --- | --- | --- |
| Elmer FEM | 三维静磁（铁芯非线性） | GPL-2.0-or-later | 外部程序 + 文件 | ✅（待核实具体版本条款） |
| Palace (AWS Labs) | RF 腔 3D 本征模 | Apache-2.0 | 外部程序 + 文件 | ✅ 与 GPL-3.0 兼容 |
| GetDP + Gmsh (ONELAB) | 备选 FEM | GPL-2.0-or-later | 外部程序 + 文件 | ✅ |
| Radia (ESRF) | 线圈/超导磁体 | GPL-2.0-or-later | 外部程序 | ✅（待核实） |
| OPAL / OPALX | 验证基准（算例数据） | GPL-2.0-or-later | 仅使用算例数据 | ✅（待核实） |

## 计划整合（库/工具链，需严格复核）

| 组件 | 用途 | 许可证 | 整合方式 | 兼容性 |
| --- | --- | --- | --- | --- |
| openPMD-api | 场图/粒子 IO | LGPL-3.0-or-later | 链接 | ⚠️ 待核实（LGPL 动态链接通常可行） |
| HDF5 (The HDF Group) | 容器 | BSD-3 | 链接（可选，R3-③ 计划） | ✅ |
| Kokkos | 性能可移植（CPU/GPU） | Apache-2.0 | 链接 | ⚠️ 仅在与 GPL-3.0 组合时可链接 |
| MPI | 多节点 | 各实现许可 | 链接 | ✅ |
| hypre | 并行 AMG | MIT / Apache-2.0 | 链接 | ✅ |
| PETSc | 线性求解 | BSD-2 | 链接 | ✅ |
| AMGx (NVIDIA) | GPU AMG | BSD-3 | 链接（可选） | ✅ |
| VTK | 3D 可视化 IO | BSD-3 | **未链接**：直接写 VTK XML 文本格式（`src/io/vtkwriter`） | ✅ |
| pyvista | Python 可视化 | MIT | 依赖 | ✅ |
| ADIOS2 | 高性能 IO | Apache-2.0 | 依赖 | ⚠️ 与 GPL-3.0 兼容 |

**待办（P0）**：逐项确认上表中标记"待核实"的许可证条款，并在引入时更新本文件。
