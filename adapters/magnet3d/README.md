# magnet3d —— 三维静磁铁磁场适配器（Elmer）

本目录把**外部 FEM 求解器**（首选 [Elmer](https://www.elmerfem.org/)）离线算出的
三维磁场转换为 IBSimu-Cycl 可读的场图格式，供 `CCylFieldMap3D` / `CFieldMap3D`
读取并驱动 3D 粒子跟踪。

> **松耦合设计**：Elmer 只是外部程序，通过**文件**交换数据；本项目不链接其代码，
> 从而避免许可证组合问题（Elmer 为 GPL，兼容本项目 GPL-3.0-or-later）。

## 流水线

```
几何 (.geo/.step)  ──ElmerGrid/gmsh──►  网格
        │
        └── ElmerSolver (WhitneyAV, 含铁芯非线性 B-H) ──►  B(x,y,z)  (VTU / ASCII)
                                                              │
                                          field_to_fieldmap3d.py  (重采样到规则网格)
                                                              │
                                                    IBSIMU-CYCL-FIELDMAP  (ASCII)
                                                              │
                                            CFieldMap3D::load() / CCylFieldMap3D
```

## 目录

| 文件 | 说明 |
| --- | --- |
| `elmer/magnetostatic.sif` | Elmer 三维静磁（WhitneyAV）求解器模板，含非线性铁芯与远场边界 |
| `elmer/sector.geo` | gmsh 几何模板：扇形磁铁（铁轭 + 磁极 + 气隙 + 空气域） |
| `field_to_fieldmap3d.py` | 把 `x y z Bx By Bz` 散点/节点数据重采样为规则网格并写出本项目场图格式 |
| `tests/test_field_to_fieldmap3d.py` | 转换器的自检（解析场 → 散点 → 重采样 → 校验） |

## 前置条件

- **Elmer**（可选，仅离线算场时需要）：`apt install elmerfem-csc` / conda `-c conda-forge elmerfem` / 源码构建。
- **gmsh**（可选）：生成 `.msh` 网格。
- **Python 3 + numpy**（转换器必需）。

## 使用步骤

```bash
# 1) 生成网格（示例几何，替换为你的实际磁铁尺寸）
gmsh -3 elmer/sector.geo -o sector.msh

# 2) 求解磁场（Elmer）
ElmerGrid 14 2 sector.msh -autoclean
ElmerSolver elmer/magnetostatic.sif          # 输出 magnet.vtu 等

# 3) 从 VTU 导出节点数据为 x y z Bx By Bz 文本
#    （ParaView: File > Save Data 选 ASCII；或使用 Elmer 的 SaveData/ASCII 输出）

# 4) 重采样为规则网格并写出本项目场图
python3 field_to_fieldmap3d.py nodes_B.txt bfield_ibs.txt \
        --nx 141 --ny 141 --nz 61 --k 8 --power 2

# 5) C++ 侧读取（示例）
#    CFieldMap3D B; B.load("bfield_ibs.txt");
```

## 结果的物理验证

务必用通用判据校验 FEM 结果（本项目在 P2 已用同类判据验证 PSI Ring 场图）：

- 中平面 `Bz(r)` 回旋频率 $\omega_c=qB/m$；
- 真空区 $\nabla\cdot\mathbf B\approx 0$、$\nabla\times\mathbf B\approx 0$；
- 铁芯饱和后 `|B|` 不随激励线性增长（非线性正确性）。

## 备注 / 待办

- 目前 `field_to_fieldmap3d.py` 支持**正交规则网格**输出；柱坐标 (r,θ,z) 输出与
  `CCylFieldMap3D` 原生格式的直转在 P3 一并补齐。
- Elmer 的 VTU → `x y z Bx By Bz` 文本一步建议用 ParaView 的 *Save Data (ASCII)*，
  避免在本项目内引入 XML/HDF5 解析依赖。
