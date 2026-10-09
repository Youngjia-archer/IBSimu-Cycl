# 算例数据来源与许可

本目录数据来自 **OPAL** 项目的回旋加速器示例（PSI Ring，590 MeV，8 扇形区），
用于 IBSimu-Cycl 的物理验证。

| 文件 | 来源 |
| --- | --- |
| `bfield.dat` | `OPALX-project/Manual-old` → `examples/Cyclotron/bfield.dat` |
| `rffield1.dat`, `rffield2.dat` | 同上 |
| `cyclotron1.in`, `cyclotron2.in` | 同上 |

- 上游：<https://github.com/OPALX-project/Manual-old>
- OPAL 以 GPL 发布；此处仅作为**验证基准数据**使用，**未使用其代码**。
- 引用与许可见仓库根 `THIRD_PARTY_NOTICES.md`。

## `bfield.dat` 格式（CERN "FIELD"）

- 第 1 个数：`rmin` [mm]；第 2 个：`delr` [mm]（**负值取倒数**）
- 第 3 个：`tetmin` [deg]；第 4 个：`dtet` [deg]（**负值取倒数**）
- 随后为标签/信息行：`LABEL=`、`CFELD=FIELD`、`NREC=`、`NPAR=`、`LPAR=`、`IENT=`、`IPAR=`，
  以及一段由 `lpar` 个双精度构成的参数块（内容忽略）
- 扫描到 `LREC=` 后进入数据区；数据按**半径**分块（共 `nrad` 块），
  每块包含 4 × `ntet` 个 `%16lE` 数值：`Bz`、`dBz/dθ`、`d²Bz/dθ²`、`d³Bz/dθ³`
  （首块无前缀信息行，其余每块前有 6 个信息记号）
- 单位：**kGauss**（×0.1 → Tesla）
- 利用 **8 折旋转对称**只存 1/8 扇区（`ntet = 135`，`dtet = 1/3°`）

参考实现：OPAL `Cyclotron::getFieldFromFile_Ring`
(`src/Classic/AbsBeamline/Cyclotron.cpp`)。

本项目的解析器：`src/cyclotron/ringmap.{hpp,cpp}` → `ibsimu_cycl::read_ring_field_map()`。
