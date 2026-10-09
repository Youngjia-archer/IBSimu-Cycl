# rfcavity3d —— 射频谐振腔适配器（Palace）

把**外部 3D 电磁本征模求解器**（首选 [Palace](https://github.com/awslabs/palace)，AWS Labs）
算出的腔体模式场与谐振频率转换为 IBSimu-Cycl 可读的场图，配合本项目
`CTimeVaryingField` 实现**时变 RF 场中的束流跟踪**。

> **松耦合**：与 `magnet3d` 相同，Palace 只作为外部程序通过文件交换数据（不链接），
> 避免许可证组合问题。

## 流水线

```
几何 (cavity.geo)  ──gmsh──►  网格 (cavity.msh)
        │
        └── Palace (Eigenmode) ──►  f [MHz], Q, R/Q, 模式场 (E, B)
                                              │
                        mode_field_to_fieldmap.py (重采样 / 直转)
                                              │
                              IBSIMU-CYCL-FIELDMAP (E 与 B 各一份)
                                              │
                CFieldMap3D + CTimeVaryingField ──► CBorisPusher (时变场跟踪)
```

## 目录

| 文件 | 说明 |
| --- | --- |
| `palace/cavity.json` | Palace 本征模仿真配置模板（`Problem.Type = Eigenmode`） |
| `palace/cavity.geo` | gmsh 圆柱腔几何模板 |
| `pillbox_tm010.py` | **解析**圆柱腔 TM010 模式生成器（频率 + E/B 场图），用作求解器的验证基准 |
| `tests/test_pillbox_tm010.py` | 解析解的 Maxwell 自检 |

## 解析验证基准：圆柱腔 TM010

半径为 $R$ 的理想圆柱腔，TM010 模：

$$k=\frac{x_{01}}{R},\qquad f=\frac{c\,x_{01}}{2\pi R},\qquad E_z(r)=E_0 J_0(kr),\qquad B_\theta(r)=\frac{E_0}{c}J_1(kr)$$

其中 $x_{01}=2.4048255577$ 为 $J_0$ 的第一个零点。

```bash
python3 pillbox_tm010.py --R 0.30 --E0 1.0e6 --nx 61 --ny 61 --nz 11 --outdir out
# 生成 out/pillbox_E.txt, out/pillbox_B.txt, out/pillbox_info.txt
python3 tests/test_pillbox_tm010.py
```

**这是求解器（Palace/Elmer）的正确性判据**：算出的本征频率应逼近 $f$，
模式场结构应满足 $E_z\propto J_0(kr)$、$B_\theta\propto J_1(kr)$，
且满足真空中 $\nabla\times\mathbf E=i\omega\mathbf B$。

## Palace 使用步骤（需自行安装 Palace）

```bash
# 1) 网格
gmsh -3 palace/cavity.geo -o cavity.msh

# 2) 本征模求解
palace palace/cavity.json          # 或在集群上: mpirun -np N palace cavity.json

# 3) 导出模式场为 x y z Ex Ey Ez / Bx By Bz 文本，再转为本项目场图
#   （与 magnet3d/field_to_fieldmap3d.py 的输入格式一致，可复用其重采样）

# 4) C++ 侧：用 CTimeVaryingField 包装场图，配合 CBorisPusher 做时变场跟踪
```

## 待办

- `mode_field_to_fieldmap.py`（Palace/HDF5 模式场 → 本项目场图）在装好 Palace 后补齐；
  当前可先用 `magnet3d/field_to_fieldmap.py` 处理导出的散点场。
- 腔体 R/Q、Q 值读取与记录。
