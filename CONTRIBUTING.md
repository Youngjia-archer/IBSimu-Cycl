# 贡献指南

感谢参与 IBSimu-Cycl。本项目是 [IBSimu](https://ibsimu.sourceforge.net/) 的派生作品，
请在贡献前阅读以下约定。

## 开发流程

1. 从 `master` 切出特性分支，例如 `feat/fieldmap3d`、`fix/...`。
2. 遵循 `docs/DEVELOPMENT.md` 的构建与代码规范。
3. 提交前确保：
   ```bash
   ./reconf && ./configure && make -j"$(nproc)" && make check
   ```
4. 提交 Pull Request，描述中说明**动机 / 改动 / 验证方式**。

## 许可证

- 本项目整体以 **GPL-3.0-or-later** 发布。
- **不得移除**上游 IBSimu 源文件中的原始版权与许可证声明。
- 新增文件请加入 `GPL-3.0-or-later` 许可头。
- 引入新的第三方依赖前，请先在 `THIRD_PARTY_NOTICES.md` 登记并确认许可证兼容性。

## 物理正确性

涉及物理的改动（场求解、粒子推进、空间电荷等）应附带**验证**：

- 解析解对照（如均匀场圆周运动、圆柱腔 TM010 频率），或
- 与 `OPAL-cycl` 标准算例对比（见 `examples/cyclotron/`）。

## 提交信息

使用简洁的祈使句，例如：

```
Add CVectorField3D with trilinear interpolation

Implements a general 3D vector field on rectilinear/cylindrical grids,
replacing the axisymmetric-only CAxisymmetricVectorField for cyclotron
field maps.
```
