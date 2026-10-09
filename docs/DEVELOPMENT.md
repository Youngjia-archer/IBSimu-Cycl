# 开发说明（DEVELOPMENT）

## 环境依赖

| 类别 | 依赖 |
| --- | --- |
| 工具链 | `build-essential` `autoconf` `automake` `libtool` `pkg-config` |
| 必需库 | `zlib` `libpng` `fontconfig` `freetype2` `libm` `gsl` |
| 可选库 | `gtk+-3.0`（GUI）、`gtkglext-3.0`（OpenGL 渲染）、`UMFPACK`、`csg` |

Ubuntu/Debian 一键安装：

```bash
sudo apt-get install -y build-essential autoconf automake libtool pkg-config \
    zlib1g-dev libpng-dev libfontconfig1-dev libfreetype-dev libgsl-dev \
    libgtk-3-dev
```

> 无 GUI 环境（服务器/CI）请加 `--without-gtk3 --without-opengl`。

## 构建

```bash
./reconf                              # 重新生成 configure（autotools）
./configure [--without-gtk3 --without-opengl]
make -j"$(nproc)"
make check                            # 运行上游测试套件
```

产物：`src/.libs/libibsimu-1.0.6dev.{a,so}`。

## 构建系统约定

P0 阶段**沿用上游 GNU autotools**（不引入 CMake），以保证始终可构建。
上游源码是**扁平布局**（全部在 `src/`）。本项目新增模块放在子目录
`src/<module>/`，并采用 automake 支持的"子目录源文件直接列入 `_SOURCES`"方式接入：

1. 头文件包含路径：在 `src/Makefile.am` 顶部
   ```make
   AM_CPPFLAGS = -Wall -I$(srcdir)/cyclotron -I$(srcdir)/<module>
   ```
2. 源文件：在 `src/Makefile.am` 的
   `libibsimu_1_0_6dev_la_SOURCES` 列表**末尾**追加
   （注意该文件里 `xygraph.hpp` 出现两次：SOURCES 与 include_HEADERS，
   追加时要锚定 `SOURCES` 那一处，它后面紧跟 `_includedir=` 行）：
   ```make
   xygraph.hpp \
   cyclotron/cycl_info.cpp \
   cyclotron/cycl_info.hpp
   ```
3. 若为公开安装头文件，需使用 `nobase_` 前缀
   （`nobase_libibsimu_1_0_6dev_la_include_HEADERS`）以支持子目录；
   内部头文件只需列在 `_SOURCES` 中即可。

> GPU/Kokkos 阶段（P5）如需 CUDA 编译规则，再评估是否引入 CMake 或自定义 automake 规则。

## 代码风格与许可证头

- 新增源文件请在文件头注明用途，并标注项目许可 **GPL-3.0-or-later**
  （上游文件保持其原始 GPL-2.0-or-later 头不动）。
- 命名空间统一使用 `ibsimu_cycl`。

## 开发环境备忘（本仓库作者机器）

- OS：Arch Linux；Shell：**fish**（注意：写命令用 fish 语法，
  不要用 bash 的 `for .. do/done`、`[[ ]]`；`ls` 被别名成带 icons 的工具，
  脚本中请用 `command ls`）。
- GPU：NVIDIA RTX 2060 Max-Q（CUDA 可用；暂不考虑非 NVIDIA 卡）。

## 验证算例

见 `examples/cyclotron/README.md`（PSI Ring，OPAL-cycl 算例）。
