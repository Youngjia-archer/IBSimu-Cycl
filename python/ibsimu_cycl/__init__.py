"""IBSimu-Cycl Python 工具包。

现有功能：
- `vtk_io`：极简 VTK XML（.vti/.vtp）读取器（仅需 numpy），
  既是 C++ 导出器的**独立交叉验证**，也是无 VTK 环境下的后处理手段。
- `hdf5_io`：HDF5 二进制后端（ROADMAP R3-④）的读取器，面向大网格/大粒子数；
  与 `vtk_io` 返回同型的 `ImageData`/`PolyData`，两种格式可互换使用。
  h5py 为**软依赖**：未安装时 `hdf5_io.available()` 返回 False。

计划功能（P6）:
- 场图读写（openPMD 完整互操作）
- 后处理与可视化（pyvista / matplotlib）
- 参数扫描与批量运行
"""

from . import hdf5_io  # noqa: F401
from . import vtk_io  # noqa: F401

__version__ = "0.1.0"
