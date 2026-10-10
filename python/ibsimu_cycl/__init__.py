"""IBSimu-Cycl Python 工具包。

现有功能：
- `vtk_io`：极简 VTK XML（.vti/.vtp）读取器（仅需 numpy），
  既是 C++ 导出器的**独立交叉验证**，也是无 VTK 环境下的后处理手段。

计划功能（P6）:
- 场图读写（openPMD / HDF5 大规模后端）
- 后处理与可视化（pyvista / matplotlib）
- 参数扫描与批量运行
"""

from . import vtk_io  # noqa: F401

__version__ = "0.1.0"
