"""极简 VTK XML 读取器（.vti / .vtp），只依赖标准库 + numpy。

为什么单独实现一个读取器，而不直接依赖 VTK/PyVista：

1. **交叉验证**：C++ 写出器（``src/io/vtkwriter.cpp``）与这里的读取器是两套
   独立实现。用后者去校验前者产生的文件，才能证明格式真的对，而不是"自己验自己"。
2. **零依赖后处理**：不装 VTK 也能做数值检查、统计与简单绘图。

只支持本项目写出器产生的 **ASCII 内联**格式（VTK XML 的一个子集），
不追求覆盖 VTK 全特性。

典型用法::

    from ibsimu_cycl.vtk_io import read

    fd = read("cycl_vtk_field.vti")        # ImageData
    fd.dims, fd.origin, fd.spacing
    bmag = fd.scalar("Bmag")               # 一维，i 最快
    grid = fd.scalar_grid("Bmag")          # (nz, ny, nx) 便于切片/绘图

    traj = read("cycl_vtk_orbit.vtp")      # PolyData
    traj.points.shape, len(traj.lines)
    t = traj.arrays["t"]
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np


def _floats(text: str | None) -> np.ndarray:
    """把 DataArray 的文本体解析为一维 float64 数组（按空白切分）。"""
    if not text or not text.strip():
        return np.zeros(0, dtype=np.float64)
    return np.asarray(text.split(), dtype=np.float64)


def _collect(elem: ET.Element, xpath: str) -> Dict[str, np.ndarray]:
    """收集 xpath 下所有 DataArray，返回 name -> 数组。

    多分量数组返回形状 (n, ncomp)，标量返回 (n,)。
    """
    out: Dict[str, np.ndarray] = {}
    for da in elem.findall(xpath):
        name = da.get("Name")
        if name is None:
            continue
        ncomp = int(da.get("NumberOfComponents", "1"))
        vals = _floats(da.text)
        if ncomp > 1:
            if vals.size % ncomp != 0:
                raise ValueError(f"DataArray '{name}': 长度不是分量数的整数倍")
            vals = vals.reshape(-1, ncomp)
        out[name] = vals
    return out


@dataclass
class ImageData:
    """规则网格数据（VTK ImageData）。"""

    dims: Tuple[int, int, int]        # (nx, ny, nz) 节点数
    origin: np.ndarray                # (3,)
    spacing: np.ndarray               # (3,)
    arrays: Dict[str, np.ndarray]     # name -> (n,) 或 (n,3)，i 最快

    def scalar(self, name: str) -> np.ndarray:
        return self.arrays[name].reshape(-1)

    def vector(self, name: str) -> np.ndarray:
        v = self.arrays[name]
        if v.ndim != 2 or v.shape[1] != 3:
            raise ValueError(f"'{name}' 不是三分量矢量数组")
        return v

    def scalar_grid(self, name: str) -> np.ndarray:
        """重排为 (nz, ny, nx)，配合 matplotlib ``imshow`` / 切片更直观。"""
        nx, ny, nz = self.dims
        return self.scalar(name).reshape(nz, ny, nx)

    def magnitude_grid(self, name: str) -> np.ndarray:
        nx, ny, nz = self.dims
        return np.linalg.norm(self.vector(name), axis=1).reshape(nz, ny, nx)

    def axis(self, d: int) -> np.ndarray:
        """第 d 个方向的节点坐标。"""
        return self.origin[d] + self.spacing[d] * np.arange(self.dims[d])


@dataclass
class PolyData:
    """折线/点集数据（VTK PolyData）。"""

    points: np.ndarray                # (npts, 3)
    lines: List[np.ndarray]           # 每条折线的点索引
    arrays: Dict[str, np.ndarray]     # 点数据

    def __len__(self) -> int:
        return len(self.points)


def read(path: str) -> ImageData | PolyData:
    """读取 .vti 或 .vtp（自动识别类型）。"""
    root = ET.parse(path).getroot()
    typ = root.get("type")

    if typ == "ImageData":
        img = root.find("ImageData")
        if img is None:
            raise ValueError(f"{path}: 缺少 <ImageData>")
        ext = [int(v) for v in img.get("WholeExtent", "").split()]
        if len(ext) != 6:
            raise ValueError(f"{path}: WholeExtent 非法")
        dims = (ext[1] - ext[0] + 1, ext[3] - ext[2] + 1, ext[5] - ext[4] + 1)
        origin = np.asarray(img.get("Origin", "0 0 0").split(), dtype=np.float64)
        spacing = np.asarray(img.get("Spacing", "1 1 1").split(), dtype=np.float64)
        return ImageData(dims, origin, spacing,
                         _collect(img, ".//PointData/DataArray"))

    if typ == "PolyData":
        pd = root.find("PolyData")
        if pd is None:
            raise ValueError(f"{path}: 缺少 <PolyData>")
        piece = pd.find("Piece")
        if piece is None:
            raise ValueError(f"{path}: 缺少 <Piece>")

        pts = _collect(piece, ".//Points/DataArray")
        points = pts.get("Points")
        if points is None:
            # 允许无名坐标数组
            points = next(iter(pts.values())) if pts else np.zeros((0, 3))

        arrays = _collect(piece, ".//PointData/DataArray")

        conn = _collect(piece, ".//Lines/DataArray")
        c, o = conn.get("connectivity"), conn.get("offsets")
        lines: List[np.ndarray] = []
        if c is not None and o is not None:
            start = 0
            for end in o.astype(int):
                lines.append(c[start:end].astype(int))
                start = end
        return PolyData(points, lines, arrays)

    raise ValueError(f"{path}: 不支持的 VTK 类型 '{typ}'")
