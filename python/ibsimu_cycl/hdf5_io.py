"""IBSimu-Cycl HDF5 后端读取器（ROADMAP R3-④）。

本模块读取 ``src/io/hdf5writer.cpp`` 写出的文件。布局见 ``hdf5writer.hpp``：

- 网格：``/data/<iter>/meshes/<name>/``，属性沿用 openPMD 2.0 的命名
  （``gridSpacing`` / ``gridGlobalOffset`` / ``axisLabels`` / ``unitSI``…），
  数据集形状为 **(nz, ny, nx)**（轴序反转），因此 C 序连续内存即 **i 最快**，
  与 ``vtk_io`` 的 ``.vti`` 完全一致。
- 轨迹：``/data/<iter>/particles/<name>/``，扁平的 ``position(N,3)`` /
  ``time(N)`` + ``offset(M+1)`` / ``count(M)``。

本模块是 C++ 写出器的**独立实现**（只共用格式约定，不共用代码），因此
「写出 → 这里读回」构成跨语言交叉验证。

依赖：numpy + **h5py**（软依赖：未安装时只有调用读取函数才会报错，
``available()`` 可用于探测，便于在没有 h5py 的环境里跳过相关测试）。
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

from .vtk_io import ImageData, PolyData

try:  # 软依赖
    import h5py  # type: ignore
    HAVE_H5PY = True
except ImportError:  # pragma: no cover - 取决于环境
    h5py = None  # type: ignore
    HAVE_H5PY = False

H5PY_HINT = ("读取 HDF5 需要 h5py（pip install h5py 或 pacman -S python-h5py；"
             "Debian/Ubuntu: apt install python3-h5py）")

__all__ = ["HAVE_H5PY", "available", "Mesh", "read", "read_meshes", "read_tracks"]


def available() -> bool:
    """是否装有 h5py（没有时读取函数会抛 ImportError）。"""
    return HAVE_H5PY


def _require_h5py() -> None:
    if not HAVE_H5PY:
        raise ImportError(H5PY_HINT)


def _text(value) -> str:
    """h5py 的字符串属性可能是 str / bytes / vlen / 定长，统一成 str。"""
    if isinstance(value, bytes):
        return value.decode()
    if isinstance(value, np.ndarray):
        return _text(value.flat[0] if value.size else b"")
    return str(value)


@dataclass
class Mesh:
    """一个网格组（mesh）的内容。"""

    name: str
    data: Dict[str, np.ndarray]           # 数据集名 -> 原样数组（标量 (nz,ny,nx)）
    units: Dict[str, float]               # 数据集名 -> unitSI
    grid_spacing: np.ndarray              # (3,)
    grid_global_offset: np.ndarray        # (3,)
    axis_labels: List[str]
    geometry: str = "cartesian"

    @property
    def dims(self) -> Tuple[int, int, int]:
        """(nx, ny, nz)——注意与数据集的 (nz,ny,nx) 正好相反。

        只看**三维标量**数据集：矢量数据集是 (nz,ny,nx,3)，直接取末三维
        会得到 (ny,nx,3) 这种错误的尺寸（曾被测试抓到）。
        """
        for arr in self.data.values():
            if arr.ndim == 3:
                nz, ny, nx = arr.shape
                return (int(nx), int(ny), int(nz))
        raise ValueError(f"mesh '{self.name}' 里没有三维标量数据集")

    def _spatial(self, name: str) -> np.ndarray:
        """把数据集整理成 i 最快的扁平数组（标量 (n,) / 矢量 (n,3)）。"""
        arr = self.data[name]
        if arr.ndim == 3:
            return arr.reshape(-1)
        if arr.ndim == 4 and arr.shape[-1] == 3:
            # (nz,ny,nx,3)：按点交错，C 序连续 ⇒ 直接展平
            return arr.reshape(-1, 3)
        if arr.ndim == 4 and arr.shape[0] == 3:
            # 兼容分量优先的 (3,nz,ny,nx) 布局
            return np.moveaxis(arr, 0, -1).reshape(-1, 3)
        raise ValueError(f"'{name}': 不支持的形状 {arr.shape}")

    def scalar(self, name: str) -> np.ndarray:
        return self._spatial(name)

    # origin/spacing 是 ImageData 里的字段名；这里给同义属性，
    # 使得 HDF5 的 Mesh 与 VTK 的 ImageData 可以直接互换使用。
    @property
    def origin(self) -> np.ndarray:
        return self.grid_global_offset

    @property
    def spacing(self) -> np.ndarray:
        return self.grid_spacing

    def scalar_grid(self, name: str) -> np.ndarray:
        """重排为 (nz, ny, nx)，与 ``ImageData.scalar_grid`` 一致。"""
        nx, ny, nz = self.dims
        return self.scalar(name).reshape(nz, ny, nx)

    def vector(self, name: str) -> np.ndarray:
        v = self._spatial(name)
        if v.ndim != 2 or v.shape[1] != 3:
            raise ValueError(f"'{name}' 不是三分量矢量数据集")
        return v

    def magnitude_grid(self, name: str) -> np.ndarray:
        nx, ny, nz = self.dims
        return np.linalg.norm(self.vector(name), axis=1).reshape(nz, ny, nx)

    def as_image_data(self) -> ImageData:
        """转换成与 ``vtk_io`` 相同的 ``ImageData``（两套格式可互换使用）。"""
        arrays = {name: self._spatial(name) for name in self.data}
        return ImageData(self.dims, self.grid_global_offset.copy(),
                         self.grid_spacing.copy(), arrays)


def _iterations(handle) -> List[str]:
    base = "/data"
    if base not in handle:
        return []
    return sorted(str(k) for k in handle[base].keys())


def read_meshes(path: str | Path) -> Dict[str, Mesh]:
    """读取文件里所有网格组，返回 ``{name: Mesh}``。"""
    _require_h5py()
    out: Dict[str, Mesh] = {}
    with h5py.File(str(path), "r") as f:
        for it in _iterations(f):
            g = f[f"/data/{it}"]
            if "meshes" not in g:
                continue
            for name in g["meshes"].keys():
                mg = g["meshes"][name]
                attrs = mg.attrs
                data: Dict[str, np.ndarray] = {}
                units: Dict[str, float] = {}
                for dname in mg.keys():
                    ds = mg[dname]
                    if not isinstance(ds, h5py.Dataset):
                        continue
                    data[dname] = np.asarray(ds[...])
                    units[dname] = float(ds.attrs.get("unitSI", 1.0))
                out[name] = Mesh(
                    name=name,
                    data=data,
                    units=units,
                    grid_spacing=np.asarray(
                        attrs.get("gridSpacing", [1.0, 1.0, 1.0]), dtype=float),
                    grid_global_offset=np.asarray(
                        attrs.get("gridGlobalOffset", [0.0, 0.0, 0.0]),
                        dtype=float),
                    axis_labels=[_text(v) for v in
                                 np.asarray(attrs.get("axisLabels", []))],
                    geometry=_text(attrs.get("geometry", "cartesian")),
                )
    return out


def read_tracks(path: str | Path, name: str | None = None) -> PolyData:
    """读取粒子轨迹，返回 ``PolyData``（与 ``.vtp`` 的读取结果同型）。

    ``name`` 为轨迹组名；为空时取文件里第一个。
    """
    _require_h5py()
    with h5py.File(str(path), "r") as f:
        found = []
        for it in _iterations(f):
            g = f[f"/data/{it}"]
            if "particles" in g:
                found += [str(k) for k in g["particles"].keys()]
        if not found:
            raise ValueError(f"{path}: 文件里没有 particles 组")
        if name is None:
            name = found[0]
        elif name not in found:
            raise ValueError(f"{path}: 没有轨迹组 '{name}'（有 {found}）")

        base = None
        for it in _iterations(f):
            if f"/data/{it}/particles/{name}" in f:
                base = f[f"/data/{it}/particles/{name}"]
                break
        assert base is not None
        position = np.asarray(base["position"][...], dtype=float)
        time = np.asarray(base["time"][...], dtype=float)
        if "offset" in base:
            offset = np.asarray(base["offset"][...], dtype=np.int64)
            lines = [np.arange(offset[m], offset[m + 1], dtype=np.int64)
                     for m in range(len(offset) - 1)]
        else:  # pragma: no cover - 本写出器总是写 offset
            lines = [np.arange(len(position))]
    return PolyData(points=position, lines=lines, arrays={"t": time})


def read(path: str | Path):
    """读取只含单一内容的 HDF5 文件（网格→``ImageData``，轨迹→``PolyData``）。

    同时含网格与轨迹时请显式调用 :func:`read_meshes` / :func:`read_tracks`。
    """
    meshes = read_meshes(path)
    tracks = []
    try:
        _require_h5py()
        with h5py.File(str(path), "r") as f:
            for it in _iterations(f):
                g = f[f"/data/{it}"]
                if "particles" in g:
                    tracks += [str(k) for k in g["particles"].keys()]
    except (KeyError, OSError):  # pragma: no cover
        pass

    if meshes and tracks:
        raise ValueError(
            f"{path}: 同时含 {len(meshes)} 个网格与 {len(tracks)} 组轨迹，"
            "请用 read_meshes() / read_tracks() 明确指定")
    if meshes:
        if len(meshes) > 1:
            raise ValueError(f"{path}: 含多个网格 {list(meshes)}，请用 read_meshes()")
        return next(iter(meshes.values())).as_image_data()
    if tracks:
        return read_tracks(path)
    raise ValueError(f"{path}: 既没有网格也没有轨迹")
