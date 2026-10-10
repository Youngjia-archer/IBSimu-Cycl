/** @file hdf5writer.hpp
 *  @brief HDF5 二进制写出器：大网格场图与大量粒子轨迹（ROADMAP R3-④）。
 *
 *  ## 为什么需要它
 *
 *  `vtkwriter`（VTK XML）是**纯文本**的：一个 200³ 的网格就接近 1 GB，
 *  写盘时间以分钟计。HDF5 通路提供：
 *    - **二进制 + deflate 压缩**（光滑场图通常压到 1/4 以下）；
 *    - 分块存储，便于后续做部分读取；
 *    - 属性里自带**单位与网格定义**，数据自描述。
 *
 *  ## 与 openPMD 的关系（重要，避免误读）
 *
 *  属性的**命名与单位约定沿用 openPMD 2.0**（`basePath`、`meshes`、`gridSpacing`、
 *  `gridGlobalOffset`、`gridUnitSI`、`unitSI`、`timeOffset`、`axisLabels`…），
 *  便于和 openPMD 生态互操作。但本写出器**只实现了网格 + 轨迹这两个子集，
 *  并未逐条核对 openPMD 规范**，因此不要把它当作“openPMD 认证”的产出；
 *  需要严格互操作时，用 h5py/openPMD-api 转写即可（结构是刻意对齐的）。
 *
 *  粒子的轨迹布局是本项目的**简化布局**（扁平的 `position`/`time` 数组 +
 *  `offset`/`count` 前缀和），不是 openPMD 的 `particles` 记录——
 *  后者以“每个时间步一份粒子表”为模型，与“每个粒子一条轨迹”不是同一件事。
 *
 *  ## 索引与内存次序
 *
 *  数据集形状为 **(nz, ny, nx)**（即轴序反转），这样 HDF5 的 C 序连续内存
 *  正好是 **i 最快**——与本项目网格、与 `vtkwriter` 的 `.vti` 完全一致。
 *  NumPy 里读出来 `d[k, j, i]` 即对应 `(x_i, y_j, z_k)`，无需转置。
 *  矢量数据集形状为 **(nz, ny, nx, 3)**，**最后一维为分量**（按点交错，
 *  与 `.vti` 和 Python 侧 `(npts, 3)` 的习惯一致）。
 */

#ifndef IBSIMU_CYCL_HDF5WRITER_HPP
#define IBSIMU_CYCL_HDF5WRITER_HPP 1

#include <string>
#include <vector>
#include "types.hpp"
#include "vec3d.hpp"
#include "vtkwriter.hpp"     // 复用 TrajectoryPoint


namespace ibsimu_cycl {


/*! \brief 是否编译进了 HDF5 后端（由 configure 探测决定）。 */
bool hdf5_available();

/*! \brief HDF5 库版本字符串；未编译进后端时返回空串。 */
std::string hdf5_version_string();


/*! \brief 把规则网格上的标量场写成 HDF5。
 *
 *  \param filename   输出文件名，建议以 `.h5` 结尾。
 *  \param mesh_name  网格组名（对应 `meshes/<mesh_name>`）。
 *  \param dims       各方向**节点数** (nx, ny, nz)。
 *  \param origo      第一个节点 (i=j=k=0) 的坐标（写进 `gridGlobalOffset`）。
 *  \param spacing    相邻节点间距（写进 `gridSpacing`）。
 *  \param scalar_name 标量数据集名。
 *  \param scalar     标量数据，长度必须为 nx*ny*nz，**按 i 最快**排列。
 *  \param scalar_unit_si 标量的 SI 单位因子（如 V/m 写 1.0，T 写 1.0，
 *                    keV 写 1.602176634e-16）；openPMD 的 `unitSI` 语义。
 *
 *  \throws std::invalid_argument 数据长度与 \a dims 不符。
 *  \throws std::runtime_error    未编译进 HDF5 后端，或 HDF5 调用失败。
 */
void hdf5_write_image_data( const std::string &filename,
			    const std::string &mesh_name,
			    const Int3D &dims, const Vec3D &origo,
			    const Vec3D &spacing,
			    const std::string &scalar_name,
			    const std::vector<double> &scalar,
			    double scalar_unit_si = 1.0 );


/*! \brief 同前，但额外写出一个**三分量矢量场**。
 *
 *  \param scalar_unit_si 标量数据集的 SI 单位因子（见上一个重载）。
 *  \param vector_name 矢量数据集名；\a vector 长度为 3*nx*ny*nz，
 *                     按点排列，每点依次为 (vx, vy, vz)。
 *  \param vector_unit_si 矢量各分量的 SI 单位因子。
 */
void hdf5_write_image_data( const std::string &filename,
			    const std::string &mesh_name,
			    const Int3D &dims, const Vec3D &origo,
			    const Vec3D &spacing,
			    const std::string &scalar_name,
			    const std::vector<double> &scalar,
			    double scalar_unit_si,
			    const std::string &vector_name,
			    const std::vector<double> &vector,
			    double vector_unit_si = 1.0 );


/*! \brief 把若干粒子轨迹写成 HDF5（扁平数组 + 前缀和偏移量）。
 *
 *  布局（\a group_name 默认建议 "tracks"）::
 *
 *      /data/0/particles/<group_name>/
 *          @positionsUnitSI = 1.0     ; m
 *          @timeUnitSI      = 1.0     ; s
 *          /position  (N, 3) float64  ; 各轨迹依次拼接
 *          /time      (N,)   float64
 *          /offset    (M+1,) int64    ; 第 m 条轨迹占 [offset[m], offset[m+1])
 *          /count     (M,)   int64
 *
 *  \throws std::runtime_error 未编译进 HDF5 后端，或 HDF5 调用失败。
 */
void hdf5_write_polylines( const std::string &filename,
			   const std::string &group_name,
			   const std::vector<std::vector<TrajectoryPoint>> &lines );


} // namespace ibsimu_cycl


#endif // IBSIMU_CYCL_HDF5WRITER_HPP
