/** @file vtkwriter.hpp
 *  @brief VTK XML 写出器：规则网格标量/矢量场（.vti）与粒子轨迹折线（.vtp）。
 *
 *  设计目标：**零外部依赖**地生成 ParaView / PyVista / VisIt 可直接打开的文件，
 *  作为 IBSimu-Cycl 的三维可视化出口（ROADMAP §4F / R3）。
 *
 *  为什么是 VTK XML 而不是先上 HDF5：
 *    - `.vti`/`.vtp` 是纯文本 XML，不需要链接任何库（HDF5/VTK 均未进构建系统）；
 *    - ParaView 与 PyVista 原生支持，无需转换脚本；
 *    - 文件可 diff、可被测试脚本逐字节校验，适合放进 `make check`。
 *  大规模数据（>1e8 点）的二进制/追加段后端（HDF5/openPMD、`.pvti`）见 ROADMAP R3-③。
 *
 *  数据排列约定：与本项目/IBSimu 网格一致 —— **i 最快、k 最慢**，
 *  即索引 `a = i + nx*(j + ny*k)`。这与 VTK ImageData 的期望顺序相同，无需转置。
 */

#ifndef IBSIMU_CYCL_VTKWRITER_HPP
#define IBSIMU_CYCL_VTKWRITER_HPP 1

#include <string>
#include <vector>
#include "types.hpp"
#include "vec3d.hpp"


namespace ibsimu_cycl {


/*! \brief 轨迹上的一个采样点：物理时刻 + 位置。 */
struct TrajectoryPoint {
    double t;
    Vec3D  x;
};


/*! \brief 把规则网格上的标量场写成 VTK ImageData（.vti）。
 *
 *  \param filename 输出文件名，建议以 `.vti` 结尾。
 *  \param dims     各方向**节点数** (nx, ny, nz)。
 *  \param origo    第一个节点 (i=j=k=0) 的坐标。
 *  \param spacing  相邻节点间距 (hx, hy, hz)；IBSimu 网格为各向同性（三者相等）。
 *  \param scalar_name 标量数组名（ParaView 中显示为变量名）。
 *  \param scalar   标量数据，长度必须为 nx*ny*nz，按 i 最快排列。
 *
 *  \throws std::invalid_argument 数据长度与 \a dims 不符。
 *  \throws std::runtime_error    文件无法写入。
 */
void vtk_write_image_data( const std::string &filename,
			   const Int3D &dims, const Vec3D &origo, const Vec3D &spacing,
			   const std::string &scalar_name,
			   const std::vector<double> &scalar );


/*! \brief 同前，但额外写出一个**三分量矢量场**（如 B、E）。
 *
 *  \param vector_name 矢量数组名；\a vector 长度为 3*nx*ny*nz，
 *                     按点排列，每点依次为 (vx, vy, vz)。
 */
void vtk_write_image_data( const std::string &filename,
			   const Int3D &dims, const Vec3D &origo, const Vec3D &spacing,
			   const std::string &scalar_name,
			   const std::vector<double> &scalar,
			   const std::string &vector_name,
			   const std::vector<double> &vector );


/*! \brief 把若干粒子轨迹写成 VTK PolyData（.vtp），每条轨迹一根折线。
 *
 *  每个点附带标量数组 `t`（物理时刻），便于在 ParaView 中按时间着色。
 *  长度 < 2 的轨迹只会贡献点、不产生折线。
 */
void vtk_write_polylines( const std::string &filename,
			  const std::vector<std::vector<TrajectoryPoint>> &lines );


} // namespace ibsimu_cycl


#endif // IBSIMU_CYCL_VTKWRITER_HPP
