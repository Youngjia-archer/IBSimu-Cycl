/** @file ringmap.hpp
 *  @brief IBSimu-Cycl: OPAL/PSI "RING" 回旋加速器磁场图读取器。
 *
 *  读取 OPAL-cycl 使用的 CERN "FIELD" 格式磁场图（如 PSI Ring 的 bfield.dat）。
 *  该文件以极坐标存储(r, theta)中平面主磁场 Bz 及其对 theta 的各阶导数，
 *  并利用 n 折旋转对称只保存 1/n 扇形区。
 *
 *  本读取器将其转换为满 360deg 的 CCylFieldMap3D（nz=1，仅中平面 Bz），
 *  以便直接复用 IBSimu-Cycl 的柱坐标插值。离面(off-midplane)分量 Br/Btheta
 *  的重建（OPAL 用中平面展开 + 导数实现）留待后续。
 *
 *  参考实现：OPAL `Cyclotron::getFieldFromFile_Ring`
 *  (Classic/AbsBeamline/Cyclotron.cpp)，GPL。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#ifndef IBSIMU_CYCL_RINGMAP_HPP
#define IBSIMU_CYCL_RINGMAP_HPP 1

#include <string>

#include "cylfmap3d.hpp"

namespace ibsimu_cycl {

/*! \brief OPAL/PSI "RING" 磁场图的头部元数据。 */
struct RingFieldMapInfo {
    double rmin_mm;    /*!< \brief 最小半径 [mm]。 */
    double dr_mm;      /*!< \brief 径向步长 [mm]（负值取倒数）。 */
    double tetmin_deg; /*!< \brief 最小方位角 [deg]。 */
    double dtet_deg;   /*!< \brief 方位角步长 [deg]（负值取倒数）。 */
    int    nrad;       /*!< \brief 径向网格点数。 */
    int    ntet;       /*!< \brief 存储扇形区内的方位角网格点数。 */
    int    lpar;       /*!< \brief 第二数据块参数个数（文件内信息，忽略）。 */
};

/*! \brief 读取 OPAL/PSI "RING" 磁场图并转换为满 360deg 的柱坐标场图。
 *
 *  \param filename  bfield.dat 路径
 *  \param b_scale   场值单位缩放（kG -> T 为 0.1，OPAL 约定）
 *  \param nsymmetry 旋转对称重数（PSI Ring 为 8）
 *  \param info      可选，输出头部元数据
 *  \return 满 360deg 柱坐标场图（nz = 1，仅设置 Bz）
 *
 *  \throws std::runtime_error 文件无法打开、格式不符或对称性与网格不一致
 */
CCylFieldMap3D read_ring_field_map( const std::string &filename,
                                    double b_scale = 0.1,
                                    int nsymmetry = 8,
                                    RingFieldMapInfo *info = 0 );

} // namespace ibsimu_cycl

#endif // IBSIMU_CYCL_RINGMAP_HPP
