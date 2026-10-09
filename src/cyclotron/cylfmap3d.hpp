/** @file cylfmap3d.hpp
 *  @brief IBSimu-Cycl: 柱坐标 (r, theta, z) 三维矢量场图。
 *
 *  上游 MeshVectorField 只支持"均匀立方网格"的笛卡尔场或轴对称(x,r)场，
 *  无法表示随方位角 theta 变化的三维场——而回旋加速器的扇形磁铁/螺旋扇
 *  恰恰依赖方位角变化。本类填补该空缺。
 *
 *  数据存于规则网格：
 *      r[i]     = r0 + i*dr          (i = 0..nr-1)
 *      theta[j] = j*dtheta           (j = 0..nt-1, dtheta = 2*pi/nt, 周期)
 *      z[k]     = z0 + k*dz          (k = 0..nz-1)
 *  每个节点存 (Br, Btheta, Bz)。求值时输入/输出均为直角坐标：
 *      x=(x,y,z) -> (r,theta,z) -> 三线性插值 -> (Br,Btheta,Bz) -> (Bx,By,Bz)
 *
 *  r/z 方向越界按最近边界(clamp)；theta 方向周期处理。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#ifndef IBSIMU_CYCL_CYLFMAP3D_HPP
#define IBSIMU_CYCL_CYLFMAP3D_HPP 1

#include <cstddef>
#include <string>
#include <vector>
#include <iostream>

#include "vectorfield.hpp"
#include "vec3d.hpp"

namespace ibsimu_cycl {

/*! \brief 柱坐标三维矢量场图（支持方位角变化）。 */
class CCylFieldMap3D : public VectorField {

    std::size_t _nr, _nt, _nz;   /*!< \brief 各方向节点数。 */
    double      _r0, _dr;        /*!< \brief r 网格起点与步长。 */
    double      _z0, _dz;        /*!< \brief z 网格起点与步长。 */
    double      _dtheta;         /*!< \brief 方位角步长 = 2*pi/nt。 */

    std::vector<double> _br;     /*!< \brief Br   分量, 索引 (i,j,k)。 */
    std::vector<double> _bt;     /*!< \brief Btheta 分量, 索引 (i,j,k)。 */
    std::vector<double> _bz;     /*!< \brief Bz   分量, 索引 (i,j,k)。 */

    std::size_t idx( std::size_t i, std::size_t j, std::size_t k ) const
    { return( (i*_nt + j)*_nz + k ); }

    double component( int c, int i, int j, int k ) const
    {
        const std::vector<double> &v = (c == 0 ? _br : (c == 1 ? _bt : _bz));
        return( v[idx((std::size_t)i,(std::size_t)j,(std::size_t)k)] );
    }

public:

    /*! \brief 构造空场（全部置零）。
     *
     *  \param nr, ntheta, nz  各方向节点数（均 >= 1）
     *  \param r0, dr          r 网格起点与步长
     *  \param z0, dz          z 网格起点与步长
     *  方位角范围固定为 [0, 2*pi)，步长 dtheta = 2*pi/ntheta。
     */
    CCylFieldMap3D( std::size_t nr, std::size_t ntheta, std::size_t nz,
                    double r0, double dr, double z0, double dz );

    virtual ~CCylFieldMap3D();

    std::size_t size_r() const { return( _nr ); }
    std::size_t size_theta() const { return( _nt ); }
    std::size_t size_z() const { return( _nz ); }
    double dtheta() const { return( _dtheta ); }

    /*! \brief 设置节点 (i,j,k) 的 (Br,Btheta,Bz)。 */
    void set_value( std::size_t i, std::size_t j, std::size_t k,
                    double br, double btheta, double bz );

    /*! \brief 读取节点 (i,j,k) 的 (Br,Btheta,Bz)。 */
    Vec3D node_value( std::size_t i, std::size_t j, std::size_t k ) const;

    /*! \brief 场求值：输入输出均为直角坐标。 */
    virtual const Vec3D operator()( const Vec3D &x ) const;

    /*! \brief 保存为本项目原生 ASCII 格式。 */
    void save_ascii( std::ostream &os ) const;

    /*! \brief 从本项目原生 ASCII 格式载入。 */
    void load_ascii( std::istream &is );

    /*! \brief 保存到文件。 */
    void save( const std::string &filename ) const;

    /*! \brief 从文件载入。 */
    void load( const std::string &filename );
};

} // namespace ibsimu_cycl

#endif // IBSIMU_CYCL_CYLFMAP3D_HPP
