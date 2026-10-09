/** @file fieldmap3d.hpp
 *  @brief IBSimu-Cycl: 笛卡尔三维矢量场图（每轴独立步长）。
 *
 *  上游 MeshVectorField 只支持**各向同性**步长 h 的均匀立方网格。本类支持
 *  三轴**独立**起点与步长，采用三线性插值，越界按最近边界 clamp。
 *  这是外部 FEM 求解器（如 Elmer）输出转换为本项目可读格式的目标容器。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#ifndef IBSIMU_CYCL_FIELDMAP3D_HPP
#define IBSIMU_CYCL_FIELDMAP3D_HPP 1

#include <cstddef>
#include <string>
#include <vector>
#include <iostream>

#include "vectorfield.hpp"
#include "vec3d.hpp"

namespace ibsimu_cycl {

/*! \brief 笛卡尔三维矢量场图（规则网格，每轴独立步长）。 */
class CFieldMap3D : public VectorField {

    std::size_t _nx, _ny, _nz;      /*!< \brief 各轴节点数。 */
    double      _x0, _dx;           /*!< \brief x 起点与步长。 */
    double      _y0, _dy;           /*!< \brief y 起点与步长。 */
    double      _z0, _dz;           /*!< \brief z 起点与步长。 */

    std::vector<double> _fx, _fy, _fz;

    std::size_t idx( std::size_t i, std::size_t j, std::size_t k ) const
    { return( (i*_ny + j)*_nz + k ); }

    double component( int c, int i, int j, int k ) const
    {
        const std::vector<double> &v = (c == 0 ? _fx : (c == 1 ? _fy : _fz));
        return( v[idx((std::size_t)i,(std::size_t)j,(std::size_t)k)] );
    }

public:

    /*! \brief 构造空场（全部置零）。各方向节点数需 >= 1。 */
    CFieldMap3D( std::size_t nx, std::size_t ny, std::size_t nz,
                 double x0, double dx,
                 double y0, double dy,
                 double z0, double dz );

    virtual ~CFieldMap3D();

    std::size_t size_x() const { return( _nx ); }
    std::size_t size_y() const { return( _ny ); }
    std::size_t size_z() const { return( _nz ); }
    double x0() const { return( _x0 ); }
    double dx() const { return( _dx ); }
    double y0() const { return( _y0 ); }
    double dy() const { return( _dy ); }
    double z0() const { return( _z0 ); }
    double dz() const { return( _dz ); }

    /*! \brief 设置节点 (i,j,k) 的矢量值。 */
    void set_value( std::size_t i, std::size_t j, std::size_t k,
                    double fx, double fy, double fz );

    /*! \brief 读取节点 (i,j,k) 的矢量值。 */
    Vec3D node_value( std::size_t i, std::size_t j, std::size_t k ) const;

    /*! \brief 场求值（三线性插值，越界 clamp）。 */
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

#endif // IBSIMU_CYCL_FIELDMAP3D_HPP
