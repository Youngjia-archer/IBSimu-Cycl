/** @file ringfield3d.hpp
 *  @brief IBSimu-Cycl: 由中平面 Bz 重建的三维回旋加速器磁场。
 *
 *  回旋加速器磁场图通常只给出**中平面**主磁场 Bz(r, theta)。在无电流真空区
 *  (nabla.B = 0, nabla x B = 0)，且中平面对称 (Bz 关于 z 偶、Br/Btheta 关于 z 奇)，
 *  可由 Bz|z=0 = b(r,theta) 作 Laplace 展开重建三维场：
 *
 *      T b      = d2b/dr2 + (1/r) db/dr + (1/r^2) d2b/dtheta2     (横向 Laplacian)
 *      Bz       = b - (1/2) Tb z^2
 *      Br       = db/dr z - (1/6) d(Tb)/dr z^3
 *      Btheta   = (1/r) db/dtheta z - (1/(6r)) d(Tb)/dtheta z^3
 *
 *  与 OPAL-cycl `Cyclotron::getField` 的中平面展开同源。
 *
 *  输入为 `CCylFieldMap3D`（nz=1，仅 Bz，满 360deg），本类在其节点上数值求导并缓存。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#ifndef IBSIMU_CYCL_RINGFIELD3D_HPP
#define IBSIMU_CYCL_RINGFIELD3D_HPP 1

#include <string>
#include <vector>

#include "vectorfield.hpp"
#include "vec3d.hpp"
#include "cylfmap3d.hpp"
#include "ringmap.hpp"

namespace ibsimu_cycl {

/*! \brief 由中平面 Bz 重建的三维回旋加速器磁场。 */
class CRingFieldMap3D : public VectorField {

    std::size_t _nr, _nt;   /*!< \brief r / theta 节点数。 */
    double      _r0, _dr;   /*!< \brief r 起点与步长。 */
    double      _dtheta;    /*!< \brief theta 步长 (rad)。 */

    /*! \brief 节点数据的分量索引（按节点交错存储，见 \a _data）。 */
    enum { C_B = 0, C_DBR, C_DBTH, C_TB, C_DTRB, C_DTTB, C_NCOMP };

    /*! \brief 交错存储的节点数据：_data[(i*_nt + k)*C_NCOMP + c]。
     *
     *  b / dbr / dbth / tb / dtrb / dttb 六个量按**节点交错**存放。一次双线性
     *  插值只触碰 2 个节点邻域(约 4 条缓存行)；若按分量各存一个大数组，同样
     *  一次插值要跨越 6 个独立的 1.2 MB 区间(约 12 条缓存行)，在真实场图上
     *  会明显受内存延迟限制。
     */
    std::vector<double> _data;

    std::size_t node( std::size_t i, std::size_t k ) const
    { return( (i*_nt + k)*(std::size_t)C_NCOMP ); }

    double  val( std::size_t i, std::size_t k, int c ) const
    { return( _data[node(i,k) + (std::size_t)c] ); }
    double &ref( std::size_t i, std::size_t k, int c )
    { return( _data[node(i,k) + (std::size_t)c] ); }

    std::size_t idx( std::size_t i, std::size_t k ) const { return( i*_nt + k ); }
    double r_at( std::size_t i ) const { return( _r0 + _dr*(double)i ); }

    void compute_derivatives();

public:

    /*! \brief 由中平面柱坐标场图构造（取其 k=0 切片，忽略 Br/Btheta）。 */
    explicit CRingFieldMap3D( const CCylFieldMap3D &midplane );

    virtual ~CRingFieldMap3D();

    std::size_t size_r() const { return( _nr ); }
    std::size_t size_theta() const { return( _nt ); }
    double r0() const { return( _r0 ); }
    double dr() const { return( _dr ); }
    double dtheta() const { return( _dtheta ); }

    /*! \brief 中平面网格点的 Bz 值。 */
    double bz_midplane( std::size_t i, std::size_t k ) const { return( val(i,k,C_B) ); }

    /*! \brief 场求值（输入输出均为直角坐标，单位 T）。 */
    virtual const Vec3D operator()( const Vec3D &x ) const;
};

/*! \brief 便捷函数：直接读取 OPAL/PSI RING 磁场图并重建三维场。
 *
 *  \param filename  bfield.dat 路径
 *  \param b_scale   kG -> T 缩放（0.1）
 *  \param nsymmetry 旋转对称重数（PSI Ring = 8）
 *  \param info      可选输出头部元数据
 */
CRingFieldMap3D read_ring_field_map3d( const std::string &filename,
                                       double b_scale = 0.1,
                                       int nsymmetry = 8,
                                       RingFieldMapInfo *info = 0 );

} // namespace ibsimu_cycl

#endif // IBSIMU_CYCL_RINGFIELD3D_HPP
