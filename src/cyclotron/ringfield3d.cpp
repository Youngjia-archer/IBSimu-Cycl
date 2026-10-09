/** @file ringfield3d.cpp
 *  @brief IBSimu-Cycl: 由中平面 Bz 重建的三维回旋加速器磁场（实现）。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <cmath>
#include <stdexcept>

#include "ringfield3d.hpp"

namespace ibsimu_cycl {

CRingFieldMap3D::CRingFieldMap3D( const CCylFieldMap3D &midplane )
    : _nr( midplane.size_r() ), _nt( midplane.size_theta() ),
      _r0( midplane.r0() ), _dr( midplane.dr() ), _dtheta( midplane.dtheta() )
{
    if( _nr < 3 || _nt < 3 )
        throw( std::invalid_argument("CRingFieldMap3D: need at least 3x3 grid") );
    if( _r0 <= 0.0 )
        throw( std::invalid_argument("CRingFieldMap3D: r0 must be > 0") );

    // 按节点交错分配：每个节点 C_NCOMP 个 double
    _data.assign( _nr*_nt*(std::size_t)C_NCOMP, 0.0 );
    for( std::size_t i = 0; i < _nr; ++i )
        for( std::size_t k = 0; k < _nt; ++k )
            ref( i, k, C_B ) = midplane.node_value( i, k, 0 )[2];

    compute_derivatives();
}

CRingFieldMap3D::~CRingFieldMap3D()
{
}

void CRingFieldMap3D::compute_derivatives()
{
    const std::size_t nr = _nr, nt = _nt;

    const long ntl = (long)nt;
    auto wrap = [ntl]( long k ) -> std::size_t { return (std::size_t)(((k % ntl) + ntl) % ntl); };

    // --- 一阶导数 ---
    for( std::size_t i = 0; i < nr; ++i ) {
        for( std::size_t k = 0; k < nt; ++k ) {
            // d/dtheta (周期中心差分)
            double bp = val( i, wrap((long)k+1), C_B );
            double bm = val( i, wrap((long)k-1), C_B );
            ref( i, k, C_DBTH ) = (bp - bm)/(2.0*_dtheta);

            // d/dr (内部中心差分, 边界单侧)
            if( i == 0 )
                ref( i, k, C_DBR ) = ( val(1,k,C_B) - val(0,k,C_B) )/_dr;
            else if( i == nr-1 )
                ref( i, k, C_DBR ) = ( val(nr-1,k,C_B) - val(nr-2,k,C_B) )/_dr;
            else
                ref( i, k, C_DBR ) = ( val(i+1,k,C_B) - val(i-1,k,C_B) )/(2.0*_dr);
        }
    }

    // --- 二阶导数与横向 Laplacian T b = b_rr + (1/r) b_r + (1/r^2) b_thth ---
    for( std::size_t i = 0; i < nr; ++i ) {
        double r = r_at( i );
        for( std::size_t k = 0; k < nt; ++k ) {
            double bthth = ( val(i,wrap((long)k+1),C_B)
                           - 2.0*val(i,k,C_B)
                           + val(i,wrap((long)k-1),C_B) )/(_dtheta*_dtheta);

            double brr;
            if( i == 0 )
                brr = ( val(2,k,C_B) - 2.0*val(1,k,C_B) + val(0,k,C_B) )/(_dr*_dr);
            else if( i == nr-1 )
                brr = ( val(nr-1,k,C_B) - 2.0*val(nr-2,k,C_B) + val(nr-3,k,C_B) )/(_dr*_dr);
            else
                brr = ( val(i+1,k,C_B) - 2.0*val(i,k,C_B) + val(i-1,k,C_B) )/(_dr*_dr);

            ref( i, k, C_TB ) = brr + val(i,k,C_DBR)/r + bthth/(r*r);
        }
    }

    // --- Tb 的一阶导数 ---
    for( std::size_t i = 0; i < nr; ++i ) {
        for( std::size_t k = 0; k < nt; ++k ) {
            ref( i, k, C_DTTB ) = ( val(i,wrap((long)k+1),C_TB)
                                  - val(i,wrap((long)k-1),C_TB) )/(2.0*_dtheta);

            if( i == 0 )
                ref( i, k, C_DTRB ) = ( val(1,k,C_TB) - val(0,k,C_TB) )/_dr;
            else if( i == nr-1 )
                ref( i, k, C_DTRB ) = ( val(nr-1,k,C_TB) - val(nr-2,k,C_TB) )/_dr;
            else
                ref( i, k, C_DTRB ) = ( val(i+1,k,C_TB) - val(i-1,k,C_TB) )/(2.0*_dr);
        }
    }
}

const Vec3D CRingFieldMap3D::operator()( const Vec3D &x ) const
{
    const double x0 = x[0], y0 = x[1], z = x[2];

    double r = std::sqrt( x0*x0 + y0*y0 );
    if( r == 0.0 )  // 轴上：Br/Btheta 项含 1/r，由轴对称极限取 0
        return( Vec3D( 0.0, 0.0, val(0,0,C_B) ) );

    double theta = std::atan2( y0, x0 );
    if( theta < 0.0 ) theta += 2.0*M_PI;

    // ---- 单元与权重只计算一次，6 个分量共用（原先每个分量各算一遍）----
    const int nrl = (int)_nr;
    const int ntl = (int)_nt;

    double fr = (r - _r0)/_dr;
    int i0;
    double tr;
    if( fr <= 0.0 ) {            // 低于下边界：clamp 到 i=0
        i0 = 0;
        tr = 0.0;
    } else {
        i0 = (int)fr;            // fr > 0，截断等价于 floor
        if( i0 >= nrl-1 ) {      // 高于上边界：clamp 到最后一个节点
            i0 = nrl-1;
            tr = 0.0;
        } else {
            tr = fr - (double)i0;
        }
    }
    const int i1 = ( i0 < nrl-1 ) ? i0+1 : i0;

    double ft = theta/_dtheta;   // theta >= 0，截断等价于 floor
    int j = (int)ft;
    double tt = ft - (double)j;
    if( tt < 0.0 ) tt = 0.0;
    j %= ntl;
    if( j < 0 ) j += ntl;
    const int j1 = ( j+1 == ntl ) ? 0 : j+1;   // theta 方向周期

    const double w00 = (1.0-tr)*(1.0-tt);
    const double w01 = (1.0-tr)*tt;
    const double w10 = tr*(1.0-tt);
    const double w11 = tr*tt;

    const double *p00 = &_data[node((std::size_t)i0,(std::size_t)j )];
    const double *p01 = &_data[node((std::size_t)i0,(std::size_t)j1)];
    const double *p10 = &_data[node((std::size_t)i1,(std::size_t)j )];
    const double *p11 = &_data[node((std::size_t)i1,(std::size_t)j1)];

    double f[C_NCOMP];
    for( int c = 0; c < (int)C_NCOMP; ++c )
        f[c] = w00*p00[c] + w01*p01[c] + w10*p10[c] + w11*p11[c];

    const double b    = f[C_B];
    const double dbr  = f[C_DBR];
    const double dbth = f[C_DBTH];
    const double tb   = f[C_TB];
    const double dtrb = f[C_DTRB];
    const double dttb = f[C_DTTB];

    const double z2 = z*z;
    const double z3 = z2*z;

    const double Bz = b - 0.5*tb*z2;
    const double Br = dbr*z - (1.0/6.0)*dtrb*z3;
    const double Bt = (dbth/r)*z - (1.0/(6.0*r))*dttb*z3;

    // 柱 -> 直角：cos/sin(theta) 直接由坐标/半径得到，省去两次超越函数
    const double cs = x0/r;
    const double sn = y0/r;
    return( Vec3D( Br*cs - Bt*sn, Br*sn + Bt*cs, Bz ) );
}

CRingFieldMap3D read_ring_field_map3d( const std::string &filename, double b_scale,
                                      int nsymmetry, RingFieldMapInfo *info )
{
    CCylFieldMap3D mid = read_ring_field_map( filename, b_scale, nsymmetry, info );
    return( CRingFieldMap3D( mid ) );
}

} // namespace ibsimu_cycl
