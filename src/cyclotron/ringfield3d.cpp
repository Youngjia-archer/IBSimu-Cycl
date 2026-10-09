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

    _b.resize( _nr*_nt );
    for( std::size_t i = 0; i < _nr; ++i )
        for( std::size_t k = 0; k < _nt; ++k )
            _b[idx(i,k)] = midplane.node_value( i, k, 0 )[2];

    compute_derivatives();
}

CRingFieldMap3D::~CRingFieldMap3D()
{
}

void CRingFieldMap3D::compute_derivatives()
{
    const std::size_t nr = _nr, nt = _nt;
    _dbr.assign( nr*nt, 0.0 );
    _dbth.assign( nr*nt, 0.0 );
    _tb.assign( nr*nt, 0.0 );
    _dtrb.assign( nr*nt, 0.0 );
    _dttb.assign( nr*nt, 0.0 );

    const long ntl = (long)nt;
    auto wrap = [ntl]( long k ) -> std::size_t { return (std::size_t)(((k % ntl) + ntl) % ntl); };

    // --- 一阶导数 ---
    for( std::size_t i = 0; i < nr; ++i ) {
        for( std::size_t k = 0; k < nt; ++k ) {
            // d/dtheta (周期中心差分)
            double bp = _b[idx(i, wrap((long)k+1))];
            double bm = _b[idx(i, wrap((long)k-1))];
            _dbth[idx(i,k)] = (bp - bm)/(2.0*_dtheta);

            // d/dr (内部中心差分, 边界单侧)
            if( i == 0 )
                _dbr[idx(i,k)] = (_b[idx(1,k)] - _b[idx(0,k)])/_dr;
            else if( i == nr-1 )
                _dbr[idx(i,k)] = (_b[idx(nr-1,k)] - _b[idx(nr-2,k)])/_dr;
            else
                _dbr[idx(i,k)] = (_b[idx(i+1,k)] - _b[idx(i-1,k)])/(2.0*_dr);
        }
    }

    // --- 二阶导数与横向 Laplacian T b = b_rr + (1/r) b_r + (1/r^2) b_thth ---
    for( std::size_t i = 0; i < nr; ++i ) {
        double r = r_at( i );
        for( std::size_t k = 0; k < nt; ++k ) {
            double bthth = ( _b[idx(i, wrap((long)k+1))]
                           - 2.0*_b[idx(i,k)]
                           + _b[idx(i, wrap((long)k-1))] )/(_dtheta*_dtheta);

            double brr;
            if( i == 0 )
                brr = ( _b[idx(2,k)] - 2.0*_b[idx(1,k)] + _b[idx(0,k)] )/(_dr*_dr);
            else if( i == nr-1 )
                brr = ( _b[idx(nr-1,k)] - 2.0*_b[idx(nr-2,k)] + _b[idx(nr-3,k)] )/(_dr*_dr);
            else
                brr = ( _b[idx(i+1,k)] - 2.0*_b[idx(i,k)] + _b[idx(i-1,k)] )/(_dr*_dr);

            _tb[idx(i,k)] = brr + _dbr[idx(i,k)]/r + bthth/(r*r);
        }
    }

    // --- Tb 的一阶导数 ---
    for( std::size_t i = 0; i < nr; ++i ) {
        for( std::size_t k = 0; k < nt; ++k ) {
            _dttb[idx(i,k)] = ( _tb[idx(i, wrap((long)k+1))]
                              - _tb[idx(i, wrap((long)k-1))] )/(2.0*_dtheta);

            if( i == 0 )
                _dtrb[idx(i,k)] = ( _tb[idx(1,k)] - _tb[idx(0,k)] )/_dr;
            else if( i == nr-1 )
                _dtrb[idx(i,k)] = ( _tb[idx(nr-1,k)] - _tb[idx(nr-2,k)] )/_dr;
            else
                _dtrb[idx(i,k)] = ( _tb[idx(i+1,k)] - _tb[idx(i-1,k)] )/(2.0*_dr);
        }
    }
}

double CRingFieldMap3D::interp( double r, double theta, const std::vector<double> &f ) const
{
    // r 方向: clamp
    double fr = (r - _r0)/_dr;
    int i0, i1; double tr;
    if( _nr <= 1 ) { i0 = i1 = 0; tr = 0.0; }
    else {
        i0 = (int)std::floor( fr ); tr = fr - i0;
        if( i0 < 0 ) { i0 = i1 = 0; tr = 0.0; }
        else if( (std::size_t)i0 >= _nr-1 ) { i0 = i1 = (int)_nr-1; tr = 0.0; }
        else i1 = i0 + 1;
    }

    // theta 方向: 周期
    long ntl = (long)_nt;
    double ft = theta/_dtheta;
    int j = (int)std::floor( ft );
    double tt = ft - j;
    if( tt < 0.0 ) tt = 0.0;
    std::size_t j0 = (std::size_t)((((long)j % ntl) + ntl) % ntl);
    std::size_t j1 = (std::size_t)(((long)j0 + 1) % ntl);

    double f00 = f[idx((std::size_t)i0,j0)], f01 = f[idx((std::size_t)i0,j1)];
    double f10 = f[idx((std::size_t)i1,j0)], f11 = f[idx((std::size_t)i1,j1)];

    return( (1.0-tr)*((1.0-tt)*f00 + tt*f01) + tr*((1.0-tt)*f10 + tt*f11) );
}

const Vec3D CRingFieldMap3D::operator()( const Vec3D &x ) const
{
    double r     = std::sqrt( x[0]*x[0] + x[1]*x[1] );
    double theta = std::atan2( x[1], x[0] );
    if( theta < 0.0 ) theta += 2.0*M_PI;
    double z = x[2];

    double b    = interp( r, theta, _b );
    double dbr  = interp( r, theta, _dbr );
    double dbth = interp( r, theta, _dbth );
    double tb   = interp( r, theta, _tb );
    double dtrb = interp( r, theta, _dtrb );
    double dttb = interp( r, theta, _dttb );

    double z2 = z*z;
    double z3 = z2*z;

    double Bz = b - 0.5*tb*z2;
    double Br = dbr*z - (1.0/6.0)*dtrb*z3;
    double Bt = (dbth/r)*z - (1.0/(6.0*r))*dttb*z3;

    double cs = std::cos( theta ), sn = std::sin( theta );
    return( Vec3D( Br*cs - Bt*sn, Br*sn + Bt*cs, Bz ) );
}

CRingFieldMap3D read_ring_field_map3d( const std::string &filename, double b_scale,
                                      int nsymmetry, RingFieldMapInfo *info )
{
    CCylFieldMap3D mid = read_ring_field_map( filename, b_scale, nsymmetry, info );
    return( CRingFieldMap3D( mid ) );
}

} // namespace ibsimu_cycl
