/** @file cylfmap3d.cpp
 *  @brief IBSimu-Cycl: 柱坐标 (r, theta, z) 三维矢量场图（实现）。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

#include "cylfmap3d.hpp"

namespace ibsimu_cycl {

static const char *CYLMAP_MAGIC = "IBSIMU-CYCL-CYLFIELDMAP";
static const int   CYLMAP_VERSION = 1;

/*! \brief 计算单调坐标 x 所在单元及插值权重（越界 clamp 到边界）。
 *
 *  返回 i0/i1 为相邻节点索引，t 为插值权重 (0..1)。
 *  当 n<=1 时退化为常数。
 */
static void cell_weights( double x, double x0, double dx, std::size_t n,
                          int &i0, int &i1, double &t )
{
    if( n <= 1 ) {
        i0 = i1 = 0;
        t = 0.0;
        return;
    }

    double f = (x - x0)/dx;
    int    i = (int)std::floor( f );
    t = f - i;

    if( i < 0 ) {                 // 低于下边界
        i0 = i1 = 0;
        t = 0.0;
        return;
    }
    if( (std::size_t)i >= n-1 ) { // 高于上边界
        i0 = i1 = (int)n-1;
        t = 0.0;
        return;
    }
    i0 = i;
    i1 = i + 1;
}

/*! \brief 计算周期方位角 theta 所在单元及权重。 */
static void theta_weights( double theta, double dtheta, std::size_t n,
                           int &j0, int &j1, double &t )
{
    if( n <= 1 ) {
        j0 = j1 = 0;
        t = 0.0;
        return;
    }

    double f = theta/dtheta;
    if( f < 0.0 )
        f += (double)n;

    int j = (int)std::floor( f );
    t = f - j;
    if( t < 0.0 ) t = 0.0;

    j0 = ((j % (int)n) + (int)n) % (int)n;
    j1 = (j0 + 1) % (int)n;
}


CCylFieldMap3D::CCylFieldMap3D( std::size_t nr, std::size_t ntheta, std::size_t nz,
                                double r0, double dr, double z0, double dz )
    : _nr(nr), _nt(ntheta), _nz(nz), _r0(r0), _dr(dr), _z0(z0), _dz(dz)
{
    if( _nr < 1 || _nt < 1 || _nz < 1 )
        throw( std::invalid_argument("CCylFieldMap3D: grid size must be >= 1 in each direction") );

    _dtheta = 2.0*M_PI/(double)_nt;

    _data.assign( 3*_nr*_nt*_nz, 0.0 );
}

CCylFieldMap3D::~CCylFieldMap3D()
{
}

void CCylFieldMap3D::set_value( std::size_t i, std::size_t j, std::size_t k,
                                double br, double btheta, double bz )
{
    if( i >= _nr || j >= _nt || k >= _nz )
        throw( std::out_of_range("CCylFieldMap3D::set_value: index out of range") );

    std::size_t n = 3*idx( i, j, k );
    _data[n]   = br;
    _data[n+1] = btheta;
    _data[n+2] = bz;
}

Vec3D CCylFieldMap3D::node_value( std::size_t i, std::size_t j, std::size_t k ) const
{
    if( i >= _nr || j >= _nt || k >= _nz )
        throw( std::out_of_range("CCylFieldMap3D::node_value: index out of range") );

    std::size_t n = 3*idx( i, j, k );
    return( Vec3D( _data[n], _data[n+1], _data[n+2] ) );
}

const Vec3D CCylFieldMap3D::operator()( const Vec3D &x ) const
{
    const double x0 = x[0], y0 = x[1], z = x[2];

    // 直角 -> 柱坐标
    double r     = std::sqrt( x0*x0 + y0*y0 );
    double theta = std::atan2( y0, x0 );
    if( theta < 0.0 )
        theta += 2.0*M_PI;

    int ir0, ir1, iz0, iz1, jt0, jt1;
    double tr, tz, tt;
    cell_weights ( r, _r0, _dr, _nr, ir0, ir1, tr );
    cell_weights ( z, _z0, _dz, _nz, iz0, iz1, tz );
    theta_weights( theta, _dtheta, _nt, jt0, jt1, tt );

    const int I[2] = { ir0, ir1 };
    const int J[2] = { jt0, jt1 };
    const int K[2] = { iz0, iz1 };
    const double U[2] = { 1.0-tr, tr };
    const double W[2] = { 1.0-tt, tt };
    const double V[2] = { 1.0-tz, tz };

    // 直接访问交错缓冲：8 个节点 × 3 个分量，无分量选择分支
    const double *dat = _data.data();
    double br = 0.0, bt = 0.0, bz = 0.0;
    for( int a = 0; a < 2; ++a )
        for( int b = 0; b < 2; ++b )
            for( int c = 0; c < 2; ++c ) {
                double w = U[a]*W[b]*V[c];
                const double *p = dat + 3*idx( (std::size_t)I[a], (std::size_t)J[b], (std::size_t)K[c] );
                br += w*p[0];
                bt += w*p[1];
                bz += w*p[2];
            }

    // 柱坐标分量 -> 直角坐标分量：cos/sin 由坐标/半径直接得到（免超越函数）
    double cs, sn;
    if( r > 0.0 ) { cs = x0/r; sn = y0/r; }
    else          { cs = 1.0;  sn = 0.0; }   // r=0: atan2(0,0)=0 -> cos=1, sin=0

    return( Vec3D( br*cs - bt*sn, br*sn + bt*cs, bz ) );
}

void CCylFieldMap3D::save_ascii( std::ostream &os ) const
{
    os << CYLMAP_MAGIC << " " << CYLMAP_VERSION << "\n";
    os << std::setprecision(17);
    os << "# nr ntheta nz\n";
    os << _nr << " " << _nt << " " << _nz << "\n";
    os << "# r0 dr z0 dz\n";
    os << _r0 << " " << _dr << " " << _z0 << " " << _dz << "\n";
    os << "# grid in order i(r), j(theta), k(z): Br Btheta Bz\n";
    for( std::size_t i = 0; i < _nr; ++i )
        for( std::size_t j = 0; j < _nt; ++j )
            for( std::size_t k = 0; k < _nz; ++k ) {
                std::size_t n = 3*idx( i, j, k );
                os << _data[n] << " " << _data[n+1] << " " << _data[n+2] << "\n";
            }
}

void CCylFieldMap3D::load_ascii( std::istream &is )
{
    std::string magic;
    int version = 0;
    is >> magic >> version;
    if( magic != CYLMAP_MAGIC || version != CYLMAP_VERSION )
        throw( std::runtime_error("CCylFieldMap3D::load_ascii: bad magic/version") );

    // 跳过 "#" 注释行
    auto skip_comments = [&is]() {
        while( is >> std::ws, is.peek() == '#' ) {
            std::string line;
            std::getline( is, line );
        }
    };

    skip_comments();
    std::size_t nr, nt, nz;
    is >> nr >> nt >> nz;
    skip_comments();
    double r0, dr, z0, dz;
    is >> r0 >> dr >> z0 >> dz;
    if( !is )
        throw( std::runtime_error("CCylFieldMap3D::load_ascii: failed to read header") );

    *this = CCylFieldMap3D( nr, nt, nz, r0, dr, z0, dz );

    skip_comments();
    for( std::size_t i = 0; i < nr; ++i )
        for( std::size_t j = 0; j < nt; ++j )
            for( std::size_t k = 0; k < nz; ++k ) {
                double br, bt, bz;
                is >> br >> bt >> bz;
                set_value( i, j, k, br, bt, bz );
            }
    if( !is )
        throw( std::runtime_error("CCylFieldMap3D::load_ascii: truncated data") );
}

void CCylFieldMap3D::save( const std::string &filename ) const
{
    std::ofstream os( filename.c_str() );
    if( !os )
        throw( std::runtime_error("CCylFieldMap3D::save: cannot open " + filename) );
    save_ascii( os );
}

void CCylFieldMap3D::load( const std::string &filename )
{
    std::ifstream is( filename.c_str() );
    if( !is )
        throw( std::runtime_error("CCylFieldMap3D::load: cannot open " + filename) );
    load_ascii( is );
}

} // namespace ibsimu_cycl
