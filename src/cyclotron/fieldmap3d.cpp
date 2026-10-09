/** @file fieldmap3d.cpp
 *  @brief IBSimu-Cycl: 笛卡尔三维矢量场图（实现）。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

#include "fieldmap3d.hpp"

namespace ibsimu_cycl {

static const char *FMAP_MAGIC = "IBSIMU-CYCL-FIELDMAP";
static const int   FMAP_VERSION = 1;

/*! \brief 单调坐标 x 所在单元及权重（越界 clamp）。 */
static void fw_cell( double x, double x0, double dx, std::size_t n,
                     int &i0, int &i1, double &t )
{
    if( n <= 1 ) { i0 = i1 = 0; t = 0.0; return; }

    double f = (x - x0)/dx;
    int    i = (int)std::floor( f );
    t = f - i;

    if( i < 0 ) { i0 = i1 = 0; t = 0.0; return; }
    if( (std::size_t)i >= n-1 ) { i0 = i1 = (int)n-1; t = 0.0; return; }
    i0 = i; i1 = i + 1;
}

CFieldMap3D::CFieldMap3D( std::size_t nx, std::size_t ny, std::size_t nz,
                          double x0, double dx, double y0, double dy,
                          double z0, double dz )
    : _nx(nx), _ny(ny), _nz(nz),
      _x0(x0), _dx(dx), _y0(y0), _dy(dy), _z0(z0), _dz(dz)
{
    if( _nx < 1 || _ny < 1 || _nz < 1 )
        throw( std::invalid_argument("CFieldMap3D: grid size must be >= 1 in each direction") );

    std::size_t total = _nx*_ny*_nz;
    _data.assign( 3*total, 0.0 );
}

CFieldMap3D::~CFieldMap3D()
{
}

void CFieldMap3D::set_value( std::size_t i, std::size_t j, std::size_t k,
                             double fx, double fy, double fz )
{
    if( i >= _nx || j >= _ny || k >= _nz )
        throw( std::out_of_range("CFieldMap3D::set_value: index out of range") );
    std::size_t n = 3*idx( i, j, k );
    _data[n] = fx; _data[n+1] = fy; _data[n+2] = fz;
}

Vec3D CFieldMap3D::node_value( std::size_t i, std::size_t j, std::size_t k ) const
{
    if( i >= _nx || j >= _ny || k >= _nz )
        throw( std::out_of_range("CFieldMap3D::node_value: index out of range") );
    std::size_t n = 3*idx( i, j, k );
    return( Vec3D( _data[n], _data[n+1], _data[n+2] ) );
}

const Vec3D CFieldMap3D::operator()( const Vec3D &x ) const
{
    int ix0, ix1, iy0, iy1, iz0, iz1;
    double tx, ty, tz;
    fw_cell( x[0], _x0, _dx, _nx, ix0, ix1, tx );
    fw_cell( x[1], _y0, _dy, _ny, iy0, iy1, ty );
    fw_cell( x[2], _z0, _dz, _nz, iz0, iz1, tz );

    const int I[2] = { ix0, ix1 };
    const int J[2] = { iy0, iy1 };
    const int K[2] = { iz0, iz1 };
    const double U[2] = { 1.0-tx, tx };
    const double V[2] = { 1.0-ty, ty };
    const double W[2] = { 1.0-tz, tz };

    const double *dat = _data.data();
    double v[3] = { 0.0, 0.0, 0.0 };
    for( int a = 0; a < 2; ++a )
        for( int b = 0; b < 2; ++b )
            for( int c = 0; c < 2; ++c ) {
                double w = U[a]*V[b]*W[c];
                const double *p = dat + 3*idx( (std::size_t)I[a], (std::size_t)J[b], (std::size_t)K[c] );
                v[0] += w*p[0];
                v[1] += w*p[1];
                v[2] += w*p[2];
            }
    return( Vec3D( v[0], v[1], v[2] ) );
}

void CFieldMap3D::save_ascii( std::ostream &os ) const
{
    os << FMAP_MAGIC << " " << FMAP_VERSION << "\n";
    os << std::setprecision(17);
    os << "# nx ny nz\n" << _nx << " " << _ny << " " << _nz << "\n";
    os << "# x0 dx y0 dy z0 dz\n"
       << _x0 << " " << _dx << " " << _y0 << " " << _dy << " " << _z0 << " " << _dz << "\n";
    os << "# data in order i(x), j(y), k(z): Fx Fy Fz\n";
    for( std::size_t i = 0; i < _nx; ++i )
        for( std::size_t j = 0; j < _ny; ++j )
            for( std::size_t k = 0; k < _nz; ++k ) {
                std::size_t n = 3*idx( i, j, k );
                os << _data[n] << " " << _data[n+1] << " " << _data[n+2] << "\n";
            }
}

void CFieldMap3D::load_ascii( std::istream &is )
{
    std::string magic; int version = 0;
    is >> magic >> version;
    if( magic != FMAP_MAGIC || version != FMAP_VERSION )
        throw( std::runtime_error("CFieldMap3D::load_ascii: bad magic/version") );

    auto skip_comments = [&is]() {
        while( is >> std::ws, is.peek() == '#' ) {
            std::string line; std::getline( is, line );
        }
    };

    skip_comments();
    std::size_t nx, ny, nz; is >> nx >> ny >> nz;
    skip_comments();
    double x0, dx, y0, dy, z0, dz; is >> x0 >> dx >> y0 >> dy >> z0 >> dz;
    if( !is )
        throw( std::runtime_error("CFieldMap3D::load_ascii: failed to read header") );

    *this = CFieldMap3D( nx, ny, nz, x0, dx, y0, dy, z0, dz );

    skip_comments();
    for( std::size_t i = 0; i < nx; ++i )
        for( std::size_t j = 0; j < ny; ++j )
            for( std::size_t k = 0; k < nz; ++k ) {
                double fx, fy, fz;
                is >> fx >> fy >> fz;
                set_value( i, j, k, fx, fy, fz );
            }
    if( !is )
        throw( std::runtime_error("CFieldMap3D::load_ascii: truncated data") );
}

void CFieldMap3D::save( const std::string &filename ) const
{
    std::ofstream os( filename.c_str() );
    if( !os )
        throw( std::runtime_error("CFieldMap3D::save: cannot open " + filename) );
    save_ascii( os );
}

void CFieldMap3D::load( const std::string &filename )
{
    std::ifstream is( filename.c_str() );
    if( !is )
        throw( std::runtime_error("CFieldMap3D::load: cannot open " + filename) );
    load_ascii( is );
}

} // namespace ibsimu_cycl
