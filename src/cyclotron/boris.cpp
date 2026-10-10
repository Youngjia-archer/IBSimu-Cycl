/** @file boris.cpp
 *  @brief IBSimu-Cycl: 非相对论 Boris 推进器（实现）。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <cmath>
#include <stdexcept>

#include "boris.hpp"

namespace ibsimu_cycl {

static const double CLIGHT = 299792458.0;   /*!< \brief 真空光速 [m/s]。 */

static Vec3D crossp( const Vec3D &a, const Vec3D &b )
{
    return( Vec3D( a[1]*b[2] - a[2]*b[1],
                   a[2]*b[0] - a[0]*b[2],
                   a[0]*b[1] - a[1]*b[0] ) );
}

CBorisPusher::CBorisPusher( double q, double m )
    : _q(q), _m(m), _relativistic(false)
{
    if( _m <= 0.0 )
        throw( std::invalid_argument("CBorisPusher: mass must be > 0") );
}

double CBorisPusher::gamma( const Vec3D &v ) const
{
    if( !_relativistic )
        return( 1.0 );
    double b2 = (v*v)/(CLIGHT*CLIGHT);
    if( b2 >= 1.0 )
        throw( std::invalid_argument("CBorisPusher: speed >= c") );
    return( 1.0/std::sqrt(1.0 - b2) );
}

void CBorisPusher::initialize( const VectorField *E, const VectorField *B,
                               const Vec3D &x, Vec3D &v, double dt ) const
{
    const double qm = _q/_m;

    Vec3D e = E ? (*E)( x ) : Vec3D( 0.0, 0.0, 0.0 );
    Vec3D b = B ? (*B)( x ) : Vec3D( 0.0, 0.0, 0.0 );
    Vec3D a = e + crossp( v, b );                 // 未乘 q/m 的加速度方向

    if( !_relativistic ) {
        v = v - a*(qm*dt/2.0);
        return;
    }

    // 相对论：在 u = gamma*v 空间做半步反踢
    Vec3D  u = v*gamma( v ) - a*(qm*dt/2.0);
    double gn = std::sqrt( 1.0 + (u*u)/(CLIGHT*CLIGHT) );
    v = u*(1.0/gn);
}

void CBorisPusher::step( const VectorField *E, const VectorField *B,
                         Vec3D &x, Vec3D &v, double dt ) const
{
    const double qm = _q/_m;

    Vec3D e = E ? (*E)( x ) : Vec3D( 0.0, 0.0, 0.0 );
    Vec3D b = B ? (*B)( x ) : Vec3D( 0.0, 0.0, 0.0 );

    if( !_relativistic ) {
        // ---- 非相对论 Boris ----
        Vec3D vminus = v + e*(qm*dt/2.0);
        Vec3D t  = b*(qm*dt/2.0);
        double t2 = t*t;                          // Vec3D::operator* 为点积
        Vec3D vprime = vminus + crossp( vminus, t );
        Vec3D s = t*(2.0/(1.0 + t2));
        Vec3D vplus = vminus + crossp( vprime, s );
        v = vplus + e*(qm*dt/2.0);
        x = x + v*dt;
        return;
    }

    // ---- 相对论 Boris（u = gamma*v 空间）----
    double v2 = v*v;
    double g  = 1.0/std::sqrt( 1.0 - v2/(CLIGHT*CLIGHT) );
    Vec3D  u  = v*g;

    u = u + e*(qm*dt/2.0);                        // 半步电场
    double u2    = u*u;
    double ghalf = std::sqrt( 1.0 + u2/(CLIGHT*CLIGHT) );
    Vec3D  t  = b*(qm*dt/(2.0*ghalf));            // 磁场旋转
    double t2 = t*t;
    Vec3D  uprime = u + crossp( u, t );
    Vec3D  s = t*(2.0/(1.0 + t2));
    u = u + crossp( uprime, s );
    u = u + e*(qm*dt/2.0);                        // 另半步电场

    double u2n = u*u;
    double gn  = std::sqrt( 1.0 + u2n/(CLIGHT*CLIGHT) );
    v = u*(1.0/gn);
    x = x + v*dt;
}

} // namespace ibsimu_cycl
