/** @file boris.cpp
 *  @brief IBSimu-Cycl: 非相对论 Boris 推进器（实现）。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <stdexcept>

#include "boris.hpp"

namespace ibsimu_cycl {

static Vec3D crossp( const Vec3D &a, const Vec3D &b )
{
    return( Vec3D( a[1]*b[2] - a[2]*b[1],
                   a[2]*b[0] - a[0]*b[2],
                   a[0]*b[1] - a[1]*b[0] ) );
}

CBorisPusher::CBorisPusher( double q, double m )
    : _q(q), _m(m)
{
    if( _m <= 0.0 )
        throw( std::invalid_argument("CBorisPusher: mass must be > 0") );
}

void CBorisPusher::step( const VectorField *E, const VectorField *B,
                         Vec3D &x, Vec3D &v, double dt ) const
{
    const double qm = _q/_m;

    Vec3D e = E ? (*E)( x ) : Vec3D( 0.0, 0.0, 0.0 );
    Vec3D b = B ? (*B)( x ) : Vec3D( 0.0, 0.0, 0.0 );

    // 半步电场加速
    Vec3D vminus = v + e*(qm*dt/2.0);

    // 磁场旋转
    Vec3D t  = b*(qm*dt/2.0);
    double t2 = t*t;                          // Vec3D::operator* 为点积
    Vec3D vprime = vminus + crossp( vminus, t );
    Vec3D s = t*(2.0/(1.0 + t2));
    Vec3D vplus = vminus + crossp( vprime, s );

    // 另半步电场加速
    v = vplus + e*(qm*dt/2.0);

    // 位置推进
    x = x + v*dt;
}

} // namespace ibsimu_cycl
