/** @file timevaryingfield.cpp
 *  @brief IBSimu-Cycl: 时变矢量场（实现）。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <cmath>
#include "timevaryingfield.hpp"

namespace ibsimu_cycl {

CTimeVaryingField::CTimeVaryingField( const VectorField *spatial, double omega,
                                      double phi0, double amplitude, double offset )
    : _spatial(spatial), _omega(omega), _phi0(phi0),
      _amplitude(amplitude), _offset(offset), _t(0.0)
{
}

CTimeVaryingField::~CTimeVaryingField()
{
}

double CTimeVaryingField::scale_factor() const
{
    return( _offset + _amplitude*std::cos( phase() ) );
}

const Vec3D CTimeVaryingField::operator()( const Vec3D &x ) const
{
    if( _spatial == 0 )
        return( Vec3D( 0.0, 0.0, 0.0 ) );

    return( (*_spatial)( x ) * scale_factor() );
}

} // namespace ibsimu_cycl
