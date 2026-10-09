/** @file ensembletracker.cpp
 *  @brief IBSimu-Cycl: 粒子系综跟踪器（实现）。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <cmath>
#include <stdexcept>

#include "boris.hpp"
#include "ensembletracker.hpp"

namespace ibsimu_cycl {

static const double CLIGHT = 299792458.0;

CEnsembleTracker::CEnsembleTracker( const VectorField *B, CTimeVaryingField *E,
                                    double q, double m )
    : _B(B), _E(E), _q(q), _m(m), _relativistic(false)
{
    if( _m <= 0.0 )
        throw( std::invalid_argument("CEnsembleTracker: mass must be > 0") );
}

void CEnsembleTracker::track( std::vector<EnsembleParticle> &p, double dt, int nsteps ) const
{
    CBorisPusher pusher( _q, _m );
    pusher.set_relativistic( _relativistic );

    const long n = (long)p.size();
    if( n <= 0 || nsteps <= 0 )
        return;

#ifdef _OPENMP
    // 单一并行区 + 每步屏障。
    // 若对每一步都开一个 `parallel for`，则 nsteps 次 fork/join 的开销在
    // 场求值较快（小场图/解析场）时会完全吃掉并行收益。这里只建立一次
    // 线程组，步间用屏障同步；时间步的相位更新放在 single 区，屏障保证
    // 同一时间步内所有粒子看到相同的 RF 相位。
#pragma omp parallel
    {
        for( int step = 0; step < nsteps; ++step ) {

#pragma omp single
            {
                if( _E )
                    _E->set_time( (double)(step+1)*dt );
            }

#pragma omp for schedule(static)
            for( long i = 0; i < n; ++i )
                pusher.step( _E, _B, p[i].x, p[i].v, dt );
        }
    }
#else
    for( int step = 0; step < nsteps; ++step ) {
        // 时变场时间在串行区设置，保证所有粒子看到相同相位
        if( _E )
            _E->set_time( (double)(step+1)*dt );

        for( long i = 0; i < n; ++i )
            pusher.step( _E, _B, p[i].x, p[i].v, dt );
    }
#endif
}

void CEnsembleTracker::track( std::vector<EnsembleParticle> &p, double dt, int nsteps,
                              std::vector<EnsembleStats> &history ) const
{
    history.clear();
    if( nsteps > 0 )
        history.reserve( (std::size_t)nsteps );

    CBorisPusher pusher( _q, _m );
    pusher.set_relativistic( _relativistic );

    const long n = (long)p.size();
    if( n <= 0 || nsteps <= 0 )
        return;

#ifdef _OPENMP
#pragma omp parallel
    {
        for( int step = 0; step < nsteps; ++step ) {

#pragma omp single
            {
                if( _E )
                    _E->set_time( (double)(step+1)*dt );
            }

#pragma omp for schedule(static)
            for( long i = 0; i < n; ++i )
                pusher.step( _E, _B, p[i].x, p[i].v, dt );

            // 屏障之后所有粒子已推进到同一时刻，统计量由单线程汇总
#pragma omp single
            history.push_back( analyze( p, _m, _relativistic ) );
        }
    }
#else
    for( int step = 0; step < nsteps; ++step ) {
        if( _E )
            _E->set_time( (double)(step+1)*dt );

        for( long i = 0; i < n; ++i )
            pusher.step( _E, _B, p[i].x, p[i].v, dt );

        history.push_back( analyze( p, _m, _relativistic ) );
    }
#endif
}

EnsembleStats CEnsembleTracker::analyze( const std::vector<EnsembleParticle> &p,
                                         double m, bool relativistic )
{
    EnsembleStats s;
    s.n = p.size();
    s.KE_min =  1.0e300;
    s.KE_max = -1.0e300;
    s.KE_mean = 0.0;
    s.r_min =  1.0e300;
    s.r_max = -1.0e300;
    if( p.empty() )
        return( s );

    for( std::size_t i = 0; i < p.size(); ++i ) {
        const Vec3D &v = p[i].v;
        double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
        double ke;
        if( relativistic ) {
            double g = 1.0/std::sqrt( 1.0 - v2/(CLIGHT*CLIGHT) );
            ke = (g - 1.0)*m*CLIGHT*CLIGHT;
        } else {
            ke = 0.5*m*v2;
        }
        double r = std::sqrt( p[i].x[0]*p[i].x[0] + p[i].x[1]*p[i].x[1] );
        if( ke < s.KE_min ) s.KE_min = ke;
        if( ke > s.KE_max ) s.KE_max = ke;
        s.KE_mean += ke;
        if( r < s.r_min ) s.r_min = r;
        if( r > s.r_max ) s.r_max = r;
    }
    s.KE_mean /= (double)p.size();
    return( s );
}

} // namespace ibsimu_cycl
