/** @file cycl_relativistic.cpp
 *  @brief IBSimu-Cycl F1 测试：相对论 Boris 与真实场强下的回旋轨道。
 *
 *  1. 均匀场：相对论回旋频率 omega_c = qB/(gamma m)，且 gamma 精确守恒；
 *  2. 对比：同一速度下非相对论推进器给出的是 qB/m（偏差 gamma 倍，回旋半径小 gamma 倍）；
 *  3. 真实 PSI Ring 场图（**真实场强**，beta~0.6）：由 p = q<B>r 设置注入动量，
 *     验证轨道有界、能量守恒、回路频率与 q<B>/(2 pi gamma m) 一致。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "boris.hpp"
#include "fieldmap3d.hpp"
#include "ringfield3d.hpp"

using namespace ibsimu_cycl;

static const double QE = 1.602176634e-19;
static const double MP = 1.67262192369e-27;
static const double CL = 299792458.0;

static int g_failures = 0;
static void check( bool ok, const std::string &msg )
{
    std::printf( "%s  %s\n", ok ? "[ ok ]" : "[FAIL]", msg.c_str() );
    if( !ok ) ++g_failures;
}

static CFieldMap3D uniform_B( double bz )
{
    CFieldMap3D f( 2, 2, 2, -5.0, 10.0, -5.0, 10.0, -5.0, 10.0 );
    for( std::size_t i = 0; i < 2; ++i )
        for( std::size_t j = 0; j < 2; ++j )
            for( std::size_t k = 0; k < 2; ++k )
                f.set_value( i, j, k, 0.0, 0.0, bz );
    return( f );
}

static double gamma_of( const Vec3D &v )
{
    double b2 = (v[0]*v[0] + v[1]*v[1] + v[2]*v[2])/(CL*CL);
    return( 1.0/std::sqrt(1.0 - b2) );
}

/*! \brief 在均匀场中推进并测量回旋角频率与 gamma 漂移。 */
static double measure_omega( const CFieldMap3D &B, double q, double m,
                             const Vec3D &v0, double dt, int nsteps,
                             bool relativistic, double &gamma_drift )
{
    CBorisPusher pusher( q, m );
    pusher.set_relativistic( relativistic );

    Vec3D x( 0.0, 0.0, 0.0 ), v = v0;
    double g0 = gamma_of( v );
    double phi_prev = std::atan2( v[1], v[0] );
    double phi_sum = 0.0;
    double gmax = g0, gmin = g0;

    for( int n = 0; n < nsteps; ++n ) {
        pusher.step( 0, &B, x, v, dt );
        double phi = std::atan2( v[1], v[0] );
        double d = phi - phi_prev;
        while( d >  M_PI ) d -= 2.0*M_PI;
        while( d < -M_PI ) d += 2.0*M_PI;
        phi_sum += d;
        phi_prev = phi;
        double g = gamma_of( v );
        if( g > gmax ) gmax = g;
        if( g < gmin ) gmin = g;
    }
    gamma_drift = (gmax - gmin)/g0;
    return( std::fabs(phi_sum)/(nsteps*dt) );
}

int main( int argc, char **argv )
{
    const double B0 = 1.0;                        // T
    const double gamma0 = 1.2;
    const double beta0 = std::sqrt( 1.0 - 1.0/(gamma0*gamma0) );
    const double v0 = beta0*CL;

    // =====================================================================
    // 1. 均匀场：相对论回旋频率
    // =====================================================================
    {
        CFieldMap3D Bf = uniform_B( B0 );
        double omega_exp = QE*B0/(gamma0*MP);
        double T = 2.0*M_PI/omega_exp;
        double dt = T/500.0;
        double gdrift = 0.0;
        // 速度垂直于 B（在 xy 平面）
        double omega_meas = measure_omega( Bf, QE, MP, Vec3D( v0, 0.0, 0.0 ),
                                           dt, 1500, true, gdrift );
        std::printf( "gamma0=%.3f  beta0=%.4f  v0=%.4e m/s\n", gamma0, beta0, v0 );
        std::printf( "1) relativistic: omega=%.6e (expect %.6e)  rel=%.2e  gamma drift=%.2e\n",
                     omega_meas, omega_exp, std::fabs(omega_meas-omega_exp)/omega_exp, gdrift );
        check( std::fabs(omega_meas - omega_exp)/omega_exp < 1e-3,
               "relativistic gyrofrequency == qB/(gamma m)" );
        check( gdrift < 1e-12, "relativistic: gamma conserved in static B" );
    }

    // =====================================================================
    // 2. 对比：非相对论推进器（同一速度）
    // =====================================================================
    {
        CFieldMap3D Bf = uniform_B( B0 );
        double omega_rel = QE*B0/(gamma0*MP);
        double omega_nr_exp = QE*B0/MP;
        double T = 2.0*M_PI/omega_rel;
        double dt = T/500.0;
        double gdrift = 0.0;
        double omega_nr = measure_omega( Bf, QE, MP, Vec3D( v0, 0.0, 0.0 ),
                                         dt, 1500, false, gdrift );
        std::printf( "2) non-relativistic: omega=%.6e (expect %.6e)  ratio=%.4f (gamma=%.4f)\n",
                     omega_nr, omega_nr_exp, omega_nr/omega_rel, gamma0 );
        check( std::fabs(omega_nr/omega_rel - gamma0)/gamma0 < 1e-3,
               "non-relativistic omega is gamma times too high" );
    }

    // =====================================================================
    // 3. 真实 PSI Ring 场图（真实场强）
    // =====================================================================
    std::vector<std::string> cands;
    if( argc > 1 ) cands.push_back( argv[1] );
    cands.push_back( "../examples/cyclotron/data/bfield.dat" );
    cands.push_back( "examples/cyclotron/data/bfield.dat" );

    std::unique_ptr<CRingFieldMap3D> map;
    std::string file;
    for( std::size_t c = 0; c < cands.size() && !map; ++c ) {
        try {
            map.reset( new CRingFieldMap3D( read_ring_field_map3d( cands[c], 0.1, 8, 0 ) ) );
            file = cands[c];
        } catch( const std::exception & ) { }
    }
    if( !map ) {
        std::printf( "\n[SKIP] real-field test: map not found\n" );
    } else {
        const double r_ref = 3.3;
        std::size_t i_ref = (std::size_t)std::lround( (r_ref - map->r0())/map->dr() );
        if( i_ref >= map->size_r() ) i_ref = map->size_r()-1;
        double Bbar = 0.0;
        for( std::size_t k = 0; k < map->size_theta(); ++k )
            Bbar += map->bz_midplane( i_ref, k );
        Bbar /= (double)map->size_theta();

        // 由 p = q <B> r 设定注入动量（相对论）
        double p   = QE*Bbar*r_ref;
        double gb  = p/(MP*CL);                       // gamma*beta
        double gam = std::sqrt(1.0 + gb*gb);
        double bet = gb/gam;
        double v   = bet*CL;
        double omega_exp = QE*Bbar/(gam*MP);
        double T_rev = 2.0*M_PI/omega_exp;
        double KE = (gam-1.0)*MP*CL*CL;

        std::printf( "\n3) real PSI Ring (file=%s)\n", file.c_str() );
        std::printf( "   <B>(r=%.2f m)=%.4f T  -> p=%.4e kg m/s  gamma=%.4f  beta=%.4f\n",
                     r_ref, Bbar, p, gam, bet );
        std::printf( "   KE=%.2f MeV  T_rev=%.2f ns  omega_exp=%.6e rad/s\n",
                     KE/QE/1e6, T_rev/1e-9, omega_exp );

        CBorisPusher pusher( QE, MP );
        pusher.set_relativistic( true );

        Vec3D x( 0.0, r_ref, 0.0 ), vv( v, 0.0, 0.0 );
        double dt = T_rev/4000.0;
        int nsteps = (int)(2.5*T_rev/dt);

        double rmin = 1e30, rmax = -1e30, g0 = gamma_of( vv ), gmax = g0, gmin = g0;
        std::vector<double> tcross;
        double yprev = x[1];
        std::ofstream csv( "cycl_relativistic.csv" );
        csv << "t,x,y,z,gamma,r,Bmag\n";

        for( int n = 1; n <= nsteps; ++n ) {
            pusher.step( 0, map.get(), x, vv, dt );
            double t = n*dt;
            double r = std::sqrt( x[0]*x[0] + x[1]*x[1] );
            if( r < rmin ) rmin = r;
            if( r > rmax ) rmax = r;
            double g = gamma_of( vv );
            if( g > gmax ) gmax = g;
            if( g < gmin ) gmin = g;
            Vec3D bb = (*map)( x );
            csv << t << "," << x[0] << "," << x[1] << "," << x[2] << ","
                << g << "," << r << ","
                << std::sqrt( bb[0]*bb[0] + bb[1]*bb[1] + bb[2]*bb[2] ) << "\n";

            if( yprev > 0.0 && x[1] <= 0.0 && x[0] > 0.0 )
                tcross.push_back( t );
            yprev = x[1];
        }

        std::printf( "   orbit r: %.4f .. %.4f m (mean %.4f)  gamma drift=%.2e\n",
                     rmin, rmax, 0.5*(rmin+rmax), (gmax-gmin)/g0 );
        check( rmin > map->r0() && rmax < map->r0() + (map->size_r()-1)*map->dr(),
               "orbit stays inside the field map" );
        check( std::fabs(0.5*(rmin+rmax) - r_ref)/r_ref < 0.15,
               "bounded orbit around r_ref (mean within 15%)" );
        check( (gmax-gmin)/g0 < 1e-9, "real field: gamma conserved" );

        if( tcross.size() >= 2 ) {
            double T_meas = (tcross.back() - tcross.front())/(double)(tcross.size()-1);
            double omega_meas = 2.0*M_PI/T_meas;
            std::printf( "   crossings=%d  T_rev measured=%.2f ns  omega=%.6e  rel.err=%.3f%%\n",
                         (int)tcross.size(), T_meas/1e-9, omega_meas,
                         100.0*std::fabs(omega_meas-omega_exp)/omega_exp );
            check( std::fabs(omega_meas - omega_exp)/omega_exp < 0.05,
                   "revolution frequency matches q<B>/(2 pi gamma m) within 5%" );
        }
    }

    std::printf( "\n%s (%d failure%s)\n",
                 g_failures ? "FAILED" : "ALL TESTS PASSED",
                 g_failures, g_failures == 1 ? "" : "s" );
    return( g_failures ? 1 : 0 );
}
