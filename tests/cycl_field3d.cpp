/** @file cycl_field3d.cpp
 *  @brief IBSimu-Cycl P1 单元测试：柱坐标三维场图 + 时变场。
 *
 *  验证：
 *   1. 均匀场在网格内外均为精确值；
 *   2. 线性场三线性插值精确；
 *   3. 方位角(theta)周期插值在节点处精确；
 *   4. 时变场按 cos(omega*t+phi0) 调制；
 *   5. ASCII 往返读写一致。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <cstdio>
#include <cmath>
#include <sstream>
#include <string>

#include "cylfmap3d.hpp"
#include "timevaryingfield.hpp"

using namespace ibsimu_cycl;

static int g_failures = 0;

static void check( bool ok, const std::string &msg )
{
    std::printf( "%s  %s\n", ok ? "[ ok ]" : "[FAIL]", msg.c_str() );
    if( !ok )
        ++g_failures;
}

static bool approx( const Vec3D &a, const Vec3D &b, double tol )
{
    for( int i = 0; i < 3; ++i )
        if( std::fabs( a[i] - b[i] ) > tol )
            return( false );
    return( true );
}

int main()
{
    const double tol = 1e-12;

    // ---- 1. 均匀 Bz 场，处处精确 ----------------------------------
    {
        CCylFieldMap3D f( 5, 8, 3, 0.0, 0.1, -0.1, 0.1 );
        for( std::size_t i = 0; i < f.size_r(); ++i )
            for( std::size_t j = 0; j < f.size_theta(); ++j )
                for( std::size_t k = 0; k < f.size_z(); ++k )
                    f.set_value( i, j, k, 0.0, 0.0, 0.5 );

        Vec3D b = f( Vec3D( 0.123, 0.045, 0.03 ) );
        check( approx( b, Vec3D(0,0,0.5), tol ), "uniform field inside grid" );
        Vec3D b2 = f( Vec3D( 0.9, 0.9, 0.9 ) );   // 越界 clamp
        check( approx( b2, Vec3D(0,0,0.5), tol ), "uniform field outside grid (clamp)" );
    }

    // ---- 2. 线性径向场 Br = r 三线性插值精确 ----------------------
    {
        CCylFieldMap3D f( 5, 8, 3, 0.0, 0.1, -0.1, 0.1 );
        for( std::size_t i = 0; i < f.size_r(); ++i )
            for( std::size_t j = 0; j < f.size_theta(); ++j )
                for( std::size_t k = 0; k < f.size_z(); ++k )
                    f.set_value( i, j, k, 0.1*(double)i, 0.0, 0.0 );

        // r=0.23, theta=pi/4 -> Bx = r*cos(theta)
        double r = 0.23, th = M_PI/4.0;
        Vec3D x( r*std::cos(th), r*std::sin(th), 0.05 );
        Vec3D b = f( x );
        Vec3D expect( r*std::cos(th), r*std::sin(th), 0.0 );
        check( approx( b, expect, tol ), "linear radial field Br=r interpolates exactly" );
    }

    // ---- 3. 方位角依赖 Bz=cos(theta) 在节点处精确 -----------------
    {
        CCylFieldMap3D f( 3, 8, 1, 0.0, 0.1, 0.0, 1.0 );
        for( std::size_t i = 0; i < f.size_r(); ++i )
            for( std::size_t j = 0; j < f.size_theta(); ++j )
                f.set_value( i, j, 0, 0.0, 0.0, std::cos( f.dtheta()*(double)j ) );

        double r = 0.2;
        for( int j = 0; j < 8; ++j ) {
            double th = f.dtheta()*(double)j;
            Vec3D x( r*std::cos(th), r*std::sin(th), 0.0 );
            Vec3D b = f( x );
            double expect = std::cos( th );
            if( !approx( b, Vec3D(0,0,expect), 1e-11 ) ) {
                check( false, "azimuthal field Bz=cos(theta) at node j=" + std::to_string(j) );
                return( 1 );
            }
        }
        check( true, "azimuthal field Bz=cos(theta) exact at all 8 nodes" );
    }

    // ---- 4. 时变场调制 --------------------------------------------
    {
        CCylFieldMap3D spatial( 2, 4, 2, 0.0, 0.1, 0.0, 0.1 );
        for( std::size_t i = 0; i < spatial.size_r(); ++i )
            for( std::size_t j = 0; j < spatial.size_theta(); ++j )
                for( std::size_t k = 0; k < spatial.size_z(); ++k )
                    spatial.set_value( i, j, k, 0.0, 0.0, 0.4 );

        CTimeVaryingField tv( &spatial, 1.0 /*omega*/, 0.0 /*phi0*/ );
        Vec3D x( 0.05, 0.05, 0.05 );

        tv.set_time( 0.0 );
        check( approx( tv(x), Vec3D(0,0,0.4), tol ), "time-varying: phase 0 -> full amplitude" );
        tv.set_time( M_PI/2.0 );
        check( approx( tv(x), Vec3D(0,0,0.0), 1e-12 ), "time-varying: phase pi/2 -> zero" );
        tv.set_time( M_PI );
        check( approx( tv(x), Vec3D(0,0,-0.4), 1e-12 ), "time-varying: phase pi -> negative" );
    }

    // ---- 5. ASCII 往返 --------------------------------------------
    {
        CCylFieldMap3D f( 3, 4, 2, 0.0, 0.05, -0.05, 0.05 );
        for( std::size_t i = 0; i < f.size_r(); ++i )
            for( std::size_t j = 0; j < f.size_theta(); ++j )
                for( std::size_t k = 0; k < f.size_z(); ++k )
                    f.set_value( i, j, k, 0.01*(double)i, -0.02*(double)j, 0.03*(double)k );

        std::stringstream ss;
        f.save_ascii( ss );

        CCylFieldMap3D g( 1, 1, 1, 0.0, 1.0, 0.0, 1.0 );
        g.load_ascii( ss );

        bool same = ( g.size_r() == f.size_r() ) &&
                    ( g.size_theta() == f.size_theta() ) &&
                    ( g.size_z() == f.size_z() );
        for( std::size_t i = 0; same && i < f.size_r(); ++i )
            for( std::size_t j = 0; same && j < f.size_theta(); ++j )
                for( std::size_t k = 0; same && k < f.size_z(); ++k )
                    if( !approx( f.node_value(i,j,k), g.node_value(i,j,k), 0.0 ) )
                        same = false;
        check( same, "ASCII save/load round-trip preserves all nodes" );
    }

    std::printf( "\n%s (%d failure%s)\n",
                 g_failures ? "FAILED" : "ALL TESTS PASSED",
                 g_failures, g_failures == 1 ? "" : "s" );
    return( g_failures ? 1 : 0 );
}
