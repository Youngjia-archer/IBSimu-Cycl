/** @file cycl_fieldmap3d.cpp
 *  @brief IBSimu-Cycl P2-C 测试：笛卡尔三维场图（每轴独立步长）。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <cstdio>
#include <cmath>
#include <sstream>
#include <string>

#include "fieldmap3d.hpp"

using namespace ibsimu_cycl;

static int g_failures = 0;

static void check( bool ok, const std::string &msg )
{
    std::printf( "%s  %s\n", ok ? "[ ok ]" : "[FAIL]", msg.c_str() );
    if( !ok ) ++g_failures;
}

static bool approx( const Vec3D &a, const Vec3D &b, double tol )
{
    for( int i = 0; i < 3; ++i )
        if( std::fabs( a[i] - b[i] ) > tol ) return( false );
    return( true );
}

int main()
{
    const double tol = 1e-12;

    // ---- 1. 均匀场（非立方步长） ----
    {
        CFieldMap3D f( 4, 5, 3, 0.0, 0.1, -0.2, 0.2, 1.0, 0.05 );
        for( std::size_t i = 0; i < f.size_x(); ++i )
            for( std::size_t j = 0; j < f.size_y(); ++j )
                for( std::size_t k = 0; k < f.size_z(); ++k )
                    f.set_value( i, j, k, 0.0, 0.0, 1.2 );

        check( approx( f( Vec3D( 0.17, -0.05, 1.07 ) ), Vec3D(0,0,1.2), tol ),
               "uniform field, non-cubic spacing: exact inside" );
        check( approx( f( Vec3D( 9, 9, 9 ) ), Vec3D(0,0,1.2), tol ),
               "uniform field: exact outside (clamp)" );
    }

    // ---- 2. 线性场精确（三线性可复现线性函数） ----
    {
        CFieldMap3D f( 5, 5, 5, 0.0, 0.1, 0.0, 0.2, 0.0, 0.05 );
        // F = (x, 2y, 3z)
        for( std::size_t i = 0; i < f.size_x(); ++i )
            for( std::size_t j = 0; j < f.size_y(); ++j )
                for( std::size_t k = 0; k < f.size_z(); ++k )
                    f.set_value( i, j, k,
                                 f.x0()+0.1*(double)i,
                                 2.0*(f.y0()+0.2*(double)j),
                                 3.0*(f.z0()+0.05*(double)k) );

        Vec3D x( 0.23, 0.37, 0.13 );
        check( approx( f(x), Vec3D( 0.23, 0.74, 0.39 ), 1e-12 ),
               "linear field F=(x,2y,3z): trilinear reproduces exactly" );
    }

    // ---- 3. ASCII 往返 ----
    {
        CFieldMap3D f( 3, 4, 2, -0.1, 0.05, 0.0, 0.1, 0.0, 0.3 );
        for( std::size_t i = 0; i < f.size_x(); ++i )
            for( std::size_t j = 0; j < f.size_y(); ++j )
                for( std::size_t k = 0; k < f.size_z(); ++k )
                    f.set_value( i, j, k, 0.01*(double)i, 0.02*(double)j, -0.03*(double)k );

        std::stringstream ss;
        f.save_ascii( ss );

        CFieldMap3D g( 1, 1, 1, 0, 1, 0, 1, 0, 1 );
        g.load_ascii( ss );

        bool same = ( g.size_x() == f.size_x() ) && ( g.size_y() == f.size_y() ) &&
                    ( g.size_z() == f.size_z() ) && ( g.dx() == f.dx() ) &&
                    ( g.dy() == f.dy() ) && ( g.dz() == f.dz() );
        for( std::size_t i = 0; same && i < f.size_x(); ++i )
            for( std::size_t j = 0; same && j < f.size_y(); ++j )
                for( std::size_t k = 0; same && k < f.size_z(); ++k )
                    if( !approx( f.node_value(i,j,k), g.node_value(i,j,k), 0.0 ) )
                        same = false;
        check( same, "ASCII save/load round-trip preserves grid and all nodes" );
    }

    std::printf( "\n%s (%d failure%s)\n",
                 g_failures ? "FAILED" : "ALL TESTS PASSED",
                 g_failures, g_failures == 1 ? "" : "s" );
    return( g_failures ? 1 : 0 );
}
