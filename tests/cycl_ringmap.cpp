/** @file cycl_ringmap.cpp
 *  @brief IBSimu-Cycl P2 测试：读取真实 OPAL/PSI RING 磁场图 (bfield.dat)。
 *
 *  用法: cycl_ringmap [bfield.dat 路径]
 *  默认路径: ../examples/cyclotron/data/bfield.dat
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <cstdio>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "ringmap.hpp"

using namespace ibsimu_cycl;

static int g_failures = 0;

static void check( bool ok, const std::string &msg )
{
    std::printf( "%s  %s\n", ok ? "[ ok ]" : "[FAIL]", msg.c_str() );
    if( !ok )
        ++g_failures;
}

static bool near( double a, double b, double tol )
{
    return( std::fabs( a - b ) <= tol );
}

int main( int argc, char **argv )
{
    std::vector<std::string> candidates;
    if( argc > 1 )
        candidates.push_back( argv[1] );
    candidates.push_back( "../examples/cyclotron/data/bfield.dat" ); // CWD = tests/
    candidates.push_back( "examples/cyclotron/data/bfield.dat" );    // CWD = repo root

    RingFieldMapInfo info;
    std::unique_ptr<CCylFieldMap3D> mapp;
    std::string file;
    for( std::size_t c = 0; c < candidates.size() && !mapp; ++c ) {
        try {
            mapp.reset( new CCylFieldMap3D(
                            read_ring_field_map( candidates[c], 0.1 /*kG->T*/, 8 /*symmetry*/, &info ) ) );
            file = candidates[c];
        } catch( const std::exception & ) { /* 尝试下一个候选路径 */ }
    }
    if( !mapp ) {
        std::printf( "[SKIP] field map not found (tried %d candidate paths)\n",
                     (int)candidates.size() );
        return( 0 );   // 数据缺失时不判失败（例如 out-of-tree 构建）
    }
    std::printf( "reading: %s\n", file.c_str() );
    CCylFieldMap3D &map = *mapp;

    // ---- 头部元数据 ----
    std::printf( "header: nrad=%d ntet=%d rmin=%.3f mm dr=%.3f mm "
                 "tetmin=%.3f deg dtet=%.6f deg lpar=%d\n",
                 info.nrad, info.ntet, info.rmin_mm, info.dr_mm,
                 info.tetmin_deg, info.dtet_deg, info.lpar );

    check( info.nrad == 141, "nrad == 141" );
    check( info.ntet == 135, "ntet == 135" );
    check( near( info.rmin_mm, 1900.0, 1e-6 ), "rmin == 1900 mm" );
    check( near( info.dr_mm, 20.0, 1e-6 ), "dr == 20 mm" );
    check( near( info.dtet_deg, 1.0/3.0, 1e-9 ), "dtet == 1/3 deg" );

    // ---- 网格几何（毫米 -> 米） ----
    check( map.size_r() == 141, "map size_r == 141" );
    check( map.size_theta() == 1080, "map size_theta == 1080 (360/0.3333)" );
    check( map.size_z() == 1, "map size_z == 1 (midplane only)" );

    // ---- 8 折旋转复制：跨扇区一致性 ----
    {
        bool same = true;
        for( std::size_t k = 0; k < 135; ++k )
            if( map.node_value( 50, k, 0 )[2] != map.node_value( 50, k + 135, 0 )[2] ) {
                same = false;
                break;
            }
        check( same, "8-fold replication: sector values repeat every 135 columns" );
    }

    // ---- 物理合理性：中平面 |Bz| 量级 ----
    {
        double bmin = 1e30, bmax = -1e30;
        std::size_t ntot = 0, ntiny = 0, nbig = 0;
        for( std::size_t i = 0; i < map.size_r(); ++i )
            for( std::size_t k = 0; k < map.size_theta(); ++k ) {
                double b = std::fabs( map.node_value( i, k, 0 )[2] );
                if( b < bmin ) bmin = b;
                if( b > bmax ) bmax = b;
                ++ntot;
                if( b < 1e-6 ) ++ntiny;
                if( b > 0.1 ) ++nbig;
            }
        std::printf( "midplane |Bz|: min=%.6e  max=%.6f T  | <1e-6: %zu/%zu  | >0.1T: %zu\n",
                     bmin, bmax, ntiny, ntot, nbig );
        std::printf( "Bz(r=%.4f m, theta=0) = %.6f T\n",
                     map.r0() + 60.0*map.dr(), map.node_value( 60, 0, 0 )[2] );
        check( bmax > 0.5 && bmax < 6.0, "|Bz| peak in (0.5, 6) T" );
        check( nbig*2 > ntot, "majority of nodes have |Bz| > 0.1 T" );
        check( bmax > bmin, "field actually varies" );
    }

    // ---- 插值一致性：节点处应精确复现 ----
    {
        Vec3D b0 = map.node_value( 100, 300, 0 );
        double r = map.r0() + 100.0*map.dr();
        double th = map.dtheta()*300.0;
        Vec3D bx( r*std::cos(th), r*std::sin(th), 0.0 );
        Vec3D bi = map( bx );
        check( near( bi[2], b0[2], 1e-9 ), "interpolation reproduces node value" );
    }

    std::printf( "\n%s (%d failure%s)\n",
                 g_failures ? "FAILED" : "ALL TESTS PASSED",
                 g_failures, g_failures == 1 ? "" : "s" );
    return( g_failures ? 1 : 0 );
}
