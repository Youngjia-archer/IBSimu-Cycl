/** @file cycl_fieldbench.cpp
 *  @brief IBSimu-Cycl 基准：三维场图求值吞吐量（单线程）。
 *
 *  场插值是三维跟踪的主要开销，也是 R1（性能）关注的核心**单核**指标。
 *  本基准对两种三维场图各做一次大规模求值计时：
 *    - CCylFieldMap3D  柱坐标 (r,theta,z)，含方位角变化（回旋加速器扇形场）
 *    - CFieldMap3D     笛卡尔 (x,y,z)，每轴独立步长
 *
 *  同时做基本正确性检查（节点处精确恢复、结果有界），避免"快但错"。
 *  采样点取绕多圈的螺旋线，避免规则访问造成不真实的缓存友好模式。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <cmath>
#include <cstdio>
#include <ctime>
#include <string>
#include <vector>

#include "cylfmap3d.hpp"
#include "fieldmap3d.hpp"

using namespace ibsimu_cycl;

static int g_failures = 0;

static void check( bool ok, const std::string &msg )
{
    std::printf( "%s  %s\n", ok ? "[ ok ]" : "[FAIL]", msg.c_str() );
    if( !ok ) ++g_failures;
}

static double now_sec()
{
    return( (double)std::clock()/(double)CLOCKS_PER_SEC );
}

/*! \brief 计时：先预热一轮，再取 reps 轮中的最短耗时（消除瞬时抖动）。 */
template< class F >
static double time_best( F f, int reps = 3 )
{
    f();
    double best = 1.0e30;
    for( int r = 0; r < reps; ++r ) {
        double t0 = now_sec();
        f();
        double dt = now_sec() - t0;
        if( dt < best ) best = dt;
    }
    return( best );
}

/*! \brief 螺旋采样点（半径 r0->r1，绕 turns 圈，带 z 起伏）。 */
static std::vector<Vec3D> make_samples( long n, double r0, double r1,
                                        double turns, double zamp )
{
    std::vector<Vec3D> p( (std::size_t)n );
    for( long i = 0; i < n; ++i ) {
        double u  = (double)i/(double)n;
        double r  = r0 + (r1-r0)*u;
        double th = 2.0*M_PI*turns*u;
        p[(std::size_t)i] = Vec3D( r*std::cos(th), r*std::sin(th), zamp*std::sin(3.0*th) );
    }
    return( p );
}

int main()
{
    const long NS = 2000000;   /*!< \brief 采样点数（足够长以消除计时噪声）。 */

    // ---------------- 柱坐标场图：含 8 折方位角变化 + z 方向 ----------------
    const std::size_t nr = 141, ntet = 1080, nzc = 9;
    const double r0 = 1.90, dr = 0.02, z0 = -0.08, dz = 0.02;

    CCylFieldMap3D cm( nr, ntet, nzc, r0, dr, z0, dz );
    for( std::size_t i = 0; i < nr; ++i ) {
        double r  = r0 + dr*(double)i;
        for( std::size_t j = 0; j < ntet; ++j ) {
            double th = 2.0*M_PI*(double)j/(double)ntet;
            double c8 = std::cos( 8.0*th ), s8 = std::sin( 8.0*th );
            double bz = 1.5*(1.0 - 0.02*(r-3.0)*(r-3.0))*(1.0 + 0.02*c8);
            double br = -0.05*(r-3.0)*(1.0 + 0.02*c8);
            double bt =  0.02*s8;
            for( std::size_t k = 0; k < nzc; ++k )
                cm.set_value( i, j, k, br, bt, bz );
        }
    }

    // 节点处必须精确恢复（同时验证索引映射）
    {
        std::size_t i = 55, k = 4;                   // r=3.0, z=0.0 恰为节点
        Vec3D node = cm.node_value( i, 0, k );
        Vec3D got  = cm( Vec3D( 3.0, 0.0, 0.0 ) );
        check( std::fabs( got[2] - node[2] ) < 1e-12 && std::fabs( got[0] - node[0] ) < 1e-12,
               "cyl map: interpolation reproduces node value exactly" );
    }

    std::vector<Vec3D> ps = make_samples( NS, 2.8, 3.3, 20.0, 0.05 );

    double acc = 0.0;
    auto eval_cyl = [&]() {
        for( long i = 0; i < NS; ++i )
            acc += cm( ps[(std::size_t)i] )[2];
    };
    double tc = time_best( eval_cyl );
    std::printf( "CCylFieldMap3D %zux%zux%zu (%5.1f MB) : %7.3f s   %8.2f Mevals/s\n",
                 nr, ntet, nzc,
                 (double)(nr*ntet*nzc*3*sizeof(double))/1048576.0,
                 tc, (double)NS/tc/1.0e6 );
    check( std::isfinite( acc ), "cyl map: sampled values are finite" );

    // ---------------- 笛卡尔场图 ----------------
    const std::size_t nx = 41, ny = 41, nzf = 41;
    const double h = 0.025;
    CFieldMap3D fm( nx, ny, nzf, -0.5, h, -0.5, h, -0.5, h );
    for( std::size_t i = 0; i < nx; ++i ) {
        double x = -0.5 + h*(double)i;
        for( std::size_t j = 0; j < ny; ++j ) {
            double y = -0.5 + h*(double)j;
            for( std::size_t k = 0; k < nzf; ++k ) {
                double z = -0.5 + h*(double)k;
                fm.set_value( i, j, k, -0.05*y, 0.05*x, 0.1*(1.0 - 0.2*z) );
            }
        }
    }

    // 节点处精确恢复（原点恰为 (20,20,20) 节点）
    {
        Vec3D node = fm.node_value( 20, 20, 20 );
        Vec3D got  = fm( Vec3D( 0.0, 0.0, 0.0 ) );
        check( std::fabs( got[2] - node[2] ) < 1e-12,
               "cartesian map: interpolation reproduces node value exactly" );
    }

    std::vector<Vec3D> qs = make_samples( NS, 0.20, 0.35, 20.0, 0.05 );

    acc = 0.0;
    auto eval_cart = [&]() {
        for( long i = 0; i < NS; ++i )
            acc += fm( qs[(std::size_t)i] )[2];
    };
    double tf = time_best( eval_cart );
    std::printf( "CFieldMap3D    %zux%zux%zu (%5.1f MB) : %7.3f s   %8.2f Mevals/s\n",
                 nx, ny, nzf,
                 (double)(nx*ny*nzf*3*sizeof(double))/1048576.0,
                 tf, (double)NS/tf/1.0e6 );
    check( std::isfinite( acc ), "cartesian map: sampled values are finite" );

    std::printf( "\n%s (%d failure%s)\n",
                 g_failures ? "FAILED" : "ALL TESTS PASSED",
                 g_failures, g_failures == 1 ? "" : "s" );
    return( g_failures ? 1 : 0 );
}
