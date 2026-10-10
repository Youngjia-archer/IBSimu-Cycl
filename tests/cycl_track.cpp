/** @file cycl_track.cpp
 *  @brief IBSimu-Cycl D1 端到端 Demo：真实磁场中的三维单粒子跟踪。
 *
 *  流程:
 *   1. 读取真实 PSI Ring 磁场图 (bfield.dat) -> CRingFieldMap3D（三维场）;
 *   2. 用 Boris 推动器积分质子在中平面磁场中的运动;
 *   3. 校验回旋频率 omega_c = qB/m（测量 vs 解析），并检查 |v| 守恒;
 *   4. 输出**整机中场图** cycl_track_map.vti（ParaView / PyVista 可直接打开）。
 *
 *  注：本算例**不**输出任何单粒子轨迹可视化。曾经额外生成过一条 0.1c 质子的
 *  "局部拉莫尔回旋"轨迹（回旋半径仅 0.2 m）并当作轨道图展示，那不是回旋加速器
 *  轨道，已删除。机器尺度的闭合轨道见 tests/cycl_closed_orbit.cpp。
 *
 *  用法: cycl_track [bfield.dat] [输出前缀(默认 cycl_track)]
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#include "ringfield3d.hpp"
#include "vtkwriter.hpp"

using namespace ibsimu_cycl;

// 物理常数 (SI)
static const double QE  = 1.602176634e-19;   // 元电荷 [C]
static const double MP  = 1.67262192369e-27; // 质子质量 [kg]

static int g_failures = 0;

static void check( bool ok, const std::string &msg )
{
    std::printf( "%s  %s\n", ok ? "[ ok ]" : "[FAIL]", msg.c_str() );
    if( !ok ) ++g_failures;
}

static Vec3D crossp( const Vec3D &a, const Vec3D &b )
{
    return( Vec3D( a[1]*b[2] - a[2]*b[1],
                   a[2]*b[0] - a[0]*b[2],
                   a[0]*b[1] - a[1]*b[0] ) );
}

static double vabs( const Vec3D &a )
{
    return( std::sqrt( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] ) );
}

/*! \brief 一步 Boris 推进（静磁场, E=0, 非相对论）。
 *
 *  t = (q/m) B dt/2;  v' = v + v x t;  v_new = v + v' x (2t/(1+|t|^2))
 */
static void boris_step( const CRingFieldMap3D &B, double qm, Vec3D &x, Vec3D &v, double dt )
{
    Vec3D b = B( x );
    Vec3D t = b * (qm*dt/2.0);
    double t2 = t[0]*t[0] + t[1]*t[1] + t[2]*t[2];
    Vec3D vprime = v + crossp( v, t );
    Vec3D s = t * (2.0/(1.0 + t2));
    v = v + crossp( vprime, s );
    x = x + v*dt;
}

struct TrajPoint {
    Vec3D x, v;
    double bmag;
};

int main( int argc, char **argv ){
    // ---- 载入磁场图 ----
    std::vector<std::string> candidates;
    if( argc > 1 ) candidates.push_back( argv[1] );
    candidates.push_back( "../examples/cyclotron/data/bfield.dat" ); // CWD = tests/
    candidates.push_back( "examples/cyclotron/data/bfield.dat" );    // CWD = repo root

    std::unique_ptr<CRingFieldMap3D> fp;
    std::string file;
    for( std::size_t c = 0; c < candidates.size() && !fp; ++c ) {
        try {
            fp.reset( new CRingFieldMap3D( read_ring_field_map3d( candidates[c], 0.1, 8, 0 ) ) );
            file = candidates[c];
        } catch( const std::exception & ) { }
    }
    if( !fp ) {
        std::printf( "[SKIP] field map not found (tried %d candidates)\n", (int)candidates.size() );
        return( 0 );
    }
    const CRingFieldMap3D &B = *fp;
    const std::string prefix = (argc > 2) ? argv[2] : "cycl_track";

    std::printf( "field map: %s\n", file.c_str() );

    const double qm = QE/MP;   // 比荷 [C/kg]

    // =====================================================================
    // 运行 A：回旋频率精度测试（小回旋半径，场近似均匀）
    // =====================================================================
    double r0 = B.r0() + 70.0*B.dr();     // 3.10 m
    double th0 = M_PI/2.0;
    Vec3D x( r0*std::cos(th0), r0*std::sin(th0), 0.0 );

    Vec3D b0 = B( x );
    double Bmag0 = vabs( b0 );
    double omega_c = QE*Bmag0/MP;         // 解析回旋角频率
    double Tc = 2.0*M_PI/omega_c;

    // 半径方向单位矢量 phi_hat = (-sin, cos, 0) -> 纯方位角初速
    double vmag = 1.0e6;                  // 1e6 m/s -> 回旋半径约 mm 量级
    Vec3D v( -std::sin(th0)*vmag, std::cos(th0)*vmag, 0.0 );

    double dt = Tc/400.0;
    int    nsteps = 800;                  // ~2 个回旋周期

    double phi_prev = std::atan2( v[1], v[0] );
    double phi_sum  = 0.0;
    double v0mag    = vabs( v );
    double emax_rel = 0.0;

    for( int n = 0; n < nsteps; ++n ) {
        boris_step( B, qm, x, v, dt );
        double phi = std::atan2( v[1], v[0] );
        double d = phi - phi_prev;
        while( d >  M_PI ) d -= 2.0*M_PI;
        while( d < -M_PI ) d += 2.0*M_PI;
        phi_sum += d;
        phi_prev = phi;

        double rel = std::fabs( vabs(v) - v0mag )/v0mag;
        if( rel > emax_rel ) emax_rel = rel;
    }

    double omega_meas = std::fabs( phi_sum )/( nsteps*dt );
    double err = std::fabs( omega_meas - omega_c )/omega_c;

    std::printf( "B(r=%.3f m) = %.6f T\n", r0, Bmag0 );
    std::printf( "omega_c (analytic) = %.6e rad/s   (T_c = %.3e s)\n", omega_c, Tc );
    std::printf( "omega   (measured) = %.6e rad/s   rel.err = %.3e\n", omega_meas, err );
    std::printf( "|v| max relative drift = %.3e\n", emax_rel );

    check( err < 1e-2, "measured cyclotron frequency matches qB/m within 1%" );
    check( emax_rel < 1e-9, "energy conservation in static B field (|v| constant)" );
    check( std::fabs( x[2] ) < 1e-12, "particle stays on midplane (z=0)" );

    // 注：这里**不再**生成任何“局部拉莫尔回旋”可视化轨迹。
    // 0.1c 质子在该场中的回旋半径只有 0.2 m，那是一段局部回旋、不是回旋加速器轨道；
    // 把它当“轨道图”展示过是误导，已删除。真正的机器尺度轨道见
    // tests/cycl_closed_orbit.cpp（求单圈映射不动点的闭合轨道）。

    // --- VTK XML：整机中场图（不再输出任何单粒子轨迹）---
    {
        // 0.1 m 步长足以分辨 8 折扇形结构。供 plot_field_map.py 绘图。
        const int    NM  = 95;
        const double MDR = 0.1;
        const double MX0 = -0.5*(NM-1)*MDR;
        std::vector<double> mmag( (std::size_t)NM*NM );
        std::vector<double> mvec( 3*(std::size_t)NM*NM );
        for( int j = 0; j < NM; ++j )
            for( int i = 0; i < NM; ++i ) {
                Vec3D bb = B( Vec3D( MX0 + MDR*i, MX0 + MDR*j, 0.0 ) );
                std::size_t a = (std::size_t)i + (std::size_t)NM*j;
                mmag[a] = vabs( bb );
                mvec[3*a+0] = bb[0];
                mvec[3*a+1] = bb[1];
                mvec[3*a+2] = bb[2];
            }
        vtk_write_image_data( prefix + "_map.vti",
                              Int3D( NM, NM, 1 ),
                              Vec3D( MX0, MX0, 0.0 ), Vec3D( MDR, MDR, MDR ),
                              "Bmag", mmag, "B", mvec );

        std::printf( "VTK XML: %s_map.vti (%dx%d 整机中场图，不含轨迹)\n",
                     prefix.c_str(), NM, NM );
    }

    std::printf( "\n%s (%d failure%s)\n",
                 g_failures ? "FAILED" : "ALL TESTS PASSED",
                 g_failures, g_failures == 1 ? "" : "s" );
    return( g_failures ? 1 : 0 );
}
