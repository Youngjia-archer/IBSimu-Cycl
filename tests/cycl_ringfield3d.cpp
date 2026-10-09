/** @file cycl_ringfield3d.cpp
 *  @brief IBSimu-Cycl P2-A 测试：由中平面 Bz 重建的三维回旋加速器磁场。
 *
 *  验证（基于真实 PSI Ring 磁场图）：
 *   1. 中平面 z=0: Br=Btheta=0 且 Bz 等于中平面值；
 *   2. 旋转对称: |B| 在 theta 与 theta+45deg 处相等；
 *   3. 真空 Maxwell: 散度 nabla.B ~ 0、旋度 nabla x B ~ 0（离面处，内点）；
 *   4. 小 z 渐近: Br ~ (dBz/dr) z。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <cstdio>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "ringfield3d.hpp"

using namespace ibsimu_cycl;

static int g_failures = 0;

static void check( bool ok, const std::string &msg )
{
    std::printf( "%s  %s\n", ok ? "[ ok ]" : "[FAIL]", msg.c_str() );
    if( !ok )
        ++g_failures;
}

static Vec3D at( const CRingFieldMap3D &f, double r, double th, double z )
{
    return( f( Vec3D( r*std::cos(th), r*std::sin(th), z ) ) );
}

int main( int argc, char **argv )
{
    std::vector<std::string> candidates;
    if( argc > 1 ) candidates.push_back( argv[1] );
    candidates.push_back( "../examples/cyclotron/data/bfield.dat" ); // CWD = tests/
    candidates.push_back( "examples/cyclotron/data/bfield.dat" );    // CWD = repo root

    std::unique_ptr<CRingFieldMap3D> fp;
    std::string file;
    for( std::size_t c = 0; c < candidates.size() && !fp; ++c ) {
        try {
            fp.reset( new CRingFieldMap3D(
                        read_ring_field_map3d( candidates[c], 0.1, 8, 0 ) ) );
            file = candidates[c];
        } catch( const std::exception & ) { }
    }
    if( !fp ) {
        std::printf( "[SKIP] field map not found (tried %d candidates)\n", (int)candidates.size() );
        return( 0 );
    }
    std::printf( "reading: %s\n", file.c_str() );
    const CRingFieldMap3D &B = *fp;

    // 选取一个内点
    const double r  = B.r0() + 70.0*B.dr();
    const double th = M_PI/2.0;         // 90 deg
    const double z  = 0.02;             // 2 cm 离面

    std::printf( "grid: nr=%zu nt=%zu  r0=%.4f m dr=%.4f m\n",
                 B.size_r(), B.size_theta(), B.r0(), B.dr() );

    // ---- 1. 中平面 ----
    {
        Vec3D b0 = at( B, r, th, 0.0 );
        std::size_t i = (std::size_t)std::lround( (r - B.r0())/B.dr() );
        std::size_t k = (std::size_t)std::lround( th/B.dtheta() ) % B.size_theta();
        std::printf( "midplane: Bz=%.6f T  (node=%.6f T)  Br=%.3e Bt=%.3e\n",
                     b0[2], B.bz_midplane( i, k ), b0[0], b0[1] );
        check( std::fabs( b0[2] - B.bz_midplane( i, k ) ) < 1e-9, "z=0: Bz equals midplane node value" );
        check( std::fabs( b0[0] ) < 1e-12 && std::fabs( b0[1] ) < 1e-12, "z=0: Br=Btheta=0" );
    }

    // ---- 2. 旋转对称 (45 deg) ----
    {
        Vec3D b1 = at( B, r, th, z );
        Vec3D b2 = at( B, r, th + M_PI/4.0, z );
        double m1 = std::sqrt( b1[0]*b1[0] + b1[1]*b1[1] + b1[2]*b1[2] );
        double m2 = std::sqrt( b2[0]*b2[0] + b2[1]*b2[1] + b2[2]*b2[2] );
        std::printf( "|B|(theta)=%.8f  |B|(theta+45deg)=%.8f\n", m1, m2 );
        check( std::fabs( m1 - m2 ) < 1e-6*std::max(1.0, m1), "|B| invariant under 45 deg rotation" );
    }

    // ---- 3. 真空 Maxwell 方程 ----
    {
        Vec3D c( r*std::cos(th), r*std::sin(th), z );
        double h = 1.0e-3;
        Vec3D bx1 = B( Vec3D(c[0]+h, c[1], c[2]) ), bx0 = B( Vec3D(c[0]-h, c[1], c[2]) );
        Vec3D by1 = B( Vec3D(c[0], c[1]+h, c[2]) ), by0 = B( Vec3D(c[0], c[1]-h, c[2]) );
        Vec3D bz1 = B( Vec3D(c[0], c[1], c[2]+h) ), bz0 = B( Vec3D(c[0], c[1], c[2]-h) );

        double div = ( (bx1[0]-bx0[0]) + (by1[1]-by0[1]) + (bz1[2]-bz0[2]) )/(2.0*h);

        double curlx = ( (by1[2]-by0[2]) - (bz1[1]-bz0[1]) )/(2.0*h);
        double curly = ( (bz1[0]-bz0[0]) - (bx1[2]-bx0[2]) )/(2.0*h);
        double curlz = ( (bx1[1]-bx0[1]) - (by1[0]-by0[0]) )/(2.0*h);
        double curl  = std::sqrt( curlx*curlx + curly*curly + curlz*curlz );

        Vec3D b = B( c );
        double bmag = std::sqrt( b[0]*b[0] + b[1]*b[1] + b[2]*b[2] );
        double scale = bmag/r;   // 典型的场梯度量级 T/m
        std::printf( "nabla.B=%.4e  |nabla x B|=%.4e  (|B|/r=%.4f T/m)\n", div, curl, scale );
        check( std::fabs(div) < 1e-2*scale, "vacuum: |nabla.B| << |B|/r" );
        check( curl < 1e-2*scale, "vacuum: |nabla x B| << |B|/r" );
    }

    // ---- 4. 小 z 渐近: Br ~ (dBz/dr) z ----
    {
        // theta=90deg 时径向单位矢量为 (0,1)，故 By 即 Br
        Vec3D bp = at( B, r, th, z );
        double slope = bp[1]/z;                       // Br(z)/z
        double b_rm = at( B, r-B.dr(), th, 0.0 )[2];
        double b_rp = at( B, r+B.dr(), th, 0.0 )[2];
        double dbr_an = (b_rp - b_rm)/(2.0*B.dr());   // dBz/dr (中平面中心差分)
        std::printf( "Br(z)/z = %.6f T/m   dBz/dr = %.6f T/m\n", slope, dbr_an );
        check( std::fabs( slope - dbr_an ) < 0.05*std::fabs(dbr_an) + 1e-9,
               "small-z: Br/z ~= dBz/dr (within 5%)" );
    }

    std::printf( "\n%s (%d failure%s)\n",
                 g_failures ? "FAILED" : "ALL TESTS PASSED",
                 g_failures, g_failures == 1 ? "" : "s" );
    return( g_failures ? 1 : 0 );
}
