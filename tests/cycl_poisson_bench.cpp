/** @file cycl_poisson_bench.cpp
 *  @brief IBSimu-Cycl 基准：三维 Poisson/空间电荷求解器的耗时与内存。
 *
 *  为什么要做这个基准？
 *  IBSimu 自带的空间电荷求解器（Gauss-Seidel / BiCGSTAB+ILU0 / 多重网格）
 *  **全部是串行单进程**实现。是否值得引入 MPI + hypre 做区域分解，取决于两件事：
 *    1. 求解时间随网格规模的**增长趋势**（若按 N^3~N^4 爆炸，优先换算法/加速器
 *       而不是加节点）；
 *    2. **内存**是否超过单节点（这才是 MPI 唯一不可替代的场景）。
 *
 *  算例：单位立方腔（六面接地 Dirichlet）+ 中心高斯电荷团。
 *  依次运行多重网格与 BiCGSTAB+ILU0，报告耗时、迭代数/循环数与常驻内存。
 *
 *  用法: cycl_poisson_bench [size ...]
 *    不带参数时只跑 65（CI 友好）；手动可传 "129 257" 等做规模扫描。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <cmath>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "diag_precond.hpp"
#include "epot_bicgstabsolver.hpp"
#include "epot_field.hpp"
#include "epot_mgsolver.hpp"
#include "geometry.hpp"
#include "ilu0_precond.hpp"
#include "meshscalarfield.hpp"

using namespace std;

static double now_sec()
{
    struct timespec ts;
    clock_gettime( CLOCK_MONOTONIC, &ts );
    return( (double)ts.tv_sec + 1.0e-9*(double)ts.tv_nsec );
}

/*! \brief 当前进程的常驻内存峰值 [MB]（读 /proc/self/status 的 VmHWM）。 */
static double peak_rss_mb()
{
    ifstream f( "/proc/self/status" );
    string line;
    while( getline( f, line ) ) {
        if( line.compare( 0, 6, "VmHWM:" ) == 0 ) {
            double kb = 0.0;
            istringstream ss( line.substr(6) );
            ss >> kb;
            return( kb/1024.0 );
        }
    }
    return( 0.0 );
}

static double max_abs_phi( const MeshScalarField &epot, const Geometry &g )
{
    double m = 0.0;
    for( uint32_t i = 0; i < g.size(0); ++i )
        for( uint32_t j = 0; j < g.size(1); ++j )
            for( uint32_t k = 0; k < g.size(2); ++k ) {
                double v = fabs( epot(i,j,k) );
                if( v > m ) m = v;
            }
    return( m );
}

static bool g_sharp = false;   /*!< \brief 是否用阶跃（陡边界）源。 */

static void bench( int size )
{
    const double L = 1.0;                    // 立方腔边长 [m]
    const double h = L/(double)(size-1);

    Geometry g( MODE_3D, Int3D(size,size,size), Vec3D(0,0,0), h );
    for( int b = 1; b <= 6; ++b )
        g.set_boundary( b, Bound(BOUND_DIRICHLET, 0.0) );
    g.build_mesh();

    EpotField       epot( g );
    MeshScalarField scharge( g );

    // 电荷分布：高斯团（光滑）或阶跃球（陡边界，考验求解器）
    const double qmax  = 1.0e-4;             // C/m^3
    const double sigma = 0.1;
    for( uint32_t i = 0; i < g.size(0); ++i ) {
        double x = i*g.h() - 0.5;
        for( uint32_t j = 0; j < g.size(1); ++j ) {
            double y = j*g.h() - 0.5;
            for( uint32_t k = 0; k < g.size(2); ++k ) {
                double z = k*g.h() - 0.5;
                double r2 = x*x + y*y + z*z;
                if( g_sharp )
                    scharge(i,j,k) = ( r2 < 0.25*0.25 ? qmax : 0.0 );
                else
                    scharge(i,j,k) = qmax*exp( -r2/(2.0*sigma*sigma) );
            }
        }
    }

    const double nodes = (double)size*(double)size*(double)size;
    printf( "\n=== %d^3 = %.2f M nodes   h = %.4g m   source=%s ===\n",
            size, nodes/1.0e6, h, g_sharp ? "step" : "gaussian" );
    printf( "  baseline RSS: %.0f MB\n", peak_rss_mb() );

    // ---- 几何多重网格 ----
    {
        EpotMGSolver mg( g );
        mg.set_levels( 4 );
        mg.set_eps( 1.0e-6 );
        mg.set_imax( 200 );

        double t0 = now_sec();
        mg.solve( epot, scharge );
        double dt = now_sec() - t0;

        printf( "  multigrid (levels=4) : %9.3f s   mgcyc=%-6u max|phi|=%.4g V   RSS=%.0f MB\n",
                dt, mg.get_mgcyc(), max_abs_phi(epot,g), peak_rss_mb() );
        printf( "    -> %.2f us / node / solve\n", 1.0e6*dt/nodes );
    }

    // ---- BiCGSTAB + ILU0 ----
    {
        EpotBiCGSTABSolver bc( g );
        ILU0_Precond pc;
        bc.set_preconditioner( pc );
        bc.set_eps( 1.0e-6 );
        bc.set_imax( 20000 );

        double t0 = now_sec();
        bc.solve( epot, scharge );
        double dt = now_sec() - t0;

        printf( "  BiCGSTAB + ILU0      : %9.3f s   iter=%-7u max|phi|=%.4g V   RSS=%.0f MB\n",
                dt, bc.get_iter(), max_abs_phi(epot,g), peak_rss_mb() );
        printf( "    -> %.2f us / node / solve\n", 1.0e6*dt/nodes );
    }
}

int main( int argc, char **argv )
{
    vector<int> sizes;
    for( int a = 1; a < argc; ++a ) {
        if( string(argv[a]) == "--sharp" )
            g_sharp = true;
        else
            sizes.push_back( atoi( argv[a] ) );
    }
    if( sizes.empty() )
        sizes.push_back( 65 );               // CI 友好的默认值

    for( size_t i = 0; i < sizes.size(); ++i )
        bench( sizes[i] );

    printf( "\nbenchmark complete\n" );
    return( 0 );
}
