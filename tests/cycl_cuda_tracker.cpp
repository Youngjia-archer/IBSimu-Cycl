/** @file cycl_cuda_tracker.cpp
 *  @brief IBSimu-Cycl R1-③ 测试：GPU (CUDA) 粒子系综跟踪。
 *
 *  验证：
 *   1. **正确性**：GPU 与 CPU（CEnsembleTracker，同一算法）在真实 PSI Ring
 *      场图上的结果一致——短程（10 步）用紧公差，长程（400 步）用松公差；
 *   2. **γ 守恒**：相对论静磁场下 GPU 端的 γ 漂移应在机器精度量级；
 *   3. **性能**：与 CPU 单线程 / OpenMP 全线程对比，并分别报告 kernel 时间与
 *      数据传输时间（GPU 的瓶颈常在后者）。
 *
 *  无 CUDA 支持或无可用 GPU 时输出 [SKIP] 并返回成功，因此 CI（无 GPU）不受影响。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <memory>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ensembletracker.hpp"
#include "gputracker.hpp"
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

static double now_sec()
{
#ifdef _OPENMP
    return( omp_get_wtime() );
#else
    return( (double)std::clock()/(double)CLOCKS_PER_SEC );
#endif
}

/*! \brief 简单确定性 LCG（保证可复现）。 */
static double lcg( unsigned long &s )
{
    s = s*6364136223846793005UL + 1442695040888963407UL;
    return( ((double)((s >> 11) & 0x1FFFFFFFFFFFFFUL))/9007199254740992.0 );
}

/*! \brief GPU 与 CPU 结果的差异（位置绝对值 / 速度绝对值与相对值）。 */
struct DiffResult {
    double dx, dv, dv_rel;
};

static DiffResult compare( const std::vector<EnsembleParticle> &a,
                           const std::vector<EnsembleParticle> &b,
                           double v_scale )
{
    DiffResult r;
    r.dx = 0.0;
    r.dv = 0.0;
    for( std::size_t i = 0; i < a.size(); ++i )
        for( int c = 0; c < 3; ++c ) {
            r.dx = std::max( r.dx, std::fabs( a[i].x[c] - b[i].x[c] ) );
            r.dv = std::max( r.dv, std::fabs( a[i].v[c] - b[i].v[c] ) );
        }
    r.dv_rel = r.dv/v_scale;
    return( r );
}

static double gamma_of( const Vec3D &v )
{
    double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    return( 1.0/std::sqrt( 1.0 - v2/(CL*CL) ) );
}

int main( int argc, char **argv )
{
    std::vector<std::string> cands;
    if( argc > 1 ) cands.push_back( argv[1] );
    cands.push_back( "../examples/cyclotron/data/bfield.dat" );
    cands.push_back( "examples/cyclotron/data/bfield.dat" );

    std::unique_ptr<CRingFieldMap3D> map;
    for( std::size_t c = 0; c < cands.size() && !map; ++c ) {
        try {
            map.reset( new CRingFieldMap3D( read_ring_field_map3d( cands[c], 0.1, 8, 0 ) ) );
        } catch( const std::exception & ) { }
    }
    if( !map ) {
        std::printf( "[SKIP] field map not found; CUDA test needs the real map\n" );
        return( 0 );
    }

    std::printf( "CUDA compiled in : %s\n", CGpuEnsembleTracker::compiled_in() ? "yes" : "no" );
    std::string reason;
    if( !CGpuEnsembleTracker::available( &reason ) ) {
        std::printf( "[SKIP] GPU unavailable: %s\n", reason.c_str() );
        return( 0 );
    }
    std::printf( "CUDA device      : %s\n", CGpuEnsembleTracker::device_name().c_str() );

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    std::printf( "OpenMP max threads = %d\n", nthreads );

    // ---- 参考轨道（与 cycl_omp_tracker / cycl_relativistic 一致）----
    const double r_ref = 3.3;
    std::size_t i_ref = (std::size_t)std::lround( (r_ref - map->r0())/map->dr() );
    if( i_ref >= map->size_r() ) i_ref = map->size_r()-1;
    double Bbar = 0.0;
    for( std::size_t k = 0; k < map->size_theta(); ++k )
        Bbar += map->bz_midplane( i_ref, k );
    Bbar /= (double)map->size_theta();

    double p    = QE*Bbar*r_ref;
    double gb   = p/(MP*CL);
    double gam  = std::sqrt( 1.0 + gb*gb );
    double bet  = gb/gam;
    double v0   = bet*CL;
    double T_c  = 2.0*M_PI*gam*MP/(QE*Bbar);
    double dt   = T_c/1400.0;
    std::printf( "orbit: r=%.2f m  <B>=%.4f T  gamma=%.4f  beta=%.4f  dt=%.3e s\n",
                 r_ref, Bbar, gam, bet, dt );

    // ---- 粒子初始化 ----
    const std::size_t N = 100000;
    std::vector<EnsembleParticle> p0( N );
    unsigned long seed = 20261009UL;
    for( std::size_t i = 0; i < N; ++i ) {
        double dr   = (lcg(seed) - 0.5)*0.02;      // +/- 10 mm
        double dth  = (lcg(seed) - 0.5)*2.0e-3;
        double dvf  = (lcg(seed) - 0.5)*2.0e-3;
        double th   = M_PI/2.0 + dth;
        double r    = r_ref + dr;
        p0[i].x = Vec3D( r*std::cos(th), r*std::sin(th), (lcg(seed)-0.5)*0.01 );
        double v = v0*(1.0 + dvf);
        p0[i].v = Vec3D( -std::sin(th)*v, std::cos(th)*v, 0.0 );
    }

    CEnsembleTracker cpu_ref( map.get(), 0, QE, MP );
    cpu_ref.set_relativistic( true );

    CGpuEnsembleTracker gpu( *map, QE, MP );
    gpu.set_relativistic( true );

    // ================= 1. 正确性：短程 =================
    {
        const int M = 10;
        std::vector<EnsembleParticle> init( p0.begin(), p0.begin()+2000 );
        std::vector<EnsembleParticle> a = init;
        std::vector<EnsembleParticle> b = init;

#ifdef _OPENMP
        omp_set_num_threads( 1 );
#endif
        cpu_ref.track( a, dt, M );
        gpu.track( b, dt, M );

        DiffResult d = compare( a, b, v0 );
        std::printf( "short run (%d steps): dx=%.3e m  dv/v=%.3e\n", M, d.dx, d.dv_rel );
        check( d.dv_rel < 1e-14 && d.dx < 1e-8,
               "GPU matches CPU after 10 steps (dv/v < 1e-14, dx < 1e-8 m)" );
    }

    // ================= 2. 正确性：长程 + gamma 守恒 =================
    {
        const int M = 400;
        std::vector<EnsembleParticle> init( p0.begin(), p0.begin()+2000 );
        std::vector<EnsembleParticle> a = init;
        std::vector<EnsembleParticle> b = init;

#ifdef _OPENMP
        omp_set_num_threads( 1 );
#endif
        cpu_ref.track( a, dt, M );
        gpu.track( b, dt, M );

        DiffResult d = compare( a, b, v0 );
        std::printf( "long run  (%d steps): dx=%.3e m  dv/v=%.3e\n", M, d.dx, d.dv_rel );
        check( d.dv_rel < 1e-12 && d.dx < 1e-5,
               "GPU matches CPU after 400 steps (dv/v < 1e-12, dx < 1e-5 m)" );

        // gamma 守恒：逐粒子比较初始与末态（而非系综内的分布宽度）
        double dg = 0.0, dgc = 0.0;
        for( std::size_t i = 0; i < b.size(); ++i ) {
            double g0 = gamma_of( init[i].v );
            dg  = std::max( dg,  std::fabs( gamma_of( b[i].v ) - g0 )/g0 );
            dgc = std::max( dgc, std::fabs( gamma_of( a[i].v ) - g0 )/g0 );
        }
        std::printf( "gamma drift over %d steps: GPU=%.3e  CPU=%.3e (relative)\n", M, dg, dgc );
        check( dg < 1e-12, "GPU conserves gamma in static B (relative drift < 1e-12)" );
    }

    // ================= 3. 性能 =================
    const int M = 200;
    std::printf( "\nperformance: N=%zu  steps=%d\n", N, M );

    double t_cpu1 = 0.0, t_cpup = 0.0, t_gpu = 0.0;
    {
        std::vector<EnsembleParticle> a = p0;
#ifdef _OPENMP
        omp_set_num_threads( 1 );
#endif
        double t0 = now_sec(); cpu_ref.track( a, dt, M ); t_cpu1 = now_sec() - t0;

        std::vector<EnsembleParticle> b = p0;
#ifdef _OPENMP
        omp_set_num_threads( nthreads );
#endif
        t0 = now_sec(); cpu_ref.track( b, dt, M ); t_cpup = now_sec() - t0;

        std::vector<EnsembleParticle> c = p0;
        gpu.track( c, dt, M );
        t_gpu = 0.0;
        // 取 3 轮最优（首轮含 NVRTC 之后的首次启动开销）
        for( int r = 0; r < 3; ++r ) {
            std::vector<EnsembleParticle> d = p0;
            t0 = now_sec(); gpu.track( d, dt, M );
            double dtg = now_sec() - t0;
            if( r == 0 || dtg < t_gpu ) t_gpu = dtg;
        }
    }

    const double steps_total = (double)N*(double)M;
    std::printf( "  CPU  1 thread : %7.3f s   %8.2f Mev/s\n", t_cpu1, steps_total/t_cpu1/1.0e6 );
    std::printf( "  CPU %2d threads: %7.3f s   %8.2f Mev/s   (%.2fx vs 1 thread)\n",
                 nthreads, t_cpup, steps_total/t_cpup/1.0e6, t_cpu1/t_cpup );
    std::printf( "  GPU           : %7.3f s   %8.2f Mev/s   (%.2fx vs CPU-1, %.2fx vs CPU-N)\n",
                 t_gpu, steps_total/t_gpu/1.0e6, t_cpu1/t_gpu, t_cpup/t_gpu );
    std::printf( "  GPU breakdown : kernel %.3f s   transfer %.3f s (%.1f%%)\n",
                 gpu.last_kernel_time(), gpu.last_transfer_time(),
                 100.0*gpu.last_transfer_time()/std::max(t_gpu, 1e-12) );

    check( t_gpu < t_cpu1, "GPU is faster than single-thread CPU" );

    std::printf( "\n%s (%d failure%s)\n",
                 g_failures ? "FAILED" : "ALL TESTS PASSED",
                 g_failures, g_failures == 1 ? "" : "s" );
    return( g_failures ? 1 : 0 );
}
