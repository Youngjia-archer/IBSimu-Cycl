/** @file cycl_omp_tracker.cpp
 *  @brief IBSimu-Cycl G1 测试：粒子系综的 OpenMP 并行化与性能基准。
 *
 *  验证：
 *   1. **正确性**：1 线程与 N 线程跟踪结果逐位一致（粒子间独立，无竞争）；
 *   2. **加速比**：报告单线程 vs 多线程的实测加速比（> 1 时说明并行有效）；
 *   3. 附带吞吐量（粒子·步/秒）。
 *
 *  使用真实 PSI Ring 场图 + 相对论推进器，负载接近实际（场插值是主要开销）。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <memory>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ensembletracker.hpp"
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

static double now_sec()
{
#ifdef _OPENMP
    return( omp_get_wtime() );
#else
    return( (double)std::clock()/(double)CLOCKS_PER_SEC );
#endif
}

struct BenchResult {
    double t_serial, t_par, speedup, maxdiff;
};

/*! \brief 在同一系综上测单线程/多线程跟踪时间与一致性。 */
static BenchResult bench( const VectorField *B,
                          const std::vector<EnsembleParticle> &p0,
                          double dt, int M, int nthreads, int reps = 3 )
{
    BenchResult r;
    CEnsembleTracker tracker( B, 0, QE, MP );
    tracker.set_relativistic( true );

    {   // 预热
        std::vector<EnsembleParticle> warm = p0;
        tracker.track( warm, dt, 50 );
    }

    // 重复多轮取最优，消除瞬时系统负载造成的噪声
    std::vector<EnsembleParticle> a0 = p0, b0 = p0;
    r.t_serial = 1.0e30;
    r.t_par    = 1.0e30;
    for( int rep = 0; rep < reps; ++rep ) {
        std::vector<EnsembleParticle> a = p0;
#ifdef _OPENMP
        omp_set_num_threads( 1 );
#endif
        double t0 = now_sec();
        tracker.track( a, dt, M );
        double ts = now_sec() - t0;
        if( ts < r.t_serial ) { r.t_serial = ts; a0 = a; }

        std::vector<EnsembleParticle> b = p0;
#ifdef _OPENMP
        omp_set_num_threads( nthreads );
#endif
        t0 = now_sec();
        tracker.track( b, dt, M );
        double tp = now_sec() - t0;
        if( tp < r.t_par ) { r.t_par = tp; b0 = b; }
    }

    r.speedup = r.t_serial/r.t_par;
    r.maxdiff = 0.0;
    for( std::size_t i = 0; i < p0.size(); ++i )
        for( int c = 0; c < 3; ++c ) {
            r.maxdiff = std::max( r.maxdiff, std::fabs( a0[i].x[c] - b0[i].x[c] ) );
            r.maxdiff = std::max( r.maxdiff, std::fabs( a0[i].v[c] - b0[i].v[c] ) );
        }
    return( r );
}

/*! \brief 简单确定性 LCG（避免依赖随机库，保证可复现）。 */
static double lcg( unsigned long &s )
{
    s = s*6364136223846793005UL + 1442695040888963407UL;
    return( ((double)((s >> 11) & 0x1FFFFFFFFFFFFFUL))/9007199254740992.0 );  // [0,1)
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
        std::printf( "[SKIP] field map not found; OpenMP benchmark needs the real map\n" );
        return( 0 );
    }

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    std::printf( "OpenMP max threads = %d\n", nthreads );

    // ---- 参考轨道参数（与 cycl_relativistic 一致）----
    const double r_ref = 3.3;
    std::size_t i_ref = (std::size_t)std::lround( (r_ref - map->r0())/map->dr() );
    if( i_ref >= map->size_r() ) i_ref = map->size_r()-1;
    double Bbar = 0.0;
    for( std::size_t k = 0; k < map->size_theta(); ++k )
        Bbar += map->bz_midplane( i_ref, k );
    Bbar /= (double)map->size_theta();

    double p = QE*Bbar*r_ref;
    double gb = p/(MP*CL);
    double gam = std::sqrt(1.0 + gb*gb);
    double bet = gb/gam;
    double v0 = bet*CL;

    // ---- 粒子系综：参考轨道 + 小扰动 ----
    const std::size_t N = 2000;
    const int M = 400;                      // 基准步数
    const double omega = QE*Bbar/(gam*MP);
    const double T_rev = 2.0*M_PI/omega;
    const double dt = T_rev/4000.0;

    std::vector<EnsembleParticle> p0( N );
    unsigned long seed = 12345UL;
    for( std::size_t i = 0; i < N; ++i ) {
        double dr   = (lcg(seed) - 0.5)*0.02;         // +/- 10 mm
        double dth  = (lcg(seed) - 0.5)*2.0e-3;       // 角度扰动
        double dvf  = (lcg(seed) - 0.5)*2.0e-3;       // 相对速度扰动
        double th   = M_PI/2.0 + dth;
        double r    = r_ref + dr;
        p0[i].x = Vec3D( r*std::cos(th), r*std::sin(th), (lcg(seed)-0.5)*0.01 );
        double v = v0*(1.0 + dvf);
        p0[i].v = Vec3D( -std::sin(th)*v, std::cos(th)*v, 0.0 );
    }
    std::printf( "ensemble: N=%zu  steps=%d  dt=%.3e s  gamma=%.4f\n", N, M, dt, gam );

    // ---- OpenMP 运行时自检：独立的平凡并行循环 ----
    double sanity_ratio = 1.0;
    {
        const long NB = 30000000;
        auto busy = []( long N ) {
            double s = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:s) schedule(static)
#endif
            for( long i = 0; i < N; ++i ) s += std::sqrt( (double)(i % 1000) );
            return( s );
        };
#ifdef _OPENMP
        omp_set_num_threads( 1 );
#endif
        double ts = now_sec(); double r1 = busy( NB ); ts = now_sec() - ts;
#ifdef _OPENMP
        omp_set_num_threads( nthreads );
#endif
        double tp = now_sec(); double r2 = busy( NB ); tp = now_sec() - tp;
        sanity_ratio = ts/tp;
        std::printf( "openmp sanity: serial=%.3f s  parallel=%.3f s  ratio=%.2f  (%s)\n",
                     ts, tp, sanity_ratio,
                     std::fabs(r1 - r2) < 1e-6*std::fabs(r1) ? "checksums match" : "CHECKSUMS DIFFER" );
    }

    // ---- 基准 1：小场数据（均匀场）—— 隔离并行开销/竞争 ----
    CFieldMap3D Buniform( 2, 2, 2, -5.0, 10.0, -5.0, 10.0, -5.0, 10.0 );
    for( std::size_t i = 0; i < 2; ++i )
        for( std::size_t j = 0; j < 2; ++j )
            for( std::size_t k = 0; k < 2; ++k )
                Buniform.set_value( i, j, k, 0.0, 0.0, Bbar );

    BenchResult r_uni = bench( &Buniform, p0, dt, M, nthreads );
    std::printf( "uniform field : serial=%.3f s  parallel=%.3f s  speedup=%.2fx  max|diff|=%.1e\n",
                 r_uni.t_serial, r_uni.t_par, r_uni.speedup, r_uni.maxdiff );

    BenchResult r_map = bench( map.get(), p0, dt, M, nthreads );
    std::printf( "real PSI map  : serial=%.3f s  parallel=%.3f s  speedup=%.2fx  max|diff|=%.1e\n",
                 r_map.t_serial, r_map.t_par, r_map.speedup, r_map.maxdiff );
    std::printf( "single-thread rate (real map): %.2e particle-steps/s\n",
                 (double)N*(double)M/r_map.t_serial );

    check( r_uni.maxdiff == 0.0, "uniform field: 1-thread and N-thread results bitwise identical" );
    check( r_map.maxdiff == 0.0, "real field map: 1-thread and N-thread results bitwise identical" );
    if( nthreads > 1 && sanity_ratio >= 1.5 ) {
        check( r_uni.speedup > 1.5, "uniform field: parallel speedup > 1.5x" );
        std::printf( "[info] real PSI map speedup = %.2fx "
                     "(field interpolation is memory/latency bound)\n", r_map.speedup );
    } else {
        std::printf( "[info] environment cannot parallelize (threads=%d, sanity=%.2fx); "
                     "skipping speedup assertion\n", nthreads, sanity_ratio );
    }

    std::printf( "\n%s (%d failure%s)\n",
                 g_failures ? "FAILED" : "ALL TESTS PASSED",
                 g_failures, g_failures == 1 ? "" : "s" );
    return( g_failures ? 1 : 0 );
}