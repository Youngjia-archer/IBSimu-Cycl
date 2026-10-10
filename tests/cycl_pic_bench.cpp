/** @file cycl_pic_bench.cpp
 *  @brief IBSimu-Cycl 基准：自洽 PIC 迭代（trace <-> solve）的时间构成与线程扩展性。
 *
 *  稳态空间电荷求解必须迭代：追踪一代粒子 -> 沉积电荷 -> 求解 Poisson -> 用新场
 *  重新追踪，直到收敛。本基准用负离子源算例（与 tests/nsimp_plasma3d.cpp 同构）
 *  分别计时每代的 ①空间电荷求解 与 ②轨迹追踪，并扫描 IBSimu::set_thread_count()，
 *  以回答两个问题：
 *
 *    - 时间究竟花在 solve（求解）还是 trace（追踪）？
 *    - 轨迹追踪的线程扩展性如何？上游的电荷沉积把 4 次 double 累加放在一把
 *      **全局互斥锁**里、且每条轨迹片段都要加解锁，怀疑是扩展性瓶颈。
 *
 *  用法: cycl_pic_bench [generations] [threads...] [--bicg]
 *    默认 generations=2, threads=1；--bicg 改用未并行的 BiCGSTAB 作为对照。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>

#include "epot_bicgstabsolver.hpp"
#include "epot_efield.hpp"
#include "epot_field.hpp"
#include "epot_mgsolver.hpp"
#include "func_solid.hpp"
#include "geometry.hpp"
#include "ibsimu.hpp"
#include "meshscalarfield.hpp"
#include "meshvectorfield.hpp"
#include "particledatabase.hpp"

using namespace std;

static double now_sec()
{
    struct timespec ts;
    clock_gettime( CLOCK_MONOTONIC, &ts );
    return( (double)ts.tv_sec + 1.0e-9*(double)ts.tv_nsec );
}

static bool solid1( double x, double y, double z )
{
    double r = sqrt( y*y + z*z );
    if( x <= 0.0 && r <= 0.005 )
        return( false );
    return( x <= 0.004 && r >= 0.0030 && r >= 4.8*x - 0.003 );
}

static bool solid2( double x, double y, double z )
{
    double r = sqrt( y*y + z*z );
    return( (x >= 0.025 || r <= 3.25*x - 0.06925) && x >= 0.0231 && r >= 0.0055 );
}

static bool initial_plasma( double x, double, double )
{
    return( x <= 0.00055 );
}

static void run_case( int generations, int threads, bool use_bicg )
{
    Geometry geom( MODE_3D, Int3D(61,81,81), Vec3D(-0.001,-2.0e-2,-2.0e-2), 5.0e-4 );

    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    geom.set_solid( 8, s2 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, +6.0e3) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 5, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 6, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, +6.0e3) );
    geom.build_mesh();

    EpotField       epot( geom );
    MeshScalarField scharge( geom );
    MeshVectorField bfield;
    EpotEfield      efield( epot );
    field_extrpl_e  efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                                      FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                                      FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    InitialPlasma initp( AXIS_X, 0.0006 );

    EpotMGSolver       *mg  = 0;
    EpotBiCGSTABSolver *bcg = 0;
    if( use_bicg ) {
        bcg = new EpotBiCGSTABSolver( geom );
        bcg->set_nsimp_initial_plasma( &initp );
    } else {
        mg = new EpotMGSolver( geom );
        mg->set_levels( 3 );
        mg->set_nsimp_initial_plasma( &initp );
    }

    ParticleDataBase3D pdb( geom );
    bool pmirror[6] = { false, false, true, false, true, false };
    pdb.set_mirror( pmirror );
    pdb.set_polyint( true );

    if( threads > 0 )
        ibsimu.set_thread_count( threads );
    const int eff_threads = (int)ibsimu.get_thread_count();

    double t_solve = 0.0, t_trace = 0.0;
    for( int i = 0; i < generations; ++i ) {

        if( i == 1 ) {
            double rhoneg = fabs( pdb.get_rhosum() );
            vector<double> Ei, rhoi;
            Ei.push_back( 1.0 );
            rhoi.push_back( 0.5*rhoneg );
            if( bcg ) bcg->set_nsimp_plasma( 0.5*rhoneg, 10.0, rhoi, Ei );
            else      mg->set_nsimp_plasma( 0.5*rhoneg, 10.0, rhoi, Ei );
        }

        double t0 = now_sec();
        if( bcg ) bcg->solve( epot, scharge );
        else      mg->solve( epot, scharge );
        t_solve += now_sec() - t0;

        efield.recalculate();

        pdb.clear();
        // H- : 1 mA total -> 35.37 A/m2
        pdb.add_cylindrical_beam_with_energy( 50000, -35.37, -1.0, 1.0,
                                              5.0, 0.0, 0.5,
                                              Vec3D(-1e-3,0,0), Vec3D(0,1,0), Vec3D(0,0,1), 0.005 );
        // e- : 20 mA total -> 707.4 A/m2
        pdb.add_cylindrical_beam_with_energy( 50000, -707.4, -1.0, 1.0/1836.15,
                                              5.0, 0.0, 0.5,
                                              Vec3D(-1e-3,0,0), Vec3D(0,1,0), Vec3D(0,0,1), 0.005 );

        t0 = now_sec();
        pdb.iterate_trajectories( scharge, efield, bfield );
        t_trace += now_sec() - t0;
    }

    const double total = t_solve + t_trace;
    printf( "threads=%-3d %-9s : solve %7.3f s (%3.0f%%)   trace %7.3f s (%3.0f%%)   total %7.3f s\n",
            eff_threads, use_bicg ? "BiCGSTAB" : "multigrid",
            t_solve, 100.0*t_solve/total, t_trace, 100.0*t_trace/total, total );

    delete mg;
    delete bcg;
}

int main( int argc, char **argv )
{
    int generations = 2;
    bool use_bicg = false;
    vector<int> nums;

    for( int a = 1; a < argc; ++a ) {
        string s( argv[a] );
        if( s == "--bicg" )      use_bicg = true;
        else                     nums.push_back( atoi( argv[a] ) );
    }

    // 参数约定：[generations] [thread1 thread2 ...]
    //   0 个参数 -> 2 代、使用默认线程数
    //   1 个参数 -> 该值作为代数、使用默认线程数
    //   >=2 个   -> 第 1 个为代数，其余为线程数
    vector<int> threads;
    if( nums.size() >= 2 ) {
        generations = nums[0];
        threads.assign( nums.begin()+1, nums.end() );
    } else {
        if( nums.size() == 1 )
            generations = nums[0];
        threads.push_back( 0 );   // 0 = 不干预，使用 IBSimu 默认（硬件并发数）
    }

    for( size_t i = 0; i < threads.size(); ++i )
        run_case( generations, threads[i], use_bicg );

    printf( "\nbenchmark complete\n" );
    return( 0 );
}
