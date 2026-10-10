/** @file cycl_stepper_modes.cpp
 *  @brief  IBSimu 两条粒子推进（步长）方案在回旋加速器场景下的可行性验证。
 *
 *  @test  均匀 Bz 磁场中的单粒子回旋运动（有解析解），分别用：
 *           ① ParticleDataBase::iterate_trajectories()      —— GSL 自适应步长（上游默认）
 *           ② ParticleDataBase::step_particles( ..., dt )   —— Boris 固定步长（逐步推进）
 *         在「非相对论」与「相对论 gamma=1.22」两种束流条件下各追踪 2 圈，
 *         与解析解比较终点轨道误差与速率守恒，确认两条路径均可用。
 *
 *  设计要点（步长选择的两种方案）：
 *    - 方案 ① 由 GSL 依据 epsabs/epsrel 自适应选步长，终点由 max_time / 边界决定；
 *              适合「一次性把轨道追踪到底」的用法。
 *    - 方案 ② 由用户给定 dt 逐步入推进，物理时间精确可控；
 *              适合需要与 RF 周期同步、或需要固定时间栅格的加速器仿真。
 *    两者的相对论能力必须一致，否则「可选」不成立。
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "constants.hpp"
#include "epot_bicgstabsolver.hpp"
#include "epot_efield.hpp"
#include "epot_field.hpp"
#include "geometry.hpp"
#include "ibsimu.hpp"
#include "meshvectorfield.hpp"
#include "particledatabase.hpp"
#include "particles.hpp"


using namespace std;


namespace {

const double PI = 3.14159265358979323846;

/*! \brief 单个算例的解析参数。 */
struct Scenario {
    const char *name;
    bool        relativistic;
    double      B;      //!< 均匀 Bz 场 (T)
    double      r;      //!< 回旋半径 (m)
    double      v;      //!< 速率 (m/s)
    double      gamma;  //!< 洛伦兹因子
    double      omega;  //!< 解析回旋角频率 qB/(gamma m) (rad/s)
    double      T;      //!< 回旋周期 (s)
};


/*! \brief 非相对论算例：由 B 与 r 定出 v，omega = qB/m。 */
Scenario scenario_nonrel( double B, double r )
{
    Scenario s;
    s.name         = "非相对论";
    s.relativistic = false;
    s.B            = B;
    s.r            = r;
    s.omega        = CHARGE_E*B/MASS_U;
    s.v            = s.omega*r;
    double beta    = s.v/SPEED_C;
    s.gamma        = 1.0/sqrt( 1.0 - beta*beta );
    s.T            = 2.0*PI/s.omega;
    return( s );
}


/*! \brief 相对论算例：给定 gamma 与 r，定出 v、B，omega = qB/(gamma m)。 */
Scenario scenario_rel( double gamma, double r )
{
    Scenario s;
    s.name         = "相对论";
    s.relativistic = true;
    s.gamma        = gamma;
    s.r            = r;
    double beta    = sqrt( 1.0 - 1.0/(gamma*gamma) );
    s.v            = beta*SPEED_C;
    s.B            = gamma*MASS_U*s.v/(CHARGE_E*r);
    s.omega        = CHARGE_E*s.B/(gamma*MASS_U);
    s.T            = 2.0*PI/s.omega;
    return( s );
}


/*! \brief 解析解：引导中心位于原点，t=0 时位于 (r,0,0)、速度沿 -y。 */
void analytic( const Scenario &s, double t, double &x, double &y )
{
    x =  s.r*cos( s.omega*t );
    y = -s.r*sin( s.omega*t );
}


/*! \brief 一次追踪的终点结果。 */
struct Result {
    double t;          //!< 终点物理时刻 (s)
    double x, y;       //!< 终点位置 (m)
    double v;          //!< 终点速率 (m/s)
    double gamma;      //!< 终点洛伦兹因子
    double orbit_err;  //!< 相对回旋半径的轨道误差
    double drift;      //!< 速率相对漂移
};


double gamma_of( const Vec3D &v )
{
    // 注意：Vec3D::norm2() 返回的是 2-范数本身（不是范数的平方），
    // 平方和要用 ssqr()。
    double b2 = v.ssqr()/SPEED_C2;
    return( 1.0/sqrt( 1.0 - b2 ) );
}


/*! \brief 从粒子数据库读出终点状态并算出误差指标。 */
Result collect( ParticleDataBase3D &pdb, const Scenario &s )
{
    const Particle3D &p = pdb.particle(0);
    Vec3D loc = p.location();
    Vec3D vel = p.velocity();

    Result r;
    r.t     = p(0);
    r.x     = loc[0];
    r.y     = loc[1];
    r.v     = vel.norm2();
    r.gamma = gamma_of( vel );

    double xa, ya;
    analytic( s, r.t, xa, ya );
    r.orbit_err = sqrt( (r.x-xa)*(r.x-xa) + (r.y-ya)*(r.y-ya) )/s.r;
    r.drift     = fabs( r.v - s.v )/s.v;
    return( r );
}


/*! \brief 方案 ①：GSL 自适应步长，追踪到 max_time。 */
Result run_adaptive( Geometry &geom, const VectorField &efield, const VectorField &bfield,
                     const Scenario &s )
{
    MeshScalarField     scharge( geom );
    ParticleDataBase3D  pdb( geom );

    pdb.set_thread_count( 1 );
    pdb.set_relativistic( s.relativistic );
    pdb.set_max_steps( 2000000 );
    pdb.set_max_time( 2.0*s.T );
    pdb.add_particle( 0.0, 1.0, 1.0,
                      ParticleP3D( 0.0,  s.r, 0.0,  0.0, -s.v,  0.0, 0.0 ) );
    pdb.iterate_trajectories( scharge, efield, bfield );

    return( collect( pdb, s ) );
}


/*! \brief 方案 ②：Boris 固定步长，推进 nsteps 步，每步 dt。 */
Result run_fixed( Geometry &geom, const VectorField &efield, const VectorField &bfield,
                  const Scenario &s, uint32_t nsteps, double dt )
{
    MeshScalarField     scharge( geom );
    ParticleDataBase3D  pdb( geom );

    pdb.set_thread_count( 1 );
    pdb.set_relativistic( s.relativistic );
    pdb.add_particle( 0.0, 1.0, 1.0,
                      ParticleP3D( 0.0,  s.r, 0.0,  0.0, -s.v,  0.0, 0.0 ) );

    for( uint32_t n = 0; n < nsteps; n++ ) {
        // step_particles 每次调用推进一个 dt，并在内部清理/沉积 scharge
        pdb.step_particles( scharge, efield, bfield, dt );
    }

    return( collect( pdb, s ) );
}


/*! \brief 打印一个算例的结果行并做判定。 */
int report( const char *path, const Scenario &s, const Result &r,
            double tol_orbit, int *err )
{
    double xa, ya;
    analytic( s, r.t, xa, ya );
    bool bad = (r.orbit_err > tol_orbit);

    std::printf( "  %-8s %-10s  t=%.6e s (2T=%.6e)  x=% .6e y=% .6e  "
                 "轨道误差=%.3e  漂移=%.3e  gamma=%.6f  %s\n",
                 path, s.name, r.t, 2.0*s.T, r.x, r.y,
                 r.orbit_err, r.drift, r.gamma, bad ? "*** 超差 ***" : "OK" );
    if( bad )
        *err = 1;
    return( bad ? 1 : 0 );
}

} // namespace


void test( int argc, char **argv )
{
    (void)argc;
    (void)argv;

    // ---- 几何：边长 0.1 m 的立方体（回旋半径 0.03 m，轨道完全在网格内）----
    Geometry geom( MODE_3D, Int3D(41,41,41), Vec3D(-0.05,-0.05,-0.05), 0.0025 );
    for( int b = 1; b <= 6; b++ )
        geom.set_boundary( b, Bound(BOUND_DIRICHLET, 0.0) );
    geom.build_mesh();

    // ---- 零电场（全域 Dirichlet 0 且零空间电荷）----
    EpotField            epot( geom );
    MeshScalarField      scharge( geom );
    EpotBiCGSTABSolver   solver( geom );
    solver.solve( epot, scharge );
    EpotEfield           efield( epot );

    // ---- 两个算例 ----
    Scenario scenarios[2];
    scenarios[0] = scenario_nonrel( 0.05, 0.03 );     // beta ~ 4.8e-4
    scenarios[1] = scenario_rel( 1.22, 0.03 );        // 质子 gamma=1.22

    // ---- 均匀 Bz 场（2x2x2 节点覆盖全域）----
    // 注意：fout[i]=true 表示该分量「存在」，false 表示恒为零。
    bool fout[3] = { false, false, true };

    int err = 0;

    std::printf( "\n=== IBSimu 两种步长方案验证（均匀 Bz 场、单粒子回旋运动、追踪 2 圈）===\n\n" );

    for( int is = 0; is < 2; is++ ) {
        Scenario &s = scenarios[is];

        // 固定步长：每圈 1000 步，共 2000 步 = 2 圈（与自适应路径同样时长）
        const uint32_t NSTEP = 2000;
        const double   dt    = s.T/1000.0;

        std::printf( "%-10s beta=%.4f  gamma=%.5f  B=%.6f T  r=%.4f m  "
                     "omega=%.6e rad/s  T=%.6e s  dt=T/1000=%.6e s\n",
                     s.name, s.v/SPEED_C, s.gamma, s.B, s.r, s.omega, s.T, dt );

        MeshVectorField bfield( MODE_3D, fout, Int3D(2,2,2),
                                Vec3D(-0.05,-0.05,-0.05), 0.1 );
        for( int i = 0; i < 2; i++ )
            for( int j = 0; j < 2; j++ )
                for( int k = 0; k < 2; k++ )
                    bfield.set( i, j, k, Vec3D( 0.0, 0.0, s.B ) );

        Result ra = run_adaptive( geom, efield, bfield, s );
        report( "自适应", s, ra, 1.0e-4, &err );

        Result rf = run_fixed( geom, efield, bfield, s, NSTEP, dt );
        report( "固定步长", s, rf, 1.0e-3, &err );

        std::printf( "\n" );
    }

    if( err )
        std::printf( "结果：*** 有方案超差 ***\n" );
    else
        std::printf( "结果：两条路径均与解析解一致\n" );

    exit( err ? 1 : 0 );
}


int main( int argc, char **argv )
{
    try {
        ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
        ibsimu.set_thread_count( 1 );
        test( argc, argv );
    } catch( Error &e ) {
        e.print_error_message( std::cerr );
        return( 1 );
    }
    return( 0 );
}
