/** @file cycl_acceptance.cpp
 *  @brief 需求验收测试：对照 ROADMAP 三大需求与里程碑验收标准（DoD）的顶层核验。
 *
 *  本测试不引入新的物理，而是把**已交付的能力**按需求编号（R1/R2/R3）与
 *  P0 工程合规逐条复核一遍，输出一张“需求 -> 结论”的验收矩阵。它是整个
 *  测试套件（make check）的**门禁汇总**：任何一条验收判据回归都会让
 *  `make check` 失败。
 *
 *  | 编号 | 需求 | 本测试核验的内容 |
 *  | --- | --- | --- |
 *  | P0 | 立项与合规 | LICENSE/NOTICE/第三方清单/CI 配置齐备；可选后端探测结果 |
 *  | R1 | 三维仿真效率 | OpenMP 可用；默认线程数为硬件并发；并行结果与串行逐位一致；并行加速自检；CUDA 后端可选核验 |
 *  | R2 | 回旋加速器物理 | 真实 PSI Ring 场图解析；离面重建满足真空 Maxwell；相对论回旋频率；RF 渡越能量增益；真实场闭合轨道 |
 *  | R3 | 可视化与 IO | VTK XML 写出/独立回读；HDF5 可选后端写出/回读；示例与绘图资产齐备 |
 *
 *  设计约束（与 CI 兼容）：
 *    - 全部检查控制在数秒级（重型场景如 12 圈加速仍由 `cycl_accel` 等专项
 *      测试覆盖）；
 *    - CUDA/HDF5 为可选后端：缺失或不可用时输出 `[SKIP]` 并计为通过；
 *    - 性能类断言带**运行时自检**（环境本身就慢/核数少时跳过数值断言），
 *      使 2 核 CI 稳定通过，而开发机上仍能抓住真实回归。
 *
 *  用法: cycl_acceptance [bfield.dat 路径]
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "config.h"
#include "boris.hpp"
#include "closedorbit.hpp"
#include "ensembletracker.hpp"
#include "fieldmap3d.hpp"
#include "gputracker.hpp"
#include "ibsimu.hpp"
#include "ringfield3d.hpp"
#include "ringmap.hpp"
#include "timevaryingfield.hpp"
#include "vec3d.hpp"
#include "vectorfield.hpp"
#include "vtkwriter.hpp"
#include "hdf5writer.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef IBSIMU_CYCL_HAVE_HDF5
#include <hdf5.h>
#endif

using namespace ibsimu_cycl;


// ===========================================================================
// 基础工具
// ===========================================================================

static const double QE = 1.602176634e-19;
static const double MP = 1.67262192369e-27;
static const double CL = 299792458.0;

static int g_fail = 0;   /*!< \brief 失败检查数（全局）。 */
static int g_skip = 0;   /*!< \brief 跳过检查数（可选后端等）。 */
static int g_nchk = 0;   /*!< \brief 总检查数。 */

static void check( bool ok, const std::string &msg )
{
    ++g_nchk;
    std::printf( "%s  %s\n", ok ? "[ ok ]" : "[FAIL]", msg.c_str() );
    if( !ok )
        ++g_fail;
}

static void skip( const std::string &msg )
{
    ++g_skip;
    std::printf( "[SKIP]  %s\n", msg.c_str() );
}

/*! \brief 文件是否存在（用输入流打开试探）。 */
static bool file_exists( const std::string &path )
{
    std::ifstream f( path.c_str() );
    return( (bool)f );
}

/*! \brief 兼容两种运行目录：make check（CWD=tests/）与仓库根。 */
static std::string pick( const std::string &rel )
{
    if( file_exists( std::string("../") + rel ) )
        return( std::string("../") + rel );
    if( file_exists( rel ) )
        return( rel );
    return( std::string() );
}

static bool check_file( const std::string &rel, const char *what )
{
    bool ok = !pick( rel ).empty();
    check( ok, std::string(what) + " 存在 (" + rel + ")" );
    return( ok );
}

static void banner( const char *title )
{
    std::printf( "\n======================================================================\n" );
    std::printf( " %s\n", title );
    std::printf( "======================================================================\n" );
}

/*! \brief 读整个文件为字符串（失败返回空串）。 */
static std::string slurp( const std::string &path )
{
    std::ifstream f( path.c_str() );
    if( !f )
        return( std::string() );
    std::ostringstream ss;
    ss << f.rdbuf();
    return( ss.str() );
}

/*! \brief 在文本中解析形如 name="v0 v1 v2" 的三分量属性。 */
static bool parse_vec3_attr( const std::string &s, const std::string &name, Vec3D &v )
{
    const std::string key = name + "=\"";
    std::size_t p = s.find( key );
    if( p == std::string::npos )
        return( false );
    p += key.size();
    const char *c = s.c_str() + p;
    char *e = 0;
    for( int a = 0; a < 3; ++a ) {
        v[a] = std::strtod( c, &e );
        if( e == c )
            return( false );
        c = e;
    }
    return( true );
}

/*! \brief 解析 VTK 文件中指定 DataArray 的全部数值（ASCII 内联格式）。 */
static bool parse_vtk_array( const std::string &text, const std::string &name,
                             std::vector<double> &out )
{
    const std::string key = "Name=\"" + name + "\"";
    std::size_t p = text.find( key );
    if( p == std::string::npos )
        return( false );
    std::size_t gt = text.find( '>', p );
    std::size_t end = ( gt == std::string::npos )
        ? std::string::npos : text.find( "</DataArray>", gt );
    if( gt == std::string::npos || end == std::string::npos )
        return( false );
    const char *c = text.c_str() + gt + 1;
    const char *e = text.c_str() + end;
    while( c < e ) {
        char *nx = 0;
        double v = std::strtod( c, &nx );
        if( nx == c ) {          // 非数值字符（空白/标记）跳过
            ++c;
            continue;
        }
        out.push_back( v );
        c = nx;
    }
    return( true );
}

/*! \brief 确定性 LCG（不依赖随机库，保证可复现）。 */
static double lcg( unsigned long &s )
{
    s = s*6364136223846793005UL + 1442695040888963407UL;
    return( ((double)((s >> 11) & 0x1FFFFFFFFFFFFFUL))/9007199254740992.0 );
}

static double now_sec()
{
#ifdef _OPENMP
    return( omp_get_wtime() );
#else
    return( (double)std::clock()/(double)CLOCKS_PER_SEC );
#endif
}


// ===========================================================================
// 公共夹具：真实 PSI Ring 场图 + 参考轨道参数
// ===========================================================================

struct Fixture {
    std::unique_ptr<CRingFieldMap3D> map;
    std::string                      file;
    double Bbar = 0.0;    /*!< <Bz> at r_ref [T] */
    double r_ref = 3.3;   /*!< 参考半径 [m] */
    double gamma = 0.0;   /*!< 参考洛伦兹因子（p = q<B>r） */
    double v0 = 0.0;      /*!< 参考速率 [m/s] */
    double omega = 0.0;   /*!< 解析回旋角频率 [rad/s] */
    bool   ok = false;
};

static Fixture load_fixture( int argc, char **argv )
{
    Fixture fx;
    std::vector<std::string> cands;
    if( argc > 1 )
        cands.push_back( argv[1] );
    cands.push_back( "../examples/cyclotron/data/bfield.dat" );   // CWD = tests/
    cands.push_back( "examples/cyclotron/data/bfield.dat" );      // CWD = 仓库根

    for( std::size_t c = 0; c < cands.size() && !fx.map; ++c ) {
        try {
            fx.map.reset( new CRingFieldMap3D(
                            read_ring_field_map3d( cands[c], 0.1, 8, 0 ) ) );
            fx.file = cands[c];
        } catch( const std::exception & ) { }
    }
    if( !fx.map )
        return( fx );

    // <Bz>(r_ref)，方位角平均 720 点（与 find_closed_orbit 内部一致）
    const int NS = 720;
    for( int i = 0; i < NS; ++i ) {
        const double th = 2.0*M_PI*(double)i/(double)NS;
        fx.Bbar += (*fx.map)( Vec3D( fx.r_ref*std::cos(th), fx.r_ref*std::sin(th), 0.0 ) )[2];
    }
    fx.Bbar /= (double)NS;

    const double gb = QE*fx.Bbar*fx.r_ref/(MP*CL);   // gamma*beta
    fx.gamma = std::sqrt( 1.0 + gb*gb );
    fx.v0    = (gb/fx.gamma)*CL;
    fx.omega = QE*fx.Bbar/(fx.gamma*MP);
    fx.ok    = true;
    return( fx );
}


// ===========================================================================
// §0 P0 立项与合规（里程碑 DoD）
// ===========================================================================

static void section_p0( void )
{
    banner( "§0 P0 立项与工程化：许可证 / 文档 / CI 齐备" );

    check_file( "LICENSE", "许可证" );
    check_file( "NOTICE", "版权归属声明" );
    check_file( "THIRD_PARTY_NOTICES.md", "第三方组件清单" );
    check_file( "CONTRIBUTING.md", "贡献指南" );
    check_file( "README.md", "项目说明" );
    check_file( "configure.ac", "autotools 构建入口" );
    check_file( ".github/workflows/ci.yml", "CI 工作流" );
    check_file( "docs/ROADMAP.md", "路线图" );
    check_file( "docs/DEVELOPMENT.md", "开发说明" );
    check_file( "docs/UPSTREAM_ARCHITECTURE.md", "上游架构总结" );
    check_file( "docs/WORK_LOG.md", "工作日志" );

    // LICENSE 内容抽检：确为 GNU 许可证全文而非占位符
    const std::string lic = slurp( pick( "LICENSE" ) );
    check( lic.find( "GENERAL PUBLIC LICENSE" ) != std::string::npos
           && lic.size() > 20000,
           "LICENSE 为 GNU 许可证全文（非占位符）" );

    // 可选后端探测结果（信息行，不作为判据）
    int omp_max = 1;
#ifdef _OPENMP
    omp_max = omp_get_max_threads();
#endif
    std::printf( "[info] 后端探测: OpenMP=%d 线程; HDF5=%s; CUDA=%s\n",
                 omp_max,
                 hdf5_available() ? hdf5_version_string().c_str() : "未编译",
                 CGpuEnsembleTracker::compiled_in() ? "已编译" : "未编译" );
}


// ===========================================================================
// §1 R1 三维仿真效率（并行化 / CUDA）
// ===========================================================================

static void section_r1( const Fixture &fx )
{
    banner( "§1 R1 三维仿真效率：OpenMP 正确性/加速 + CUDA 可选核验" );

    // ---- 1.1 OpenMP 编译可用 ----
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
    check( nthreads >= 1, "OpenMP 已启用（omp_get_max_threads>=1）" );
#else
    check( false, "OpenMP 已启用（当前构建未定义 _OPENMP）" );
#endif
    const unsigned hw = std::thread::hardware_concurrency();
    if( hw >= 2 ) {
#ifdef _OPENMP
        check( nthreads >= 2, "本机多核（>=2 线程）可用" );
#endif
    }
    std::printf( "[info] 硬件并发 = %u, OpenMP 最大线程 = %d\n", hw, nthreads );

    // ---- 1.2 默认线程数 = 硬件并发（IBSIMU_THREADS 未设置时）----
    {
        const uint32_t def = ibsimu.get_thread_count();
        const char *env = std::getenv( "IBSIMU_THREADS" );
        if( env ) {
            std::printf( "[info] IBSIMU_THREADS=%s -> get_thread_count()=%u（用户覆盖，跳过默认值断言）\n",
                         env, (unsigned)def );
        } else {
            const uint32_t expect = hw > 0 ? (uint32_t)hw : 1u;
            check( def == expect,
                   "IBSimu 默认线程数 == 硬件并发数（R1 PIC 修正）" );
        }
    }

    if( !fx.ok ) {
        skip( "真实场图不可用；跳过 R1 并行一致性/加速与 CUDA 核验" );
        return;
    }

    // ---- 参考系综（与 cycl_omp_tracker 相同的构造）----
    const double T_rev = 2.0*M_PI/fx.omega;
    const double dt    = T_rev/2000.0;

    const std::size_t N = 512;
    std::vector<EnsembleParticle> p0( N );
    unsigned long seed = 20261010UL;
    for( std::size_t i = 0; i < N; ++i ) {
        const double dr  = (lcg(seed) - 0.5)*0.02;      // ±10 mm
        const double dth = (lcg(seed) - 0.5)*2.0e-3;
        const double dvf = (lcg(seed) - 0.5)*2.0e-3;
        const double th  = M_PI/2.0 + dth;
        const double r   = fx.r_ref + dr;
        p0[i].x = Vec3D( r*std::cos(th), r*std::sin(th), (lcg(seed)-0.5)*0.01 );
        const double v = fx.v0*(1.0 + dvf);
        p0[i].v = Vec3D( -std::sin(th)*v, std::cos(th)*v, 0.0 );
    }

    // ---- 1.3 并行与串行结果逐位一致（真实场图）----
    {
        CEnsembleTracker tracker( fx.map.get(), 0, QE, MP );
        tracker.set_relativistic( true );

        std::vector<EnsembleParticle> a = p0, b = p0;
#ifdef _OPENMP
        omp_set_num_threads( 1 );
#endif
        tracker.track( a, dt, 40 );
#ifdef _OPENMP
        omp_set_num_threads( nthreads );
#endif
        tracker.track( b, dt, 40 );

        double maxdiff = 0.0;
        for( std::size_t i = 0; i < N; ++i )
            for( int c = 0; c < 3; ++c ) {
                maxdiff = std::max( maxdiff, std::fabs( a[i].x[c] - b[i].x[c] ) );
                maxdiff = std::max( maxdiff, std::fabs( a[i].v[c] - b[i].v[c] ) );
            }
        check( maxdiff == 0.0,
               "并行与串行跟踪结果逐位一致（max|diff| == 0）" );
    }

    // ---- 1.4 并行加速自检（均匀场，运行时条件断言）----
    // 测量参数沿用 cycl_omp_tracker 的稳健配置（大系综 x 长窗口 x 3 轮取最优）：
    // 测量窗口太短时，一次瞬时调度抖动就会让 12 线程的屏障同步集体停顿，
    // 造成 0.3x 级别的假性回退（本测试曾在套件长跑中出现过一次）。
    {
        CFieldMap3D Buniform( 2, 2, 2, -5.0, 10.0, -5.0, 10.0, -5.0, 10.0 );
        for( std::size_t i = 0; i < 2; ++i )
            for( std::size_t j = 0; j < 2; ++j )
                for( std::size_t k = 0; k < 2; ++k )
                    Buniform.set_value( i, j, k, 0.0, 0.0, fx.Bbar );

        std::vector<EnsembleParticle> p_big;
        p_big.reserve( p0.size()*40 );
        for( int rep = 0; rep < 40; ++rep )
            p_big.insert( p_big.end(), p0.begin(), p0.end() );

        CEnsembleTracker tracker( &Buniform, 0, QE, MP );
        tracker.set_relativistic( true );

        {   // 预热
            std::vector<EnsembleParticle> w = p_big;
            tracker.track( w, dt, 50 );
        }

        // 独立平凡并行循环做环境自检（负载/调度异常时不冤枉实现）
        double sanity = 1.0;
        {
            const long NB = 8000000;
            auto busy = []( long n ) {
                double s = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:s) schedule(static)
#endif
                for( long i = 0; i < n; ++i )
                    s += std::sqrt( (double)(i % 1000) );
                return( s );
            };
#ifdef _OPENMP
            omp_set_num_threads( 1 );
#endif
            double t0 = now_sec();
            double r1 = busy( NB );
            t0 = now_sec() - t0;
            const double ts = t0;
#ifdef _OPENMP
            omp_set_num_threads( nthreads );
#endif
            t0 = now_sec();
            double r2 = busy( NB );
            t0 = now_sec() - t0;
            const double tp = t0;
            sanity = ts/tp;
            if( std::fabs( r1 - r2 ) > 1e-6*std::fabs(r1) )
                sanity = 0.0;
        }

        double t_ser = 1.0e30, t_par = 1.0e30;
        double speedup = 0.0;
        int    rounds = 0;
        for( int round = 0; round < 2; ++round ) {
            ++rounds;
            for( int rep = 0; rep < 3; ++rep ) {
                std::vector<EnsembleParticle> a = p_big, b = p_big;
#ifdef _OPENMP
                omp_set_num_threads( 1 );
#endif
                double t0 = now_sec();
                tracker.track( a, dt, 100 );
                t_ser = std::min( t_ser, now_sec() - t0 );
#ifdef _OPENMP
                omp_set_num_threads( nthreads );
#endif
                t0 = now_sec();
                tracker.track( b, dt, 100 );
                t_par = std::min( t_par, now_sec() - t0 );
            }
            speedup = std::max( speedup, t_ser/t_par );
            if( speedup > 1.5 )
                break;   // 达标即止；未达标时自动再测一轮（抵御瞬时抖动）
        }
        std::printf( "[info] 均匀场 N=%zu x100 步: serial=%.3f s parallel=%.3f s -> %.2fx"
                     " (环境自检 %.2fx, %d 轮 x3 次取最优)\n",
                     p_big.size(), t_ser, t_par, speedup, sanity, rounds );

        if( nthreads >= 4 && sanity >= 2.0 )
            check( speedup > 1.5, "并行加速比 > 1.5x（R1 KPI 方向正确）" );
        else
            std::printf( "        （仅 %d 线程/自检 %.2fx，跳过数值加速断言——与 CI 一致的降级路径）\n",
                         nthreads, sanity );
    }

    // ---- 1.5 CUDA 后端（可选）----
    {
        if( !CGpuEnsembleTracker::compiled_in() ) {
            skip( "CUDA 后端未编译（--without-cuda）" );
        } else {
            std::string reason;
            if( !CGpuEnsembleTracker::available( &reason ) ) {
                skip( "GPU 运行期不可用: " + reason );
            } else {
                CEnsembleTracker cpu( fx.map.get(), 0, QE, MP );
                cpu.set_relativistic( true );
                CGpuEnsembleTracker gpu( *fx.map, QE, MP );
                gpu.set_relativistic( true );

                std::vector<EnsembleParticle> a = p0, b = p0;
#ifdef _OPENMP
                omp_set_num_threads( std::min( nthreads, 4 ) );
#endif
                cpu.track( a, dt, 100 );
                gpu.track( b, dt, 100 );

                double dx = 0.0, dv = 0.0;
                for( std::size_t i = 0; i < N; ++i )
                    for( int c = 0; c < 3; ++c ) {
                        dx = std::max( dx, std::fabs( a[i].x[c] - b[i].x[c] ) );
                        dv = std::max( dv, std::fabs( a[i].v[c] - b[i].v[c] ) );
                    }
                std::printf( "[info] GPU=%s: max|dx|=%.2e m, max|dv|/v0=%.2e\n",
                             CGpuEnsembleTracker::device_name().c_str(),
                             dx, dv/fx.v0 );
                check( dx < 1e-9 && dv/fx.v0 < 1e-9,
                       "GPU 与 CPU 结果一致（相对机器精度级）" );
            }
        }
    }
}


// ===========================================================================
// §2 R2 回旋加速器物理（三维磁场 / RF / 时变场 / 闭合轨道）
// ===========================================================================

/*! \brief 平滑 RF 间隙电场剖面（仅 x 分量，复刻 cycl_rfgap 的解析模型）。 */
class CGapEField : public VectorField {
    double _e0, _h;
public:
    CGapEField( double e0, double h ) : _e0(e0), _h(h) {}
    virtual const Vec3D operator()( const Vec3D &x ) const {
        const double u = x[0]/_h;
        if( std::fabs(u) >= 1.0 )
            return( Vec3D( 0.0, 0.0, 0.0 ) );
        const double f = 0.5*(1.0 + std::cos( M_PI*u ));
        return( Vec3D( _e0*f, 0.0, 0.0 ) );
    }
};

static double gap_profile( double x, double h )
{
    if( std::fabs(x) >= h )
        return( 0.0 );
    return( 0.5*(1.0 + std::cos( M_PI*x/h )) );
}

/*! \brief 精确渡越积分（Simpson）。 */
static double exact_dW( double e0, double h, double q, double v,
                        double omega, double phi0, double x_entry )
{
    const int n = 20000;
    const double a = -h, b = h, dx = (b-a)/n;
    double sum = 0.0;
    for( int i = 0; i <= n; ++i ) {
        const double x = a + i*dx;
        const double t = (x - x_entry)/v;
        const double f = gap_profile( x, h )*std::cos( omega*t + phi0 );
        const double w = (i == 0 || i == n) ? 1.0 : ((i % 2) ? 4.0 : 2.0);
        sum += w*f;
    }
    return( q*e0*sum*dx/3.0 );
}

/*! \brief 用 Boris 跟踪穿越间隙，返回 dW。 */
static double track_dW( double e0, double h, double v0, double omega, double phi0,
                        double x_entry, double x_exit, double dt )
{
    CGapEField gap( e0, h );
    CTimeVaryingField efield( &gap, omega, phi0 );
    CBorisPusher pusher( QE, MP );

    Vec3D x( x_entry, 0.0, 0.0 ), v( v0, 0.0, 0.0 );
    const int nsteps = (int)std::ceil( (x_exit - x_entry)/v0/dt ) + 10;
    for( int n = 0; n < nsteps; ++n ) {
        efield.set_time( n*dt );
        pusher.step( &efield, 0, x, v, dt );
        if( x[0] >= x_exit )
            break;
    }
    return( 0.5*MP*(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) - 0.5*MP*v0*v0 );
}

static CFieldMap3D uniform_B( double bz )
{
    CFieldMap3D f( 2, 2, 2, -5.0, 10.0, -5.0, 10.0, -5.0, 10.0 );
    for( std::size_t i = 0; i < 2; ++i )
        for( std::size_t j = 0; j < 2; ++j )
            for( std::size_t k = 0; k < 2; ++k )
                f.set_value( i, j, k, 0.0, 0.0, bz );
    return( f );
}

static double gamma_of( const Vec3D &v )
{
    const double b2 = (v[0]*v[0] + v[1]*v[1] + v[2]*v[2])/(CL*CL);
    return( 1.0/std::sqrt(1.0 - b2) );
}

/*! \brief 均匀场中推进并测量回旋角频率与 gamma 漂移（复刻 cycl_relativistic）。 */
static double measure_omega( const CFieldMap3D &B, const Vec3D &v0, double dt,
                             int nsteps, double &gamma_drift )
{
    CBorisPusher pusher( QE, MP );
    pusher.set_relativistic( true );

    Vec3D x( 0.0, 0.0, 0.0 ), v = v0;
    const double g0 = gamma_of( v );
    double phi_prev = std::atan2( v[1], v[0] );
    double phi_sum = 0.0;
    double gmax = g0, gmin = g0;
    for( int n = 0; n < nsteps; ++n ) {
        pusher.step( 0, &B, x, v, dt );
        const double phi = std::atan2( v[1], v[0] );
        double d = phi - phi_prev;
        while( d >  M_PI ) d -= 2.0*M_PI;
        while( d < -M_PI ) d += 2.0*M_PI;
        phi_sum += d;
        phi_prev = phi;
        const double g = gamma_of( v );
        gmax = std::max( gmax, g );
        gmin = std::min( gmin, g );
    }
    gamma_drift = (gmax - gmin)/g0;
    return( std::fabs(phi_sum)/(nsteps*dt) );
}

static void section_r2( int argc, char **argv, const Fixture &fx )
{
    banner( "§2 R2 回旋加速器物理：真实场图 / 离面重建 / 相对论 / RF / 闭合轨道" );

    if( !fx.ok ) {
        skip( "真实场图不可用；跳过 R2 全部核验（检查数据文件是否随仓库提供）" );
        return;
    }
    std::printf( "[info] 场图: %s  <Bz>(r=%.2f m)=%.4f T  gamma=%.4f  v0=%.3e m/s\n",
                 fx.file.c_str(), fx.r_ref, fx.Bbar, fx.gamma, fx.v0 );
    const CRingFieldMap3D &B = *fx.map;

    // ---- 2.1 真实场图解析（OPAL/PSI RING 格式）----
    {
        RingFieldMapInfo info;
        try {
            read_ring_field_map3d( fx.file, 0.1, 8, &info );
        } catch( const std::exception &e ) {
            check( false, std::string("场图解析头部失败: ") + e.what() );
        }
        check( info.nrad == 141 && info.ntet == 135,
               "场图网格 141 x 135（PSI Ring 原数据）" );
        check( std::fabs( info.rmin_mm - 1900.0 ) < 1e-6
               && std::fabs( info.dr_mm - 20.0 ) < 1e-6
               && std::fabs( info.dtet_deg - 1.0/3.0 ) < 1e-9,
               "网格步长 rmin=1900 mm, dr=20 mm, dtheta=1/3 deg" );
        check( B.size_r() == 141 && B.size_theta() == 1080,
               "8 折对称复制到满 360 度（1080 列）" );

        const double bz31 = B.bz_midplane( 60, 0 );   // r=3.1 m 为节点 i=60
        std::printf( "[info] Bz(r=3.1 m, theta=0) = %.6f T\n", bz31 );
        check( std::fabs( bz31 - 1.526363 ) < 2.0e-4,
               "Bz(3.1 m) 与参考值一致（1.5264 T）" );
    }

    // ---- 2.2 离面重建满足真空 Maxwell 方程 ----
    {
        const double r  = B.r0() + 70.0*B.dr();   // 3.3 m
        const double th = M_PI/2.0;
        const double z  = 0.02;
        Vec3D c( r*std::cos(th), r*std::sin(th), z );
        const double h = 1.0e-3;
        const Vec3D bx1 = B( Vec3D(c[0]+h, c[1], c[2]) ), bx0 = B( Vec3D(c[0]-h, c[1], c[2]) );
        const Vec3D by1 = B( Vec3D(c[0], c[1]+h, c[2]) ), by0 = B( Vec3D(c[0], c[1]-h, c[2]) );
        const Vec3D bz1 = B( Vec3D(c[0], c[1], c[2]+h) ), bz0 = B( Vec3D(c[0], c[1], c[2]-h) );

        const double div = ( (bx1[0]-bx0[0]) + (by1[1]-by0[1]) + (bz1[2]-bz0[2]) )/(2.0*h);
        const double curlx = ( (by1[2]-by0[2]) - (bz1[1]-bz0[1]) )/(2.0*h);
        const double curly = ( (bz1[0]-bz0[0]) - (bx1[2]-bx0[2]) )/(2.0*h);
        const double curlz = ( (bx1[1]-bx0[1]) - (by1[0]-by0[0]) )/(2.0*h);
        const double curl  = std::sqrt( curlx*curlx + curly*curly + curlz*curlz );
        const Vec3D b = B( c );
        const double bmag = std::sqrt( b[0]*b[0] + b[1]*b[1] + b[2]*b[2] );
        std::printf( "[info] 真空检查: |div B|=%.3e  |curl B|=%.3e  (|B|/r=%.4f T/m)\n",
                     std::fabs(div), curl, bmag/r );
        check( std::fabs(div) < 1e-2*(bmag/r), "真空: |div B| << |B|/r" );
        check( curl < 1e-2*(bmag/r), "真空: |curl B| << |B|/r" );

        // 45 度旋转不变性
        const Vec3D b2 = B( Vec3D( r*std::cos(th + M_PI/4.0),
                                   r*std::sin(th + M_PI/4.0), z ) );
        const double m2 = std::sqrt( b2[0]*b2[0] + b2[1]*b2[1] + b2[2]*b2[2] );
        check( std::fabs( bmag - m2 ) < 1e-6*std::max( 1.0, bmag ),
               "|B| 在 45 度旋转下不变（8 折对称）" );

        // 小 z 渐近: Br/z ~ dBz/dr
        const Vec3D bp = B( c );
        const double slope = bp[1]/z;   // theta=90° 时 Br 即 By
        const double b_rm = B( Vec3D( (r-B.dr())*std::cos(th), (r-B.dr())*std::sin(th), 0.0 ) )[2];
        const double b_rp = B( Vec3D( (r+B.dr())*std::cos(th), (r+B.dr())*std::sin(th), 0.0 ) )[2];
        const double dbr_an = (b_rp - b_rm)/(2.0*B.dr());
        check( std::fabs( slope - dbr_an ) < 0.05*std::fabs(dbr_an) + 1e-9,
               "小 z 渐近: Br/z ≈ dBz/dr（5% 内）" );
    }

    // ---- 2.3 相对论回旋频率 ω_c = qB/(γm) ----
    {
        const double B0 = 1.0, gamma0 = 1.2;
        const double beta0 = std::sqrt( 1.0 - 1.0/(gamma0*gamma0) );
        const double v0 = beta0*CL;
        CFieldMap3D Bf = uniform_B( B0 );

        const double omega_exp = QE*B0/(gamma0*MP);
        const double dt = (2.0*M_PI/omega_exp)/500.0;
        double gdrift = 0.0;
        const double omega_meas = measure_omega( Bf, Vec3D( v0, 0.0, 0.0 ), dt, 1500, gdrift );
        const double rel = std::fabs( omega_meas - omega_exp )/omega_exp;
        std::printf( "[info] 相对论回旋: omega=%.6e (期望 %.6e) rel=%.2e gamma漂移=%.2e\n",
                     omega_meas, omega_exp, rel, gdrift );
        check( rel < 1e-3, "相对论回旋频率 == qB/(gamma m)（1e-3 内）" );
        check( gdrift < 1e-12, "静磁场中 gamma 精确守恒" );
    }

    // ---- 2.4 RF 渡越能量增益（1D 平滑间隙）----
    {
        const double h = 0.02;            // 间隙半宽 [m]
        const double e0 = 1000.0/h;       // V0 = 1000 V
        const double v0 = 5.0e6;          // ~130 keV 质子
        const double omega = 2.0*M_PI*10.0e6;
        const double x_entry = -0.06, x_exit = 0.06;
        const double dt = 2.0e-11;

        const double dW_num = track_dW( e0, h, v0, omega, 0.0, x_entry, x_exit, dt );
        const double dW_ana = exact_dW( e0, h, QE, v0, omega, 0.0, x_entry );
        const double rel = std::fabs( dW_num - dW_ana )/std::fabs( dW_ana );
        std::printf( "[info] RF 渡越: dW_num=%.3f eV  dW_ana=%.3f eV  rel=%.2e\n",
                     dW_num/QE, dW_ana/QE, rel );
        check( rel < 0.02, "RF 能量增益与精确渡越积分一致（2% 内）" );

        const double dW_pi = track_dW( e0, h, v0, omega, M_PI, x_entry, x_exit, dt );
        check( dW_num > 0.0 && dW_pi < 0.0, "RF 相位: +V 加速 / -V 减速" );
    }

    // ---- 2.5 真实场闭合轨道（R2-D 的起步条件）----
    {
        ClosedOrbit co;
        bool ok = true;
        try {
            co = find_closed_orbit( B, QE, MP, fx.gamma, M_PI/2.0, fx.r_ref, true );
        } catch( const std::exception &e ) {
            ok = false;
            check( false, std::string("闭合轨道求解失败: ") + e.what() );
        }
        if( ok ) {
            const double scallop = co.rmax - co.rmin;
            const double vr_frac = std::fabs(co.vr)/co.vmag;
            const double omega_rev = 2.0*M_PI/co.T_rev;
            const double omega_ana = QE*co.Bbar/(fx.gamma*MP);
            const double rel = std::fabs( omega_rev - omega_ana )/omega_ana;
            std::printf( "[info] 闭合轨道: r=%.4f m  v_r=%+.3e m/s (%.2f%% of |v|)  "
                         "scallop=%.1f mm  T_rev=%.3f ns  omega 偏差=%.2e\n",
                         co.r, co.vr, 100.0*vr_frac, 1000.0*scallop,
                         co.T_rev/1e-9, rel );
            check( std::fabs( co.r - fx.r_ref ) < 0.15,
                   "闭合轨道半径落在参考半径 15 cm 内" );
            check( co.vr < 0.0 && vr_frac > 0.02 && vr_frac < 0.15,
                   "不动点带非零径向速度（非切向发射，2%~15% |v|）" );
            check( co.rmax > co.rmin && scallop < 0.35,
                   "扇贝幅度 < 0.35 m（远小于朴素初值的 0.62 m）" );
            check( rel < 0.035,
                   "回路频率与 q<B>/(gamma m) 一致（3.5% 内）" );
        }
    }
}


// ===========================================================================
// §3 R3 可视化与 IO
// ===========================================================================

static void section_r3( void )
{
    banner( "§3 R3 可视化与 IO：VTK XML / HDF5 / 示例资产" );

    // ---- 3.1 VTK XML 写出 + 独立回读 ----
    {
        const std::string fname = "cycl_acceptance_field.vti";
        const Int3D dims( 3, 2, 2 );
        const Vec3D origo( 0.0, 0.0, 0.0 );
        const Vec3D spacing( 0.5, 1.0, 1.5 );

        const std::size_t n = 3*2*2;
        std::vector<double> scalar( n ), vector( 3*n );
        for( std::size_t k = 0; k < 2; ++k )
            for( std::size_t j = 0; j < 2; ++j )
                for( std::size_t i = 0; i < 3; ++i ) {
                    const std::size_t a = i + 3*(j + 2*k);
                    const double f = (double)i + 10.0*(double)j + 100.0*(double)k + 0.5;
                    scalar[a] = f;
                    vector[3*a+0] =  f;
                    vector[3*a+1] = 2.0*f;
                    vector[3*a+2] = -f;
                }
        vtk_write_image_data( fname, dims, origo, spacing,
                              "acc_scalar", scalar, "acc_vector", vector );

        const std::string text = slurp( fname );
        check( text.find( "<VTKFile type=\"ImageData\"" ) != std::string::npos
               && text.find( "WholeExtent=\"0 2 0 1 0 1\"" ) != std::string::npos,
               "VTK .vti 头部与网格范围正确" );

        Vec3D io( 0.0, 0.0, 0.0 ), sp( 0.0, 0.0, 0.0 );
        check( parse_vec3_attr( text, "Origin", io ) && parse_vec3_attr( text, "Spacing", sp )
               && io[0] == 0.0 && sp[0] == 0.5 && sp[1] == 1.0 && sp[2] == 1.5,
               "VTK Origin/Spacing 属性回读一致" );

        std::vector<double> rs, rv;
        bool ok = parse_vtk_array( text, "acc_scalar", rs ) && rs.size() == n
                  && parse_vtk_array( text, "acc_vector", rv ) && rv.size() == 3*n;
        if( ok )
            for( std::size_t a = 0; a < n && ok; ++a )
                ok = ( rs[a] == scalar[a] && rv[3*a] == vector[3*a]
                       && rv[3*a+1] == vector[3*a+1] && rv[3*a+2] == vector[3*a+2] );
        check( ok, "VTK 标量/矢量数据回读逐值一致" );

        // 轨迹折线（时刻取 0.125 的倍数：二进制/十进制均可精确表示，
        // 保证 12 位有效数字写出后能逐位回读——1e-9 的普通倍数做不到）
        const std::string pname = "cycl_acceptance_orbit.vtp";
        std::vector<TrajectoryPoint> line( 4 );
        for( int a = 0; a < 4; ++a ) {
            line[a].t = 0.125*(double)a;
            line[a].x = Vec3D( 0.25*(double)a, 1.5, -2.25 );
        }
        vtk_write_polylines( pname, std::vector<std::vector<TrajectoryPoint>>( 1, line ) );
        const std::string ptext = slurp( pname );
        std::vector<double> tt, pts;
        ok = ptext.find( "<VTKFile type=\"PolyData\"" ) != std::string::npos
             && ptext.find( "NumberOfPoints=\"4\"" ) != std::string::npos
             && parse_vtk_array( ptext, "t", tt ) && tt.size() == 4
             && parse_vtk_array( ptext, "Points", pts ) && pts.size() == 12;
        if( ok )
            for( int a = 0; a < 4 && ok; ++a ) {
                ok = ( pts[3*a] == 0.25*(double)a ) && tt[a] == 0.125*(double)a;
            }
        if( !ok )
            std::printf( "        (t.size=%zu pts.size=%zu)\n", tt.size(), pts.size() );
        check( ok, "VTK .vtp 轨迹点与时刻回读一致" );
    }

    // ---- 3.2 HDF5 二进制后端（可选）----
    if( !hdf5_available() ) {
        skip( "HDF5 后端未编译（--without-hdf5）" );
    } else {
#ifdef IBSIMU_CYCL_HAVE_HDF5
        const std::string fname = "cycl_acceptance_io.h5";
        const Int3D dims( 4, 3, 2 );
        const std::size_t n = 4*3*2;
        std::vector<double> scalar( n ), vector( 3*n );
        for( std::size_t k = 0; k < 2; ++k )
            for( std::size_t j = 0; j < 3; ++j )
                for( std::size_t i = 0; i < 4; ++i ) {
                    const std::size_t a = i + 4*(j + 3*k);
                    const double f = (double)i + 10.0*(double)j + 100.0*(double)k;
                    scalar[a] = f;
                    vector[3*a+0] = f;
                    vector[3*a+1] = 2.0*f;
                    vector[3*a+2] = -f;
                }
        hdf5_write_image_data( fname, "acctest", dims, Vec3D(0,0,0), Vec3D(0.5,1.0,1.5),
                               "acc_s", scalar, 1.0, "acc_v", vector, 1.0 );

        hid_t fid = H5Fopen( fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
        check( fid >= 0, "HDF5 文件可重新打开" );
        if( fid >= 0 ) {
            // 独立回读（HDF5 C API，不复用写出器封装）
            auto read_ds = []( hid_t file, const char *path,
                               std::vector<double> &out, std::vector<hsize_t> &dims_out ) {
                hid_t ds = H5Dopen2( file, path, H5P_DEFAULT );
                if( ds < 0 )
                    return( false );
                hid_t sp = H5Dget_space( ds );
                const int rank = (sp >= 0) ? H5Sget_simple_extent_ndims( sp ) : 0;
                hsize_t d[4] = {0,0,0,0};
                if( rank >= 1 )
                    H5Sget_simple_extent_dims( sp, d, 0 );
                std::size_t nn = 1;
                for( int a = 0; a < rank; ++a )
                    nn *= (std::size_t)d[a];
                out.assign( nn, 0.0 );
                const bool okread = ( rank >= 1 )
                    && H5Dread( ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                                H5P_DEFAULT, &out[0] ) >= 0;
                dims_out.assign( d, d + std::max( rank, 0 ) );
                if( sp >= 0 ) H5Sclose( sp );
                H5Dclose( ds );
                return( okread );
            };

            std::vector<double> rs, rv;
            std::vector<hsize_t> sd, vd;
            const bool okr = read_ds( fid, "/data/0/meshes/acctest/acc_s", rs, sd )
                             && read_ds( fid, "/data/0/meshes/acctest/acc_v", rv, vd );
            check( okr && rs.size() == n && sd.size() == 3
                   && sd[0] == 2 && sd[1] == 3 && sd[2] == 4
                   && vd.size() == 4 && vd[3] == 3,
                   "HDF5 数据集形状 (nz,ny,nx) 与分量序正确" );
            bool same = okr && rs.size() == scalar.size() && rv.size() == vector.size();
            if( same )
                for( std::size_t a = 0; a < scalar.size() && same; ++a )
                    same = ( rs[a] == scalar[a] ) && ( rv[a] == vector[a] );
            check( same, "HDF5 标量/矢量数据回读逐位一致" );
            H5Fclose( fid );
        }
#else
        skip( "HDF5 可用但缺失编译宏（不应当发生）" );
#endif
    }

    // ---- 3.3 示例与绘图资产齐备（R3 DoD：示例可复现）----
    {
        const char *assets[] = {
            "examples/cyclotron/README.md",
            "examples/cyclotron/plot_acceleration.py",
            "examples/cyclotron/plot_closed_orbit.py",
            "examples/cyclotron/plot_field_map.py",
            "examples/cyclotron/view_3d.py",
            "python/ibsimu_cycl/vtk_io.py",
            "python/ibsimu_cycl/hdf5_io.py",
            "docs/img/cyclotron_field_map.png",
            "docs/img/cyclotron_acceleration.png",
            "docs/img/cyclotron_real_orbit.png",
        };
        int have = 0;
        const int total = (int)(sizeof(assets)/sizeof(assets[0]));
        for( int a = 0; a < total; ++a )
            if( !pick( assets[a] ).empty() )
                ++have;
        check( have == total,
               "可视化示例与成图资产齐备（10/10）" );
        if( have != total )
            for( int a = 0; a < total; ++a )
                if( pick( assets[a] ).empty() )
                    std::printf( "        缺失: %s\n", assets[a] );
    }
}


// ===========================================================================
// main
// ===========================================================================

int main( int argc, char **argv )
{
    std::printf( "IBSimu-Cycl 需求验收测试（R1/R2/R3 + P0 合规）\n" );
    std::printf( "ROADMAP: docs/ROADMAP.md  §5 里程碑验收标准\n" );

    const Fixture fx = load_fixture( argc, argv );
    if( fx.ok )
        std::printf( "场图: %s\n", fx.file.c_str() );
    else
        std::printf( "场图: [不可用——将跳过依赖真实场图的检查]\n" );

    const int f0 = g_fail, s0 = g_skip;

    section_p0();
    const int p0_fail = g_fail - f0;

    const int f1 = g_fail, s1 = g_skip;
    section_r1( fx );
    const int r1_fail = g_fail - f1, r1_skip = g_skip - s1;

    const int f2 = g_fail, s2 = g_skip;
    section_r2( argc, argv, fx );
    const int r2_fail = g_fail - f2, r2_skip = g_skip - s2;

    const int f3 = g_fail, s3 = g_skip;
    section_r3();
    const int r3_fail = g_fail - f3, r3_skip = g_skip - s3;

    banner( "验收汇总（需求 -> 结论）" );
    auto row = []( const char *name, int fails, int skips ) {
        std::printf( " %-28s : %s", name, fails ? "[FAIL]" : "[PASS]" );
        if( skips )
            std::printf( "（跳过 %d 项）", skips );
        std::printf( "\n" );
    };
    row( "P0 立项与工程化", p0_fail, g_skip - s0 );
    row( "R1 三维仿真效率（并行/GPU）", r1_fail, r1_skip );
    row( "R2 回旋加速器物理", r2_fail, r2_skip );
    row( "R3 可视化与 IO", r3_fail, r3_skip );

    std::printf( "\n总计 %d 项检查, %d 失败, %d 跳过\n", g_nchk, g_fail, g_skip );
    std::printf( "\n%s\n", g_fail ? "ACCEPTANCE FAILED" : "ACCEPTANCE PASSED" );
    return( g_fail ? 1 : 0 );
}
