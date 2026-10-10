/** @file cycl_sc_beam.cpp
 *  @brief IBSimu-Cycl：注入束流的多粒子空间电荷自洽模拟（PSI Ring 实场，72 MeV）。
 *
 *  ## 场景（用户规格）
 *
 *  以 `cyclotron2.in` 的单粒子信息（72 MeV 质子、r0 = 2045 mm）为基准：
 *    - 发射面扩展为 **1 mm 半径的标准圆盘均匀分布**（初始速度=参考粒子，零发射度）；
 *    - 宏粒子总数 **1e6**；
 *    - 束流强度 **1 mA**（h=6、f_RF=50.65 MHz → 单束团电荷 Qb = I/f = 19.74 pC，
 *      每个宏粒子 = 123 个质子）；
 *    - **空间电荷通过迭代自洽求解**（冻结场不动点迭代 / Vlasov 式迭代）。
 *
 *  ## 方法
 *
 *  1. 用 `find_closed_orbit()` 在真实场图中求 72 MeV 闭合轨道；沿轨道逐波前
 *     建立随动坐标系 (ŝ, ê1, ẑ)：粒子局部坐标 (ζ, ξ, η)。
 *  2. 束团系 3D PIC：一个**随束团平动的立方盒**（纵向 ζ 沿束团，横向 ξ,η），
 *     用 IBSimu 原生沉积（`scharge_add_step_pic`，CIC）+ 多重网格求解
 *     （`EpotMGSolver`，Dirichlet 盒壁）得到自场 E(ζ,ξ,η)。轨迹按弧长分段
 *     （nseg 段），每段一张场图（由上一轮迭代中参考粒子过该段中心时刻的
 *     分布沉积求出）。
 *  3. 推进用相对论 Boris（`CBorisPusher`，新增显式 E 矢量重载——空间电荷是
 *     逐粒子的场）。横向空间电荷力乘 (1-β²)（运动电荷之间的磁场抵消因子）。
 *  4. 迭代：pass k 用图 k-1 推进 → 由结果重建图 k → 直到场图相对变化
 *     r_E < tol（自洽）。最终用收敛图冻结推进（含多圈外推）。
 *  5. 验证：
 *     - 场链（沉积+求解+插值）与**直接两两求和**对比（小样本）；
 *     - 与无限长线电荷解析场对比（量级/单位检查）；
 *     - 静磁场下 |v| 守恒（Boris 精确性）；
 *     - 时间步长 T/4000 与 T/2000 对比（手动命令行复核）。
 *  6. 分析：单圈线性映射拟合 → 分数工作点 ν_r/ν_z（SC 与无 SC 对比 = 空间电荷
 *     工作点漂移）；RMS 发射度沿圈演化；包络 σ(s)；束团纵向拉伸；SC 场剖面。
 *
 *  ## 范围的说明（重要）
 *
 *  - RF 加速关闭（研究注入束团在整机磁场中的自洽演化；190 圈 × 1e6 × 迭代
 *    在单机不可行，加速段的空间电荷留待后续分级接入）；
 *  - 只模拟 h=6 中的**一个束团**（相邻束团距离 2.14 m，其场贡献 ~0.04 V/m，
 *    比自场低 5 个量级，可忽略——代码中打印该估计）；
 *  - 「1 mA」通过 Qb = I/f_RF 折算到单个束团的电荷。束团长度 L 未在源数据中
 *    给出，作为参数（默认 10 cm ≈ 16° @50.65 MHz），SC 强度 ∝ 1/L；
 *  - 盒壁 Dirichlet 与有限束团长度带来的截断误差由「直接求和」验证量化。
 *
 *  用法（CI 为快速默认；全尺寸用 --full）:
 *    cycl_sc_beam                          # CI: N=4000, 1/4 圈, 2 迭代
 *    cycl_sc_beam --full                   # N=1e6, 3 圈, 6 迭代, nseg=16
 *    cycl_sc_beam --n 200000 --nsteps 4000 --turns 3 --iters 6 --nseg 16
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "boris.hpp"
#include "closedorbit.hpp"
#include "constants.hpp"
#include "epot_field.hpp"
#include "epot_mgsolver.hpp"
#include "geometry.hpp"
#include "hdf5writer.hpp"
#include "meshscalarfield.hpp"
#include "ringfield3d.hpp"
#include "scharge.hpp"
#include "vec3d.hpp"
#include "vectorfield.hpp"
#include "vtkwriter.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace ibsimu_cycl;


// ===========================================================================
// 常量与工具
// ===========================================================================

static const double QE = 1.602176634e-19;      /*!< \brief 元电荷 [C]。 */
static const double MP = 1.67262192369e-27;    /*!< \brief 质子质量 [kg]。 */
static const double CL = 299792458.0;          /*!< \brief 光速 [m/s]。 */
static const double EPS0 = 8.854187817e-12;    /*!< \brief 真空介电常数 [F/m]。 */

static int g_fail = 0;

static void check( bool ok, const std::string &msg )
{
    std::printf( "%s  %s\n", ok ? "[ ok ]" : "[FAIL]", msg.c_str() );
    if( !ok )
        ++g_fail;
}

static double now_sec()
{
#ifdef _OPENMP
    return( omp_get_wtime() );
#else
    return( (double)std::clock()/(double)CLOCKS_PER_SEC );
#endif
}

static double lcg( unsigned long &s )
{
    s = s*6364136223846793005UL + 1442695040888963407UL;
    return( ((double)((s >> 11) & 0x1FFFFFFFFFFFFFUL))/9007199254740992.0 );
}

static double norm( const Vec3D &a )
{
    return( std::sqrt( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] ) );
}

static Vec3D crossp( const Vec3D &a, const Vec3D &b )
{
    return( Vec3D( a[1]*b[2] - a[2]*b[1],
                   a[2]*b[0] - a[0]*b[2],
                   a[0]*b[1] - a[1]*b[0] ) );
}


// ===========================================================================
// 参数
// ===========================================================================

struct Params {
    int    nparts   = 4000;       /*!< \brief 宏粒子数（CI 默认）。 */
    int    turns    = 1;          /*!< \brief 模拟圈数。 */
    int    maxit    = 2;          /*!< \brief SC 迭代上限。 */
    int    nseg     = 4;          /*!< \brief 弧长分段数（每段一张场图）。 */
    int    nsteps   = 1000;       /*!< \brief 每圈步数（CI=1/4 圈）。 */
    double bunch_len = 0.10;      /*!< \brief 束团长度 [m]（沿轨道）。 */
    double beam_r    = 0.001;     /*!< \brief 发射面半径 [m]。 */
    double current   = 1.0e-3;    /*!< \brief 束流强度 [A]。 */
    double f_bunch   = 50.65e6;   /*!< \brief 束团频率 [Hz]（f_RF）。 */
    double ke_mev    = 72.0;      /*!< \brief 动能 [MeV]。 */
    double r_guess   = 2.045;     /*!< \brief 闭合轨道半径初值 [m]。 */
    double div       = 0.0;       /*!< \brief 发射角散半径 [rad]（默认 0 = 严格按规格；
                                       工作点由独立微分探针测量，不依赖此值）。 */
    int    ksc       = 4;         /*!< \brief SC 场采样间隔（步）。 */
    double sc_tol    = 2.0e-3;    /*!< \brief 场图收敛判据（相对峰值）。 */
    int    env_stride = 5;        /*!< \brief 包络记录步长间隔。 */
    int    env_npick  = 20000;    /*!< \brief 包络子样本数。 */
    int    traj_n     = 2000;     /*!< \brief 轨迹导出子集数。 */
    int    traj_stride = 20;      /*!< \brief 轨迹记录步长间隔。 */
    int    nthreads   = 0;        /*!< \brief OpenMP 线程数（0=全部）。 */
    bool   full       = false;
};

static void print_usage()
{
    std::printf(
        "用法: cycl_sc_beam [--full] [--n N] [--turns T] [--iters K] [--nseg M]\n"
        "                   [--nsteps S] [--len L(m)] [--nr R(mm)] [--current I(mA)]\n"
        "                   [--div θ(mrad)]\n"
        "       --full 等价于 --n 1000000 --turns 3 --iters 6 --nseg 16 --nsteps 4000\n" );
}

static Params parse_args( int argc, char **argv )
{
    Params p;
    for( int a = 1; a < argc; ++a ) {
        std::string s( argv[a] );
        auto next = [&]() -> const char * {
            if( a+1 >= argc ) {
                std::printf( "缺少参数值: %s\n", s.c_str() );
                std::exit( 1 );
            }
            return( argv[++a] );
        };
        if( s == "--full" ) {
            p.full = true;
            p.nparts = 1000000; p.turns = 3; p.maxit = 6; p.nseg = 16; p.nsteps = 4000;
        } else if( s == "--n" )       p.nparts = std::atoi( next() );
        else if( s == "--turns" )     p.turns = std::atoi( next() );
        else if( s == "--iters" )     p.maxit = std::atoi( next() );
        else if( s == "--nseg" )      p.nseg = std::atoi( next() );
        else if( s == "--nsteps" )    p.nsteps = std::atoi( next() );
        else if( s == "--len" )       p.bunch_len = std::atof( next() );
        else if( s == "--nr" )        p.beam_r = 0.001*std::atof( next() );
        else if( s == "--current" )   p.current = 1.0e-3*std::atof( next() );
        else if( s == "--div" )       p.div = 1.0e-3*std::atof( next() );
        else if( s == "--threads" )   p.nthreads = std::atoi( next() );
        else if( s == "--help" )      { print_usage(); std::exit(0); }
        else {
            std::printf( "未知参数: %s\n", s.c_str() );
            print_usage();
            std::exit( 1 );
        }
    }
    return( p );
}


// ===========================================================================
// 随动坐标系（参考轨道逐波前系）
// ===========================================================================

struct Frame {
    int    n1 = 0;                 /*!< \brief 每圈步数。 */
    int    ntot = 0;               /*!< \brief 总步数（含 t=0）。 */
    double dt = 0.0;               /*!< \brief 时间步 [s]。 */
    double T1 = 0.0;               /*!< \brief 一圈参考时长 [s]。 */
    double vref = 0.0;             /*!< \brief 参考速率 [m/s]。 */
    double bet2 = 0.0;             /*!< \brief β²（横向 SC 力 1-β² 修正用）。 */
    std::vector<Vec3D> R;          /*!< \brief 参考位置（索引 0..ntot）。 */
    std::vector<Vec3D> S;          /*!< \brief 单位切向 ŝ。 */
    std::vector<Vec3D> E1;         /*!< \brief 单位横向 ê1 = normalize(ŝ × ẑ)。 */
};

/*! \brief 由闭合轨道构建逐波前系（参考粒子连续推进 turns 圈）。 */
static Frame build_frame( const CRingFieldMap3D &B, double gamma, double r_guess,
                          int n1, int nturn, double &r_co, double &vr_co,
                          double &T_rev )
{
    const ClosedOrbit co = find_closed_orbit( B, QE, MP, gamma, M_PI/2.0, r_guess, true );
    r_co  = co.r;
    vr_co = co.vr;
    T_rev = co.T_rev;

    const double th0 = M_PI/2.0;
    const double er[2]  = { std::cos(th0), std::sin(th0) };
    const double eth[2] = { -std::sin(th0), std::cos(th0) };
    const double vth = -std::sqrt( co.vmag*co.vmag - co.vr*co.vr );
    Vec3D x( er[0]*co.r, er[1]*co.r, 0.0 );
    Vec3D v( er[0]*co.vr + eth[0]*vth, er[1]*co.vr + eth[1]*vth, 0.0 );

    Frame F;
    F.n1  = n1;
    F.ntot = n1*nturn;
    F.dt  = T_rev/(double)n1;
    F.T1  = F.dt*(double)n1;
    F.vref = norm( v );
    F.bet2 = (F.vref*F.vref)/(CL*CL);

    F.R.resize( F.ntot+1 );
    F.S.resize( F.ntot+1 );
    F.E1.resize( F.ntot+1 );

    CBorisPusher pusher( QE, MP );
    pusher.set_relativistic( true );
    pusher.initialize( Vec3D(0.0,0.0,0.0), &B, x, v, F.dt );

    const Vec3D ez( 0.0, 0.0, 1.0 );
    for( int n = 0; n <= F.ntot; ++n ) {
        F.R[n] = x;
        F.S[n] = v*(1.0/norm(v));
        F.E1[n] = crossp( F.S[n], ez );
        F.E1[n] = F.E1[n]*(1.0/norm(F.E1[n]));
        if( n < F.ntot )
            pusher.step( Vec3D(0.0,0.0,0.0), &B, x, v, F.dt );
    }
    return( F );
}

/*! \brief 参考轨道在 θ0 附近的细采样（发射用）。
 *
 *  步长 dt/16（弧长约 0.8 mm），覆盖 ±span。沿轨道发射时每个粒子必须取
 *  **其自身弧长位置处**的切向；否则 ζ=±50 mm 的粒子会带 ζ/r 量级的假角度
 *  误差（曾实测一圈后放大为 ±100 mm 的假发散，corr(Δr,ζ0)=0.9994）。
 */
struct FinePath {
    std::vector<Vec3D> P;   /*!< \brief 位置。 */
    std::vector<Vec3D> T;   /*!< \brief 单位切向。 */
    double ds = 0.0;        /*!< \brief 相邻样本弧长。 */
    int    ic = 0;          /*!< \brief ζ=0 中心索引。 */
};

static FinePath build_fine_path( const CRingFieldMap3D &B, const Frame &F, double span )
{
    const double dtf = F.dt/16.0;
    const int nside = (int)std::ceil( span/(F.vref*dtf) ) + 2;
    FinePath fp;
    fp.ds = F.vref*dtf;
    fp.ic = nside;
    const int N = 2*nside + 1;
    fp.P.resize( N );
    fp.T.resize( N );

    CBorisPusher pusher( QE, MP );
    pusher.set_relativistic( true );

    // 正向
    {
        Vec3D x = F.R[0], v = F.S[0]*F.vref;
        fp.P[fp.ic] = x;
        fp.T[fp.ic] = F.S[0];
        for( int j = fp.ic+1; j < N; ++j ) {
            pusher.step( Vec3D(0.0,0.0,0.0), &B, x, v, dtf );
            fp.P[j] = x;
            fp.T[j] = v*(1.0/norm(v));
        }
    }
    // 反向（负 dt 的 leapfrog，轨迹沿时间反演）
    {
        Vec3D x = F.R[0], v = F.S[0]*F.vref;
        for( int j = fp.ic-1; j >= 0; --j ) {
            pusher.step( Vec3D(0.0,0.0,0.0), &B, x, v, -dtf );
            fp.P[j] = x;
            fp.T[j] = v*(1.0/norm(v));
        }
    }
    return( fp );
}


// ===========================================================================
// 束团系 PIC 盒（IBSimu 原生沉积 + 多重网格）
// ===========================================================================

struct EMap {
    std::vector<float> ex, ey, ez;   /*!< \brief 节点上的 (Eζ, Eξ, Eη)。 */
};

struct ScBox {
    int    nx, ny, nz;      /*!< \brief 轴 0=ζ（纵向）, 1=ξ, 2=η。 */
    double h;               /*!< \brief 网格步长 [m]。 */
    Vec3D  origo;
    Geometry       *geom;
    MeshScalarField *scharge;
    EpotField      *epot;
    EpotMGSolver   *mg;
    double         us_per_node;

    ScBox( int nx_, int ny_, int nz_, double h_, const Vec3D &origo_ )
        : nx(nx_), ny(ny_), nz(nz_), h(h_), origo(origo_), us_per_node(0.0)
    {
        geom = new Geometry( MODE_3D, Int3D(nx,ny,nz), origo, h );
        for( int b = 1; b <= 6; ++b )
            geom->set_boundary( b, Bound( BOUND_DIRICHLET, 0.0 ) );
        geom->build_mesh();
        scharge = new MeshScalarField( *geom );
        epot    = new EpotField( *geom );
        mg      = new EpotMGSolver( *geom );
        mg->set_levels( 4 );
        mg->set_eps( 1.0e-6 );
        mg->set_imax( 500 );
    }

    ~ScBox() {
        delete mg;
        delete epot;
        delete scharge;
        delete geom;
    }

    /*! \brief 沉积（CIC）荷密度（不求解）。 */
    void deposit( const std::vector<Vec3D> &loc, double q )
    {
        scharge->clear();
        for( std::size_t a = 0; a < loc.size(); ++a )
            scharge_add_step_pic( *scharge, q, loc[a] );
        scharge_finalize_step_pic( *scharge );
    }

    /*! \brief 求解当前源，返回耗时 [s]。 */
    double solve()
    {
        const double t0 = now_sec();
        mg->solve( *epot, *scharge );
        return( now_sec() - t0 );
    }

    /*! \brief 沉积 + 求解一次（验证用）。 */
    double solve_for( const std::vector<Vec3D> &loc, double q )
    {
        deposit( loc, q );
        return( solve() );
    }

    /*! \brief 提取 E = -grad(phi) 至 map（中心差分，边界单侧）。 */
    void extract( EMap &map ) const
    {
        map.ex.assign( (std::size_t)nx*ny*nz, 0.0f );
        map.ey.assign( (std::size_t)nx*ny*nz, 0.0f );
        map.ez.assign( (std::size_t)nx*ny*nz, 0.0f );
        auto idx = [&]( int i, int j, int k ) { return( (std::size_t)i + (std::size_t)nx*((std::size_t)j + (std::size_t)ny*(std::size_t)k) ); };
        auto d = [&]( int i, int j, int k, int comp ) -> double {
            int im = i > 0 ? i-1 : i, ip = i < nx-1 ? i+1 : i;
            int jm = j > 0 ? j-1 : j, jp = j < ny-1 ? j+1 : j;
            int km = k > 0 ? k-1 : k, kp = k < nz-1 ? k+1 : k;
            double step = h*(double)( ( (comp==0 ? ip-im : comp==1 ? jp-jm : kp-km) ) );
            switch( comp ) {
            case 0: return( -((*epot)(ip,j,k) - (*epot)(im,j,k))/step );
            case 1: return( -((*epot)(i,jp,k) - (*epot)(i,jm,k))/step );
            default: return( -((*epot)(i,j,kp) - (*epot)(i,j,km))/step );
            }
        };
        for( int k = 0; k < nz; ++k )
            for( int j = 0; j < ny; ++j )
                for( int i = 0; i < nx; ++i ) {
                    const std::size_t a = idx(i,j,k);
                    map.ex[a] = (float)d(i,j,k,0);
                    map.ey[a] = (float)d(i,j,k,1);
                    map.ez[a] = (float)d(i,j,k,2);
                }
    }

    /*! \brief 束团系 (ζ,ξ,η) 处三线性插值取场（返回 (Eζ,Eξ,Eη)）。 */
    Vec3D sample_local( const EMap &map, double zeta, double xi, double eta ) const
    {
        double p[3] = { zeta, xi, eta };
        int    n[3] = { nx, ny, nz };
        const std::vector<float> *arr[3] = { &map.ex, &map.ey, &map.ez };

        double w[3][2];
        int    i0[3];
        for( int a = 0; a < 3; ++a ) {
            double t = (p[a] - origo[a])/h;
            if( t < 0.0 ) t = 0.0;
            double tmax = (double)(n[a]-1);
            if( t > tmax ) t = tmax;
            int i = (int)std::floor( t );
            if( i > n[a]-2 ) i = n[a]-2;
            i0[a] = i;
            w[a][0] = 1.0 - (t - (double)i);
            w[a][1] = t - (double)i;
        }
        Vec3D e( 0.0, 0.0, 0.0 );
        for( int c = 0; c < 3; ++c ) {
            double s = 0.0;
            for( int dk = 0; dk < 2; ++dk )
                for( int dj = 0; dj < 2; ++dj )
                    for( int di = 0; di < 2; ++di ) {
                        const std::size_t a = (std::size_t)(i0[0]+di)
                            + (std::size_t)nx*((std::size_t)(i0[1]+dj)
                            + (std::size_t)ny*(std::size_t)(i0[2]+dk));
                        s += w[0][di]*w[1][dj]*w[2][dk]*(*arr[c])[a];
                    }
            e[c] = s;
        }
        return( e );
    }

    /*! \brief 场图最大/均方值（用于收敛判据）。 */
    static void map_stats( const EMap &m, double &peak, double &rms )
    {
        double p = 0.0, s = 0.0;
        for( std::size_t a = 0; a < m.ex.size(); ++a ) {
            const double mm = (double)m.ex[a]*m.ex[a] + (double)m.ey[a]*m.ey[a] + (double)m.ez[a]*m.ez[a];
            p = std::max( p, mm );
            s += mm;
        }
        peak = std::sqrt( p );
        rms  = std::sqrt( s/(double)m.ex.size() );
    }

    /*! \brief 两张场图的最大/均方差（收敛判据）。 */
    static void map_diff( const EMap &a, const EMap &b, double &dpeak, double &drms )
    {
        double p = 0.0, s = 0.0, norm = 0.0;
        for( std::size_t c = 0; c < a.ex.size(); ++c ) {
            double dx = (double)a.ex[c] - (double)b.ex[c];
            double dy = (double)a.ey[c] - (double)b.ey[c];
            double dz = (double)a.ez[c] - (double)b.ez[c];
            double dd = dx*dx + dy*dy + dz*dz;
            p = std::max( p, dd );
            s += dd;
            double bb = (double)b.ex[c]*b.ex[c] + (double)b.ey[c]*b.ey[c] + (double)b.ez[c]*b.ez[c];
            norm += bb;
        }
        dpeak = std::sqrt( p );
        drms  = std::sqrt( s/std::max(norm,1e-300) );   // 相对均方
    }
};


// ===========================================================================
// 线性映射拟合（最小二乘，正规方程 6x6）
// ===========================================================================

struct LinMapFit {
    bool   ok = false;
    double M[6][6];        /*!< \brief w_i = M_ij u_j。 */
    double nu_r = 0.0, nu_z = 0.0;   /*!< \brief 分数工作点（acos 法）。 */
    double det_r = 0.0, det_z = 0.0;
    double resid = 0.0;
};

/*! \brief 解 6x6 线性方程组（Gauss 消元，部分主元），右端 nrhs 组。 */
static bool solve6( double A[6][6], double B[6][6], double X[6][6] )
{
    int piv[6];
    for( int i = 0; i < 6; ++i ) piv[i] = i;
    for( int c = 0; c < 6; ++c ) {
        int best = c;
        for( int r = c+1; r < 6; ++r )
            if( std::fabs(A[r][c]) > std::fabs(A[best][c]) ) best = r;
        if( std::fabs(A[best][c]) < 1.0e-300 )
            return( false );
        if( best != c ) {
            for( int j = 0; j < 6; ++j ) std::swap( A[c][j], A[best][j] );
            for( int j = 0; j < 6; ++j ) std::swap( B[c][j], B[best][j] );
            std::swap( piv[c], piv[best] );
        }
        const double d = A[c][c];
        for( int r = c+1; r < 6; ++r ) {
            const double f = A[r][c]/d;
            if( f == 0.0 ) continue;
            for( int j = c; j < 6; ++j ) A[r][j] -= f*A[c][j];
            for( int j = 0; j < 6; ++j ) B[r][j] -= f*B[c][j];
        }
    }
    for( int c = 5; c >= 0; --c ) {
        for( int j = 0; j < 6; ++j ) {
            double s = B[c][j];
            for( int t = c+1; t < 6; ++t ) s -= A[c][t]*X[t][j];
            X[c][j] = s/A[c][c];
        }
    }
    return( true );
}

/*! \brief 取第 step 步的局部相空间坐标 (ξ, ξ', η, η', ζ, ζ')。 */
static void local_coords( const Frame &F, int step, const Vec3D &x, const Vec3D &v,
                          double u[6] )
{
    const Vec3D &R = F.R[step], &S = F.S[step], &E1 = F.E1[step];
    Vec3D d = x - R;
    const double zeta = d[0]*S[0] + d[1]*S[1] + d[2]*S[2];
    const double xi   = d[0]*E1[0] + d[1]*E1[1] + d[2]*E1[2];
    const double eta  = d[2];
    const double vs   = v[0]*S[0] + v[1]*S[1] + v[2]*S[2];
    const double v1   = v[0]*E1[0] + v[1]*E1[1] + v[2]*E1[2];
    u[0] = xi;
    u[1] = v1/vs;
    u[2] = eta;
    u[3] = v[2]/vs;
    u[4] = zeta;
    u[5] = vs/F.vref - 1.0;
}

static LinMapFit fit_map( const Frame &F, const std::vector<Vec3D> &x0,
                          const std::vector<Vec3D> &v0, int step0,
                          const std::vector<Vec3D> &x1, const std::vector<Vec3D> &v1,
                          int step1 )
{
    LinMapFit fit;
    const std::size_t N = x0.size();
    // 第一遍：估计各坐标 RMS 尺度。输入输出用同一套尺度归一（D⁻¹MD 的
    // 迹与行列式不变 → 工作点提取不受缩放影响），条件数大幅改善。
    double su2[6] = {0}, sw2[6] = {0};
    for( std::size_t ip = 0; ip < N; ++ip ) {
        double u[6], w[6];
        local_coords( F, step0, x0[ip], v0[ip], u );
        local_coords( F, step1, x1[ip], v1[ip], w );
        for( int i = 0; i < 6; ++i ) { su2[i] += u[i]*u[i]; sw2[i] += w[i]*w[i]; }
    }
    double s[6];
    for( int i = 0; i < 6; ++i ) {
        s[i] = std::sqrt( std::max( su2[i], sw2[i] )/(double)N );
        if( s[i] < 1.0e-30 ) s[i] = 1.0;
    }
    // 第二遍：正规方程（缩放变量）
    double A[6][6] = {{0}}, B[6][6] = {{0}};
    for( std::size_t ip = 0; ip < N; ++ip ) {
        double u[6], w[6];
        local_coords( F, step0, x0[ip], v0[ip], u );
        local_coords( F, step1, x1[ip], v1[ip], w );
        for( int i = 0; i < 6; ++i ) u[i] /= s[i];
        for( int j = 0; j < 6; ++j ) w[j] /= s[j];
        for( int i = 0; i < 6; ++i ) {
            for( int j = 0; j < 6; ++j ) A[i][j] += u[i]*u[j];
            for( int j = 0; j < 6; ++j ) B[i][j] += u[i]*w[j];
        }
    }
    double X[6][6];
    // 岭正则化（相对 1e-10）：处理零方差坐标（如 --div 0 时的 η'）导致的奇异。
    {
        double tr = 0.0;
        for( int i = 0; i < 6; ++i ) tr += A[i][i];
        const double lam = 1.0e-10*tr/6.0;
        for( int i = 0; i < 6; ++i ) A[i][i] += lam;
    }
    if( !solve6( A, B, X ) )
        return( fit );
    for( int j = 0; j < 6; ++j )
        for( int i = 0; i < 6; ++i )
            fit.M[j][i] = X[i][j];     // w_j = Σ_i X[i][j] u_i

    auto tune2x2 = []( double a, double b, double c, double d,
                       double &nu, double &det ) {
        det = a*d - b*c;
        double tr = a + d;
        double disc = 4.0*det - tr*tr;
        if( det <= 0.0 || disc <= 0.0 ) { nu = -1.0; return; }
        double mu = std::atan2( 0.5*std::sqrt(disc), 0.5*tr );
        nu = mu/(2.0*M_PI);
    };
    tune2x2( fit.M[0][0], fit.M[0][1], fit.M[1][0], fit.M[1][1], fit.nu_r, fit.det_r );
    tune2x2( fit.M[2][2], fit.M[2][3], fit.M[3][2], fit.M[3][3], fit.nu_z, fit.det_z );
    fit.ok = true;
    return( fit );
}


// ===========================================================================
// 记录器与推进内核
// ===========================================================================

struct TurnRecorder {
    /*! \brief [nsample][15]（第一圈）：n, Σξ, Σξ², Ση, Ση², Σζ, Σζ², Σdx, Σdy,
     *         Σξ', Σξ'², Σξξ', Ση', Ση'², Σηη'。 */
    std::vector<double> env;
    std::vector<double> traj;    /*!< \brief [ntraj*npts*3] 轨迹子集位置。 */
    int nsamples = 0;
    int ntraj = 0, npts = 0;
    int env_g0 = 0;              /*!< \brief 包络记录所在圈的起点（绝对步）。 */
    int traj_g0 = 0;             /*!< \brief 轨迹记录的起点（绝对步）。 */

    void alloc( int n1, int total_steps, const Params &p, int N )
    {
        nsamples = n1/p.env_stride + 1;
        env.assign( (std::size_t)nsamples*15, 0.0 );
        ntraj = std::min( p.traj_n, N/10 );
        if( ntraj < 1 ) ntraj = 0;
        npts  = total_steps/p.traj_stride + 2;
        if( ntraj ) traj.assign( (std::size_t)ntraj*npts*3, 0.0 );
        env_g0 = 0;
        traj_g0 = 0;
    }
};

static Vec3D sample_sc( const ScBox &box, const std::vector<EMap> &maps,
                        const Frame &F, int g, int n1, int nseg,
                        const Vec3D &x, double fcorr )
{
    const int m = (g % n1)*nseg/n1;
    const Vec3D &R = F.R[g], &S = F.S[g], &E1 = F.E1[g];
    Vec3D d = x - R;
    const double zeta = d[0]*S[0] + d[1]*S[1] + d[2]*S[2];
    const double xi   = d[0]*E1[0] + d[1]*E1[1] + d[2]*E1[2];
    const double eta  = d[2];
    const Vec3D e = box.sample_local( maps[m], zeta, xi, eta );
    // 回到实验室系；横向乘 (1-β²)（运动电荷对的磁场抵消）
    const Vec3D ez( 0.0, 0.0, 1.0 );
    return( S*e[0] + E1*(e[1]*fcorr) + ez*(e[2]*fcorr) );
}

/*! \brief 推进 [step0, step0+nsteps) 步（OpenMP 并行；无跨块同步）。
 *
 *  maps_apply 非空时按 ksc 间隔采样 SC 场（步间冻结）；包络/轨迹写入 rec
 *  （按线程号累加，调用者最后用 reduce_env 归并）。
 */
static void track_chunk( const CRingFieldMap3D &B, const Frame &F,
                         int step0, int nsteps,
                         std::vector<Vec3D> &x, std::vector<Vec3D> &v,
                         const std::vector<EMap> *maps_apply, const ScBox *box,
                         const Params &P, TurnRecorder *rec, bool env_active,
                         std::vector<double> &acc, int nthreads )
{
    const double dt = F.dt;
    const double fcorr = 1.0 - F.bet2;
    const int n1 = F.n1;
    const int N  = (int)x.size();
    const int T  = (env_active && rec) ? rec->nsamples : 0;
    const int pick_stride = std::max( 1, N/P.env_npick );
    const int tstride = (rec && rec->ntraj > 0) ? std::max( 1, N/rec->ntraj ) : 0;

    CBorisPusher pusher( QE, MP );
    pusher.set_relativistic( true );

#pragma omp parallel num_threads(nthreads)
    {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        double *myacc = (T > 0) ? &acc[(std::size_t)tid*T*15] : nullptr;

#pragma omp for schedule(static)
        for( int ip = 0; ip < N; ++ip ) {
            Vec3D xp = x[ip], vp = v[ip];
            Vec3D ep( 0.0, 0.0, 0.0 );
            const bool is_traj = (tstride > 0) && (ip % tstride == 0)
                && (ip/tstride < rec->ntraj);
            const int tline = is_traj ? ip/tstride : -1;

            for( int n = 0; n < nsteps; ++n ) {
                const int g = step0 + n;

                if( maps_apply && (n % P.ksc == 0) )
                    ep = sample_sc( *box, *maps_apply, F, g, n1, P.nseg, xp, fcorr );

                if( is_traj && (g - rec->traj_g0) % P.traj_stride == 0 ) {
                    const int j = (g - rec->traj_g0)/P.traj_stride;
                    if( j < rec->npts ) {
                        double *tp = &rec->traj[((std::size_t)tline*rec->npts + j)*3];
                        tp[0] = xp[0]; tp[1] = xp[1]; tp[2] = xp[2];
                    }
                }

                pusher.step( ep, &B, xp, vp, dt );

                if( myacc && (g >= rec->env_g0)
                    && ((g - rec->env_g0) % P.env_stride == 0)
                    && (ip % pick_stride == 0) ) {
                    const int s = (g - rec->env_g0)/P.env_stride;
                    if( s < T ) {
                        const Vec3D &R = F.R[g+1], &S = F.S[g+1], &E1 = F.E1[g+1];
                        Vec3D d = xp - R;
                        const double zeta = d[0]*S[0] + d[1]*S[1] + d[2]*S[2];
                        const double xi   = d[0]*E1[0] + d[1]*E1[1] + d[2]*E1[2];
                        const double eta  = d[2];
                        const double vs = vp[0]*S[0] + vp[1]*S[1] + vp[2]*S[2];
                        const double v1 = vp[0]*E1[0] + vp[1]*E1[1] + vp[2]*E1[2];
                        const double xip = v1/vs, etap = vp[2]/vs;
                        double *a = myacc + (std::size_t)s*15;
                        a[0]  += 1.0;
                        a[1]  += xi;        a[2]  += xi*xi;
                        a[3]  += eta;       a[4]  += eta*eta;
                        a[5]  += zeta;      a[6]  += zeta*zeta;
                        a[7]  += d[0];      a[8]  += d[1];
                        a[9]  += xip;       a[10] += xip*xip;  a[11] += xi*xip;
                        a[12] += etap;      a[13] += etap*etap; a[14] += eta*etap;
                    }
                }
            }
            x[ip] = xp;
            v[ip] = vp;
        }
    }
}

static void reduce_env( const std::vector<double> &acc, TurnRecorder &rec, int nthreads )
{
    const int T = rec.nsamples;
    for( int t = 0; t < T; ++t )
        for( int a = 0; a < 15; ++a ) {
            double s = 0.0;
            for( int th = 0; th < nthreads; ++th )
                s += acc[((std::size_t)th*T + t)*15 + a];
            rec.env[(std::size_t)t*15 + a] = s;
        }
}

/*! \brief 单扇区（Δθ=2π/8）截面微分映射的工作点测量。
 *
 *  在 θ0 截面取 (r,vr) 与 (z,vz) 两组微扰（能量固定、按**方位角**而非时间
 *  推进 → 纵向自由度被截面消去，无 ζ 滑移耦合污染），差商得 2×2 映射，
 *  ν = 8·acos(tr/(2√det))/(2π)。与 tests/cycl_closed_orbit.cpp §4 同一
 *  已标定方法（该处实测 ν_r=1.263）。SC 用传入冻结图，直接测空间电荷漂移。
 */


struct TuneProbe {
    bool   ok = false;
    double nu_r = -1.0, nu_z = -1.0;
    double det_r = 0.0, det_z = 0.0;   /*!< \brief 复用：相位估计相干度 ∈[0,1]。 */
    double arg_r = 0.0, arg_z = 0.0;   /*!< \brief 保留未用。 */
};

/*! \brief 细采样序列的 FFT 峰值频率（含去趋势 + Hann 窗 + 抛物线细化）。
 *
 *  输入为每 krec 步记录一次的坐标序列；返回真实连续分数频率（每周圈）
 *  与谱显著性（峰值幅 / 带内平均幅）。频率分辨率 ~ n1/(N·krec)，细化后
 *  ~0.1 像素。高频/低频都不混叠（每圈 n1/krec 个样本 ≫ 2·ν）。
 */
static void fft_peak( const std::vector<double> &sig, int n1, int krec,
                      double &nu, double &prom )
{
    const int M = (int)sig.size();
    nu = 0.0; prom = 0.0;
    if( M < 32 ) return;
    // 去趋势：减 2 圈滑动均值（压 ζ-滑移类低频漂移，保 betatron 峰）
    std::vector<double> x( M );
    {
        const int W = std::min( M/4, std::max( 4, 2*n1/krec ) );
        double acc = 0.0;
        for( int n = 0; n < M; ++n ) {
            acc += sig[n];
            if( n >= W ) acc -= sig[n-W];
            x[n] = acc/(double)std::min( n+1, W );
        }
        for( int n = 0; n < W/2; ++n ) x[n] = x[W/2];   // 头部填充
        for( int n = 0; n < M; ++n ) x[n] = sig[n] - x[n];
    }
    int N = 1; while( N < M ) N <<= 1;
    std::vector<double> re( N, 0.0 ), im( N, 0.0 );
    for( int n = 0; n < M; ++n ) {
        const double w = 0.5 - 0.5*std::cos( 2.0*M_PI*(double)n/(double)(M-1) );
        re[n] = x[n]*w;
    }
    // 迭代 Cooley-Tukey
    for( int i = 1, j = 0; i < N; ++i ) {
        int bit = N >> 1;
        for( ; j & bit; bit >>= 1 ) j ^= bit;
        j |= bit;
        if( i < j ) { std::swap( re[i], re[j] ); std::swap( im[i], im[j] ); }
    }
    for( int len = 2; len <= N; len <<= 1 ) {
        const double ang = -2.0*M_PI/(double)len;
        const double wr = std::cos( ang ), wi = std::sin( ang );
        for( int i = 0; i < N; i += len ) {
            double cwr = 1.0, cwi = 0.0;
            for( int k = 0; k < len/2; ++k ) {
                const double ur = re[i+k], ui = im[i+k];
                const double vr = re[i+k+len/2]*cwr - im[i+k+len/2]*cwi;
                const double vi = re[i+k+len/2]*cwi + im[i+k+len/2]*cwr;
                re[i+k] = ur+vr;        im[i+k] = ui+vi;
                re[i+k+len/2] = ur-vr;  im[i+k+len/2] = ui-vi;
                const double nwr = cwr*wr - cwi*wi;
                cwi = cwr*wi + cwi*wr;  cwr = nwr;
            }
        }
    }
    auto amp2 = [&]( int b ) { return re[b]*re[b] + im[b]*im[b]; };
    double best = -1.0;
    int bi = 1;
    for( int b = 2; b < N/2; ++b )
        if( amp2( b ) > best ) { best = amp2( b ); bi = b; }
    double d = 0.0;
    if( bi > 2 && bi < N/2-1 ) {
        const double a1 = std::log( amp2(bi-1) + 1e-300 );
        const double a2 = std::log( amp2(bi  ) + 1e-300 );
        const double a3 = std::log( amp2(bi+1) + 1e-300 );
        const double den = a1 - 2.0*a2 + a3;
        if( std::fabs( den ) > 1.0e-300 ) d = 0.5*(a1 - a3)/den;
    }
    nu = ((double)bi + d)*(double)n1/((double)N*(double)krec);
    double s = 0.0;
    int c = 0;
    for( int b = 2; b < N/2; ++b ) { s += std::sqrt( amp2( b ) ); ++c; }
    prom = std::sqrt( best )/(s/std::max( c, 1 ) + 1.0e-300);
}

static TuneProbe probe_tunes( const CRingFieldMap3D &B, const Frame &F,
                              const std::vector<EMap> *maps, const ScBox *box,
                              const Params &P, int npturns )
{
    // 细采样 FFT 探针（无混叠）：两个 0.1mm 偏置单粒子各跟踪 npturns 圈，
    // 每 krec 步记录一次 (ξ, η) 时间序列 → FFT 峰得**真实连续分数频率**。
    // 与逐圈采样不同，ν 与 1−ν 可区分；对幅度增长、局部近双曲映射、
    // ζ-滑移耦合均免疫；两平面用同一估计器，SC 前后可直接作差。
    if( npturns < 8 ) npturns = 64;
    const int krec = 8;
    const double dx = 1.0e-4;     // 0.1 mm 偏置
    const double fcorr = 1.0 - F.bet2;
    const Vec3D ez( 0.0, 0.0, 1.0 );
    const int mod = (F.ntot > 0) ? F.ntot : F.n1;

    CBorisPusher pusher( QE, MP );
    pusher.set_relativistic( true );

    std::vector<double> sig_r, sig_z;
    sig_r.reserve( (std::size_t)npturns*(F.n1/krec+1) );
    sig_z.reserve( (std::size_t)npturns*(F.n1/krec+1) );
    for( int pl = 0; pl < 2; ++pl ) {
        Vec3D xp = F.R[0] + (pl == 0 ? F.E1[0]*dx : ez*dx);
        Vec3D vp = F.S[0]*F.vref;
        Vec3D ep( 0.0, 0.0, 0.0 );
        for( int t = 0; t < npturns; ++t ) {
            for( int n = 0; n < F.n1; ++n ) {
                const int g = t*F.n1 + n;
                if( maps && (n % P.ksc == 0) )
                    ep = sample_sc( *box, *maps, F, g % mod, F.n1, P.nseg, xp, fcorr );
                pusher.step( ep, &B, xp, vp, F.dt );
                if( (n+1) % krec == 0 ) {
                    double u[6];
                    local_coords( F, (g+1) % mod, xp, vp, u );
                    if( pl == 0 ) sig_r.push_back( u[0] );
                    else          sig_z.push_back( u[2] );
                }
            }
        }
    }

    TuneProbe tp;
    fft_peak( sig_r, F.n1, krec, tp.nu_r, tp.det_r );
    fft_peak( sig_z, F.n1, krec, tp.nu_z, tp.det_z );
    tp.arg_r = tp.arg_z = 0.0;
    tp.ok = tp.nu_r > 1.0e-4 && tp.nu_z > 1.0e-4
        && tp.det_r > 4.0 && tp.det_z > 4.0;
    if( getenv( "CYCL_SC_DEBUG" ) )
        std::printf( "   [probe] ν=(%.5f, %.5f) 谱显著性=(%.1f, %.1f)（%d 圈细采样 FFT）\n",
                     tp.nu_r, tp.nu_z, tp.det_r, tp.det_z, npturns );
    return( tp );
}

/*! \brief 在某步（段起点）用当前粒子分布重建一张段场图。 */
static void build_map_at( ScBox &box, const Frame &F, int g,
                          const std::vector<Vec3D> &x, const std::vector<Vec3D> &v,
                          double qm, std::vector<Vec3D> &locbuf, EMap &out )
{
    (void)v;
    const Vec3D &R = F.R[g], &S = F.S[g], &E1 = F.E1[g];
    const int N = (int)x.size();
    for( int a = 0; a < N; ++a ) {
        Vec3D d = x[a] - R;
        locbuf[a][0] = d[0]*S[0] + d[1]*S[1] + d[2]*S[2];
        locbuf[a][1] = d[0]*E1[0] + d[1]*E1[1] + d[2]*E1[2];
        locbuf[a][2] = d[2];
    }
    box.deposit( locbuf, qm );
    box.solve();
    box.extract( out );
}

/*! \brief 推进一圈（或指定步数）；maps_build 非空时按段分块并在段起点重建图。 */
static void drive_pass( const CRingFieldMap3D &B, const Frame &F,
                        int step0, int nsteps,
                        std::vector<Vec3D> &x, std::vector<Vec3D> &v,
                        const std::vector<EMap> *maps_apply,
                        std::vector<EMap> *maps_build, ScBox *box, double qm,
                        const Params &P, TurnRecorder *rec, bool env_active,
                        int nthreads )
{
    std::vector<double> acc( (std::size_t)nthreads
                             * (size_t)std::max(1,(rec?rec->nsamples:1)) * 15, 0.0 );
    if( !maps_build ) {
        track_chunk( B, F, step0, nsteps, x, v, maps_apply, box, P, rec,
                     env_active, acc, nthreads );
    } else {
        std::vector<Vec3D> locbuf( x.size() );
        const int seglen = F.n1/P.nseg;
        build_map_at( *box, F, step0, x, v, qm, locbuf, (*maps_build)[0] );
        int done = 0;
        for( int m = 0; m < P.nseg; ++m ) {
            const int len = (m == P.nseg-1) ? nsteps-done : seglen;
            track_chunk( B, F, step0+done, len, x, v, maps_apply, box, P, rec,
                         env_active, acc, nthreads );
            done += len;
            if( m+1 < P.nseg )
                build_map_at( *box, F, step0+done, x, v, qm, locbuf, (*maps_build)[m+1] );
        }
    }
    if( env_active && rec )
        reduce_env( acc, *rec, nthreads );
}

/*! \brief 束流矩（相对参考轨迹）。 */
struct BeamMoments {
    double sig_xi = 0.0, sig_eta = 0.0, sig_zeta = 0.0;
    double eps_r = 0.0, eps_z = 0.0;      /*!< \brief RMS 几何发射度 [m·rad]。 */
    double cen_xi = 0.0, cen_eta = 0.0;   /*!< \brief 束团中心横移 [m]。 */
    double p_spread = 0.0;                /*!< \brief |v|/vref-1 的 RMS。 */
};

static BeamMoments beam_moments( const Frame &F, int step,
                                 const std::vector<Vec3D> &x, const std::vector<Vec3D> &v )
{
    BeamMoments b;
    const int N = (int)x.size();
    double s_xi=0, s_xi2=0, s_xip=0, s_xip2=0, s_xixip=0;
    double s_eta=0, s_eta2=0, s_etap=0, s_etap2=0, s_etaetap=0;
    double s_zeta=0, s_zeta2=0, s_psp=0, s_psp2=0;
    for( int a = 0; a < N; ++a ) {
        double u[6];
        local_coords( F, step, x[a], v[a], u );
        s_xi += u[0];  s_xi2 += u[0]*u[0];
        s_xip += u[1]; s_xip2 += u[1]*u[1]; s_xixip += u[0]*u[1];
        s_eta += u[2]; s_eta2 += u[2]*u[2];
        s_etap += u[3]; s_etap2 += u[3]*u[3]; s_etaetap += u[2]*u[3];
        s_zeta += u[4]; s_zeta2 += u[4]*u[4];
        s_psp += u[5]; s_psp2 += u[5]*u[5];
    }
    const double n = (double)N;
    auto var = []( double s1, double s2, double n_ ) {
        double m = s1/n_;
        double vv = s2/n_ - m*m;
        return( vv > 0.0 ? vv : 0.0 );
    };
    b.cen_xi  = s_xi/n;
    b.cen_eta = s_eta/n;
    b.sig_xi  = std::sqrt( var( s_xi, s_xi2, n ) );
    b.sig_eta = std::sqrt( var( s_eta, s_eta2, n ) );
    b.sig_zeta = std::sqrt( var( s_zeta, s_zeta2, n ) );
    {
        const double vx2 = var( s_xi, s_xi2, n ), vxp2 = var( s_xip, s_xip2, n );
        const double cv = s_xixip/n - (s_xi/n)*(s_xip/n);
        b.eps_r = std::sqrt( std::max( 0.0, vx2*vxp2 - cv*cv ) );
    }
    {
        const double vy2 = var( s_eta, s_eta2, n ), vyp2 = var( s_etap, s_etap2, n );
        const double cv = s_etaetap/n - (s_eta/n)*(s_etap/n);
        b.eps_z = std::sqrt( std::max( 0.0, vy2*vyp2 - cv*cv ) );
    }
    b.p_spread = std::sqrt( var( s_psp, s_psp2, n ) );
    return( b );
}

// ===========================================================================
// 导出
// ===========================================================================

/*! \brief 完整状态快照（float32 二进制 + 首行自描述头）与子样本 CSV（局部坐标）。 */
static void write_states( const std::string &tag, int turn, const Params &P,
                          const Frame &F, int step,
                          const std::vector<Vec3D> &x, const std::vector<Vec3D> &v )
{
    (void)P;
    {
        const std::string fn = "cycl_sc_beam_states_" + tag + "_t"
            + std::to_string(turn) + ".bin";
        std::ofstream f( fn.c_str(), std::ios::binary );
        f << "# IBSIMU-CYCL-BEAM-STATE v1 tag=" << tag << " turn=" << turn
          << " step=" << step << " N=" << x.size()
          << " cols=x,y,z,vx,vy,vz(f32,SI)\n";
        std::vector<float> buf( (std::size_t)x.size()*6 );
        for( std::size_t a = 0; a < x.size(); ++a ) {
            buf[6*a+0] = (float)x[a][0]; buf[6*a+1] = (float)x[a][1]; buf[6*a+2] = (float)x[a][2];
            buf[6*a+3] = (float)v[a][0]; buf[6*a+4] = (float)v[a][1]; buf[6*a+5] = (float)v[a][2];
        }
        f.write( (const char*)( buf.empty() ? 0 : &buf[0] ),
                 (std::streamsize)( buf.size()*sizeof(float) ) );
    }
    {
        const std::string fn = "cycl_sc_beam_sample_" + tag + "_t"
            + std::to_string(turn) + ".csv";
        std::ofstream f( fn.c_str() );
        f.precision( 10 );
        f << "xi_m,xi_rad,eta_m,eta_rad,zeta_m\n";
        const int stride = std::max( 1, (int)x.size()/50000 );
        for( int a = 0; a < (int)x.size(); a += stride ) {
            double u[6];
            local_coords( F, step, x[a], v[a], u );
            f << u[0] << "," << u[1] << "," << u[2] << "," << u[3] << "," << u[4] << "\n";
        }
    }
}

/*! \brief 包络/发射度沿圈曲线（第一圈）。 */
static void write_env_csv( const std::string &tag, const TurnRecorder &rec,
                           const Params &P, const Frame &F )
{
    std::ofstream f( ("cycl_sc_beam_envelope_" + tag + ".csv").c_str() );
    f.precision( 10 );
    f << "step,t_ns,sigma_xi_mm,sigma_eta_mm,sigma_zeta_mm,eps_r_m,eps_z_m,cen_xi_mm,cen_eta_mm\n";
    for( int s = 0; s < rec.nsamples; ++s ) {
        const double *a = &rec.env[(std::size_t)s*15];
        if( a[0] < 1.0 )
            continue;
        const int step = rec.env_g0 + s*P.env_stride;
        const double n = a[0];
        auto var = []( double s1, double s2, double n_ ) {
            double m = s1/n_;
            double vv = s2/n_ - m*m;
            return( vv > 0.0 ? vv : 0.0 );
        };
        const double vx2  = var( a[1], a[2], n ),  vxp2 = var( a[9], a[10], n );
        const double vy2  = var( a[3], a[4], n ),  vyp2 = var( a[12], a[13], n );
        const double cvx  = a[11]/n - (a[1]/n)*(a[9]/n);
        const double cvy  = a[14]/n - (a[3]/n)*(a[12]/n);
        const double er   = std::sqrt( std::max(0.0, vx2*vxp2 - cvx*cvx) );
        const double ez   = std::sqrt( std::max(0.0, vy2*vyp2 - cvy*cvy) );
        f << step << "," << (double)(step-rec.env_g0)*F.dt/1e-9 << ","
          << std::sqrt(vx2)*1e3 << "," << std::sqrt(vy2)*1e3 << ","
          << std::sqrt(var(a[5],a[6],n))*1e3 << ","
          << er << "," << ez << ","
          << (a[1]/n)*1e3 << "," << (a[3]/n)*1e3 << "\n";
    }
}

/*! \brief 轨迹子集写盘（HDF5 可用则 .h5，否则 VTK .vtp）。 */
static void write_tracks( const std::string &tag, const TurnRecorder &rec,
                          const Params &P, const Frame &F, int total_steps )
{
    if( rec.ntraj <= 0 )
        return;
    const int jmax = std::min( rec.npts-1, total_steps/P.traj_stride );
    std::vector<std::vector<TrajectoryPoint>> lines( rec.ntraj );
    for( int i = 0; i < rec.ntraj; ++i ) {
        lines[i].resize( jmax+1 );
        for( int j = 0; j <= jmax; ++j ) {
            const double *p = &rec.traj[((std::size_t)i*rec.npts + j)*3];
            lines[i][j].t = (double)(rec.traj_g0 + j*P.traj_stride)*F.dt;
            lines[i][j].x = Vec3D( p[0], p[1], p[2] );
        }
    }
    if( hdf5_available() )
        hdf5_write_polylines( "cycl_sc_beam_tracks_" + tag + ".h5", "tracks", lines );
    else
        vtk_write_polylines( "cycl_sc_beam_tracks_" + tag + ".vtp", lines );
}

/*! \brief 收敛场剖面：过束轴的 ξ / ζ 方向自场。 */
static void write_scprofile( const ScBox &box, const std::vector<EMap> &maps )
{
    std::ofstream f( "cycl_sc_beam_scprofile.csv" );
    f.precision( 10 );
    f << "axis,coord_mm,E_Vpm\n";
    for( int a = -120; a <= 120; ++a ) {          // 沿 ξ：±3 mm
        const double xi = 3.0e-3*(double)a/120.0;
        f << "xi," << xi*1e3 << "," << box.sample_local( maps[0], 0.0, xi, 0.0 )[1] << "\n";
    }
    for( int a = -100; a <= 100; ++a ) {          // 沿 ζ：±束团半长（用第 0 图采样至盒内）
        const double zeta = box.h*(double)(box.nx-1)*0.5*(double)a/100.0;
        f << "zeta," << zeta*1e3 << "," << box.sample_local( maps[0], zeta, 0.0, 0.0 )[0] << "\n";
    }
}

static void write_kv_csv( const std::string &fn,
                          const std::vector<std::pair<std::string,double>> &kv )
{
    std::ofstream f( fn.c_str() );
    f.precision( 12 );
    f << "key,value\n";
    for( std::size_t a = 0; a < kv.size(); ++a )
        f << kv[a].first << "," << kv[a].second << "\n";
}


// ===========================================================================
// 场链验证：直接两两求和 vs 沉积+MG+插值
// ===========================================================================

static void validate_field_chain( ScBox &box, const Params &P, double Qb )
{
    std::printf( "\n--- 场链验证：直接求和 vs (CIC + MG + 插值) ---\n" );

    const int S = 20000;
    const double qs = Qb/(double)S;      // 每宏粒子电荷（与生产同量级）
    std::vector<Vec3D> loc( S );
    unsigned long seed = 4321UL;
    for( int a = 0; a < S; ++a ) {
        const double rr  = P.beam_r*std::sqrt( lcg(seed) );
        const double ph  = 2.0*M_PI*lcg(seed);
        const double zt  = (lcg(seed) - 0.5)*P.bunch_len;
        loc[a] = Vec3D( zt, rr*std::cos(ph), rr*std::sin(ph) );
    }

    const double tsolve = box.solve_for( loc, qs );
    EMap map;
    box.extract( map );
    std::printf( "  盒 %dx%dx%d h=%.3g mm  求解 %.2f s\n",
                 box.nx, box.ny, box.nz, box.h*1e3, tsolve );

    // 探针点（束团系）：(ζ, ξ, η)
    const double probes[][3] = {
        { 0.0,      0.5e-3, 0.0 },
        { 0.0,     -0.5e-3, 0.0 },
        { 0.0,      1.0e-3, 0.0 },
        { 0.0,     -1.0e-3, 0.0 },
        { 0.025,    1.0e-3, 0.0 },
        { 0.050,    1.0e-3, 0.0 },
        { 0.0,      0.0,    1.0e-3 },
        { 0.0,      1.5e-3, 0.0 },
    };
    const int NP = (int)(sizeof(probes)/sizeof(probes[0]));

    double worst = 0.0;
    std::printf( "  probe(ζ,ξ,η)[mm]      E_direct(V/m)   E_mesh(V/m)   相对差\n" );
    for( int a = 0; a < NP; ++a ) {
        Vec3D p( probes[a][0], probes[a][1], probes[a][2] );
        Vec3D ed( 0.0, 0.0, 0.0 );
        // 软化（ε=h）两两求和：与 CIC 沉积的平滑尺度一致，比较的是平滑场。
        // 未软化时 r ≲ 0.5 mm 的探针会被最近邻离散奇点污染（曾误报 45% 差异）。
        const double eps2 = box.h*box.h;
        for( int b = 0; b < S; ++b ) {
            Vec3D d = p - loc[b];
            const double r2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2] + eps2;
            const double r3 = r2*std::sqrt(r2);
            ed = ed + d*(qs/(4.0*M_PI*EPS0*r3));
        }
        Vec3D em = box.sample_local( map, p[0], p[1], p[2] );
        const double md = norm( ed ), mm = norm( em );
        const double rel = (md > 1.0) ? std::fabs( mm-md )/md : std::fabs( mm-md )/std::max( mm, 1.0 );
        if( md > 10.0 )
            worst = std::max( worst, rel );
        std::printf( "  (%+6.2f,%+5.2f,%+5.2f)   %12.4g  %12.4g    %6.2f%%\n",
                     p[0]*1e3, p[1]*1e3, p[2]*1e3, md, mm, 100.0*rel );
    }

    // 无限长线电荷参考值（单位/量级检查）
    const double lambda = Qb/P.bunch_len;
    const double e_line = lambda/(2.0*M_PI*EPS0*1.0e-3);
    std::printf( "  无限长线电荷参考：λ=Qb/L=%.4g C/m,  E_r(1mm)=%.4g V/m\n", lambda, e_line );

    const double tol = 0.25;
    check( worst < tol, "场链与直接求和一致（束内探针 < 25%，含盒壁截断与 CIC 误差）" );
}


// ===========================================================================
// 主程序
// ===========================================================================

int main( int argc, char **argv )
{
    const double t_begin = now_sec();
    Params P = parse_args( argc, argv );

    int nthreads = P.nthreads;
#ifdef _OPENMP
    if( nthreads <= 0 ) nthreads = omp_get_max_threads();
#else
    nthreads = 1;
#endif

    std::printf( "IBSimu-Cycl 注入束流空间电荷自洽模拟（PSI Ring, 72 MeV, 实场）\n" );
    std::printf( "N=%d 圈数=%d 迭代上限=%d nseg=%d 每圈步数=%d 线程=%d\n",
                 P.nparts, P.turns, P.maxit, P.nseg, P.nsteps, nthreads );

    // ---- 场图 ----
    std::vector<std::string> cands;
    cands.push_back( "../examples/cyclotron/data/bfield.dat" );
    cands.push_back( "examples/cyclotron/data/bfield.dat" );
    std::unique_ptr<CRingFieldMap3D> map;
    std::string mapfile;
    for( std::size_t c = 0; c < cands.size() && !map; ++c ) {
        try {
            map.reset( new CRingFieldMap3D( read_ring_field_map3d( cands[c], 0.1, 8, 0 ) ) );
            mapfile = cands[c];
        } catch( const std::exception & ) { }
    }
    if( !map ) {
        std::printf( "[SKIP] 场图不可用（需要 examples/cyclotron/data/bfield.dat）\n" );
        return( 0 );
    }
    std::printf( "场图: %s\n", mapfile.c_str() );

    // ---- 束流参数 ----
    const double gamma = 1.0 + P.ke_mev*1.0e6*QE/(MP*CL*CL);
    const double Qb    = P.current/P.f_bunch;          // 单束团电荷 [C]
    const double qm    = Qb/(double)P.nparts;          // 宏粒子电荷 [C]
    const double protons_per_macro = qm/QE;
    std::printf( "γ=%.6f  Qb=%.6g C (%.2f pC)  宏粒子=%.4g C (%.1f 质子)\n",
                 gamma, Qb, Qb*1e12, qm, protons_per_macro );

    // 相邻束团贡献估计（可忽略性检查）
    {
        const double dsep = 2.0*M_PI*P.r_guess/6.0;
        const double e_far = Qb/(4.0*M_PI*EPS0*dsep*dsep);
        std::printf( "相邻束团距离 %.2f m，其场贡献 %.3g V/m（自场为 kV/m 量级）\n",
                     dsep, e_far );
    }

    // ---- 参考轨道与随动系 ----
    const double t_frame0 = now_sec();
    double r_co, vr_co, T_rev;
    Frame F = build_frame( *map, gamma, P.r_guess, P.nsteps, P.turns, r_co, vr_co, T_rev );
    std::printf( "闭合轨道: r=%.4f m  vr=%+.3e m/s  T_rev=%.3f ns  dt=%.3f ps\n",
                 r_co, vr_co, T_rev/1e-9, F.dt/1e-12 );
    std::printf( "参考系构建耗时 %.1f s（%d 步）\n", now_sec()-t_frame0, F.ntot );
    {
        const double closure = norm( F.R[F.n1] - F.R[0] );
        std::printf( "参考轨道整圈闭合残差 = %.3e m\n", closure );
    }

    // ---- 场图剖面探针（诊断：径向梯度 / 方位波纹 / 纵向结构）----
    if( getenv( "CYCL_SC_BPROFILE" ) ) {
        std::printf( "\n--- 场图剖面（θ=0, z=0） ---\n" );
        for( double r = 1.70; r <= 2.551; r += 0.05 ) {
            const Vec3D b = (*map)( Vec3D( r, 0.0, 0.0 ) );
            std::printf( "  r=%.3f m  B=(%+.4f, %+.4f, %+.4f) T  |B|=%.4f\n",
                         r, b[0], b[1], b[2], norm(b) );
        }
        std::printf( "--- 方位角波纹（r=%.3f, z=0） ---\n", r_co );
        for( int k = 0; k < 16; ++k ) {
            const double th = k*(2.0*M_PI/16.0);
            const Vec3D b = (*map)( Vec3D( r_co*std::cos(th), r_co*std::sin(th), 0.0 ) );
            std::printf( "  θ=%.3f  |B|=%.5f  Bz=%.5f\n", th, norm(b), b[2] );
        }
        std::printf( "--- z 结构（r=%.3f, θ=0） ---\n", r_co );
        for( int k = 0; k <= 6; ++k ) {
            const double z = 0.001*(double)k;
            const Vec3D b = (*map)( Vec3D( r_co, 0.0, z ) );
            std::printf( "  z=%+.3f mm  |B|=%.5f  Bz=%.5f\n", z*1e3, norm(b), b[2] );
        }
    }

    // ---- 随动盒 ----
    // 纵向：束团 ±L/2 + 边距；横向：±4 mm。h=0.25 mm。
    const double margin = 0.050;
    const double half   = P.bunch_len/2.0 + margin;
    const int    nx = 2*(int)std::lround( half/0.00025 ) + 1;
    const int    nyt = 2*(int)std::lround( 0.004/0.00025 ) + 1;
    ScBox box( nx, nyt, nyt, 0.00025, Vec3D( -half, -0.004, -0.004 ) );
    std::printf( "SC 盒: %dx%dx%d 节点 (%.1f MB/图)\n", nx, nyt, nyt,
                 (double)nx*nyt*nyt*3*sizeof(float)/1e6 );

    // ---- 场链验证 ----
    validate_field_chain( box, P, Qb );

    // ---- 束流初始化 ----
    // 关键：沿参考轨道发射——每个粒子的初速度必须取**其自身弧长位置处**的
    // 轨道切向（细采样插值；见 build_fine_path 说明）。
    std::vector<Vec3D> x0( P.nparts ), v0( P.nparts );
    {
        const FinePath fp = build_fine_path( *map, F, 0.5*P.bunch_len + 0.020 );
        unsigned long seed = 20261010UL;
        const Vec3D ez( 0.0, 0.0, 1.0 );
        for( int a = 0; a < P.nparts; ++a ) {
            const double rr = P.beam_r*std::sqrt( lcg(seed) );
            const double ph = 2.0*M_PI*lcg(seed);
            const double xi  = rr*std::cos(ph);
            const double eta = rr*std::sin(ph);
            const double zt  = (lcg(seed) - 0.5)*P.bunch_len;
            const double rr2 = P.div*std::sqrt( lcg(seed) );
            const double ph2 = 2.0*M_PI*lcg(seed);
            const double th_xi  = rr2*std::cos(ph2);
            const double th_eta = rr2*std::sin(ph2);

            double t = (double)fp.ic + zt/fp.ds;
            int j0 = (int)std::floor( t );
            if( j0 < 0 ) j0 = 0;
            if( j0 > (int)fp.P.size()-2 ) j0 = (int)fp.P.size()-2;
            const double fr = t - (double)j0;
            const Vec3D base = fp.P[j0] + (fp.P[j0+1]-fp.P[j0])*fr;
            Vec3D tan = fp.T[j0] + (fp.T[j0+1]-fp.T[j0])*fr;
            tan = tan*(1.0/norm(tan));
            Vec3D e1 = crossp( tan, ez );
            e1 = e1*(1.0/norm(e1));
            x0[a] = base + e1*xi + ez*eta;
            Vec3D vt = tan + e1*th_xi + ez*th_eta;
            v0[a] = vt*(F.vref/norm(vt));
        }
    }
    {
        const BeamMoments bmi = beam_moments( F, 0, x0, v0 );
        std::printf( "发射状态: σξ=%.4f mm  ση=%.4f mm  σζ=%.4f mm（零能散/零发射度）\n",
                     bmi.sig_xi*1e3, bmi.sig_eta*1e3, bmi.sig_zeta*1e3 );
    }

    // ---- 段长整除性 & 记录器 ----
    while( P.nseg > 1 && (F.n1 % P.nseg) != 0 )
        --P.nseg;                       // 段长必须整除，否则调整
    TurnRecorder rec_nosc, rec_sc;
    rec_nosc.alloc( F.n1, F.n1*P.turns, P, P.nparts );
    rec_sc.alloc(   F.n1, F.n1*P.turns, P, P.nparts );

    // 发射初始状态快照（两种运行共用）
    write_states( "init", 0, P, F, 0, x0, v0 );

    // ---- 无 SC 参考运行（第 1 圈顺带构建 maps_1）----
    std::printf( "\n=== 无空间电荷参考运行（%d 圈） ===\n", P.turns );
    std::vector<EMap> maps( P.nseg );
    std::vector<Vec3D> xA = x0, vA = v0;
    std::vector<BeamMoments> mom_nosc;
    std::vector<double> vmag_drift( P.turns, 0.0 );
    for( int t = 0; t < P.turns; ++t ) {
        const double t0 = now_sec();
        drive_pass( *map, F, t*F.n1, F.n1, xA, vA, nullptr,
                    (t == 0 ? &maps : nullptr), &box, qm, P,
                    &rec_nosc, (t == 0), nthreads );
        double dvmax = 0.0;
        for( int a = 0; a < P.nparts; ++a )
            dvmax = std::max( dvmax, std::fabs( norm(vA[a]) - F.vref )/F.vref );
        vmag_drift[t] = dvmax;
        const BeamMoments bm = beam_moments( F, (t+1)*F.n1, xA, vA );
        mom_nosc.push_back( bm );
        std::printf( "  turn %d: %.1f s  σξ=%.4f mm  ε_r=%.3e  |v|漂移=%.1e\n",
                     t+1, now_sec()-t0, bm.sig_xi*1e3, bm.eps_r, dvmax );
        write_states( "nosc", t+1, P, F, (t+1)*F.n1, xA, vA );
    }
    check( vmag_drift[0] < 1.0e-11, "静磁场中 |v| 守恒（Boris 精确性，<1e-11）" );
    const TuneProbe dm_nosc = probe_tunes( *map, F, nullptr, nullptr, P, 0 );
    std::printf( "无 SC 工作点（64 圈细采样 FFT）: ν_r=%.5f（显著性=%.1f） ν_z=%.5f（显著性=%.1f）\n",
                 dm_nosc.nu_r, dm_nosc.det_r, dm_nosc.nu_z, dm_nosc.det_z );
    check( dm_nosc.ok, "无 SC 工作点探针（64 圈细采样 FFT，无混叠）" );
    write_tracks( "nosc", rec_nosc, P, F, F.n1*P.turns );
    write_env_csv( "nosc", rec_nosc, P, F );

    // ---- 空间电荷自洽迭代（冻结场/不动点迭代）----
    std::printf( "\n=== 空间电荷自洽迭代（maps_1 由无 SC 分布构建）===\n" );
    struct IterRow { int it; double rpeak, rrms, sxi, seta, nu_r, nu_z, sec; };
    std::vector<IterRow> iters;
    bool converged = false;
    int K = 1;                                  // 当前图编号（maps_K）
    for( int it = 1; it <= P.maxit && !converged; ++it ) {
        std::vector<EMap> maps_new( P.nseg );
        std::vector<Vec3D> xC = x0, vC = v0;
        const double t0 = now_sec();
        drive_pass( *map, F, 0, F.n1, xC, vC, &maps, &maps_new, &box, qm, P,
                    nullptr, false, nthreads );
        const double dsec = now_sec()-t0;

        double dpeak = 0.0, drrms = 0.0, peak = 0.0;
        for( int m = 0; m < P.nseg; ++m ) {
            double dp, dr, pk, rm;
            ScBox::map_diff( maps_new[m], maps[m], dp, dr );
            ScBox::map_stats( maps_new[m], pk, rm );
            dpeak = std::max( dpeak, dp );
            drrms = std::max( drrms, dr );
            peak  = std::max( peak, pk );
        }
        const double rpeak = dpeak/std::max( peak, 1.0e-30 );
        const BeamMoments bm = beam_moments( F, F.n1, xC, vC );
        const TuneProbe dm = probe_tunes( *map, F, &maps_new, &box, P, 0 );
        iters.push_back( IterRow{ it, rpeak, drrms, bm.sig_xi, bm.sig_eta, dm.nu_r, dm.nu_z, dsec } );
        std::printf( "  iter %d: |ΔE|峰值/|E|峰值=%.3e  RMS=%.3e  ν=(%+.5f, %+.5f)  σξ=%.4f ση=%.4f mm  (%.0f s)\n",
                     it, rpeak, drrms, dm.nu_r, dm.nu_z, bm.sig_xi*1e3, bm.sig_eta*1e3, dsec );
        maps = maps_new;
        K = it+1;
        if( rpeak < P.sc_tol )
            converged = true;
    }
    {
        double pk, rm;
        ScBox::map_stats( maps[0], pk, rm );
        std::printf( "SC 场峰值（段 0） = %.4g V/m；迭代 %s（maps_%d）\n",
                     pk, converged ? "已收敛" : "达上限未完全收敛", K );
    }

    // ---- 收敛场下的最终运行 ----
    std::printf( "\n=== 收敛场下最终运行（%d 圈；第 1 圈自洽验证，其余圈冻结图外推）===\n",
                 P.turns );
    std::vector<Vec3D> xS = x0, vS = v0;
    std::vector<BeamMoments> mom_sc;
    for( int t = 0; t < P.turns; ++t ) {
        const double t0 = now_sec();
        drive_pass( *map, F, t*F.n1, F.n1, xS, vS, &maps, nullptr, &box, qm, P,
                    &rec_sc, (t == 0), nthreads );
        const BeamMoments bm = beam_moments( F, (t+1)*F.n1, xS, vS );
        mom_sc.push_back( bm );
        std::printf( "  turn %d: %.1f s  σξ=%.4f mm  ση=%.4f mm  ε_r=%.3e m·rad\n",
                     t+1, now_sec()-t0, bm.sig_xi*1e3, bm.sig_eta*1e3, bm.eps_r );
        write_states( "sc", t+1, P, F, (t+1)*F.n1, xS, vS );
    }
    const TuneProbe dm_sc = probe_tunes( *map, F, &maps, &box, P, 0 );
    std::printf( "收敛场工作点（64 圈细采样 FFT）: ν_r=%.5f（显著性=%.1f） ν_z=%.5f（显著性=%.1f）\n",
                 dm_sc.nu_r, dm_sc.det_r, dm_sc.nu_z, dm_sc.det_z );
    write_tracks( "sc", rec_sc, P, F, F.n1*P.turns );
    write_env_csv( "sc", rec_sc, P, F );

    // ---- 分析汇总 ----
    std::printf( "\n================ 束流动力学分析（第 1 圈）================\n" );
    const BeamMoments &b0 = mom_nosc[0];
    const BeamMoments &bs = mom_sc[0];
    const double dnu_r = dm_sc.nu_r - dm_nosc.nu_r;
    const double dnu_z = dm_sc.nu_z - dm_nosc.nu_z;
    const bool ok_r = dm_nosc.ok && dm_sc.ok;
    const bool ok_z = dm_nosc.ok && dm_sc.ok;
    std::printf( "  工作点（探针）无 SC: ν_r=%+.5f ν_z=%+.5f\n", dm_nosc.nu_r, dm_nosc.nu_z );
    std::printf( "                有 SC: ν_r=%+.5f ν_z=%+.5f\n", dm_sc.nu_r, dm_sc.nu_z );
    if( ok_r && ok_z )
        std::printf( "  空间电荷工作点漂移 Δν_r=%+.5f (%.2f%%)  Δν_z=%+.5f (%.2f%%)\n",
                     dnu_r, 100.0*dnu_r/dm_nosc.nu_r, dnu_z, 100.0*dnu_z/dm_nosc.nu_z );
    else
        std::printf( "  空间电荷工作点漂移：不可测\n" );
    std::printf( "  包络(σ)     无 SC: σξ=%.4f mm ση=%.4f mm\n", b0.sig_xi*1e3, b0.sig_eta*1e3 );
    std::printf( "              有 SC: σξ=%.4f mm ση=%.4f mm（变化 %+.2f%% / %+.2f%%）\n",
                 bs.sig_xi*1e3, bs.sig_eta*1e3,
                 100.0*(bs.sig_xi-b0.sig_xi)/b0.sig_xi, 100.0*(bs.sig_eta-b0.sig_eta)/b0.sig_eta );
    std::printf( "  RMS 发射度  无 SC: ε_r=%.3e ε_z=%.3e m·rad\n", b0.eps_r, b0.eps_z );
    std::printf( "              有 SC: ε_r=%.3e ε_z=%.3e m·rad\n", bs.eps_r, bs.eps_z );
    std::printf( "  束团纵向    无 SC: σζ=%.4f mm  有 SC: σζ=%.4f mm\n", b0.sig_zeta*1e3, bs.sig_zeta*1e3 );
    std::printf( "  动量散度(RMS)：无 SC %.2e → 有 SC %.2e\n", b0.p_spread, bs.p_spread );
    {
        // 交叉检查（仅供参考）：系综最小二乘拟合——横向块受束团内 ζ 关联
        // 污染（两者耦合不在 4D 子空间内），定量工作点以上面的探针为准。
        const LinMapFit lfit = fit_map( F, x0, v0, 0, xS, vS, F.n1 );
        std::printf( "  [参考] 系综拟合 ν_r=%.5f ν_z=%.5f（仅供对照）\n",
                     lfit.nu_r, lfit.nu_z );
    }
    {
        const double d1 = iters.empty() ? 0.0
            : std::fabs( iters.back().sxi - bs.sig_xi )/bs.sig_xi;
        std::printf( "  自洽性：最终运行 σξ 与末次迭代差 %.2e（相对）\n", d1 );
        check( d1 < 0.05, "自洽迭代稳定（最终运行与末次迭代统计一致 <5%）" );
    }
    check( (!ok_r || std::fabs(dnu_r) < 0.25) && (!ok_z || std::fabs(dnu_z) < 0.25),
           "空间电荷工作点漂移量级合理（可测平面 |Δν|<0.25；"
           "1mm/100mm 束团解析预期 0.05-0.2）" );

    // ---- 输出文件 ----
    {
        std::ofstream f( "cycl_sc_beam_convergence.csv" );
        f.precision( 10 );
        f << "iter,rpeak_rel,rrms_rel,sig_xi_mm,sig_eta_mm,nu_r,nu_z,sec\n";
        for( std::size_t a = 0; a < iters.size(); ++a )
            f << iters[a].it << "," << iters[a].rpeak << "," << iters[a].rrms << ","
              << iters[a].sxi*1e3 << "," << iters[a].seta*1e3 << ","
              << iters[a].nu_r << "," << iters[a].nu_z << "," << iters[a].sec << "\n";
    }
    {
        std::ofstream f( "cycl_sc_beam_map_probe.csv" );
        f.precision( 10 );
        f << "run,nu_r,nu_z,prom_r,prom_z,ok\n";
        f << "nosc," << dm_nosc.nu_r << "," << dm_nosc.nu_z << ","
          << dm_nosc.det_r << "," << dm_nosc.det_z << "," << (dm_nosc.ok?1:0) << "\n";
        f << "sc," << dm_sc.nu_r << "," << dm_sc.nu_z << ","
          << dm_sc.det_r << "," << dm_sc.det_z << "," << (dm_sc.ok?1:0) << "\n";
    }
    write_scprofile( box, maps );
    {
        std::vector<std::pair<std::string,double>> kv;
        kv.push_back( std::make_pair( "N", (double)P.nparts ) );
        kv.push_back( std::make_pair( "turns", (double)P.turns ) );
        kv.push_back( std::make_pair( "nseg", (double)P.nseg ) );
        kv.push_back( std::make_pair( "dt_ps", F.dt/1e-12 ) );
        kv.push_back( std::make_pair( "gamma", gamma ) );
        kv.push_back( std::make_pair( "r_co_m", r_co ) );
        kv.push_back( std::make_pair( "vr_co", vr_co ) );
        kv.push_back( std::make_pair( "T_rev_ns", T_rev/1e-9 ) );
        kv.push_back( std::make_pair( "Qb_C", Qb ) );
        kv.push_back( std::make_pair( "qm_C", qm ) );
        kv.push_back( std::make_pair( "bunch_len_m", P.bunch_len ) );
        kv.push_back( std::make_pair( "beam_r_m", P.beam_r ) );
        kv.push_back( std::make_pair( "nu_r_nosc", dm_nosc.nu_r ) );
        kv.push_back( std::make_pair( "nu_z_nosc", dm_nosc.nu_z ) );
        kv.push_back( std::make_pair( "nu_r_sc", dm_sc.nu_r ) );
        kv.push_back( std::make_pair( "nu_z_sc", dm_sc.nu_z ) );
        kv.push_back( std::make_pair( "dnu_r", dnu_r ) );
        kv.push_back( std::make_pair( "dnu_z", dnu_z ) );
        kv.push_back( std::make_pair( "eps_r_nosc", b0.eps_r ) );
        kv.push_back( std::make_pair( "eps_r_sc", bs.eps_r ) );
        kv.push_back( std::make_pair( "eps_z_nosc", b0.eps_z ) );
        kv.push_back( std::make_pair( "eps_z_sc", bs.eps_z ) );
        kv.push_back( std::make_pair( "sig_xi_nosc_mm", b0.sig_xi*1e3 ) );
        kv.push_back( std::make_pair( "sig_xi_sc_mm", bs.sig_xi*1e3 ) );
        kv.push_back( std::make_pair( "sig_eta_nosc_mm", b0.sig_eta*1e3 ) );
        kv.push_back( std::make_pair( "sig_eta_sc_mm", bs.sig_eta*1e3 ) );
        kv.push_back( std::make_pair( "sig_zeta_nosc_mm", b0.sig_zeta*1e3 ) );
        kv.push_back( std::make_pair( "sig_zeta_sc_mm", bs.sig_zeta*1e3 ) );
        kv.push_back( std::make_pair( "iters_done", (double)iters.size() ) );
        kv.push_back( std::make_pair( "converged", converged ? 1.0 : 0.0 ) );
        write_kv_csv( "cycl_sc_beam_summary.csv", kv );
    }

    std::printf( "\n完成（%.1f s）。输出: cycl_sc_beam_*.{csv,bin,h5,vtp}\n", now_sec()-t_begin );
    std::printf( "%s (%d 失败)\n", g_fail ? "FAILED" : "ALL TESTS PASSED", g_fail );
    return( g_fail ? 1 : 0 );
}
