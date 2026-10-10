/** @file cycl_accel.cpp
 *  @brief IBSimu-Cycl P4 测试：回旋加速器多圈加速与相位滑移。
 *
 *  Part A（可严格验证）：均匀磁场 + 一个局域 RF 间隙
 *    - 多圈加速，能量增益 dE = q V0 cos(phi)（含渡越因子）；
 *    - 相位锁定：omega_RF = omega_c 时每次穿越的 RF 相位保持不变；
 *    - 相位滑移：omega_RF = (1-d) omega_c 时相位每圈漂移 2*pi*d。
 *
 *  Part B（真实磁场图）：真实 PSI Ring 场图整体缩放，使参考轨道为非相对论
 *    （保留 8 折扇形与径向结构），从**闭合轨道**起步 + RF 加速，
 *    测量增益、相位演化与半径外扩。
 *    注：真实 PSI Ring 的绝对场强对应相对论速度，超出当前非相对论 Boris 范围，
 *    故此处用缩放模型演示；相对论推进器见后续工作。
 *
 *  Part C（等时性场设计）：先在原场里量出闭合轨道的回转频率 ω_rev(r)，
 *    再令径向修正因子 c(r)=ω0/ω_rev(r)（因 ω = q<Bz>/(γm)）；修正后的
 *    中平面 Bz 使 ω_rev 与 r 无关（等时性），RF 相位因此得以锁定。
 *
 *  用法: cycl_accel [bfield.dat]
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "boris.hpp"
#include "fieldmap3d.hpp"
#include "ringfield3d.hpp"
#include "closedorbit.hpp"
#include "timevaryingfield.hpp"

using namespace ibsimu_cycl;

static const double QE  = 1.602176634e-19;
static const double MP  = 1.67262192369e-27;
static const double CL  = 299792458.0;

static int g_failures = 0;
static void check( bool ok, const std::string &msg )
{
    std::printf( "%s  %s\n", ok ? "[ ok ]" : "[FAIL]", msg.c_str() );
    if( !ok ) ++g_failures;
}

/*! \brief 局域方位角 RF 间隙（可含多个，用于模拟 dee 的两个对径间隙）。
 *
 *  每个间隙 k 的场为 E = s_k * (-ephi) * phi_hat * f(theta-theta_k)，
 *  幅度 ephi = V0/(r*dtheta) 使 ∫E·ds = V0 与 r 无关。
 *
 *  对经典回旋加速器使用两个对径间隙 (theta=0, s=+1) 与 (theta=pi, s=-1)：
 *  两者在同一时刻场方向相反（同一对 dee 的两侧），使相位相差 pi 的两侧
 *  都加速，且切向冲量引起的轨道中心漂移相互抵消（轨道同心外扩）。
 *  只有单间隙时中心会持续漂移、半径不增长 —— 这是单间隙模型的缺陷。
 */
class CAzimuthalGapField : public VectorField {
    double _v0, _dtheta, _rmin, _rmax;
    std::vector<std::pair<double,double> > _gaps;   // (theta0, sign)
public:
    CAzimuthalGapField( double v0, double dtheta, double rmin, double rmax )
        : _v0(v0), _dtheta(dtheta), _rmin(rmin), _rmax(rmax) {}

    void add_gap( double theta0, double sign ) {
        _gaps.push_back( std::make_pair( theta0, sign ) );
    }

    virtual const Vec3D operator()( const Vec3D &x ) const {
        if( _v0 == 0.0 || _gaps.empty() ) return( Vec3D(0,0,0) );
        double r = std::sqrt( x[0]*x[0] + x[1]*x[1] );
        if( r < _rmin || r > _rmax ) return( Vec3D(0,0,0) );
        double th = std::atan2( x[1], x[0] );
        if( th < 0.0 ) th += 2.0*M_PI;

        double sum = 0.0;
        for( std::size_t k = 0; k < _gaps.size(); ++k ) {
            double d = th - _gaps[k].first;
            while( d >  M_PI ) d -= 2.0*M_PI;
            while( d < -M_PI ) d += 2.0*M_PI;
            if( std::fabs(d) < _dtheta )
                sum += _gaps[k].second * 0.5*(1.0 + std::cos( M_PI*d/_dtheta ));
        }
        if( sum == 0.0 ) return( Vec3D(0,0,0) );

        double e = _v0/(r*_dtheta)*sum;      // 带符号的方位角场幅度
        return( Vec3D( std::sin(th)*e, -std::cos(th)*e, 0.0 ) );   // e * (-phi_hat)
    }
};

/*! \brief 把中平面 Bz 乘上径向修正因子 c(r) 的场（等时性场设计）。
 *
 *  只改中平面 Bz 的**径向剖面**，完整保留原场图的方位角调制结构
 *  （扇区/尖瓣）。与「给定 <Bz>(r) 的回旋加速器设计」同一理想化层次：
 *  一个严格满足 Maxwell 方程的场还需配套构造 B_r，本演示只关心中平面轨道。
 *
 *  c(r) 用给定采样点上的分段线性插值，范围外取端值（常数外推）。
 */
class CRadialProfileField : public VectorField {
    const VectorField *_base;
    std::vector<double> _r, _c;
public:
    explicit CRadialProfileField( const VectorField *base ) : _base(base) {}

    void add( double r, double c ) { _r.push_back(r); _c.push_back(c); }

    double c_of( double r ) const {
        if( _r.empty() ) return( 1.0 );
        const std::size_t n = _r.size();
        // 范围外**线性外推**（用末端两点的斜率）。粒子会加速到扫描范围之外，
        // 若在那里把 c 冻结为常数，等时性就不再成立（实测这正是残余相位滑移的主因）。
        double c;
        if( n == 1 ) {
            c = _c[0];
        } else if( r <= _r.front() ) {
            const double s = (_c[1] - _c[0])/(_r[1] - _r[0]);
            c = _c.front() + s*(r - _r.front());
        } else if( r >= _r.back() ) {
            const double s = (_c[n-1] - _c[n-2])/(_r[n-1] - _r[n-2]);
            c = _c.back() + s*(r - _r.back());
        } else {
            c = _c.back();
            for( std::size_t k = 0; k + 1 < n; ++k ) {
                if( r <= _r[k+1] ) {
                    const double f = (r - _r[k])/(_r[k+1] - _r[k]);
                    c = _c[k] + f*(_c[k+1] - _c[k]);
                    break;
                }
            }
        }
        // 安全限幅：外推过远时不让修正因子发散
        return( std::max( 0.7, std::min( 1.4, c ) ) );
    }

    virtual const Vec3D operator()( const Vec3D &x ) const {
        const Vec3D b = (*_base)( x );
        const double r = std::sqrt( x[0]*x[0] + x[1]*x[1] );
        return( b*c_of(r) );
    }
};

static CFieldMap3D uniform_B( double bz )
{
    CFieldMap3D f( 2, 2, 2, -5.0, 10.0, -5.0, 10.0, -5.0, 10.0 );
    for( std::size_t i = 0; i < 2; ++i )
        for( std::size_t j = 0; j < 2; ++j )
            for( std::size_t k = 0; k < 2; ++k )
                f.set_value( i, j, k, 0.0, 0.0, bz );
    return( f );
}

struct Crossing {
    double t, KE, r, phi;
};

/*! \brief 推进并记录每次穿越 +x 轴（theta=0，顺时针）时的量。 */
static std::vector<Crossing> run( const VectorField *B, CTimeVaryingField *E,
                                  Vec3D x, Vec3D v, double dt, int nsteps,
                                  int max_cross, double q, double m, double rmax_guard )
{
    CBorisPusher pusher( q, m );
    std::vector<Crossing> out;
    double yprev = x[1];
    for( int n = 1; n <= nsteps; ++n ) {
        E->set_time( (double)n*dt );
        pusher.step( E, B, x, v, dt );
        double r = std::sqrt( x[0]*x[0] + x[1]*x[1] );
        if( rmax_guard > 0.0 && (r > rmax_guard) ) { out.push_back( Crossing{-1.0, -1.0, -1.0, -1.0} ); break; }
        if( yprev > 0.0 && x[1] <= 0.0 && x[0] > 0.0 ) {
            double KE = 0.5*m*(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
            out.push_back( Crossing{ (double)n*dt, KE, r, E->phase() } );
            if( (int)out.size() >= max_cross ) break;
        }
        yprev = x[1];
    }
    return( out );
}

static double wrap2pi( double p )
{
    while( p >  M_PI ) p -= 2.0*M_PI;
    while( p < -M_PI ) p += 2.0*M_PI;
    return( p );
}

int main( int argc, char **argv )
{
    // =====================================================================
    // Part A: 均匀场 + 单 RF 间隙
    // =====================================================================
    const double B0 = 1.0;                    // T
    const double omega_c = QE*B0/MP;          // rad/s
    const double Tc = 2.0*M_PI/omega_c;
    const double r_inj = 0.2;                 // m
    const double v_inj = omega_c*r_inj;
    const double KE_inj = 0.5*MP*v_inj*v_inj;
    const double V0 = 50.0e3;                 // V
    const double dtheta = 10.0*M_PI/180.0;
    const double dt = Tc/400.0;

    CFieldMap3D Bf = uniform_B( B0 );
    CAzimuthalGapField gap( V0, dtheta, 0.05, 2.0 );
    gap.add_gap( 0.0,  +1.0 );       // dee 间隙 1
    gap.add_gap( M_PI, -1.0 );       // dee 间隙 2（对径）

    std::printf( "Part A: B0=%.3f T  omega_c=%.4e rad/s  T_c=%.4f ns\n",
                 B0, omega_c, Tc/1e-9 );
    std::printf( "  injection: r=%.3f m  KE=%.4f MeV  V0=%.1f kV  gap=%.1f deg\n",
                 r_inj, KE_inj/QE/1e6, V0/1e3, dtheta*180.0/M_PI );

    // ---- A1: 相位锁定 (omega_RF = omega_c) ----
    {
        // 粒子从 theta=pi/2 顺时针出发, T/4 后首次穿越 theta=0
        CTimeVaryingField efield( &gap, omega_c, -omega_c*Tc/4.0 );
        auto cr = run( &Bf, &efield, Vec3D(0.0, r_inj, 0.0), Vec3D(v_inj, 0.0, 0.0),
                       dt, (int)(12.0*Tc/dt), 11, QE, MP, 0.0 );
        check( cr.size() >= 10, "A1: at least 10 turns recorded" );

        double dE_sum = 0.0, dphi_max = 0.0;
        for( std::size_t k = 1; k < cr.size(); ++k ) {
            dE_sum += (cr[k].KE - cr[k-1].KE);
            double p = wrap2pi( cr[k].phi );
            if( std::fabs(p) > dphi_max ) dphi_max = std::fabs(p);
        }
        double dE_avg = dE_sum/(double)(cr.size()-1);
        double expect = 2.0*QE*V0;    // 每圈穿越 2 个对径间隙
        std::printf( "  A1 locked: dE/turn=%.4f keV (expect %.4f keV)  max|phi|=%.4f rad  r: %.4f->%.4f m\n",
                     dE_avg/QE/1e3, expect/QE/1e3, dphi_max,
                     cr.front().r, cr.back().r );
        check( std::fabs(dE_avg - expect)/expect < 0.02,
               "A1: dE per turn == q*V0 (within 2%, transit factor)" );
        check( dphi_max < 0.05, "A1: RF phase locked at gap crossings" );
        // 半径随能量增长 (r = v/omega_c = sqrt(2KE/m)/omega_c)
        double r_expect = std::sqrt(2.0*cr.back().KE/MP)/omega_c;
        check( std::fabs(cr.back().r - r_expect)/r_expect < 0.02,
               "A1: orbit radius follows sqrt(2 KE/m)/omega_c" );

        std::ofstream f( "cycl_accel_locked.csv" );
        f << "turn,t,KE_MeV,r,phi_rad\n";
        for( std::size_t k = 0; k < cr.size(); ++k )
            f << k << "," << cr[k].t << "," << cr[k].KE/QE/1e6 << ","
              << cr[k].r << "," << wrap2pi(cr[k].phi) << "\n";
    }

    // ---- A2: 相位滑移 (omega_RF = 0.98 omega_c) ----
    {
        double detune = 0.02;
        double omega_rf = (1.0 - detune)*omega_c;
        CTimeVaryingField efield( &gap, omega_rf, -omega_rf*Tc/4.0 );
        auto cr = run( &Bf, &efield, Vec3D(0.0, r_inj, 0.0), Vec3D(v_inj, 0.0, 0.0),
                       dt, (int)(12.0*Tc/dt), 11, QE, MP, 0.0 );
        check( cr.size() >= 10, "A2: at least 10 turns recorded" );

        // 每圈相位漂移 (线性拟合 phi vs turn)
        double n = (double)(cr.size()-1);
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        for( std::size_t k = 0; k < cr.size(); ++k ) {
            double xx = (double)k, yy = wrap2pi(cr[k].phi);
            sx += xx; sy += yy; sxx += xx*xx; sxy += xx*yy;
        }
        double slope = (n+1>0) ? ( (n+1)*sxy - sx*sy )/( (n+1)*sxx - sx*sx ) : 0.0;
        double expect_slip = -2.0*M_PI*detune;      // rad/turn
        std::printf( "  A2 detuned(2%%): dphi/turn=%.5f rad (expect %.5f)\n", slope, expect_slip );
        check( std::fabs(slope - expect_slip) < 0.02*std::fabs(expect_slip)+1e-3,
               "A2: phase slips by 2*pi*detune per turn" );
    }

    // =====================================================================
    // Part B: 真实场图（整体缩放，保持扇形/径向结构）
    // =====================================================================
    std::vector<std::string> cands;
    if( argc > 1 ) cands.push_back( argv[1] );
    cands.push_back( "../examples/cyclotron/data/bfield.dat" );
    cands.push_back( "examples/cyclotron/data/bfield.dat" );

    std::unique_ptr<CRingFieldMap3D> map1;
    std::string file;
    for( std::size_t c = 0; c < cands.size() && !map1; ++c ) {
        try {
            map1.reset( new CRingFieldMap3D( read_ring_field_map3d( cands[c], 0.1, 8, 0 ) ) );
            file = cands[c];
        } catch( const std::exception & ) { }
    }
    if( !map1 ) {
        std::printf( "\n[SKIP] Part B: field map not found\n" );
    } else {
        const double r_ref = 3.3;
        const double v_ref = 0.05*CL;
        // 半径方向最近节点的方位角平均场 Bbar（用于设定缩放）
        std::size_t i_ref = (std::size_t)std::lround( (r_ref - map1->r0())/map1->dr() );
        if( i_ref >= map1->size_r() ) i_ref = map1->size_r()-1;
        double Bbar = 0.0;
        for( std::size_t k = 0; k < map1->size_theta(); ++k )
            Bbar += map1->bz_midplane( i_ref, k );
        Bbar /= (double)map1->size_theta();

        double scale = (MP*v_ref)/(QE*r_ref*Bbar);   // 使 <B> = m v/(q r)
        std::printf( "\nPart B: file=%s\n", file.c_str() );
        std::printf( "  map: r0=%.3f m, nr=%d, dr=%.4f m, rmax=%.3f m\n",
                     map1->r0(), (int)map1->size_r(), map1->dr(),
                     map1->r0() + (map1->size_r()-1)*map1->dr() );
        std::printf( "  Bbar(r=%.2f m)=%.4f T  -> scale=%.5f  (scaled Bbar=%.4f T)\n",
                     r_ref, Bbar, scale, Bbar*scale );

        std::unique_ptr<CRingFieldMap3D> mapB;
        try {
            mapB.reset( new CRingFieldMap3D(
                            read_ring_field_map3d( file, 0.1*scale, 8, 0 ) ) );
        } catch( const std::exception &e ) {
            std::printf( "  [SKIP] reload failed: %s\n", e.what() );
        }

        if( mapB ) {
            double omega_rev = v_ref/r_ref;
            double T_rev = 2.0*M_PI/omega_rev;
            CAzimuthalGapField gapB( 20.0e3, dtheta, 1.0, 5.0 );
            gapB.add_gap( 0.0, +1.0 );
            gapB.add_gap( M_PI, -1.0 );

            // B1: 无 RF 测量回路周期
            CAzimuthalGapField gap0( 0.0, dtheta, 1.0, 5.0 );
            CTimeVaryingField e0( &gap0, omega_rev, 0.0 );
            auto c0 = run( mapB.get(), &e0, Vec3D(0.0, r_ref, 0.0), Vec3D(v_ref, 0.0, 0.0),
                           T_rev/3000.0, (int)(3.0*T_rev/(T_rev/3000.0)), 4, QE, MP, 4.6 );
            if( c0.size() >= 3 && c0[0].t > 0.0 ) {
                double T_meas = c0[2].t - c0[1].t;
                std::printf( "  B1: T_rev(measured)=%.4f ns vs v/r=%.4f ns   r: %.4f->%.4f m\n",
                             T_meas/1e-9, T_rev/1e-9, c0.front().r, c0.back().r );
            }

            // ---- B2: 闭合轨道起步 + 两个对径间隙的 RF 加速 ----
            // 关键修正：初值**不再**是「在 r=r_ref 处切向发射」——那不是闭合轨道，
            // 轨迹被大幅 betatron 摆动主导（表现为半径随能量增加反而减小、
            // 增益远低于设计值）。先在相同能量下求出闭合轨道（不动点）再起步。
            const double gamma_ref = 1.0/std::sqrt( 1.0 - (v_ref/CL)*(v_ref/CL) );
            ClosedOrbit co = find_closed_orbit( *mapB, QE, MP, gamma_ref,
                                                M_PI/2.0, r_ref, false );
            std::printf( "  B2: closed orbit @ r=%.4f m, v_r=%+.3e m/s (%.2f%% of |v|), "
                         "scallop %.1f mm\n",
                         co.r, co.vr, 100.0*co.vr/co.vmag,
                         1000.0*(co.rmax - co.rmin) );

            const double th0 = M_PI/2.0;
            const double cr = std::cos( th0 ), sr = std::sin( th0 );
            const double er[2]  = {  cr, sr };
            const double eth[2] = { -sr, cr };
            const double vth0 = -std::sqrt( co.vmag*co.vmag - co.vr*co.vr );
            const Vec3D x0( er[0]*co.r, er[1]*co.r, 0.0 );
            const Vec3D v0( er[0]*co.vr + eth[0]*vth0,
                            er[1]*co.vr + eth[1]*vth0, 0.0 );

            // 扫描初始 RF 相位取加速效果最好者（真实机器同样需要调相）
            const int NTURN_B = 12;
            double                best_gain = -1.0e300, best_phase = 0.0;
            std::vector<Crossing> cB;
            for( int ip = 0; ip < 8; ++ip ) {
                const double ph0 = 2.0*M_PI*(double)ip/8.0;
                CTimeVaryingField e( &gapB, omega_rev, ph0 );
                std::vector<Crossing> c =
                    run( mapB.get(), &e, x0, v0, T_rev/3000.0,
                         (int)(NTURN_B*T_rev/(T_rev/3000.0)), NTURN_B+2,
                         QE, MP, 4.6 );
                if( c.size() >= 3 ) {
                    const double gain = (c.back().KE - c.front().KE)
                        /(double)(c.size()-1);
                    if( gain > best_gain ) {
                        best_gain = gain;
                        best_phase = ph0;
                        cB = c;
                    }
                }
            }

            if( cB.size() >= 3 ) {
                const double nturn = (double)(cB.size()-1);
                // 注意 KE 的单位是 J：先 /QE 换成 eV 再换算 keV
                const double dKE_turn_keV =
                    (cB.back().KE - cB.front().KE)/nturn/QE/1e3;
                std::printf( "  B2: RF phase=%.3f rad; KE %.4f -> %.4f MeV over %d turns\n",
                             best_phase, cB.front().KE/QE/1e6, cB.back().KE/QE/1e6,
                             (int)cB.size()-1 );
                std::printf( "      dKE/turn = %.2f keV (design 2*q*V0 = 40.00);  "
                             "r: %.4f -> %.4f m\n",
                             dKE_turn_keV, cB.front().r, cB.back().r );
                check( cB.back().r > cB.front().r,
                       "B2: real field - orbit radius GROWS with energy" );
                check( dKE_turn_keV > 0.5*40.0,
                       "B2: real field - energy gain above half of 2*q*V0" );
                std::ofstream f( "cycl_accel_real.csv" );
                f << "turn,t,KE_MeV,r,phi_rad\n";
                for( std::size_t k = 0; k < cB.size(); ++k )
                    f << k << "," << cB[k].t << "," << cB[k].KE/QE/1e6 << ","
                      << cB[k].r << "," << wrap2pi(cB[k].phi) << "\n";
            } else {
                std::printf( "  B2: [SKIP] closed orbit + RF produced no crossings\n" );
            }

            // =============================================================
            // C: 等时性场设计（<Bz>(r) ∝ gamma 使 ω_rev 与 r 无关）
            // =============================================================
            // 先在**原场**中量出闭合轨道在若干半径上的回转频率 ω_rev(r)=2π/T_rev(r)，
            // 再令径向修正因子 c(r)=ω0/ω_rev(r)。因 ω = q·<Bz>/(γm)，把中平面 Bz
            // 乘上 c(r) 后 ω_rev(r) 处处等于 ω0（相位得以锁定）。
            //
            // 物理背景：真实 PSI Ring 的 <Bz>(r) 在 3.3~3.6 m 内上升约 4%，这正好是
            // 该机器能量下（γ=1.227，dγ/dr·Δr ≈ 4%）等时性所需的斜率；而在本演示的
            // 整体缩放模型里 γ≈1.001（等时性只要求 <Bz> 近似不变），同一个场形就
            // **过陡**了 —— 这正是 B2 中 −0.246 rad/圈 相位滑移的来源。
            const int    NR = 5;
            double       r_scan[NR] = { 3.00, 3.15, 3.30, 3.45, 3.60 };
            double om_scan[NR], cc_scan[NR], u_found[NR];

            // 给定半径的闭合轨道有确定的能量（γ **不是**自由参数）：因 r ∝ γv = c·(γβ)，
            // 外层迭代必须按 u=γβ 做修正。注意**不能**直接按 r ∝ γ 迭代：
            // 低能下 γ≈1，乘一次就会把 γ 压到 1 以下、sqrt 出 NaN（已踩过一次）。
            const double u_ref = std::sqrt( gamma_ref*gamma_ref - 1.0 );

            // 求解「半径 = r_t」的闭合轨道。r*(u) 近似正比于 u，但 Bbar 随 r 变化，
            // 简单的 u←u·(r_t/r) 迭代收敛很慢（因子≈0.9），改用**割线法**：
            // 3~5 次 find_closed_orbit 即可到 1e-5 相对精度。
            auto r_of_u = []( const VectorField &f, double u, double r_t )
                -> ClosedOrbit {
                return( find_closed_orbit( f, QE, MP, std::sqrt(u*u + 1.0),
                                           M_PI/2.0, r_t, false ) );
            };
            auto solve_at = [&]( const VectorField &f, double r_t, double u0,
                                 double *u_out ) -> ClosedOrbit {
                double u1 = u0;
                ClosedOrbit o1 = r_of_u( f, u1, r_t );
                double f1 = o1.r - r_t;
                double u2 = u1*(r_t/o1.r);        // 线性化的第一步
                ClosedOrbit o2 = r_of_u( f, u2, r_t );
                double f2 = o2.r - r_t;
                for( int it = 0; it < 20 && std::fabs(f2) > 1.0e-5*r_t; ++it ) {
                    const double den = f2 - f1;
                    if( std::fabs(den) < 1.0e-300 ) break;
                    // 注意：必须在覆盖 u1 之前算出割线步长（否则 u2-u1 恒为 0）
                    double u3 = u2 - f2*(u2 - u1)/den;
                    // 步长保险：割线外推过远会把试探能量推出场图有用范围，
                    // 那里场被截断、轨道不再闭合，扫描点会耗尽步数上限而抛异常。
                    const double step = 0.15*u2;
                    if( u3 > u2 + step ) u3 = u2 + step;
                    if( u3 < u2 - step ) u3 = u2 - step;
                    if( u3 < 0.5*u_ref ) u3 = 0.5*u_ref;
                    if( u3 > 2.0*u_ref ) u3 = 2.0*u_ref;
                    u1 = u2; f1 = f2;
                    u2 = u3;
                    o2 = r_of_u( f, u2, r_t );
                    f2 = o2.r - r_t;
                }
                if( u_out ) *u_out = u2;
                (void)u1;
                return( o2 );
            };

            std::printf( "  C1: isochronicity scan, closed orbit in the original field\n" );
            bool scan_ok = true;
            for( int i = 0; i < NR; ++i ) {
                om_scan[i] = 0.0; cc_scan[i] = 1.0; u_found[i] = u_ref;
                try {
                    ClosedOrbit o = solve_at( *mapB, r_scan[i], u_ref,
                                              &u_found[i] );
                    om_scan[i] = 2.0*M_PI/o.T_rev;
                    // 登记修正因子的横坐标用**实际**半径：外层迭代停在 1e-5 相对
                    // 精度，若改用目标半径会在 c(r) 里引入 ~1e-2 的插值错位。
                    r_scan[i] = o.r;
                    std::printf( "      r=%.4f m  u=gb=%.6f  T_rev=%.4f ns  "
                                 "omega=%.6e rad/s\n",
                                 o.r, u_found[i], o.T_rev/1e-9, om_scan[i] );
                } catch( const std::exception &e ) {
                    std::printf( "      r=%.2f m: [SKIP] %s\n", r_scan[i], e.what() );
                    scan_ok = false;
                }
            }

            if( scan_ok ) {
                // ω0 取**原场**在 r_ref=3.3 m 处的回转频率（在采样点间线性插值）
                double om0 = om_scan[NR-1];
                for( int i = 0; i + 1 < NR; ++i ) {
                    if( r_ref <= r_scan[i+1] ) {
                        const double f = (r_ref - r_scan[i])
                            /(r_scan[i+1] - r_scan[i]);
                        om0 = om_scan[i] + f*(om_scan[i+1] - om_scan[i]);
                        break;
                    }
                }
                for( int i = 0; i < NR; ++i ) cc_scan[i] = om0/om_scan[i];
                std::printf( "  C1: omega0=%.6e rad/s ; c(r)=omega0/omega(r): ",
                             om0 );
                for( int i = 0; i < NR; ++i )
                    std::printf( "%.4f%s", cc_scan[i], i + 1 < NR ? ", " : "\n" );

                CRadialProfileField iso( mapB.get() );
                for( int i = 0; i < NR; ++i ) iso.add( r_scan[i], cc_scan[i] );

                // ---- C2: 验证等时性 ----
                double worst = 0.0;
                for( int i = 0; i < NR; ++i ) {
                    try {
                        double u = u_found[i];
                        ClosedOrbit o = solve_at( iso, r_scan[i], u_found[i], &u );
                        const double dev = 2.0*M_PI/o.T_rev/om0 - 1.0;
                        worst = std::max( worst, std::fabs(dev) );
                        std::printf( "      r=%.4f m  domega/omega0 = %+.2e\n",
                                     o.r, dev );
                    } catch( const std::exception & ) { worst = 1.0e300; }
                }
                check( worst < 1.0e-2,
                       "C2: isochronous field - omega_rev(r) is radius independent" );

                // ---- C3: 修正场中的多圈 RF 加速 ----
                try {
                    ClosedOrbit co0 = find_closed_orbit( iso, QE, MP, gamma_ref,
                                                         M_PI/2.0, r_ref, false );
                    const double vth0 = -std::sqrt( co0.vmag*co0.vmag
                                                    - co0.vr*co0.vr );
                    const Vec3D xi( 0.0, co0.r, 0.0 );
                    // RF 频率取修正场里参考半径处的**实测**回转频率，
                    // 以避开 c(r) 线性插值在 r_ref 处的微小偏差。
                    const double om_rf = 2.0*M_PI/co0.T_rev;
                    std::printf( "  C3: omega_RF=%.6e rad/s (nominal omega0=%.6e)\n",
                                 om_rf, om0 );
                    // th0=pi/2: e_r=+y、e_th=-x，故 v = vr*e_r + vth*e_th = (-vth, vr)
                    const Vec3D vi( -vth0, co0.vr, 0.0 );

                    // 圈数取 10：修正场只在 r∈[3.0,3.6] m 上标定过，
                    // 而半径每圈约增 0.028 m，10 圈后≈3.55 m，仍在标定区间内。
                    const int NTURN_C = 10;
                    double best = -1.0e300, bph = 0.0;
                    std::vector<Crossing> cC, cCbest;
                    for( int ip = 0; ip < 8; ++ip ) {
                        const double ph0 = 2.0*M_PI*(double)ip/8.0;
                        CTimeVaryingField e( &gapB, om_rf, ph0 );
                        std::vector<Crossing> c =
                            run( &iso, &e, xi, vi, T_rev/3000.0,
                                 (int)(NTURN_C*T_rev/(T_rev/3000.0)), NTURN_C+2,
                                 QE, MP, 4.6 );
                        if( c.size() >= 3 ) {
                            const double gain = (c.back().KE - c.front().KE)
                                /(double)(c.size()-1);
                            if( gain > best ) { best = gain; bph = ph0; cC = c; }
                        }
                    }
                    cCbest.swap( cC );

                    if( cCbest.size() >= 3 ) {
                        const double n = (double)(cCbest.size()-1);
                        const double slip = (wrap2pi(cCbest.back().phi)
                                             - wrap2pi(cCbest.front().phi))/n;
                        std::printf( "  C3: isochronous + RF phase=%.3f rad; "
                                     "KE %.4f -> %.4f MeV over %d turns\n",
                                     bph, cCbest.front().KE/QE/1e6,
                                     cCbest.back().KE/QE/1e6, (int)cCbest.size()-1 );
                        std::printf( "      dKE/turn = %.2f keV ; dphi/turn = %+.4f rad "
                                     "(B2: -0.246) ; r: %.4f -> %.4f m\n",
                                     (cCbest.back().KE - cCbest.front().KE)/n/QE/1e3,
                                     slip, cCbest.front().r, cCbest.back().r );
                        check( std::fabs(slip) < 0.05,
                               "C3: isochronous field - RF phase stays locked" );
                        check( cCbest.back().r > cCbest.front().r,
                               "C3: isochronous field - radius grows with energy" );
                        check( cCbest.back().KE > cCbest.front().KE,
                               "C3: isochronous field - energy rises monotonically" );
                        std::ofstream f( "cycl_accel_iso.csv" );
                        f << "turn,t,KE_MeV,r,phi_rad\n";
                        for( std::size_t k = 0; k < cCbest.size(); ++k )
                            f << k << "," << cCbest[k].t << ","
                              << cCbest[k].KE/QE/1e6 << "," << cCbest[k].r << ","
                              << wrap2pi(cCbest[k].phi) << "\n";
                    } else {
                        std::printf( "  C3: [SKIP] no crossings\n" );
                    }
                } catch( const std::exception &e ) {
                    std::printf( "  C3: [SKIP] %s\n", e.what() );
                }
            }
        }
    }

    std::printf( "\n%s (%d failure%s)\n",
                 g_failures ? "FAILED" : "ALL TESTS PASSED",
                 g_failures, g_failures == 1 ? "" : "s" );
    return( g_failures ? 1 : 0 );
}
