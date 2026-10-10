/** @file cycl_accel.cpp
 *  @brief IBSimu-Cycl P4 测试：回旋加速器多圈加速与相位滑移。
 *
 *  Part A（可严格验证）：均匀磁场 + 一个局域 RF 间隙
 *    - 多圈加速，能量增益 dE = q V0 cos(phi)（含渡越因子）；
 *    - 相位锁定：omega_RF = omega_c 时每次穿越的 RF 相位保持不变；
 *    - 相位滑移：omega_RF = (1-d) omega_c 时相位每圈漂移 2*pi*d。
 *
 *  Part B（真实磁场图）：真实 PSI Ring 场图整体缩放，使参考轨道为非相对论
 *    （保留 8 折扇形与径向结构），测量回路频率、能量增益与相位演化。
 *    注：真实 PSI Ring 的绝对场强对应相对论速度，超出当前非相对论 Boris 范围，
 *    故此处用缩放模型演示；相对论推进器见后续工作。
 *
 *  用法: cycl_accel [bfield.dat]
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

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
        }
    }

    std::printf( "\n%s (%d failure%s)\n",
                 g_failures ? "FAILED" : "ALL TESTS PASSED",
                 g_failures, g_failures == 1 ? "" : "s" );
    return( g_failures ? 1 : 0 );
}
