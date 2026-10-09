/** @file cycl_rfgap.cpp
 *  @brief IBSimu-Cycl P3-A 测试：时变场中的 RF 间隙加速。
 *
 *  模型：一维 RF 间隙，电场 E_x(x,t) = E0 * f(x) * cos(omega t + phi0)，
 *  其中 f(x) = 0.5(1+cos(pi x/h)) 于 |x|<h（平滑、C^1），外部为 0，
 *  且 int f dx = h，故有效间隙电压 V0 = E0 * h。
 *
 *  验证：
 *   1. 跟踪得到的能量增益 dW 与**精确渡越积分**
 *      dW = q E0 ∫ f(x) cos(omega t(x) + phi0) dx  一致；
 *   2. RF 相位依赖：phi0=pi 时为减速（dW<0）；
 *   3. 横向无受力（y=z=0 保持）；
 *   4. 场关闭时速率守恒。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <cmath>
#include <cstdio>
#include <string>

#include "boris.hpp"
#include "timevaryingfield.hpp"

using namespace ibsimu_cycl;

static const double QE = 1.602176634e-19;
static const double MP = 1.67262192369e-27;

static int g_failures = 0;
static void check( bool ok, const std::string &msg )
{
    std::printf( "%s  %s\n", ok ? "[ ok ]" : "[FAIL]", msg.c_str() );
    if( !ok ) ++g_failures;
}

/*! \brief 平滑 RF 间隙电场剖面（仅 x 分量）。 */
class CGapEField : public VectorField {
    double _e0, _h;
public:
    CGapEField( double e0, double h ) : _e0(e0), _h(h) {}
    virtual const Vec3D operator()( const Vec3D &x ) const {
        double u = x[0]/_h;
        if( std::fabs(u) >= 1.0 ) return( Vec3D( 0.0, 0.0, 0.0 ) );
        double f = 0.5*(1.0 + std::cos( M_PI*u ));
        return( Vec3D( _e0*f, 0.0, 0.0 ) );
    }
};

static double gap_profile( double x, double h )
{
    if( std::fabs(x) >= h ) return( 0.0 );
    return( 0.5*(1.0 + std::cos( M_PI*x/h )) );
}

/*! \brief 精确渡越积分（Simpson），假设 v 近似常量。 */
static double exact_dW( double e0, double h, double q, double v,
                        double omega, double phi0, double x_entry )
{
    int n = 20000;
    double a = -h, b = h, dx = (b-a)/n, sum = 0.0;
    for( int i = 0; i <= n; ++i ) {
        double x = a + i*dx;
        double t = (x - x_entry)/v;                    // 进入点 x_entry, t=0
        double f = gap_profile( x, h )*std::cos( omega*t + phi0 );
        double w = (i == 0 || i == n) ? 1.0 : ((i % 2) ? 4.0 : 2.0);
        sum += w*f;
    }
    return( q*e0*sum*dx/3.0 );
}

/*! \brief 用 Boris 跟踪穿越间隙，返回 dW。 */
static double track_dW( double e0, double h, double v0, double omega, double phi0,
                        double x_entry, double x_exit, double dt, Vec3D &vout )
{
    CGapEField gap( e0, h );
    CTimeVaryingField efield( &gap, omega, phi0 );
    CBorisPusher pusher( QE, MP );

    Vec3D x( x_entry, 0.0, 0.0 ), v( v0, 0.0, 0.0 );
    int nsteps = (int)std::ceil( (x_exit - x_entry)/v0/dt ) + 10;
    for( int n = 0; n < nsteps; ++n ) {
        efield.set_time( n*dt );
        pusher.step( &efield, 0, x, v, dt );
        if( x[0] >= x_exit ) break;
    }
    vout = v;
    return( 0.5*MP*(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) - 0.5*MP*v0*v0 );
}

int main()
{
    const double h   = 0.02;        // 间隙半宽 [m]
    const double e0  = 1000.0/h;    // V0 = 1000 V -> E0 = 50 kV/m
    const double v0  = 5.0e6;       // 5e6 m/s -> KE ~ 130 keV
    const double fRF = 10.0e6;      // 10 MHz
    const double omega = 2.0*M_PI*fRF;
    const double x_entry = -0.06, x_exit = 0.06;
    const double dt = 2.0e-11;

    double KE = 0.5*MP*v0*v0;
    std::printf( "KE = %.1f keV,  V0 = %.1f V,  f_RF = %.1f MHz,  transit = %.2f ns\n",
                 KE/QE/1e3, e0*h, fRF/1e6, 2.0*h/v0/1e-9 );

    // ---- 1. 能量增益 vs 精确渡越积分 ----
    for( double phi0 : { 0.0, M_PI/2.0, M_PI } ) {
        Vec3D vout;
        double dW_num = track_dW( e0, h, v0, omega, phi0, x_entry, x_exit, dt, vout );
        double dW_ana = exact_dW( e0, h, QE, v0, omega, phi0, x_entry );
        double rel = std::fabs( dW_num - dW_ana )/std::fabs(dW_ana);
        std::printf( "phi0=%.3f rad: dW_num=%.6e eV   dW_analytic=%.6e eV   rel=%.2e\n",
                     phi0, dW_num/QE, dW_ana/QE, rel );
        check( rel < 0.02, "RF gain matches exact transit-time integral (2%)" );
        check( std::fabs( vout[1] ) < 1e-12 && std::fabs( vout[2] ) < 1e-12,
               "no transverse deflection (y=z=0)" );
    }

    // ---- 2. 相位依赖：pi 相位为减速 ----
    {
        Vec3D vout;
        double dW_pi = track_dW( e0, h, v0, omega, M_PI, x_entry, x_exit, dt, vout );
        double dW_0  = track_dW( e0, h, v0, omega, 0.0, x_entry, x_exit, dt, vout );
        std::printf( "dW(0)=%.4f eV  dW(pi)=%.4f eV\n", dW_0/QE, dW_pi/QE );
        check( dW_0 > 0.0 && dW_pi < 0.0, "RF phase: +V accelerates, -V decelerates" );
    }

    // ---- 3. 无场时速率守恒 ----
    {
        CBorisPusher pusher( QE, MP );
        Vec3D x( 0.0, 0.0, 0.0 ), v( 5.0e6, 1.0e6, -2.0e6 );
        double v0m = std::sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
        for( int n = 0; n < 500; ++n )
            pusher.step( 0, 0, x, v, dt );
        double vm = std::sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
        check( std::fabs( vm - v0m )/v0m < 1e-14, "field-free: speed conserved" );
    }

    std::printf( "\n%s (%d failure%s)\n",
                 g_failures ? "FAILED" : "ALL TESTS PASSED",
                 g_failures, g_failures == 1 ? "" : "s" );
    return( g_failures ? 1 : 0 );
}
