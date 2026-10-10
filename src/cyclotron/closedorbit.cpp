/** @file closedorbit.cpp
 *  @brief 闭合轨道求解的实现（见 closedorbit.hpp 的方法与精度限制说明）。
 *
 *  本实现是从 `tests/cycl_closed_orbit.cpp` 中验证过的算法提取而来，
 *  数值细节（尤其是半步速度插值与 Jacobian 步长量纲）见 WORK_LOG §8.6。
 */

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "closedorbit.hpp"
#include "boris.hpp"
#include "constants.hpp"


namespace ibsimu_cycl {

namespace {

/*! \brief 中平面轨道状态：方位角处的半径与径向速度。 */
struct State {
    double r;
    double vr;
};

/*! \brief 由洛伦兹因子求速率。 */
double v_of_gamma( double g )
{
    return( SPEED_C*std::sqrt( 1.0 - 1.0/(g*g) ) );
}

/*! \brief 在 th0 处构造状态 (r,v_r) 对应的位置与速度（顺时针）。 */
void state_to_xv( double th0, double vmag, const State &s, Vec3D &x, Vec3D &v )
{
    const double c = std::cos( th0 ), sn = std::sin( th0 );
    const double er[2]  = {  c, sn };
    const double eth[2] = { -sn,  c };

    const double vth2 = vmag*vmag - s.vr*s.vr;
    if( vth2 <= 0.0 )
	throw( std::runtime_error( "closedorbit: |v_r| exceeds |v|" ) );
    const double vth = -std::sqrt( vth2 );

    x = Vec3D( er[0]*s.r, er[1]*s.r, 0.0 );
    v = Vec3D( er[0]*s.vr + eth[0]*vth,
	       er[1]*s.vr + eth[1]*vth, 0.0 );
}

/*! \brief 从 th0 处的状态积分到方位角变化 dtheta，返回终点状态。 */
State one_sector( const VectorField &B, const CBorisPusher &pusher,
		  double th0, double vmag, double dt, const State &s0,
		  double dtheta,
		  double *rmin_out = 0, double *rmax_out = 0 )
{
    Vec3D x, v;
    state_to_xv( th0, vmag, s0, x, v );
    // 半步初始化：leapfrog 的标准启动，消除 O(omega*dt) 初始瞬态
    pusher.initialize( 0, &B, x, v, dt );

    double acc = 0.0, th_prev = th0;
    double rmin = s0.r, rmax = s0.r;

    for( int n = 0; n < 400000; ++n ) {
	Vec3D x_prev = x, v_prev = v;
	pusher.step( 0, &B, x, v, dt );

	double th = std::atan2( x[1], x[0] );
	double d = th - th_prev;
	while( d >  M_PI ) d -= 2.0*M_PI;
	while( d < -M_PI ) d += 2.0*M_PI;

	const double rr = std::sqrt( x[0]*x[0] + x[1]*x[1] );
	rmin = std::min( rmin, rr );
	rmax = std::max( rmax, rr );

	const double acc_before = acc;
	acc += d;

	if( (dtheta < 0.0 && acc <= dtheta) || (dtheta > 0.0 && acc >= dtheta) ) {
	    const double f = (dtheta - acc_before)/d;      // in (0,1]
	    const Vec3D xf = x_prev + (x - x_prev)*f;
	    // 速度是**半步量**：v_prev 在 t_n-dt/2、v 在 t_n+dt/2，
	    // 所以要取 t_n+f*dt 处的速度，插值系数必须是 (f+0.5)。
	    // 用 f 会引入 ~a*dt/2 的系统偏差（约 1e5 m/s），
	    // 那会变成 Newton 的残差地板。
	    const Vec3D vf = v_prev + (v - v_prev)*(f + 0.5);
	    const double thf = th0 + dtheta;
	    const double erf[2] = { std::cos(thf), std::sin(thf) };
	    State out;
	    out.r  = std::sqrt( xf[0]*xf[0] + xf[1]*xf[1] );
	    out.vr = vf[0]*erf[0] + vf[1]*erf[1];
	    if( rmin_out ) *rmin_out = rmin;
	    if( rmax_out ) *rmax_out = rmax;
	    return( out );
	}
	th_prev = th;
    }
    throw( std::runtime_error( "closedorbit: did not complete the sector" ) );
}

} // namespace


ClosedOrbit find_closed_orbit( const VectorField &B, double q, double m,
			       double gamma, double th0, double r_guess,
			       bool relativistic )
{
    const double vmag = v_of_gamma( gamma );
    const double DTHETA = -M_PI/4.0;        // 45°：PSI Ring 的 8 折周期

    // 方位平均场 -> 估计回旋周期，用于选 dt
    double Bbar = 0.0;
    const int NS = 720;
    for( int i = 0; i < NS; ++i ) {
	const double th = 2.0*M_PI*(double)i/(double)NS;
	Bbar += B( Vec3D( r_guess*std::cos(th), r_guess*std::sin(th), 0.0 ) )[2];
    }
    Bbar /= (double)NS;
    if( std::fabs( Bbar ) < 1.0e-12 )
	throw( std::runtime_error( "closedorbit: zero average field" ) );

    const double omega = std::fabs(q)*std::fabs(Bbar)/(gamma*m);
    const double dt = (2.0*M_PI/omega)/3000.0;

    CBorisPusher pusher( q, m );
    pusher.set_relativistic( relativistic );

    // ---- 1. 粗扫描定位不动点附近 ----
    // 残差 ~0.2 m 时已超出单扇区映射的线性区，直接 Newton 会在强非线性区振荡；
    // 先扫描拿到好初值，再让 Newton 快速收敛。
    const double R0 = 1.0e-3;      // r 的归一化尺度 [m]
    const double V0 = 5.0e4;       // v_r 的归一化尺度 [m/s]（≈ omega*R0）
    State x { r_guess, 0.0 };
    double best = 1.0e300;

    for( int ir = 0; ir <= 60; ++ir ) {
	const double rr = r_guess*(1.0 - 0.10) + r_guess*0.20*((double)ir/60.0);
	for( int iv = 0; iv <= 60; ++iv ) {
	    const double vv = vmag*(-0.25 + 0.40*((double)iv/60.0));
	    const State s { rr, vv };
	    try {
		const State f = one_sector( B, pusher, th0, vmag, dt, s, DTHETA );
		const double a = (f.r - rr)/R0, b = (f.vr - vv)/V0;
		const double n = std::sqrt( a*a + b*b );
		if( n < best ) {
		    best = n;
		    x = s;
		}
	    } catch( const std::exception & ) {
	    }
	}
    }

    // ---- 2. 信赖域阻尼 Newton ----
    const double DR = 1.0e-4;      // [m]
    const double DV = 2.0e3;       // [m/s]——必须与 1/omega 的响应尺度匹配
    const double MAX_DR = 0.01;    // 信赖域：单步半径上限 [m]
    const double MAX_DV = 1.0e5;   // [m/s]

    int iters = 0;
    for( int it = 0; it < 400; ++it ) {
	iters = it;
	const State f0 = one_sector( B, pusher, th0, vmag, dt, x, DTHETA );
	const double Fr = f0.r - x.r;
	const double Fv = f0.vr - x.vr;
	const double norm = std::sqrt( (Fr/R0)*(Fr/R0) + (Fv/V0)*(Fv/V0) );
	if( std::fabs(Fr) < 1.0e-9 && std::fabs(Fv) < 1.0e-3 )
	    break;

	State xa = x; xa.r  += DR;
	State xb = x; xb.vr += DV;
	const State fa = one_sector( B, pusher, th0, vmag, dt, xa, DTHETA );
	const State fb = one_sector( B, pusher, th0, vmag, dt, xb, DTHETA );

	const double j11 = (fa.r  - f0.r)/DR,  j21 = (fa.vr - f0.vr)/DR;
	const double j12 = (fb.r  - f0.r)/DV,  j22 = (fb.vr - f0.vr)/DV;
	const double det = j11*j22 - j12*j21;
	if( std::fabs(det) < 1.0e-30 )
	    break;

	const double d1 = (-Fr*j22 + Fv*j12)/det;
	const double d2 = (-j11*Fv + j21*Fr)/det;

	double lambda = 1.0;
	if( std::fabs(d1) > MAX_DR )
	    lambda = std::min( lambda, MAX_DR/std::fabs(d1) );
	if( std::fabs(d2) > MAX_DV )
	    lambda = std::min( lambda, MAX_DV/std::fabs(d2) );

	bool improved = false;
	State xn = x;
	for( int ls = 0; ls < 30; ++ls ) {
	    xn.r  = x.r  + lambda*d1;
	    xn.vr = x.vr + lambda*d2;
	    if( xn.r > 0.1*r_guess && xn.r < 3.0*r_guess ) {
		const State fn = one_sector( B, pusher, th0, vmag, dt, xn, DTHETA );
		const double a = (fn.r - xn.r)/R0, b = (fn.vr - xn.vr)/V0;
		if( std::sqrt( a*a + b*b ) < norm ) {
		    improved = true;
		    break;
		}
	    }
	    lambda *= 0.5;
	}
	if( !improved )
	    break;
	x = xn;
    }

    // ---- 3. 一整圈的形状 ----
    double rmn = x.r, rmx = x.r;
    one_sector( B, pusher, th0, vmag, dt, x, -2.0*M_PI, &rmn, &rmx );

    ClosedOrbit co;
    co.r = x.r;
    co.vr = x.vr;
    co.vmag = vmag;
    co.rmin = rmn;
    co.rmax = rmx;
    co.Bbar = Bbar;
    co.iterations = iters;
    return( co );
}


} // namespace ibsimu_cycl
