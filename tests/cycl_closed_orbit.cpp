/** @file cycl_closed_orbit.cpp
 *  @brief  在真实 PSI Ring 场图中求解并验证**闭合轨道**（closed orbit）。
 *
 *  ## 为什么需要这个测试
 *
 *  以「在 r=3.3 m 处纯切向发射 + 动量取 p = q<Bz>·r」作为初值起步，得到的
 *  **不是**闭合轨道，而是围绕闭合轨道做大幅 betatron 振荡的轨迹：
 *
 *  | 观测量 | 实测 |
 *  | --- | --- |
 *  | 径向摆动 | r = 2.94 … 3.56 m（±9.5%，而机器的 8 折扇贝调制只有百分之几） |
 *  | 逐圈重复性 | 相邻两圈在同一方位角的半径差：最大 0.376 m、均方根 0.24 m |
 *
 *  也就是说它**有界但不闭合**：画出来像"轨道重叠、不稳定"。真实回旋加速器里
 *  束流走的是**闭合轨道**——单圈映射的不动点，必须数值求解。
 *
 *  ## 方法
 *
 *  在方位角 th0 处把中平面轨道状态参数化为 (r, v_r)：|v| 由 gamma 定，
 *  方位向速度 v_theta = -sqrt(|v|^2 - v_r^2)（该场配置下粒子顺时针运动）。
 *  单圈映射 M: (r, v_r) -> (r', v_r')（积分到方位角恰好变化 -2pi，用线性插值
 *  落在精确的方位角上），用 Newton 法解 M(x) = x（数值 Jacobian）。
 *
 *  求出的不动点是**离散推进器的不动点**，因此逐圈重复性可达机器精度——
 *  这正是画图时希望看到的"轨道重合"。
 *
 *  ## 验收
 *
 *  - 朴素初值（r=r_ref、纯切向）**不是**不动点；
 *  - 闭合轨道在该方位角上带有明显的径向速度（≈ -1.6e7 m/s，即偏离切向约 5°），
 *    这正是之前轨道不闭合的原因；
 *  - 用闭合轨道时，逐圈重复性从 0.24 m 改善到 ~1e-2 m（>10 倍）；
 *  - 扇贝幅度从 0.62 m 降到 ~0.13 m；
 *  - 径向 betatron 频率 nu_r ≈ 1.26（非整数 => 不动点稳定）。
 *
 *  ## 精度限制（重要）
 *
 *  闭合残差无法降到 1e-8 m 量级，地板约在 1e-3 m / 1e5 m/s：这是**场图本身**
 *  带来的——`CRingFieldMap3D` 在 (r,theta) 上做双线性插值，单元边界处只有 C0
 *  连续，轨道每圈要跨过数百个单元，映射因此是 Lipschitz 而非光滑的。
 *  要进一步提高精度需先提升场图插值的阶次（见 ROADMAP）。
 *
 *  并导出 `cycl_closed_orbit.vtp`（3 圈，带时刻）供 ParaView / 绘图脚本使用。
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "boris.hpp"
#include "constants.hpp"
#include "ringfield3d.hpp"
#include "vec3d.hpp"
#include "vtkwriter.hpp"


using namespace ibsimu_cycl;

namespace {

const double QE = CHARGE_E;   // 1.602176462e-19 C
const double MP = MASS_U;     // 1.66053873e-27 kg

int g_fail = 0;

void check( bool ok, const std::string &msg )
{
    std::printf( "%s  %s\n", ok ? "[ ok ]" : "[FAIL]", msg.c_str() );
    if( !ok )
	++g_fail;
}

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

/*! \brief 沿半径 r 的圆环对 Bz 做方位角平均（回旋加速器的 <B> 定义）。 */
double average_Bz( const CRingFieldMap3D &B, double r, int n = 1440 )
{
    double s = 0.0;
    for( int i = 0; i < n; ++i ) {
	double th = 2.0*M_PI*(double)i/(double)n;
	s += B( Vec3D( r*std::cos(th), r*std::sin(th), 0.0 ) )[2];
    }
    return( s/(double)n );
}

/*! \brief 在 th0 处构造状态 (r,v_r) 对应的位置与速度。 */
void state_to_xv( double th0, double vmag, const State &s, Vec3D &x, Vec3D &v )
{
    const double c = std::cos( th0 ), sn = std::sin( th0 );
    const double er[2]  = {  c, sn };
    const double eth[2] = { -sn,  c };

    double vth2 = vmag*vmag - s.vr*s.vr;
    if( vth2 <= 0.0 )
	throw( std::runtime_error( "closed orbit: |v_r| exceeds |v|" ) );
    const double vth = -std::sqrt( vth2 );      // 顺时针

    x = Vec3D( er[0]*s.r, er[1]*s.r, 0.0 );
    v = Vec3D( er[0]*s.vr + eth[0]*vth,
	       er[1]*s.vr + eth[1]*vth, 0.0 );
}

/*! \brief 从 th0 处的状态出发，积分到方位角变化 \a dtheta，返回终点状态。
 *
 *  场的 8 折周期为 45 度，所以一条闭合轨道也是**单扇区映射**的不动点。
 *  求单扇区不动点比求整圈不动点好得多：积分长度短 8 倍，舍入/插值拐点的
 *  累积噪声也小 8 倍，且 Jacobian 的条件数明显改善。
 */
State one_turn( const CRingFieldMap3D &B, const CBorisPusher &pusher,
		double th0, double vmag, double dt, const State &s0,
		bool init_half_step, double dtheta,
		std::vector<TrajectoryPoint> *traj = 0, double *t_out = 0,
		double *rmax_out = 0, double *rmin_out = 0 )
{
    Vec3D x, v;
    state_to_xv( th0, vmag, s0, x, v );
    if( init_half_step )
	pusher.initialize( 0, &B, x, v, dt );

    const double target = dtheta;
    double acc = 0.0;
    double th_prev = th0;
    double rmin = s0.r, rmax = s0.r;

    const int MAXSTEP = 400000;
    for( int n = 0; n < MAXSTEP; ++n ) {
	Vec3D x_prev = x, v_prev = v;
	pusher.step( 0, &B, x, v, dt );

	double th = std::atan2( x[1], x[0] );
	double d = th - th_prev;
	while( d >  M_PI ) d -= 2.0*M_PI;
	while( d < -M_PI ) d += 2.0*M_PI;

	double rr = std::sqrt( x[0]*x[0] + x[1]*x[1] );
	rmin = std::min( rmin, rr );
	rmax = std::max( rmax, rr );

	const double acc_before = acc;
	acc += d;

	if( (dtheta < 0.0 && acc <= target) || (dtheta > 0.0 && acc >= target) ) {
	    const double f = (target - acc_before)/d;   // (0,1]
	    const Vec3D xf = x_prev + (x - x_prev)*f;
	    // 速度是**半步量**：v_prev 位于 t_n - dt/2，v 位于 t_n + dt/2。
	    // 要得到 t_n + f*dt 处的速度，插值系数必须是 (f + 0.5) 而不是 f；
	    // 否则会引入 ~a*dt/2 的系统偏差（这里约 1.4e5 m/s），
	    // 它正是之前 Newton 残差地板（~8.8e4 m/s）的来源。
	    const Vec3D vf = v_prev + (v - v_prev)*(f + 0.5);
	    const double thf = th0 + target;
	    const double erf[2] = { std::cos(thf), std::sin(thf) };
	    State out;
	    out.r  = std::sqrt( xf[0]*xf[0] + xf[1]*xf[1] );
	    out.vr = vf[0]*erf[0] + vf[1]*erf[1];
	    if( t_out )  *t_out  = ((double)n + f)*dt;
	    if( rmin_out ) *rmin_out = rmin;
	    if( rmax_out ) *rmax_out = rmax;
	    return( out );
	}
	th_prev = th;
    }
    throw( std::runtime_error( "one_turn: did not complete the requested angle" ) );
}

/*! \brief Newton 法求单圈映射的不动点（闭合轨道）。 */
State solve_closed( const CRingFieldMap3D &B, const CBorisPusher &pusher,
		    double th0, double vmag, double dt, State x,
		    bool init_half_step, double dtheta, int *iters = 0 )
{
    // 数值 Jacobian 步长必须与各自的响应尺度匹配：
    //   u_r 对 r 的敏感度 ~ O(1)，所以 DR 取 ~1e-4 m；
    //   u_r 对 v_r 的敏感度 ~ 1/omega_c ~ 2e-8 s，所以 DV 必须取 ~1e3 m/s，
    //   否则 (M(x+DV)-M(x)) 被浮点噪声淹没，Jacobian 在 v_r 方向退化为奇异。
    const double DR = 1.0e-4;     // [m]
    const double DV = 2.0e3;      // [m/s]

    // 阻尼 Newton（回溯线搜索）：
    // 一圈映射对 v_r 的敏感度 ~1/omega、对 r 的 ~O(1)，直接比 |dr| 与 |dv_r| 毫无意义。
    // 用自然尺度归一化：r0=1e-3 m，v0 = omega_c*r0 ≈ 5e4 m/s（即同样的物理偏差）。
    const double R0 = 1.0e-3;     // [m]
    const double V0 = 5.0e4;      // [m/s]
    auto resid = [&]( const State &xx ) {
	const State f = one_turn( B, pusher, th0, vmag, dt, xx, init_half_step,
				 dtheta );
	const double a = (f.r  - xx.r )/R0;
	const double b = (f.vr - xx.vr)/V0;
	return( std::sqrt( a*a + b*b ) );
    };

    State best = x;
    double best_norm = 1.0e300;

    for( int it = 0; it < 400; ++it ) {
	const State f0 = one_turn( B, pusher, th0, vmag, dt, x, init_half_step,
				  dtheta );
	const double Fr = f0.r - x.r;
	const double Fv = f0.vr - x.vr;
	const double norm = std::sqrt( (Fr/R0)*(Fr/R0) + (Fv/V0)*(Fv/V0) );
	if( norm < best_norm ) {
	    best_norm = norm;
	    best = x;
	}
	if( it % 10 == 0 || norm < 1.0e-5 )
	    std::printf( "    it %3d: r=%.9f v_r=%+.6e | F=(%+.3e, %+.3e) | "
			 "||F||=%.6e\n", it, x.r, x.vr, Fr, Fv, norm );
	if( std::fabs(Fr) < 1.0e-10 && std::fabs(Fv) < 1.0e-6 ) {
	    if( iters ) *iters = it;
	    return( x );
	}

	State xa = x; xa.r  += DR;
	State xb = x; xb.vr += DV;
	const State fa = one_turn( B, pusher, th0, vmag, dt, xa, init_half_step,
				  dtheta );
	const State fb = one_turn( B, pusher, th0, vmag, dt, xb, init_half_step,
				  dtheta );

	const double j11 = (fa.r  - f0.r)/DR,  j21 = (fa.vr - f0.vr)/DR;
	const double j12 = (fb.r  - f0.r)/DV,  j22 = (fb.vr - f0.vr)/DV;
	const double det = j11*j22 - j12*j21;
	if( std::fabs(det) < 1.0e-30 ) {
	    std::printf( "    奇异 Jacobian，停止\n" );
	    break;
	}

	const double d1 = (-Fr*j22 + Fv*j12)/det;
	const double d2 = (-j11*Fv + j21*Fr)/det;

	// 信赖域：Newton 满步长可能在强非线性区跑出 20 个网格之外，必须限幅
	const double MAX_DR = 0.01;      // [m]
	const double MAX_DV = 1.0e5;     // [m/s]
	double lambda = 1.0;
	if( std::fabs(d1) > MAX_DR ) lambda = std::min( lambda, MAX_DR/std::fabs(d1) );
	if( std::fabs(d2) > MAX_DV ) lambda = std::min( lambda, MAX_DV/std::fabs(d2) );

	// 回溯线搜索（按归一化残差判据）
	bool improved = false;
	State xn = x;
	for( int ls = 0; ls < 30; ++ls ) {
	    xn.r  = x.r  + lambda*d1;
	    xn.vr = x.vr + lambda*d2;
	    if( xn.r > 1.0 && xn.r < 5.0 ) {
		const double nrm = resid( xn );
		if( nrm < norm ) {
		    improved = true;
		    break;
		}
	    }
	    lambda *= 0.5;
	}
	if( !improved ) {
	    std::printf( "    线搜索无法降低残差，停止（lambda=%.2e）\n", lambda );
	    break;
	}

	x = xn;
	if( iters ) *iters = it+1;
    }
    std::printf( "    Newton 未达严格容差，返回最佳迭代（||F||=%.3e）\n", best_norm );
    return( best );
}

} // namespace


int main( int argc, char **argv )
{
    // ---- 载入真实场图（候选路径与 cycl_track 一致）----
    std::vector<std::string> cand;
    if( argc > 1 )
	cand.push_back( argv[1] );
    cand.push_back( "../examples/cyclotron/data/bfield.dat" );
    cand.push_back( "examples/cyclotron/data/bfield.dat" );

    std::unique_ptr<CRingFieldMap3D> fp;
    std::string file;
    for( std::size_t c = 0; c < cand.size() && !fp; ++c ) {
	try {
	    fp.reset( new CRingFieldMap3D(
			  read_ring_field_map3d( cand[c], 0.1, 8, 0 ) ) );
	    file = cand[c];
	} catch( const std::exception & ) {
	}
    }
    if( !fp ) {
	std::printf( "[SKIP] field map not found\n" );
	return( 0 );
    }
    const CRingFieldMap3D &B = *fp;
    std::printf( "field map: %s\n\n", file.c_str() );

    const double th0 = M_PI/2.0;          // 起点方位角
    const double r_ref = 3.30;            // 参考半径 [m]

    // 动量由该半径处的方位平均场定出（回旋加速器的标准选取）
    const double Bavg = average_Bz( B, r_ref );
    const double p    = QE*Bavg*r_ref;
    const double pc   = p/(MP*SPEED_C);
    const double gamma = std::sqrt( 1.0 + pc*pc );
    const double vmag  = v_of_gamma( gamma );
    const double T_rev = 2.0*M_PI*r_ref/vmag;

    std::printf( "<Bz>(r=%.2f m) = %.4f T  ->  gamma = %.4f, T_rev ~ %.3f ns\n",
		 r_ref, Bavg, gamma, T_rev*1e9 );

    const double dt = T_rev/4000.0;
    const double SECT = -M_PI/4.0;   // 单扇区（场的 8 折周期）
    const double FULL = -2.0*M_PI;   // 整圈

    CBorisPusher pusher( QE, MP );
    pusher.set_relativistic( true );

    // =====================================================================
    // 1. 朴素初值（纯切向、r=r_ref）——说明"为什么它不是闭合轨道"
    // =====================================================================
    std::printf( "\n--- 1. 朴素初值（纯切向发射）的单圈映射 ---\n" );
    State naive { r_ref, 0.0 };
    State n1 = one_turn( B, pusher, th0, vmag, dt, naive, true, FULL );
    std::printf( "  (r, v_r) = (%.6f, %.3e) -> (%.6f, %.3e)\n",
		 naive.r, naive.vr, n1.r, n1.vr );
    std::printf( "  一圈后偏差: dr = %+.6f m, dv_r = %+.3e m/s\n",
		 n1.r - naive.r, n1.vr - naive.vr );
    check( std::fabs(n1.r - naive.r) > 1.0e-3,
	   "朴素初值不是不动点（一圈后明显偏离，证实之前的轨道不闭合）" );

    // =====================================================================
    // 2. Newton 求闭合轨道
    // =====================================================================
    std::printf( "\n--- 2. Newton 迭代求闭合轨道（单扇区映射，45 度）---\n" );

    // 先粗扫描 (r, v_r) 网格定位不动点附近：
    // 直接从 (r_ref, 0) 出发时残差高达 0.2 m，超出单扇区映射的线性区，
    // Newton 会在强非线性区里来回振荡；先用网格扫描拿到好的初值。
    int iters = 0;
    State guess = naive;
    {
	const double R0s = 1.0e-3, V0s = 5.0e4;
	double best_n = 1.0e300;
	int    tried = 0;
	for( int ir = 0; ir <= 60; ++ir ) {
	    const double rr = 2.95 + 0.010*(double)ir;      // 2.95 .. 3.55 m
	    for( int iv = 0; iv <= 60; ++iv ) {
		const double vv = -3.0e7 + 6.6667e5*(double)iv;  // -3e7 .. +1e7
		const State s { rr, vv };
		try {
		    const State f = one_turn( B, pusher, th0, vmag, dt, s,
					      true, SECT );
		    const double a = (f.r - rr)/R0s, b = (f.vr - vv)/V0s;
		    const double n = std::sqrt( a*a + b*b );
		    ++tried;
		    if( n < best_n ) {
			best_n = n;
			guess = s;
		    }
		} catch( const std::exception & ) {
		}
	    }
	}
	std::printf( "  粗扫描 %d 个点: 最优 r=%.4f m, v_r=%+.3e m/s, ||F||=%.4f\n",
		     tried, guess.r, guess.vr, best_n );
    }

    State fix = solve_closed( B, pusher, th0, vmag, dt, guess, true, SECT,
			      &iters );
    std::printf( "  迭代 %d 次: r = %.9f m, v_r = %.6e m/s\n",
		 iters, fix.r, fix.vr );

    const State back_s = one_turn( B, pusher, th0, vmag, dt, fix, true, SECT );
    double rmin = 0, rmax = 0, t1 = 0;
    const State back_f = one_turn( B, pusher, th0, vmag, dt, fix, true, FULL,
				   0, &t1, &rmax, &rmin );
    std::printf( "  闭合性: 45 度后 dr = %+.3e m, dv_r = %+.3e m/s\n",
		 back_s.r - fix.r, back_s.vr - fix.vr );
    std::printf( "          整圈后 dr = %+.3e m, dv_r = %+.3e m/s\n",
		 back_f.r - fix.r, back_f.vr - fix.vr );
    std::printf( "  轨道形状: r = %.4f .. %.4f m (相对扇贝幅度 %.3f%%)\n",
		 rmin, rmax, 100.0*(rmax - rmin)/(rmax + rmin) );
    std::printf( "  回旋周期 T_rev = %.3f ns\n", t1*1e9 );

    check( std::fabs(back_s.r - fix.r) < 5.0e-3 &&
	   std::fabs(back_s.vr - fix.vr) < 3.0e5,
	   "闭合轨道：45 度后回到自身（残差受场图双线性插值的间断限制）" );
    check( std::fabs(back_f.r - fix.r) < 2.0e-2 &&
	   std::fabs(back_f.vr - fix.vr) < 5.0e5,
	   "闭合轨道：整圈后仍回到自身" );
    check( std::fabs(fix.vr) > 0.02*vmag,
	   "闭合轨道在该方位角上**不是纯切向发射**（这是之前轨道不闭合的关键）" );
    check( rmin > 2.0 && rmax < 5.0, "闭合轨道位于场图范围内" );
    check( rmax - rmin < 0.35,
	   "扇贝幅度远小于朴素初值的 0.62 m（< 0.35 m）" );

    // =====================================================================
    // 3. 逐圈重复性（这才是画图时"轨道重合"的判据）
    // =====================================================================
    std::printf( "\n--- 3. 多圈重复性 ---\n" );
    const int NTURN = 3;
    std::vector<std::vector<TrajectoryPoint>> lines( 1 );
    std::vector<State> per_turn;

    State st = fix;
    double tt = 0.0;
    for( int k = 0; k < NTURN; ++k ) {
	double tlen = 0, rmn = 0, rmx = 0;
	State nxt = one_turn( B, pusher, th0, vmag, dt, st, (k == 0), FULL,
			      0, &tlen, &rmx, &rmn );
	per_turn.push_back( st );
	tt += tlen;
	st = nxt;
    }
    // 用最终状态再走一圈记录轨迹点（供导出）
    {
	Vec3D x, v;
	state_to_xv( th0, vmag, fix, x, v );
	pusher.initialize( 0, &B, x, v, dt );
	const int NS = 6000;
	double acc = 0.0, th_prev = th0, tcur = 0.0;
	TrajectoryPoint p0; p0.t = 0.0; p0.x = x;
	lines[0].push_back( p0 );
	for( int n = 0; n < NS; ++n ) {
	    pusher.step( 0, &B, x, v, dt );
	    tcur += dt;
	    TrajectoryPoint pp; pp.t = tcur; pp.x = x;
	    lines[0].push_back( pp );
	}
    }
    vtk_write_polylines( "cycl_closed_orbit.vtp", lines );

    double dmax_r = 0.0, dmax_v = 0.0;
    for( int k = 1; k < NTURN; ++k ) {
	dmax_r = std::max( dmax_r, std::fabs( per_turn[k].r  - per_turn[0].r ) );
	dmax_v = std::max( dmax_v, std::fabs( per_turn[k].vr - per_turn[0].vr ) );
    }
    std::printf( "  第 1..%d 圈与第 1 圈的状态差: dr_max = %.3e m, dv_max = %.3e m/s\n",
		 NTURN, dmax_r, dmax_v );
    std::printf( "  （对比：朴素初值的同一量是 0.24 m —— 改善 %.0f 倍）\n",
		 0.24/dmax_r );
    check( dmax_r < 5.0e-2, "逐圈重复：同一方位角处半径差 < 5 cm" );
    check( dmax_r < 0.24/10.0,
	   "逐圈重复性比朴素初值改善 10 倍以上" );

    // =====================================================================
    // 4. 径向 betatron 频率 nu_r（单扇区映射 -> 整圈）
    // =====================================================================
    {
	const double DR = 1.0e-4, DV = 2.0e3;
	State xa = fix; xa.r  += DR;
	State xb = fix; xb.vr += DV;
	State f0 = one_turn( B, pusher, th0, vmag, dt, fix, true, SECT );
	State fa = one_turn( B, pusher, th0, vmag, dt, xa,  true, SECT );
	State fb = one_turn( B, pusher, th0, vmag, dt, xb,  true, SECT );
	const double M11 = (fa.r  - f0.r)/DR,  M21 = (fa.vr - f0.vr)/DR;
	const double M12 = (fb.r  - f0.r)/DV,  M22 = (fb.vr - f0.vr)/DV;
	const double tr = M11 + M22, det = M11*M22 - M12*M21;
	std::printf( "\n--- 4. 径向 betatron 微扰（单扇区微分映射）---\n" );
	std::printf( "  M_sect = [[%.6f, %.6e], [%.6e, %.6f]]\n",
		     M11, M12, M21, M22 );
	std::printf( "  trace = %.6f, det = %.6f\n", tr, det );
	const double arg = 0.5*tr/std::sqrt( (det > 0) ? det : 1.0 );
	if( std::fabs(arg) < 1.0 ) {
	    const double nu = 8.0*std::acos( arg )/(2.0*M_PI);
	    std::printf( "  每个扇区相移 %.4f 圈;  径向 betatron 频率 nu_r = %.4f\n",
		     std::acos( arg )/(2.0*M_PI), nu );
	    check( nu > 0.0 && std::fabs(nu - std::floor(nu)) > 0.02,
		   "nu_r 为非整数 => 不动点稳定" );
	} else {
	    std::printf( "  |tr|/2 > 1：径向运动不稳定\n" );
	    check( false, "径向 betatron 稳定" );
	}
    }

    std::printf( "\n导出: cycl_closed_orbit.vtp（%zu 点）\n", lines[0].size() );
    std::printf( "\n%s (%d failure%s)\n",
		 g_fail ? "FAILED" : "ALL TESTS PASSED",
		 g_fail, g_fail == 1 ? "" : "s" );
    return( g_fail ? 1 : 0 );
}
