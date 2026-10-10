/** @file cycl_vtk_export.cpp
 *  @brief  R3-①：VTK XML 导出器验证（场图 `.vti` + 轨迹 `.vtp`）。
 *
 *  @test  三部分：
 *            A. 规则网格场图导出 → 逐值解析回读校验（属性 + 标量 + 矢量）；
 *            B. 解析折线导出 → 点数/折线数/连接表/坐标校验（含多折线）；
 *            C. **端到端**：用 `CFieldMap3D` 均匀场 + `CBorisPusher` 积分一条
 *               螺旋轨道，导出为 `.vtp`，并与解析解（r、ω、z 线性推进）对比。
 *
 *  产出文件可直接用 ParaView / PyVista 打开：
 *    paraview cycl_vtk_field.vti cycl_vtk_orbit.vtp
 *    python3 examples/cyclotron/view_3d.py cycl_vtk_field.vti cycl_vtk_orbit.vtp
 *
 *  说明：导出器按 ASCII 内联写法（无外部依赖），本测试用**独立实现**的极简
 *  VTK XML 解析器回读，避免"自己验自己"。
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "constants.hpp"
#include "boris.hpp"
#include "fieldmap3d.hpp"
#include "vec3d.hpp"
#include "vtkwriter.hpp"


using namespace std;
using namespace ibsimu_cycl;


namespace {

const double PI = 3.14159265358979323846;

int err_count = 0;


void check( bool ok, const string &what )
{
    printf( "  [%s] %s\n", ok ? "OK  " : "FAIL", what.c_str() );
    if( !ok )
	err_count++;
}


void check_close( double a, double b, double tol, const string &what )
{
    double scale = std::max( fabs(a), fabs(b) );
    bool ok = fabs(a-b) <= tol*scale + 1e-300;
    if( !ok )
	printf( "       got %.12g, expected %.12g (tol %.2g rel)\n", a, b, tol );
    check( ok, what );
}


/*! \brief 绝对误差判据（用于期望值接近 0 的场合）。 */
void check_abs( double v, double tol, const string &what )
{
    bool ok = fabs(v) <= tol;
    if( !ok )
	printf( "       got %.12g, |value| should be <= %.2g\n", v, tol );
    check( ok, what );
}


/* ---------------- 极简 VTK XML 读取（独立于写出器） ---------------- */

string slurp( const string &fn )
{
    ifstream is( fn.c_str() );
    if( !is )
	throw( runtime_error( "cannot open " + fn ) );
    ostringstream ss;
    ss << is.rdbuf();
    return( ss.str() );
}


/*! \brief 取属性值：attr="value"。 */
bool get_attr( const string &s, const string &attr, string &val )
{
    const string key = attr + "=\"";
    size_t p = s.find( key );
    if( p == string::npos )
	return( false );
    p += key.size();
    size_t q = s.find( '"', p );
    if( q == string::npos )
	return( false );
    val = s.substr( p, q-p );
    return( true );
}


/*! \brief 取名为 name 的 DataArray 的数值体（按空白切分）。 */
bool get_array( const string &s, const string &name, vector<double> &vals,
		int *ncomp_out = 0 )
{
    const string key = "Name=\"" + name + "\"";
    size_t p = s.find( key );
    if( p == string::npos )
	return( false );

    size_t tag_end = s.find( '>', p );
    if( tag_end == string::npos )
	return( false );

    if( ncomp_out ) {
	string ncs;
	*ncomp_out = get_attr( s.substr( p, tag_end-p ), "NumberOfComponents", ncs )
	    ? atoi( ncs.c_str() ) : 1;
    }

    size_t body_end = s.find( "</DataArray>", tag_end );
    if( body_end == string::npos )
	return( false );

    istringstream is( s.substr( tag_end+1, body_end-tag_end-1 ) );
    vals.clear();
    double v;
    while( is >> v )
	vals.push_back( v );
    return( true );
}


/*! \brief 解析 "0 8 0 6 0 4" → {0,8,0,6,0,4}。 */
vector<long> parse_ints( const string &s )
{
    istringstream is( s );
    vector<long> v;
    long x;
    while( is >> x )
	v.push_back( x );
    return( v );
}


/*! \brief 解析 "a b c" → {a,b,c}。 */
vector<double> parse_reals( const string &s )
{
    istringstream is( s );
    vector<double> v;
    double x;
    while( is >> x )
	v.push_back( x );
    return( v );
}


/* ---------------- 数值比较若干数组 ---------------- */

bool arrays_close( const vector<double> &a, const vector<double> &b, double tol )
{
    if( a.size() != b.size() )
	return( false );
    for( size_t i = 0; i < a.size(); i++ ) {
	double scale = std::max( fabs(a[i]), fabs(b[i]) );
	if( fabs(a[i]-b[i]) > tol*scale + 1e-300 )
	    return( false );
    }
    return( true );
}

} // namespace


void test( int argc, char **argv )
{
    (void)argc;
    (void)argv;

    printf( "\n=== R3-1 VTK XML 导出器验证 ===\n" );

    /* ================= A. 规则网格场图 (.vti) ================= */
    printf( "\n[A] ImageData 场图导出 / 回读\n" );

    const int    NX = 9, NY = 7, NZ = 5;
    const double X0 = -0.04, Y0 = -0.03, Z0 = -0.02, H = 0.01;
    const double B0 = 0.5, KX = 30.0, KY = 20.0;

    vector<double> bmag( (size_t)NX*NY*NZ, 0.0 );
    vector<double> bvec( 3*(size_t)NX*NY*NZ, 0.0 );

    for( int k = 0; k < NZ; k++ ) {
	for( int j = 0; j < NY; j++ ) {
	    for( int i = 0; i < NX; i++ ) {
		const size_t a = (size_t)i + NX*((size_t)j + (size_t)NY*k);
		const double x = X0 + i*H;
		const double y = Y0 + j*H;
		const double bz = B0*cos( KX*x )*cos( KY*y );
		bmag[a]     = fabs( bz );
		bvec[3*a+0] = 0.0;
		bvec[3*a+1] = 0.0;
		bvec[3*a+2] = bz;
	    }
	}
    }

    vtk_write_image_data( "cycl_vtk_field.vti",
			  Int3D( NX, NY, NZ ), Vec3D( X0, Y0, Z0 ), Vec3D( H, H, H ),
			  "Bmag", bmag, "B", bvec );

    const string vti = slurp( "cycl_vtk_field.vti" );

    string val;
    check( get_attr( vti, "type", val ) && val == "ImageData", "文件类型 = ImageData" );
    check( get_attr( vti, "WholeExtent", val ) &&
	   parse_ints( val ) == vector<long>( { 0, NX-1, 0, NY-1, 0, NZ-1 } ),
	   "WholeExtent = 节点范围" );
    check( get_attr( vti, "Origin", val ) &&
	   arrays_close( parse_reals( val ), vector<double>{ X0, Y0, Z0 }, 1e-11 ),
	   "Origin 回读一致" );
    check( get_attr( vti, "Spacing", val ) &&
	   arrays_close( parse_reals( val ), vector<double>{ H, H, H }, 1e-11 ),
	   "Spacing 回读一致" );

    {
	vector<double> read_mag, read_vec;
	int nc_mag = 0, nc_vec = 0;
	bool ok = get_array( vti, "Bmag", read_mag, &nc_mag );
	ok = ok && get_array( vti, "B", read_vec, &nc_vec );
	check( ok, "标量与矢量数组均存在" );
	check( ok && read_mag.size() == bmag.size() && nc_mag == 1,
	       "标量数组长度/分量数正确" );
	check( ok && read_vec.size() == bvec.size() && nc_vec == 3,
	       "矢量数组长度/分量数正确" );
	check( ok && arrays_close( read_mag, bmag, 1e-11 ), "标量逐值一致（1e-11）" );
	check( ok && arrays_close( read_vec, bvec, 1e-11 ), "矢量逐值一致（1e-11）" );
    }

    /* ================= B. 解析折线 (.vtp) ================= */
    printf( "\n[B] PolyData 折线导出 / 回读\n" );

    const double R_HELIX = 0.0622, W_HELIX = 4.8243e7, VZ_HELIX = 0.3e6;
    const double T_HELIX = 2.0*PI/W_HELIX;

    vector<vector<TrajectoryPoint>> lines( 2 );
    const int NP = 200;
    for( int a = 0; a <= NP; a++ ) {
	const double t = a*T_HELIX/NP;
	TrajectoryPoint p;
	p.t = t;
	p.x = Vec3D( R_HELIX*cos(W_HELIX*t), -R_HELIX*sin(W_HELIX*t), VZ_HELIX*t );
	lines[0].push_back( p );
    }
    // 第二条：短直线，检验多折线连接表
    for( int a = 0; a < 3; a++ ) {
	TrajectoryPoint p;
	p.t = a*1e-9;
	p.x = Vec3D( 0.2 + 0.001*a, -0.2, 0.0 );
	lines[1].push_back( p );
    }

    vtk_write_polylines( "cycl_vtk_orbit.vtp", lines );

    const string vtp = slurp( "cycl_vtk_orbit.vtp" );
    check( get_attr( vtp, "type", val ) && val == "PolyData", "文件类型 = PolyData" );

    string s_npts, s_nlines;
    const int N_TOT = (NP+1) + 3;              // 长折线 201 点 + 短折线 3 点
    get_attr( vtp, "NumberOfPoints", s_npts );
    get_attr( vtp, "NumberOfLines", s_nlines );
    check( atoi(s_npts.c_str()) == N_TOT, "NumberOfPoints = 204" );
    check( atoi(s_nlines.c_str()) == 2, "NumberOfLines = 2" );

    {
	vector<double> pts, tim, conn, offs;
	int nc_pts = 0, nc_tim = 0;
	bool ok = get_array( vtp, "Points", pts, &nc_pts );
	ok = ok && get_array( vtp, "t", tim, &nc_tim );
	ok = ok && get_array( vtp, "connectivity", conn );
	ok = ok && get_array( vtp, "offsets", offs );
	check( ok, "Points / t / connectivity / offsets 均存在" );
	check( ok && nc_pts == 3 && pts.size() == 3*(size_t)N_TOT, "Points 为 3 分量 × 204 点" );
	check( ok && nc_tim == 1 && tim.size() == (size_t)N_TOT, "时刻数组长度 204" );
	check( ok && conn.size() == (size_t)N_TOT && offs.size() == 2,
	       "连接表 204 项、偏移表 2 项" );
	check( ok && offs.size() == 2 && (long)offs[0] == NP+1 && (long)offs[1] == N_TOT,
	       "offsets = {201, 204}" );

	// 抽查首/末点坐标
	bool geo_ok = ok;
	if( ok ) {
	    double tol = 1e-11;
	    geo_ok = fabs( pts[0] - R_HELIX ) <= tol*R_HELIX
		&& fabs( pts[1] - 0.0 )      <= tol*R_HELIX
		&& fabs( pts[2] - 0.0 )      <= tol*R_HELIX;
	    const size_t last = 3*200;
	    double t_end = T_HELIX;   // 整两圈后回到起点
	    geo_ok = geo_ok
		&& fabs( pts[last+0] - R_HELIX*cos(W_HELIX*t_end) ) <= 1e-9
		&& fabs( pts[last+2] - VZ_HELIX*t_end )             <= 1e-11*VZ_HELIX*t_end;
	}
	check( geo_ok, "首/末点坐标与解析解一致" );
    }

    /* ================= C. 端到端：真实积分 + 导出 ================= */
    printf( "\n[C] CBorisPusher 积分螺旋轨道 → .vtp（端到端）\n" );

    const double QE = CHARGE_E;
    const double MP = MASS_U;
    const double BZ = 0.5;                       // T
    const double VPERP = 3.0e6;                  // m/s
    const double VPARA = 0.3e6;                  // m/s
    const double R_EXP = MP*VPERP/(QE*BZ);       // 解析回旋半径
    const double W_EXP = QE*BZ/MP;               // 解析角频率
    const double T_REV = 2.0*PI/W_EXP;
    const int    NREV = 2;
    const int    NSTEP = 800;
    const double DT = NREV*T_REV/NSTEP;

    // 覆盖 [-0.10, 0.10] 的均匀 Bz 场图
    CFieldMap3D Bf( 21, 21, 21, -0.10, 0.01, -0.10, 0.01, -0.10, 0.01 );
    for( size_t i = 0; i < 21; i++ )
	for( size_t j = 0; j < 21; j++ )
	    for( size_t k = 0; k < 21; k++ )
		Bf.set_value( i, j, k, 0.0, 0.0, BZ );

    CBorisPusher pusher( QE, MP );               // 非相对论（β≈0.01）

    /*! \brief 半径统计：最小/最大/均值/振幅。 */
    struct RadiusStat { double rmin, rmax, rmean, amp; };

    auto radius_stat = [&]( const vector<TrajectoryPoint> &o ) {
	RadiusStat s { 1e30, -1e30, 0.0, 0.0 };
	for( size_t n = 0; n < o.size(); n++ ) {
	    double rr = sqrt( o[n].x[0]*o[n].x[0] + o[n].x[1]*o[n].x[1] );
	    s.rmin = std::min( s.rmin, rr );
	    s.rmax = std::max( s.rmax, rr );
	    s.rmean += rr;
	}
	s.rmean /= o.size();
	s.amp = 0.5*(s.rmax - s.rmin);
	return( s );
    };

    // 起步位置/速度：引导中心在原点，方向为 +z 的螺旋。
    // init=true 时先做半步初始化（leapfrog 启动）。
    auto trace_orbit = [&]( bool init ) {
	vector<TrajectoryPoint> o;
	Vec3D x( R_EXP, 0.0, 0.0 );
	Vec3D v( 0.0, -VPERP, VPARA );
	if( init )
	    pusher.initialize( nullptr, &Bf, x, v, DT );
	for( int n = 0; n <= NSTEP; n++ ) {
	    TrajectoryPoint p;
	    p.t = n*DT;
	    p.x = x;
	    o.push_back( p );
	    pusher.step( nullptr, &Bf, x, v, DT );
	}
	return( o );
    };

    const vector<TrajectoryPoint> orb_plain = trace_orbit( false );
    const vector<TrajectoryPoint> orb_init  = trace_orbit( true );

    // 导出**已初始化**的轨道（更接近真实回旋运动）供 ParaView/PyVista 查看
    vector<vector<TrajectoryPoint>> lines2( 1 );
    lines2[0] = orb_init;
    vtk_write_polylines( "cycl_vtk_orbit.vtp", lines2 );

    const string vtp2 = slurp( "cycl_vtk_orbit.vtp" );
    get_attr( vtp2, "NumberOfPoints", s_npts );
    check( atoi(s_npts.c_str()) == NSTEP+1, "端到端：点数 = 步数+1" );

    const double t_end    = NSTEP*DT;
    const double omega_dt = W_EXP*DT;
    const RadiusStat sp = radius_stat( orb_plain );
    const RadiusStat si = radius_stat( orb_init );

    printf( "       omega*dt = %.6e；半径振幅：未初始化 %.4e，已初始化 %.4e（理论 %.4e）\n",
	    omega_dt, sp.amp, si.amp, 0.5*omega_dt*R_EXP );

    // 未初始化：半径振荡幅度的理论值就是 (omega*dt/2)*r
    check_close( sp.amp, 0.5*omega_dt*R_EXP, 0.05,
		 "未初始化：半径振幅 = (omega*dt/2)*r（理论值）" );
    // 已初始化：瞬态被显著抑制
    check( si.amp < 0.1*sp.amp, "已初始化：半径振幅被抑制 10x 以上" );
    check_close( si.rmean, R_EXP, 1e-4, "已初始化：平均回旋半径 = 解析值" );
    check_close( sp.rmean, R_EXP, 1e-4, "未初始化：平均回旋半径仍 = 解析值" );

    // 终点解析解对比（误差按 R_EXP 归一，因为整圈处分量接近零）
    const TrajectoryPoint &pend = orb_init.back();
    check_abs( (pend.x[0] - R_EXP*cos( W_EXP*t_end ))/R_EXP, 1e-3,
	       "终点 x 与解析解一致（归一 R_EXP）" );
    check_abs( (pend.x[1] + R_EXP*sin( W_EXP*t_end ))/R_EXP, 1e-3,
	       "终点 y 与解析解一致（归一 R_EXP）" );
    check_abs( (pend.x[2] - VPARA*t_end)/(VPARA*t_end), 1e-9,
	       "终点 z = v_parallel·t（精确）" );

    printf( "\n解析量：r = %.6e m, omega = %.6e rad/s, T = %.6e s, 追踪 %d 圈\n",
	    R_EXP, W_EXP, T_REV, NREV );
    printf( "\n产出：cycl_vtk_field.vti（%d 节点）, cycl_vtk_orbit.vtp（%d 点折线）\n",
	    NX*NY*NZ, NSTEP+1 );

    if( err_count ) {
	printf( "\n结果：*** %d 项检查未通过 ***\n", err_count );
	exit( 1 );
    }
    printf( "\n结果：全部检查通过\n" );
    exit( 0 );
}


int main( int argc, char **argv )
{
    try {
	test( argc, argv );
    } catch( std::exception &e ) {
	printf( "异常：%s\n", e.what() );
	return( 1 );
    }
    return( 0 );
}
