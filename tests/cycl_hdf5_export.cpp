/** @file cycl_hdf5_export.cpp
 *  @brief HDF5 二进制 IO 后端测试（ROADMAP R3-④）。
 *
 *  验证四件事：
 *    1. **写入 → 回读逐位一致**：回读刻意用 HDF5 的 C API 直接读（不复用写出器
 *       内部的任何封装），做到两条独立通路交叉验证——与 vtk_io 的做法一致。
 *    2. **索引次序**：网格取 (nx,ny,nz)=(24,20,16) 这种三维互不相同的尺寸，
 *       数值编码成 f = i + 1000*j + 1e6*k，任何转置/翻转都会立刻被抓到。
 *       同时验证数据集形状确实是 (nz,ny,nx)（i 最快，与 .vti 一致）。
 *    3. **属性自描述**：openPMD 头部、gridSpacing/gridGlobalOffset/axisLabels、
 *       数据集 unitSI。
 *    4. **大数据的实际收益**：同一份 128³ 光滑场，比较 ASCII `.vti` 与 HDF5
 *       的字节数，并演示 HDF5 的**分块读取**（只读一个平面）。
 *
 *  未编译进 HDF5 后端时输出 [SKIP] 并返回成功（与 cycl_cuda_tracker 一致），
 *  因此没有 HDF5 的机器上 CI 不受影响。
 *
 *  用法: cycl_hdf5_export
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "config.h"
#include "hdf5writer.hpp"
#include "vtkwriter.hpp"

#ifdef IBSIMU_CYCL_HAVE_HDF5
#include <hdf5.h>
#endif

using namespace ibsimu_cycl;

static int g_failures = 0;

static void check( bool ok, const std::string &msg )
{
    std::printf( "%s  %s\n", ok ? "[ ok ]" : "[FAIL]", msg.c_str() );
    if( !ok ) ++g_failures;
}

/*! \brief 文件字节数（不存在返回 0）。 */
static double file_bytes( const std::string &path )
{
    std::ifstream f( path.c_str(), std::ios::binary );
    if( !f )
        return( 0.0 );
    f.seekg( 0, std::ios::end );
    return( (double)f.tellg() );
}


#ifdef IBSIMU_CYCL_HAVE_HDF5

// ===========================================================================
// 独立回读（直接用 HDF5 C API）
// ===========================================================================

/*! \brief 用 RAII 包住一个句柄，避免失败路径泄漏。 */
class Id {
    hid_t _id;
    herr_t (*_closer)(hid_t);
public:
    Id( hid_t id, herr_t (*closer)(hid_t) ) : _id(id), _closer(closer) {}
    ~Id() { if( _id >= 0 ) _closer( _id ); }
    hid_t id() const { return( _id ); }
    bool ok() const { return( _id >= 0 ); }
private:
    Id( const Id & );
    Id &operator=( const Id & );
};

static double attr_double( hid_t obj, const char *name )
{
    Id a( H5Aopen( obj, name, H5P_DEFAULT ), H5Aclose );
    if( !a.ok() )
        throw( std::runtime_error( std::string("no attribute ") + name ) );
    double v = -1.0e300;
    if( H5Aread( a.id(), H5T_NATIVE_DOUBLE, &v ) < 0 )
        throw( std::runtime_error( std::string("cannot read ") + name ) );
    return( v );
}

static std::vector<double> attr_double_list( hid_t obj, const char *name )
{
    Id a( H5Aopen( obj, name, H5P_DEFAULT ), H5Aclose );
    if( !a.ok() )
        throw( std::runtime_error( std::string("no attribute ") + name ) );
    Id sp( H5Aget_space( a.id() ), H5Sclose );
    const int rank = H5Sget_simple_extent_ndims( sp.id() );
    std::size_t len = 1;
    if( rank >= 1 ) {
        std::vector<hsize_t> d( rank );
        H5Sget_simple_extent_dims( sp.id(), &d[0], 0 );
        len = (std::size_t)d[0];
    }
    std::vector<double> v( len );
    if( H5Aread( a.id(), H5T_NATIVE_DOUBLE, &v[0] ) < 0 )
        throw( std::runtime_error( std::string("cannot read ") + name ) );
    return( v );
}

static std::string attr_string( hid_t obj, const char *name )
{
    Id a( H5Aopen( obj, name, H5P_DEFAULT ), H5Aclose );
    if( !a.ok() )
        throw( std::runtime_error( std::string("no attribute ") + name ) );
    Id t( H5Aget_type( a.id() ), H5Tclose );
    const std::size_t sz = H5Tget_size( t.id() );
    std::vector<char> buf( sz + 1, '\0' );
    if( H5Aread( a.id(), t.id(), &buf[0] ) < 0 )
        throw( std::runtime_error( std::string("cannot read ") + name ) );
    return( std::string( &buf[0] ) );
}

/*! \brief 字符串数组属性的第 k 个元素（如 axisLabels）。 */
static std::string attr_string_element( hid_t obj, const char *name,
					std::size_t k )
{
    Id a( H5Aopen( obj, name, H5P_DEFAULT ), H5Aclose );
    if( !a.ok() )
        throw( std::runtime_error( std::string("no attribute ") + name ) );
    Id t( H5Aget_type( a.id() ), H5Tclose );
    Id sp( H5Aget_space( a.id() ), H5Sclose );
    const std::size_t sz = H5Tget_size( t.id() );
    hsize_t n = 0;
    H5Sget_simple_extent_dims( sp.id(), &n, 0 );
    std::vector<char> buf( (std::size_t)n*sz, '\0' );
    if( H5Aread( a.id(), t.id(), &buf[0] ) < 0 )
        throw( std::runtime_error( std::string("cannot read ") + name ) );
    if( k >= (std::size_t)n )
        throw( std::runtime_error( std::string("index out of range: ") + name ) );
    return( std::string( &buf[k*sz] ) );
}

struct Dataset {
    int rank;
    std::vector<hsize_t> dims;
    std::vector<double> data;      // 仅支持 double / int64 转 double
    std::vector<long long> idata;
    double unit_si;
    bool has_unit_si;
};

static Dataset read_double_dataset( hid_t file, const std::string &path )
{
    Dataset ds;
    ds.rank = 0;
    ds.has_unit_si = false;
    Id d( H5Dopen2( file, path.c_str(), H5P_DEFAULT ), H5Dclose );
    if( !d.ok() )
        throw( std::runtime_error( "cannot open dataset " + path ) );
    Id sp( H5Dget_space( d.id() ), H5Sclose );
    ds.rank = H5Sget_simple_extent_ndims( sp.id() );
    ds.dims.resize( ds.rank );
    H5Sget_simple_extent_dims( sp.id(), &ds.dims[0], 0 );
    std::size_t n = 1;
    for( int i = 0; i < ds.rank; ++i )
        n *= (std::size_t)ds.dims[i];
    ds.data.resize( n );
    if( H5Dread( d.id(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		 H5P_DEFAULT, &ds.data[0] ) < 0 )
        throw( std::runtime_error( "cannot read dataset " + path ) );

    Id a( H5Aopen( d.id(), "unitSI", H5P_DEFAULT ), H5Aclose );
    if( a.ok() ) {
        double u = 0.0;
        H5Aread( a.id(), H5T_NATIVE_DOUBLE, &u );
        ds.has_unit_si = true;
        ds.unit_si = u;
    }
    return( ds );
}

static std::vector<long long> read_int64_dataset( hid_t file,
						  const std::string &path )
{
    Id d( H5Dopen2( file, path.c_str(), H5P_DEFAULT ), H5Dclose );
    if( !d.ok() )
        throw( std::runtime_error( "cannot open dataset " + path ) );
    Id sp( H5Dget_space( d.id() ), H5Sclose );
    const int rank = H5Sget_simple_extent_ndims( sp.id() );
    std::vector<hsize_t> dims( rank );
    H5Sget_simple_extent_dims( sp.id(), &dims[0], 0 );
    std::size_t n = 1;
    for( int i = 0; i < rank; ++i )
        n *= (std::size_t)dims[i];
    std::vector<long long> v( n );
    if( H5Dread( d.id(), H5T_NATIVE_INT64, H5S_ALL, H5S_ALL,
		 H5P_DEFAULT, &v[0] ) < 0 )
        throw( std::runtime_error( "cannot read dataset " + path ) );
    return( v );
}

#endif // IBSIMU_CYCL_HAVE_HDF5


int main( int argc, char **argv )
{
    (void)argc; (void)argv;

#ifndef IBSIMU_CYCL_HAVE_HDF5
    std::printf( "[SKIP] HDF5 backend not compiled in "
		 "(configure with the HDF5 library available)\n" );
    std::printf( "\nALL TESTS PASSED (0 failures, 0 checks)\n" );
    return( 0 );
#else
    std::printf( "HDF5 backend: available, library version %s\n",
		 hdf5_version_string().c_str() );
    check( hdf5_available(), "hdf5_available() reports a compiled-in backend" );

    try {
	// -------------------------------------------------------------------
	// 1. 网格：写入 → 独立回读
	// -------------------------------------------------------------------
	const int nx = 24, ny = 20, nz = 16;      // 三维互不相同，便于抓转置
	const Int3D dims( nx, ny, nz );
	const Vec3D origo( -1.5, -2.0, -0.5 );
	const Vec3D spacing( 0.05, 0.05, 0.05 );

	std::vector<double> scalar( (std::size_t)nx*ny*nz );
	std::vector<double> vector( 3*(std::size_t)nx*ny*nz );
	for( int k = 0; k < nz; ++k )
	    for( int j = 0; j < ny; ++j )
		for( int i = 0; i < nx; ++i ) {
		    // 数值里编码 (i,j,k)：三位十进制互不重叠
		    const double f = (double)i + 1000.0*(double)j
			+ 1.0e6*(double)k;
		    const std::size_t a = (std::size_t)i + (std::size_t)nx*
			((std::size_t)j + (std::size_t)ny*(std::size_t)k);
		    scalar[a] = f;
		    vector[3*a+0] =  f;
		    vector[3*a+1] = 2.0*f;
		    vector[3*a+2] = -f;
		}

	const std::string fname = "cycl_hdf5_export.h5";
	// 用一个非平凡的单位因子，确保 unitSI 真的在传递而不是恒为 1
	const double KE_UNIT = 1.602176634e-16;      // keV -> J
	hdf5_write_image_data( fname, "Bfield", dims, origo, spacing,
			       "Bz", scalar, KE_UNIT,
			       "Bvec", vector, 1.0 );
	std::printf( "  wrote %s (%.1f kB)\n", fname.c_str(),
		     file_bytes(fname)/1024.0 );

	Id fid( H5Fopen( fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT ), H5Fclose );
	check( fid.ok(), "HDF5 file can be reopened" );

	// ---- 1a. openPMD 头部属性 ----
	check( attr_string( fid.id(), "openPMD" ) == "2.0.0",
	       "root attribute openPMD == 2.0.0" );
	check( attr_string( fid.id(), "iterationEncoding" ) == "fileBased",
	       "root attribute iterationEncoding" );
	check( attr_string( fid.id(), "basePath" ) == "/data/%T/",
	       "root attribute basePath" );

	// ---- 1b. 网格定义与单位 ----
	Id mesh( H5Gopen2( fid.id(), "/data/0/meshes/Bfield", H5P_DEFAULT ),
		 H5Gclose );
	check( mesh.ok(), "mesh group exists at /data/0/meshes/Bfield" );
	std::vector<double> gs = attr_double_list( mesh.id(), "gridSpacing" );
	check( gs.size() == 3 && gs[0] == spacing[0] && gs[1] == spacing[1]
	       && gs[2] == spacing[2],
	       "gridSpacing attribute round-trips" );
	std::vector<double> go = attr_double_list( mesh.id(), "gridGlobalOffset" );
	check( go.size() == 3 && go[0] == origo[0] && go[1] == origo[1]
	       && go[2] == origo[2],
	       "gridGlobalOffset attribute round-trips" );
	check( attr_string( mesh.id(), "geometry" ) == "cartesian",
	       "geometry == cartesian" );
	check( attr_string_element( mesh.id(), "axisLabels", 0 ) == "x"
	       && attr_string_element( mesh.id(), "axisLabels", 2 ) == "z",
	       "axisLabels == [x,y,z]" );

	// ---- 1c. 标量数据：形状 (nz,ny,nx) + 逐位一致 ----
	Dataset bz = read_double_dataset( fid.id(), "/data/0/meshes/Bfield/Bz" );
	check( bz.rank == 3 && bz.dims[0] == (hsize_t)nz
	       && bz.dims[1] == (hsize_t)ny && bz.dims[2] == (hsize_t)nx,
	       "scalar dataset shape is (nz,ny,nx) - i fastest" );

	bool same = ( bz.data.size() == scalar.size() );
	if( same ) {
	    // 回读是 C 序的 (nz,ny,nx)，与写入缓冲区应逐位相同
	    for( std::size_t a = 0; a < scalar.size() && same; ++a )
		same = ( bz.data[a] == scalar[a] );
	}
	check( same, "scalar data round-trips bit-exactly (no transpose)" );

	// 反证：若被当成 (nx,ny,nz) 解释，数值编码会立刻不合法
	bool transposed_ok = true;
	for( int k = 0; k < nz && transposed_ok; ++k )
	    for( int j = 0; j < ny && transposed_ok; ++j )
		for( int i = 0; i < nx && transposed_ok; ++i ) {
		    const std::size_t a = (std::size_t)k + (std::size_t)nz*
			((std::size_t)j + (std::size_t)ny*(std::size_t)i);
		    if( a < bz.data.size() ) {
			const double f = (double)i + 1000.0*(double)j
			    + 1.0e6*(double)k;
			if( bz.data[a] != f )
			    transposed_ok = false;
		    }
		}
	check( !transposed_ok,
	       "the same bytes do NOT satisfy a transposed (nx,ny,nz) reading" );
	check( bz.has_unit_si && bz.unit_si == KE_UNIT,
	       "scalar dataset carries the exact unitSI it was given" );

	// ---- 1d. 矢量数据：按点交错的 (nz,ny,nx,3) ----
	Dataset bv = read_double_dataset( fid.id(), "/data/0/meshes/Bfield/Bvec" );
	check( bv.rank == 4 && bv.dims[0] == (hsize_t)nz
	       && bv.dims[1] == (hsize_t)ny && bv.dims[2] == (hsize_t)nx
	       && bv.dims[3] == 3,
	       "vector dataset shape is (nz,ny,nx,3) - component last" );
	bool vsame = ( bv.data.size() == vector.size() );
	if( vsame ) {
	    for( std::size_t a = 0; a < vector.size() && vsame; ++a )
		vsame = ( bv.data[a] == vector[a] );
	}
	check( vsame, "vector data round-trips bit-exactly" );

	// 但“扁平字节相同”对**布局**不敏感：必须按声明的形状再验证一次分量结构，
	// 否则「声明 (3,nz,ny,nx) 而实际写的是按点交错数据」这种错会被漏掉
	// （这个 bug 正是 Python 侧交叉验证发现的）。
	{
	    bool comp_ok = ( bv.rank == 4 && bv.dims[3] == 3 );
	    for( int k = 0; k < nz && comp_ok; ++k )
		for( int j = 0; j < ny && comp_ok; ++j )
		    for( int i = 0; i < nx && comp_ok; ++i ) {
			const std::size_t a =
			    (((std::size_t)k*ny + j)*nx + i)*3;
			if( a + 2 >= bv.data.size() ) { comp_ok = false; break; }
			const double f = (double)i + 1000.0*(double)j
			    + 1.0e6*(double)k;
			comp_ok = ( bv.data[a+0] == f ) && ( bv.data[a+1] == 2.0*f )
			    && ( bv.data[a+2] == -f );
		    }
	    check( comp_ok,
		   "vector components are (f,2f,-f) under the declared "
		   "(nz,ny,nx,3) shape" );
	}

	// -------------------------------------------------------------------
	// 2. 轨迹：写入 → 独立回读
	// -------------------------------------------------------------------
	std::vector<std::vector<TrajectoryPoint>> lines( 3 );
	const int npts[3] = { 5, 7, 2 };
	for( int m = 0; m < 3; ++m )
	    for( int i = 0; i < npts[m]; ++i ) {
		TrajectoryPoint p;
		p.t = 1.0e-9*(double)(10*m + i);
		p.x = Vec3D( 0.1*m + 0.01*i, -0.2*m + 0.02*i, 0.3*m - 0.03*i );
		lines[m].push_back( p );
	    }

	const std::string pfname = "cycl_hdf5_export_tracks.h5";
	hdf5_write_polylines( pfname, "protons", lines );

	Id pf( H5Fopen( pfname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT ), H5Fclose );
	const std::string base = "/data/0/particles/protons";
	Dataset pos = read_double_dataset( pf.id(), base + "/position" );
	Dataset tim = read_double_dataset( pf.id(), base + "/time" );
	std::vector<long long> off = read_int64_dataset( pf.id(), base + "/offset" );
	std::vector<long long> cnt = read_int64_dataset( pf.id(), base + "/count" );

	const std::size_t ntot = 5 + 7 + 2;
	check( pos.rank == 2 && pos.dims[0] == (hsize_t)ntot && pos.dims[1] == 3,
	       "particle position dataset is (N,3)" );
	check( tim.data.size() == ntot, "particle time dataset is (N,)" );
	check( off.size() == 4 && off[0] == 0 && off[1] == 5 && off[2] == 12
	       && off[3] == 14,
	       "offset prefix sums are correct" );
	check( cnt.size() == 3 && cnt[0] == 5 && cnt[1] == 7 && cnt[2] == 2,
	       "per-track counts are correct" );

	bool psame = true;
	std::size_t k2 = 0;
	for( int m = 0; m < 3; ++m )
	    for( int i = 0; i < npts[m]; ++i, ++k2 ) {
		const TrajectoryPoint &p = lines[m][i];
		psame = psame && pos.data[3*k2+0] == p.x[0]
		    && pos.data[3*k2+1] == p.x[1]
		    && pos.data[3*k2+2] == p.x[2]
		    && tim.data[k2] == p.t;
	    }
	check( psame, "particle positions and times round-trip bit-exactly" );

	// -------------------------------------------------------------------
	// 3. 大数据：体积对比 + 分块读取
	// -------------------------------------------------------------------
	const int big = 128;
	const Int3D bdims( big, big, big );
	std::vector<double> field( (std::size_t)big*big*big );
	for( int k = 0; k < big; ++k )
	    for( int j = 0; j < big; ++j )
		for( int i = 0; i < big; ++i ) {
		    const std::size_t a = (std::size_t)i + (std::size_t)big*
			((std::size_t)j + (std::size_t)big*(std::size_t)k);
		    field[a] = 1.5*std::sin( 0.03*((double)i + 0.5*(double)j) )
			*std::cos( 0.05*((double)k + 0.3*(double)j) );
		}

	const std::string sfx = "cycl_hdf5_bench";
	const std::clock_t t0 = std::clock();
	vtk_write_image_data( sfx + ".vti", bdims, Vec3D(0,0,0),
			      Vec3D(0.01,0.01,0.01), "B", field );
	const double dt_vtk = (double)(std::clock() - t0)/CLOCKS_PER_SEC;

	const std::clock_t t1 = std::clock();
	hdf5_write_image_data( sfx + ".h5", "Bfield", bdims, Vec3D(0,0,0),
			       Vec3D(0.01,0.01,0.01), "B", field, 1.0 );
	const double dt_h5 = (double)(std::clock() - t1)/CLOCKS_PER_SEC;

	const double nbytes_vtk = file_bytes( sfx + ".vti" );
	const double nbytes_h5  = file_bytes( sfx + ".h5" );
	const double npoint = (double)big*big*big;
	std::printf( "  grid %d^3 = %.2f M points\n", big, npoint/1.0e6 );
	std::printf( "    ASCII .vti : %8.1f MB  %.2f B/point  (%.2f s)\n",
		     nbytes_vtk/1048576.0, nbytes_vtk/npoint, dt_vtk );
	std::printf( "    HDF5       : %8.1f MB  %.2f B/point  (%.2f s)\n",
		     nbytes_h5/1048576.0, nbytes_h5/npoint, dt_h5 );
	std::printf( "    ratio      : %.1fx smaller\n", nbytes_vtk/nbytes_h5 );
	// 注意：double 光滑数据的无损压缩比远不如整数/单精度数据（尾数低位近似随机），
	// 所以这里的阈值取得保守；HDF5 通路真正的收益在分块随机读取与自描述元数据。
	check( nbytes_h5 > 0.0 && nbytes_vtk > 2.0*nbytes_h5,
	       "HDF5 output is more than 2x smaller than ASCII VTK" );

	// 分块读取：只取 z=64 的一个平面（ny*nx 个值）
	{
	    Id bf( H5Fopen( (sfx + ".h5").c_str(), H5F_ACC_RDONLY, H5P_DEFAULT ),
		   H5Fclose );
	    Id d( H5Dopen2( bf.id(), "/data/0/meshes/Bfield/B", H5P_DEFAULT ),
		  H5Dclose );
	    const hsize_t start[3] = { 64, 0, 0 };
	    const hsize_t count[3] = { 1, (hsize_t)big, (hsize_t)big };
	    Id mspace( H5Screate_simple( 3, count, 0 ), H5Sclose );
	    Id fspace( H5Dget_space( d.id() ), H5Sclose );
	    if( H5Sselect_hyperslab( fspace.id(), H5S_SELECT_SET,
				     start, 0, count, 0 ) >= 0 ) {
		std::vector<double> plane( (std::size_t)big*big );
		const std::clock_t t2 = std::clock();
		const herr_t r = H5Dread( d.id(), H5T_NATIVE_DOUBLE, mspace.id(),
					  fspace.id(), H5P_DEFAULT, &plane[0] );
		const double dt_plane =
		    (double)(std::clock() - t2)/CLOCKS_PER_SEC;
		bool plane_ok = ( r >= 0 );
		for( int j = 0; j < big && plane_ok; ++j )
		    for( int i = 0; i < big && plane_ok; ++i ) {
			const std::size_t a = (std::size_t)i
			    + (std::size_t)big*(std::size_t)j;
			const std::size_t s = (std::size_t)i + (std::size_t)big*
			    ((std::size_t)j + (std::size_t)big*64);
			plane_ok = ( plane[a] == field[s] );
		    }
		std::printf( "    hyperslab read of one plane: %.0f kB in %.4f s\n",
			     (double)big*big*8.0/1024.0, dt_plane );
		check( plane_ok,
		       "single-plane hyperslab read returns the right values" );
	    } else {
		check( false, "hyperslab selection failed" );
	    }
	}

	// 大文件不留在测试目录里
	std::remove( (sfx + ".vti").c_str() );
	std::remove( (sfx + ".h5").c_str() );

    } catch( const std::exception &e ) {
	std::printf( "[FAIL] exception: %s\n", e.what() );
	++g_failures;
    }

    std::printf( "\n%s (%d failure%s)\n",
		 g_failures ? "FAILED" : "ALL TESTS PASSED",
		 g_failures, g_failures == 1 ? "" : "s" );
    return( g_failures ? 1 : 0 );
#endif // IBSIMU_CYCL_HAVE_HDF5
}
