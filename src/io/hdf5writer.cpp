/** @file hdf5writer.cpp
 *  @brief HDF5 二进制写出器实现（可选后端，见 hdf5writer.hpp 的布局说明）。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include "hdf5writer.hpp"

#include <algorithm>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "config.h"

#ifdef IBSIMU_CYCL_HAVE_HDF5
#include <hdf5.h>
#endif


namespace ibsimu_cycl {


#ifdef IBSIMU_CYCL_HAVE_HDF5


namespace {


/*! \brief HDF5 对象句柄的 RAII 包装（id < 0 表示空）。
 *
 *  HDF5 的 C API 全靠手工 H5*close，中途抛异常就会泄漏句柄；包装后
 *  所有失败路径都由析构函数兜住。
 *
 *  可移动但不可拷贝（句柄只能有一个所有者）。注意：一旦显式声明了拷贝
 *  构造函数，隐式移动构造函数就不会生成——必须自己写出移动语义，否则
 *  `return(H5Id(...))` 这种按值返回会编译失败。
 */
class H5Id {
    hid_t _id;
    herr_t (*_closer)(hid_t);
public:
    H5Id() : _id(-1), _closer(0) {}
    H5Id( hid_t id, herr_t (*closer)(hid_t) ) : _id(id), _closer(closer) {
	if( id < 0 )
	    throw( std::runtime_error( "HDF5: could not create object" ) );
    }
    H5Id( H5Id &&o ) : _id(o._id), _closer(o._closer) {
	o._id = -1;
	o._closer = 0;
    }
    H5Id &operator=( H5Id &&o ) {
	if( this != &o ) {
	    if( _id >= 0 && _closer ) _closer( _id );
	    _id = o._id; _closer = o._closer;
	    o._id = -1; o._closer = 0;
	}
	return( *this );
    }
    ~H5Id() { if( _id >= 0 && _closer ) _closer( _id ); }
    hid_t id() const { return( _id ); }
private:
    H5Id( const H5Id & );
    H5Id &operator=( const H5Id & );
};


/*! \brief 允许自动创建中间组的链接属性表。 */
H5Id link_plist()
{
    H5Id lcpl( H5Pcreate( H5P_LINK_CREATE ), H5Pclose );
    if( H5Pset_create_intermediate_group( lcpl.id(), 1 ) < 0 )
	throw( std::runtime_error( "HDF5: cannot enable intermediate groups" ) );
    return( lcpl );
}


void attr_string( hid_t obj, const char *name, const std::string &value )
{
    H5Id type( H5Tcopy( H5T_C_S1 ), H5Tclose );
    if( H5Tset_size( type.id(), value.size() + 1 ) < 0 ||
	H5Tset_strpad( type.id(), H5T_STR_NULLTERM ) < 0 )
	throw( std::runtime_error( "HDF5: cannot build string type" ) );
    H5Id space( H5Screate( H5S_SCALAR ), H5Sclose );
    H5Id attr( H5Acreate2( obj, name, type.id(), space.id(),
			   H5P_DEFAULT, H5P_DEFAULT ), H5Aclose );
    if( H5Awrite( attr.id(), type.id(), value.c_str() ) < 0 )
	throw( std::runtime_error( std::string("HDF5: cannot write attribute ")
				   + name ) );
}


/*! \brief 定长字符串数组属性（如 openPMD 的 `axisLabels = ["x","y","z"]`）。 */
void attr_string_array( hid_t obj, const char *name,
			const std::vector<std::string> &values )
{
    if( values.empty() )
	return;
    std::size_t width = 1;
    for( std::size_t i = 0; i < values.size(); ++i )
	width = std::max( width, values[i].size() + 1 );

    H5Id type( H5Tcopy( H5T_C_S1 ), H5Tclose );
    if( H5Tset_size( type.id(), width ) < 0 ||
	H5Tset_strpad( type.id(), H5T_STR_NULLTERM ) < 0 )
	throw( std::runtime_error( "HDF5: cannot build string array type" ) );

    const hsize_t n = (hsize_t)values.size();
    H5Id space( H5Screate_simple( 1, &n, 0 ), H5Sclose );

    std::vector<char> buf( (std::size_t)n*width, '\0' );
    for( std::size_t i = 0; i < values.size(); ++i )
	std::memcpy( &buf[i*width], values[i].c_str(), values[i].size() );

    H5Id attr( H5Acreate2( obj, name, type.id(), space.id(),
			   H5P_DEFAULT, H5P_DEFAULT ), H5Aclose );
    if( H5Awrite( attr.id(), type.id(), &buf[0] ) < 0 )
	throw( std::runtime_error( std::string("HDF5: cannot write attribute ")
				   + name ) );
}


void attr_double( hid_t obj, const char *name, double value )
{
    H5Id space( H5Screate( H5S_SCALAR ), H5Sclose );
    H5Id attr( H5Acreate2( obj, name, H5T_NATIVE_DOUBLE, space.id(),
			   H5P_DEFAULT, H5P_DEFAULT ), H5Aclose );
    if( H5Awrite( attr.id(), H5T_NATIVE_DOUBLE, &value ) < 0 )
	throw( std::runtime_error( std::string("HDF5: cannot write attribute ")
				   + name ) );
}


void attr_double3( hid_t obj, const char *name, const Vec3D &v )
{
    const hsize_t n = 3;
    H5Id space( H5Screate_simple( 1, &n, 0 ), H5Sclose );
    H5Id attr( H5Acreate2( obj, name, H5T_NATIVE_DOUBLE, space.id(),
			   H5P_DEFAULT, H5P_DEFAULT ), H5Aclose );
    const double buf[3] = { v[0], v[1], v[2] };
    if( H5Awrite( attr.id(), H5T_NATIVE_DOUBLE, buf ) < 0 )
	throw( std::runtime_error( std::string("HDF5: cannot write attribute ")
				   + name ) );
}


/*! \brief 写一个分块 + deflate 压缩的数据集。
 *
 *  \param loc    父组（或文件，路径里可以带中间组）。
 *  \param path   数据集路径。
 *  \param type   H5T_NATIVE_* 之一。
 *  \param rank   维数。
 *  \param dims   各维长度（C 序）。
 *  \param chunk  分块尺寸（与 \a dims 同长度）。
 *  \param data   数据缓冲区。
 */
void write_dataset( hid_t loc, const H5Id &lcpl, const char *path, hid_t type,
		    int rank, const hsize_t *dims, const hsize_t *chunk,
		    const void *data )
{
    H5Id space( H5Screate_simple( rank, dims, 0 ), H5Sclose );
    H5Id dcpl( H5Pcreate( H5P_DATASET_CREATE ), H5Pclose );
    if( H5Pset_chunk( dcpl.id(), rank, chunk ) < 0 )
	throw( std::runtime_error( "HDF5: cannot set chunking" ) );
    // shuffle 把相同量级的字节聚到一起（double 的指数/尾数字节分开），
    // 对浮点数据的 deflate 效果提升明显；先 shuffle 再 deflate。
    if( H5Pset_shuffle( dcpl.id() ) < 0 )
	throw( std::runtime_error( "HDF5: cannot enable shuffle" ) );
    // deflate 级别 6：对光滑场图通常有数倍压缩，代价很小
    if( H5Pset_deflate( dcpl.id(), 6 ) < 0 )
	throw( std::runtime_error( "HDF5: cannot enable deflate" ) );

    H5Id dset( H5Dcreate2( loc, path, type, space.id(), lcpl.id(),
			   dcpl.id(), H5P_DEFAULT ), H5Dclose );
    if( H5Dwrite( dset.id(), type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data ) < 0 )
	throw( std::runtime_error( std::string("HDF5: cannot write dataset ")
				   + path ) );
}


void set_unit( hid_t loc, const char *path, double unit_si )
{
    hid_t d = H5Dopen2( loc, path, H5P_DEFAULT );
    if( d < 0 )
	throw( std::runtime_error( std::string("HDF5: cannot reopen ") + path ) );
    try {
	attr_double( d, "unitSI", unit_si );
    } catch( ... ) {
	H5Dclose( d );
	throw;
    }
    H5Dclose( d );
}


void check_size( const Int3D &dims, std::size_t got, std::size_t per_point,
		 const char *what )
{
    const std::size_t want = (std::size_t)dims[0]*dims[1]*dims[2]*per_point;
    if( got != want ) {
	std::ostringstream os;
	os << "hdf5_write_image_data: " << what << " size mismatch ("
	   << got << " != " << want << " = " << dims[0] << "x" << dims[1]
	   << "x" << dims[2] << "x" << per_point << ")";
	throw( std::invalid_argument( os.str() ) );
    }
}


} // namespace


bool hdf5_available()
{
    return( true );
}


std::string hdf5_version_string()
{
    unsigned maj = 0, min = 0, rel = 0;
    if( H5get_libversion( &maj, &min, &rel ) < 0 )
	return( std::string() );
    std::ostringstream os;
    os << maj << "." << min << "." << rel;
    return( os.str() );
}


namespace {

/*! \brief 建好 / 及其 openPMD 头部属性、/data/0 迭代组，返回文件句柄。 */
H5Id make_file_with_openpmd_header( const std::string &filename )
{
    H5Id fid( H5Fcreate( filename.c_str(), H5F_ACC_TRUNC,
			 H5P_DEFAULT, H5P_DEFAULT ), H5Fclose );

    // 沿用 openPMD 2.0 的属性名与语义（见头文件里关于“不宣称通过校验”的说明）
    attr_string( fid.id(), "openPMD", "2.0.0" );
    attr_double( fid.id(), "openPMDextension", 0.0 );
    attr_string( fid.id(), "basePath", "/data/%T/" );
    attr_string( fid.id(), "iterationEncoding", "fileBased" );
    attr_string( fid.id(), "iterationFormat", "/data/%T/" );

    H5Id lcpl = link_plist();
    H5Id it( H5Gcreate2( fid.id(), "/data/0", lcpl.id(),
			 H5P_DEFAULT, H5P_DEFAULT ), H5Gclose );
    attr_double( it.id(), "time", 0.0 );
    attr_double( it.id(), "dt", 0.0 );
    attr_double( it.id(), "timeUnitSI", 1.0 );
    return( fid );
}

} // namespace


void hdf5_write_image_data( const std::string &filename,
			    const std::string &mesh_name,
			    const Int3D &dims, const Vec3D &origo,
			    const Vec3D &spacing,
			    const std::string &scalar_name,
			    const std::vector<double> &scalar,
			    double scalar_unit_si )
{
    check_size( dims, scalar.size(), 1, "scalar" );

    H5Id fid = make_file_with_openpmd_header( filename );
    H5Id lcpl = link_plist();

    const std::string mesh_path = "/data/0/meshes/" + mesh_name;
    {
	H5Id mesh( H5Gcreate2( fid.id(), mesh_path.c_str(), lcpl.id(),
			       H5P_DEFAULT, H5P_DEFAULT ), H5Gclose );
	attr_string( mesh.id(), "geometry", "cartesian" );
	attr_string( mesh.id(), "dataOrder", "C" );
	attr_string_array( mesh.id(), "axisLabels",
			   std::vector<std::string>{ "x", "y", "z" } );
	attr_double3( mesh.id(), "gridSpacing", spacing );
	attr_double3( mesh.id(), "gridGlobalOffset", origo );
	attr_double( mesh.id(), "gridUnitSI", 1.0 );
	attr_double( mesh.id(), "timeOffset", 0.0 );
	attr_double( mesh.id(), "timeOffsetUnitSI", 1.0 );
    }

    // 形状反转成 (nz,ny,nx)：C 序连续内存＝i 最快，与 VTK/网格一致
    const hsize_t d3[3] = { (hsize_t)dims[2], (hsize_t)dims[1],
			    (hsize_t)dims[0] };
    const hsize_t c3[3] = { 1, (hsize_t)dims[1], (hsize_t)dims[0] };
    const std::string spath = mesh_path + "/" + scalar_name;
    write_dataset( fid.id(), lcpl, spath.c_str(), H5T_NATIVE_DOUBLE, 3,
		   d3, c3, &scalar[0] );
    set_unit( fid.id(), spath.c_str(), scalar_unit_si );
}


void hdf5_write_image_data( const std::string &filename,
			    const std::string &mesh_name,
			    const Int3D &dims, const Vec3D &origo,
			    const Vec3D &spacing,
			    const std::string &scalar_name,
			    const std::vector<double> &scalar,
			    double scalar_unit_si,
			    const std::string &vector_name,
			    const std::vector<double> &vector,
			    double vector_unit_si )
{
    check_size( dims, scalar.size(), 1, "scalar" );
    check_size( dims, vector.size(), 3, "vector" );

    hdf5_write_image_data( filename, mesh_name, dims, origo, spacing,
			   scalar_name, scalar, scalar_unit_si );

    H5Id fid( H5Fopen( filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT ), H5Fclose );
    H5Id lcpl = link_plist();

    // 矢量**按点交错**：形状 (nz,ny,nx,3)，最后一维是分量。
    // 调用者给的缓冲区就是 (vx,vy,vz),(vx,vy,vz),… 这样的交错布局，
    // 与 VTK/ParaView 以及 Python 侧 (npts,3) 的习惯一致。
    // （早期版本声明成 (3,nz,ny,nx) 但写入的仍是交错缓冲——声明与
    //   数据不符；C++ 回读只比扁平字节所以没露馅，是 Python 侧交叉验证发现的。）
    const hsize_t d4[4] = { (hsize_t)dims[2], (hsize_t)dims[1],
			    (hsize_t)dims[0], 3 };
    const hsize_t c4[4] = { 1, (hsize_t)dims[1], (hsize_t)dims[0], 3 };
    const std::string vpath = "/data/0/meshes/" + mesh_name + "/" + vector_name;
    write_dataset( fid.id(), lcpl, vpath.c_str(), H5T_NATIVE_DOUBLE, 4,
		   d4, c4, &vector[0] );
    set_unit( fid.id(), vpath.c_str(), vector_unit_si );
}


void hdf5_write_polylines( const std::string &filename,
			   const std::string &group_name,
			   const std::vector<std::vector<TrajectoryPoint>> &lines )
{
    std::size_t npts = 0;
    for( std::size_t m = 0; m < lines.size(); ++m )
	npts += lines[m].size();

    H5Id fid = make_file_with_openpmd_header( filename );
    H5Id lcpl = link_plist();

    const std::string grp_path = "/data/0/particles/" + group_name;
    {
	H5Id grp( H5Gcreate2( fid.id(), grp_path.c_str(), lcpl.id(),
			      H5P_DEFAULT, H5P_DEFAULT ), H5Gclose );
	attr_double( grp.id(), "positionsUnitSI", 1.0 );   // m
	attr_double( grp.id(), "timeUnitSI", 1.0 );        // s
	attr_double( grp.id(), "nTracks", (double)lines.size() );
    }

    // 扁平化：position(N,3)、time(N)、offset(M+1)、count(M)
    std::vector<double>     pos( npts*3, 0.0 );
    std::vector<double>     tim( npts, 0.0 );
    std::vector<long long>  off( lines.size() + 1, 0 );
    std::vector<long long>  cnt( lines.size(), 0 );

    std::size_t k = 0;
    for( std::size_t m = 0; m < lines.size(); ++m ) {
	off[m] = (long long)k;
	cnt[m] = (long long)lines[m].size();
	for( std::size_t i = 0; i < lines[m].size(); ++i ) {
	    pos[3*k+0] = lines[m][i].x[0];
	    pos[3*k+1] = lines[m][i].x[1];
	    pos[3*k+2] = lines[m][i].x[2];
	    tim[k]     = lines[m][i].t;
	    ++k;
	}
    }
    off[lines.size()] = (long long)k;

    if( npts > 0 ) {
	const hsize_t dp[2] = { (hsize_t)npts, 3 };
	const hsize_t cp[2] = { (hsize_t)npts > 1024 ? 1024 : (hsize_t)npts, 3 };
	write_dataset( fid.id(), lcpl, (grp_path + "/position").c_str(),
		       H5T_NATIVE_DOUBLE, 2, dp, cp, &pos[0] );

	const hsize_t dt[1] = { (hsize_t)npts };
	const hsize_t ct[1] = { (hsize_t)(npts > 4096 ? 4096 : npts) };
	write_dataset( fid.id(), lcpl, (grp_path + "/time").c_str(),
		       H5T_NATIVE_DOUBLE, 1, dt, ct, &tim[0] );
    }

    if( !lines.empty() ) {
	const hsize_t dm[1] = { (hsize_t)(lines.size() + 1) };
	const hsize_t cm[1] = { (hsize_t)(lines.size() + 1) };
	write_dataset( fid.id(), lcpl, (grp_path + "/offset").c_str(),
		       H5T_NATIVE_INT64, 1, dm, cm, &off[0] );

	const hsize_t dn[1] = { (hsize_t)lines.size() };
	const hsize_t cn[1] = { (hsize_t)lines.size() };
	write_dataset( fid.id(), lcpl, (grp_path + "/count").c_str(),
		       H5T_NATIVE_INT64, 1, dn, cn, &cnt[0] );
    }
}


#else // IBSIMU_CYCL_HAVE_HDF5


// ---------------------------------------------------------------------------
// 未编译进 HDF5 后端时的桩实现（与 gputracker.cpp 的处理方式一致）
// ---------------------------------------------------------------------------

bool hdf5_available()          { return( false ); }
std::string hdf5_version_string() { return( std::string() ); }

static const char *NO_HDF5 =
    "HDF5 backend not built (rerun configure with the HDF5 library available)";


void hdf5_write_image_data( const std::string &, const std::string &,
			    const Int3D &, const Vec3D &, const Vec3D &,
			    const std::string &, const std::vector<double> &,
			    double )
{
    throw( std::runtime_error( std::string("hdf5_write_image_data: ")
			       + NO_HDF5 ) );
}


void hdf5_write_image_data( const std::string &, const std::string &,
			    const Int3D &, const Vec3D &, const Vec3D &,
			    const std::string &, const std::vector<double> &,
			    double,
			    const std::string &, const std::vector<double> &,
			    double )
{
    throw( std::runtime_error( std::string("hdf5_write_image_data: ")
			       + NO_HDF5 ) );
}


void hdf5_write_polylines( const std::string &, const std::string &,
			   const std::vector<std::vector<TrajectoryPoint>> & )
{
    throw( std::runtime_error( std::string("hdf5_write_polylines: ")
			       + NO_HDF5 ) );
}


#endif // IBSIMU_CYCL_HAVE_HDF5


} // namespace ibsimu_cycl
