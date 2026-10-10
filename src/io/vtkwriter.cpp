/** @file vtkwriter.cpp
 *  @brief VTK XML 写出器的实现（见 vtkwriter.hpp 的设计说明）。
 */

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>

#include "vtkwriter.hpp"


namespace ibsimu_cycl {

namespace {

/*! \brief 输出精度：12 位有效数字（约 1e-12 相对精度，足够可视化且文件不过大）。 */
const int PRECISION = 12;

/*! \brief VTK WholeExtent/Extent 字符串 "0 nx-1 0 ny-1 0 nz-1"。 */
std::string extent_string( const Int3D &dims )
{
    std::ostringstream ss;
    ss.imbue( std::locale::classic() );
    ss << "0 " << dims[0]-1 << " 0 " << dims[1]-1 << " 0 " << dims[2]-1;
    return( ss.str() );
}


/*! \brief 三分量字符串 "x y z"（用于 Origin / Spacing 属性）。 */
std::string vec3_string( const Vec3D &v )
{
    std::ostringstream ss;
    ss.imbue( std::locale::classic() );
    ss << std::setprecision( PRECISION )
       << v[0] << " " << v[1] << " " << v[2];
    return( ss.str() );
}


/*! \brief 按 VTK ASCII 内联格式写出一个 DataArray。
 *
 *  \param ncomp 每点的分量数（1 = 标量，3 = 矢量）
 *  \param nper  每行写多少个数值（仅影响可读性，VTK 按空白分割解析）
 */
template<class T>
void write_ascii_array( std::ostream &os, const std::string &type,
			const std::string &name, int ncomp,
			const std::vector<T> &data, int nper )
{
    os << "        <DataArray type=\"" << type << "\" Name=\"" << name
       << "\" NumberOfComponents=\"" << ncomp << "\" format=\"ascii\">\n";
    for( size_t a = 0; a < data.size(); a++ ) {
	if( a % nper == 0 )
	    os << "          ";
	else
	    os << " ";
	os << data[a];
	if( a % nper == nper-1 || a+1 == data.size() )
	    os << "\n";
    }
    os << "        </DataArray>\n";
}


/*! \brief 打开输出文件并设定数值格式，失败抛异常。 */
void open_out( std::ofstream &os, const std::string &filename, const char *who )
{
    os.open( filename.c_str() );
    if( !os )
	throw( std::runtime_error( std::string(who) + ": cannot open output file '" +
				   filename + "'" ) );
    os.imbue( std::locale::classic() );
    os << std::setprecision( PRECISION );
}

} // namespace


void vtk_write_image_data( const std::string &filename,
			   const Int3D &dims, const Vec3D &origo, const Vec3D &spacing,
			   const std::string &scalar_name,
			   const std::vector<double> &scalar )
{
    vtk_write_image_data( filename, dims, origo, spacing, scalar_name, scalar,
			  std::string(), std::vector<double>() );
}


void vtk_write_image_data( const std::string &filename,
			   const Int3D &dims, const Vec3D &origo, const Vec3D &spacing,
			   const std::string &scalar_name,
			   const std::vector<double> &scalar,
			   const std::string &vector_name,
			   const std::vector<double> &vector )
{
    const size_t n = (size_t)dims[0]*(size_t)dims[1]*(size_t)dims[2];
    if( scalar.size() != n )
	throw( std::invalid_argument( "vtk_write_image_data: scalar array size does not "
				      "match (nx*ny*nz)" ) );
    if( !vector.empty() && vector.size() != 3*n )
	throw( std::invalid_argument( "vtk_write_image_data: vector array size does not "
				      "match (3*nx*ny*nz)" ) );

    std::ofstream os;
    open_out( os, filename, "vtk_write_image_data" );

    const std::string extent = extent_string( dims );

    os << "<?xml version=\"1.0\"?>\n";
    os << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    os << "  <ImageData WholeExtent=\"" << extent
       << "\" Origin=\"" << vec3_string( origo )
       << "\" Spacing=\"" << vec3_string( spacing ) << "\">\n";
    os << "    <Piece Extent=\"" << extent << "\">\n";
    os << "      <PointData";
    if( !vector.empty() )
	os << " Vectors=\"" << vector_name << "\"";
    os << " Scalars=\"" << scalar_name << "\">\n";

    write_ascii_array( os, "Float64", scalar_name, 1, scalar, 6 );
    if( !vector.empty() )
	write_ascii_array( os, "Float64", vector_name, 3, vector, 3 );

    os << "      </PointData>\n";
    os << "    </Piece>\n";
    os << "  </ImageData>\n";
    os << "</VTKFile>\n";
}


void vtk_write_polylines( const std::string &filename,
			  const std::vector<std::vector<TrajectoryPoint>> &lines )
{
    size_t npts = 0;
    size_t nlines = 0;
    for( size_t a = 0; a < lines.size(); a++ ) {
	npts += lines[a].size();
	if( lines[a].size() >= 2 )
	    nlines++;
    }

    std::vector<double>  points;
    std::vector<double>  times;
    std::vector<int64_t> connectivity;
    std::vector<int64_t> offsets;
    points.reserve( 3*npts );
    times.reserve( npts );

    for( size_t a = 0; a < lines.size(); a++ ) {
	const std::vector<TrajectoryPoint> &ln = lines[a];
	const int64_t base = (int64_t)( points.size()/3 );

	for( size_t b = 0; b < ln.size(); b++ ) {
	    points.push_back( ln[b].x[0] );
	    points.push_back( ln[b].x[1] );
	    points.push_back( ln[b].x[2] );
	    times.push_back( ln[b].t );
	}

	if( ln.size() >= 2 ) {
	    for( size_t b = 0; b < ln.size(); b++ )
		connectivity.push_back( base + (int64_t)b );
	    offsets.push_back( (int64_t)connectivity.size() );
	}
    }

    std::ofstream os;
    open_out( os, filename, "vtk_write_polylines" );

    os << "<?xml version=\"1.0\"?>\n";
    os << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    os << "  <PolyData>\n";
    os << "    <Piece NumberOfPoints=\"" << npts
       << "\" NumberOfVerts=\"0\" NumberOfLines=\"" << nlines
       << "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";

    os << "      <PointData Scalars=\"t\">\n";
    write_ascii_array( os, "Float64", "t", 1, times, 6 );
    os << "      </PointData>\n";

    os << "      <Points>\n";
    write_ascii_array( os, "Float64", "Points", 3, points, 3 );
    os << "      </Points>\n";

    os << "      <Lines>\n";
    write_ascii_array( os, "Int64", "connectivity", 1, connectivity, 12 );
    write_ascii_array( os, "Int64", "offsets", 1, offsets, 12 );
    os << "      </Lines>\n";

    os << "    </Piece>\n";
    os << "  </PolyData>\n";
    os << "</VTKFile>\n";
}


} // namespace ibsimu_cycl
