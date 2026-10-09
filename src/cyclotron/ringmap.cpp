/** @file ringmap.cpp
 *  @brief IBSimu-Cycl: OPAL/PSI "RING" 磁场图读取器（实现）。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "ringmap.hpp"

namespace ibsimu_cycl {

/*! \brief 跳过 \a n 个以空白分隔的记号。 */
static void skip_tokens( std::istream &is, int n, const char *what )
{
    std::string tok;
    for( int i = 0; i < n; ++i ) {
        if( !(is >> tok) )
            throw( std::runtime_error( std::string("read_ring_field_map: unexpected EOF while skipping ") + what ) );
    }
}

static void read_doubles( std::istream &is, std::vector<double> &out,
                          std::size_t count, const char *what )
{
    out.resize( count );
    for( std::size_t i = 0; i < count; ++i ) {
        if( !(is >> out[i]) )
            throw( std::runtime_error( std::string("read_ring_field_map: failed to read data (") + what + ")" ) );
    }
}

CCylFieldMap3D read_ring_field_map( const std::string &filename, double b_scale,
                                    int nsymmetry, RingFieldMapInfo *info )
{
    std::ifstream is( filename.c_str() );
    if( !is )
        throw( std::runtime_error("read_ring_field_map: cannot open " + filename) );

    // ---- 1. 头部四个数: rmin[mm] delr[mm] tetmin[deg] dtet[deg] ----
    double rmin_mm, dr_mm, tetmin_deg, dtet_deg;
    if( !(is >> rmin_mm >> dr_mm >> tetmin_deg >> dtet_deg) )
        throw( std::runtime_error("read_ring_field_map: failed to read header numbers") );

    if( dr_mm   < 0.0 ) dr_mm   = 1.0/(-dr_mm);    // 负值取倒数（OPAL 约定）
    if( dtet_deg < 0.0 ) dtet_deg = 1.0/(-dtet_deg);

    // ---- 2. 跳过 13 个记号, 读 nrad, ntet ----
    skip_tokens( is, 13, "label tokens" );
    int nrad = 0, ntet = 0;
    if( !(is >> nrad >> ntet) )
        throw( std::runtime_error("read_ring_field_map: failed to read nrad/ntet") );
    if( nrad < 2 || ntet < 2 )
        throw( std::runtime_error("read_ring_field_map: invalid grid size") );

    // ---- 3. 跳过 5 个记号, 读 lpar ----
    skip_tokens( is, 5, "tokens before lpar" );
    int lpar = 0;
    if( !(is >> lpar) )
        throw( std::runtime_error("read_ring_field_map: failed to read lpar") );

    // ---- 4. 跳过 4 个记号, 读 lpar 个双精度（信息块, 忽略） ----
    skip_tokens( is, 4, "tokens before info block" );
    {
        std::vector<double> info_block;
        read_doubles( is, info_block, (std::size_t)lpar, "info block" );
    }

    // ---- 5. 跳过 6 个记号, 扫描到 "LREC=" ----
    skip_tokens( is, 6, "tokens before LREC" );
    {
        std::string tok;
        bool found = false;
        for( int i = 0; i < 10000; ++i ) {
            if( !(is >> tok) )
                break;
            if( tok == "LREC=" ) { found = true; break; }
        }
        if( !found )
            throw( std::runtime_error("read_ring_field_map: 'LREC=' marker not found") );
    }
    skip_tokens( is, 5, "tokens after LREC" );

    // ---- 6. 数据块: 每半径块 4 个数组(各 ntet 个), 首块无信息行 ----
    const std::size_t ntot = (std::size_t)nrad * (std::size_t)ntet;
    std::vector<double> bfld( ntot ), dbt( ntot ), dbtt( ntot ), dbttt( ntot );
    for( int i = 0; i < nrad; ++i ) {
        if( i > 0 )
            skip_tokens( is, 6, "per-block info line" );
        for( int k = 0; k < ntet; ++k )
            if( !(is >> bfld[(std::size_t)i*ntet + k]) )
                throw( std::runtime_error("read_ring_field_map: truncated bfld data") );
        for( int k = 0; k < ntet; ++k )
            if( !(is >> dbt[(std::size_t)i*ntet + k]) )
                throw( std::runtime_error("read_ring_field_map: truncated dbt data") );
        for( int k = 0; k < ntet; ++k )
            if( !(is >> dbtt[(std::size_t)i*ntet + k]) )
                throw( std::runtime_error("read_ring_field_map: truncated dbtt data") );
        for( int k = 0; k < ntet; ++k )
            if( !(is >> dbttt[(std::size_t)i*ntet + k]) )
                throw( std::runtime_error("read_ring_field_map: truncated dbttt data") );
    }

    // ---- 7. 构造满 360deg 场图（通过 nsymmetry 折旋转复制） ----
    int ntheta_full = (int)std::lround( 360.0/dtet_deg );
    if( ntheta_full != ntet*nsymmetry )
        throw( std::runtime_error("read_ring_field_map: ntet*nsymmetry != 360/dtet "
                                  "(grid/symmetry mismatch)") );

    const double mm2m = 1.0e-3;
    CCylFieldMap3D map( (std::size_t)nrad, (std::size_t)ntheta_full, 1,
                        rmin_mm*mm2m, dr_mm*mm2m, 0.0, 1.0 );

    // 扇区起始角 tetmin 对应的整体索引偏移
    const double dtheta_deg = 360.0/(double)ntheta_full;
    const int offset = (int)std::lround( tetmin_deg/dtheta_deg );

    for( int i = 0; i < nrad; ++i )
        for( int j = 0; j < nsymmetry; ++j )
            for( int k = 0; k < ntet; ++k ) {
                int col = (offset + j*ntet + k) % ntheta_full;
                double bz = bfld[(std::size_t)i*ntet + k] * b_scale;
                map.set_value( (std::size_t)i, (std::size_t)col, 0,
                               0.0, 0.0, bz );
            }

    if( info ) {
        info->rmin_mm   = rmin_mm;
        info->dr_mm     = dr_mm;
        info->tetmin_deg = tetmin_deg;
        info->dtet_deg  = dtet_deg;
        info->nrad      = nrad;
        info->ntet      = ntet;
        info->lpar      = lpar;
    }

    return( map );
}

} // namespace ibsimu_cycl
