/** @file gputracker.cpp
 *  @brief IBSimu-Cycl: GPU (CUDA) 粒子系综跟踪器（实现）。
 *
 *  实现基于 CUDA **Driver API + NVRTC**：内核源码以字符串内嵌
 *  （见 gpu_kernels.cu.h），在首次使用时由 NVRTC 编译成 PTX，再用 Driver API
 *  加载执行。这样构建系统只需要 CUDA 头文件与 libnvrtc/libcuda，无需 nvcc，
 *  也不必让 automake/libtool 处理 .cu 文件（详见 gpu_kernels.cu.h 的说明）。
 *
 *  未编译进 CUDA 支持时，本文件退化为"不可用"的桩实现，所有接口都能安全调用。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include "gputracker.hpp"
#include "gpu_kernels.cu.h"

namespace ibsimu_cycl {

#ifdef IBSIMU_CYCL_HAVE_CUDA

#include <cuda.h>
#include <nvrtc.h>

namespace {

const int IEC_BLOCK = 256;    /*!< \brief 每块线程数。 */

std::string cu_str( CUresult r )
{
    const char *s = 0;
    if( cuGetErrorString( r, &s ) != CUDA_SUCCESS || !s )
        return( "CUDA error" );
    return( std::string(s) );
}

std::string nvrtc_str( nvrtcResult r )
{
    const char *s = nvrtcGetErrorString( r );
    return( std::string( s ? s : "NVRTC error" ) );
}

#define IEC_CU(expr)                                                              \
    do {                                                                          \
        CUresult _r = (expr);                                                     \
        if( _r != CUDA_SUCCESS )                                                  \
            throw( std::runtime_error( std::string("CUDA: ") + cu_str(_r) ) );    \
    } while( 0 )

#define IEC_NVRTC(expr)                                                           \
    do {                                                                          \
        nvrtcResult _r = (expr);                                                  \
        if( _r != NVRTC_SUCCESS )                                                 \
            throw( std::runtime_error( std::string("NVRTC: ") + nvrtc_str(_r) ) );\
    } while( 0 )

/*! \brief 全局 CUDA 状态（进程内只初始化一次）。
 *
 *  故意用 new 且不释放：避免静态析构顺序导致"上下文先于使用者销毁"。
 */
struct CudaGlobal {
    bool         ok;
    std::string  err;
    CUdevice     dev;
    CUcontext    ctx;
    CUmodule     mod;
    CUfunction   fn;
    std::string  name;
    int          cc_major, cc_minor;
};

CudaGlobal *g_cuda = 0;

/*! \brief 初始化 CUDA 与内核；失败时把原因记在 \a g.err 而不抛异常。 */
CudaGlobal &cuda_global()
{
    if( g_cuda )
        return( *g_cuda );

    CudaGlobal *g      = new CudaGlobal();
    g->ok              = false;
    g->dev             = 0;
    g->ctx             = 0;
    g->mod             = 0;
    g->fn              = 0;
    g->cc_major        = 0;
    g->cc_minor        = 0;

    nvrtcProgram prog = 0;
    try {
        CUresult r = cuInit( 0 );
        if( r != CUDA_SUCCESS ) {
            g->err = "cuInit failed: " + cu_str( r );
            g_cuda = g;
            return( *g );
        }

        IEC_CU( cuDeviceGet( &g->dev, 0 ) );

        char name[256] = { 0 };
        IEC_CU( cuDeviceGetName( name, sizeof(name), g->dev ) );
        g->name = name;

        IEC_CU( cuDeviceGetAttribute( &g->cc_major,
                                      CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, g->dev ) );
        IEC_CU( cuDeviceGetAttribute( &g->cc_minor,
                                      CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, g->dev ) );

        IEC_CU( cuDevicePrimaryCtxRetain( &g->ctx, g->dev ) );
        IEC_CU( cuCtxSetCurrent( g->ctx ) );

        // ---- NVRTC：把内嵌源码编译成 PTX ----
        IEC_NVRTC( nvrtcCreateProgram( &prog, GPU_KERNEL_SOURCE, "ibsimu_cycl.cu", 0, 0, 0 ) );

        const std::string arch = "--gpu-architecture=compute_"
                               + std::to_string( g->cc_major*10 + g->cc_minor );
        // --fmad=false: 禁止 FMA 合并，使 GPU 结果与 CPU 端更接近（便于交叉验证）
        const char *opts[] = { arch.c_str(), "--std=c++14", "--fmad=false" };

        nvrtcResult nr = nvrtcCompileProgram( prog, 3, opts );
        if( nr != NVRTC_SUCCESS ) {
            std::size_t lsz = 0;
            nvrtcGetProgramLogSize( prog, &lsz );
            std::string log( lsz ? lsz : 1, '\0' );
            if( lsz > 0 )
                nvrtcGetProgramLog( prog, &log[0] );
            g->err = "NVRTC compile failed (" + arch + "): " + log;
            nvrtcDestroyProgram( &prog );
            g_cuda = g;
            return( *g );
        }

        std::size_t psz = 0;
        IEC_NVRTC( nvrtcGetPTXSize( prog, &psz ) );
        std::vector<char> ptx( psz );
        IEC_NVRTC( nvrtcGetPTX( prog, ptx.data() ) );
        nvrtcDestroyProgram( &prog );
        prog = 0;

        IEC_CU( cuModuleLoadData( &g->mod, ptx.data() ) );
        IEC_CU( cuModuleGetFunction( &g->fn, g->mod, "iec_track" ) );

        g->ok = true;
    } catch( const std::exception &e ) {
        g->err = e.what();
        if( prog )
            nvrtcDestroyProgram( &prog );
    }

    g_cuda = g;
    return( *g );
}

} // namespace

struct CGpuEnsembleTracker::Impl {
    CUdeviceptr d_fld;            /*!< \brief 显存中的场数据（交错 6 double/节点）。 */
    std::size_t fld_doubles;

    CUdeviceptr dp[6];            /*!< \brief 显存中的粒子 SoA: x,y,z,vx,vy,vz。 */
    std::size_t cap;              /*!< \brief 已分配的粒子容量。 */

    int         nr, nt;
    double      r0, dr, dtheta;
    double      qm;

    CUevent     ev0, ev1;
    double      t_kernel;
    double      t_transfer;
};

CGpuEnsembleTracker::CGpuEnsembleTracker( const CRingFieldMap3D &B, double q, double m )
    : _relativistic(false)
{
    if( m <= 0.0 )
        throw( std::invalid_argument("CGpuEnsembleTracker: mass must be > 0") );

    Impl *im = new Impl();
    im->d_fld       = 0;
    im->fld_doubles = 0;
    for( int c = 0; c < 6; ++c )
        im->dp[c] = 0;
    im->cap        = 0;
    im->nr         = (int)B.size_r();
    im->nt         = (int)B.size_theta();
    im->r0         = B.r0();
    im->dr         = B.dr();
    im->dtheta     = B.dtheta();
    im->qm         = q/m;
    im->ev0        = 0;
    im->ev1        = 0;
    im->t_kernel   = 0.0;
    im->t_transfer = 0.0;
    _impl = im;

    CudaGlobal &g = cuda_global();
    if( !g.ok )
        throw( std::runtime_error( "CGpuEnsembleTracker: " + g.err ) );

    IEC_CU( cuCtxSetCurrent( g.ctx ) );

    // ---- 场数据一次性上传（之后每步不再传场）----
    const std::size_t nc = CRingFieldMap3D::component_count();
    im->fld_doubles = (std::size_t)im->nr*(std::size_t)im->nt*nc;
    IEC_CU( cuMemAlloc( &im->d_fld, im->fld_doubles*sizeof(double) ) );
    IEC_CU( cuMemcpyHtoD( im->d_fld, B.raw_data(), im->fld_doubles*sizeof(double) ) );

    IEC_CU( cuEventCreate( &im->ev0, CU_EVENT_DEFAULT ) );
    IEC_CU( cuEventCreate( &im->ev1, CU_EVENT_DEFAULT ) );
}

CGpuEnsembleTracker::~CGpuEnsembleTracker()
{
    if( !_impl )
        return;

    CudaGlobal &g = cuda_global();
    if( g.ok )
        cuCtxSetCurrent( g.ctx );

    if( _impl->ev1 ) cuEventDestroy( _impl->ev1 );
    if( _impl->ev0 ) cuEventDestroy( _impl->ev0 );
    for( int c = 0; c < 6; ++c )
        if( _impl->dp[c] ) cuMemFree( _impl->dp[c] );
    if( _impl->d_fld ) cuMemFree( _impl->d_fld );

    delete _impl;
    _impl = 0;
}

void CGpuEnsembleTracker::track( std::vector<EnsembleParticle> &p, double dt, int nsteps )
{
    if( p.empty() || nsteps <= 0 )
        return;

    CudaGlobal &g = cuda_global();
    if( !g.ok )
        throw( std::runtime_error( "CGpuEnsembleTracker: " + g.err ) );

    IEC_CU( cuCtxSetCurrent( g.ctx ) );

    const std::size_t n = p.size();

    // ---- 按需扩容粒子缓冲 ----
    if( n > _impl->cap ) {
        for( int c = 0; c < 6; ++c ) {
            if( _impl->dp[c] )
                IEC_CU( cuMemFree( _impl->dp[c] ) );
            IEC_CU( cuMemAlloc( &_impl->dp[c], n*sizeof(double) ) );
        }
        _impl->cap = n;
    }

    // ---- 打包 SoA 并上传 ----
    std::vector<double> hx(n), hy(n), hz(n), hvx(n), hvy(n), hvz(n);
    for( std::size_t i = 0; i < n; ++i ) {
        hx[i]  = p[i].x[0]; hy[i]  = p[i].x[1]; hz[i]  = p[i].x[2];
        hvx[i] = p[i].v[0]; hvy[i] = p[i].v[1]; hvz[i] = p[i].v[2];
    }

    CUevent ev_t0 = 0, ev_t1 = 0;
    IEC_CU( cuEventCreate( &ev_t0, CU_EVENT_DEFAULT ) );
    IEC_CU( cuEventCreate( &ev_t1, CU_EVENT_DEFAULT ) );
    IEC_CU( cuEventRecord( ev_t0, 0 ) );

    IEC_CU( cuMemcpyHtoD( _impl->dp[0], hx.data(),  n*sizeof(double) ) );
    IEC_CU( cuMemcpyHtoD( _impl->dp[1], hy.data(),  n*sizeof(double) ) );
    IEC_CU( cuMemcpyHtoD( _impl->dp[2], hz.data(),  n*sizeof(double) ) );
    IEC_CU( cuMemcpyHtoD( _impl->dp[3], hvx.data(), n*sizeof(double) ) );
    IEC_CU( cuMemcpyHtoD( _impl->dp[4], hvy.data(), n*sizeof(double) ) );
    IEC_CU( cuMemcpyHtoD( _impl->dp[5], hvz.data(), n*sizeof(double) ) );

    // ---- 启动内核（每个粒子一个线程，一次跑完 nsteps 步）----
    const double *a_fld = (const double *)(uintptr_t)_impl->d_fld;
    int    a_nr = _impl->nr, a_nt = _impl->nt;
    double a_r0 = _impl->r0, a_dr = _impl->dr, a_dtheta = _impl->dtheta;
    double a_qm = _impl->qm, a_dt = dt;
    int    a_rel = _relativistic ? 1 : 0;
    int    a_nsteps = nsteps;
    double *a_px = (double *)(uintptr_t)_impl->dp[0];
    double *a_py = (double *)(uintptr_t)_impl->dp[1];
    double *a_pz = (double *)(uintptr_t)_impl->dp[2];
    double *a_vx = (double *)(uintptr_t)_impl->dp[3];
    double *a_vy = (double *)(uintptr_t)_impl->dp[4];
    double *a_vz = (double *)(uintptr_t)_impl->dp[5];
    long   a_n  = (long)n;

    void *args[] = { &a_fld, &a_nr, &a_nt, &a_r0, &a_dr, &a_dtheta,
                     &a_qm, &a_dt, &a_rel, &a_nsteps,
                     &a_px, &a_py, &a_pz, &a_vx, &a_vy, &a_vz, &a_n };

    const int blocks = (int)((n + IEC_BLOCK - 1)/IEC_BLOCK);

    IEC_CU( cuEventRecord( _impl->ev0, 0 ) );
    IEC_CU( cuLaunchKernel( g.fn, (unsigned)blocks, 1, 1, IEC_BLOCK, 1, 1,
                            0, 0, args, 0 ) );
    IEC_CU( cuEventRecord( _impl->ev1, 0 ) );
    IEC_CU( cuEventSynchronize( _impl->ev1 ) );

    float ms = 0.0f;
    IEC_CU( cuEventElapsedTime( &ms, _impl->ev0, _impl->ev1 ) );
    _impl->t_kernel = (double)ms/1000.0;

    // ---- 回传并解包 ----
    IEC_CU( cuMemcpyDtoH( hx.data(),  _impl->dp[0], n*sizeof(double) ) );
    IEC_CU( cuMemcpyDtoH( hy.data(),  _impl->dp[1], n*sizeof(double) ) );
    IEC_CU( cuMemcpyDtoH( hz.data(),  _impl->dp[2], n*sizeof(double) ) );
    IEC_CU( cuMemcpyDtoH( hvx.data(), _impl->dp[3], n*sizeof(double) ) );
    IEC_CU( cuMemcpyDtoH( hvy.data(), _impl->dp[4], n*sizeof(double) ) );
    IEC_CU( cuMemcpyDtoH( hvz.data(), _impl->dp[5], n*sizeof(double) ) );

    IEC_CU( cuEventRecord( ev_t1, 0 ) );
    IEC_CU( cuEventSynchronize( ev_t1 ) );
    IEC_CU( cuEventElapsedTime( &ms, ev_t0, ev_t1 ) );
    _impl->t_transfer = (double)ms/1000.0 - _impl->t_kernel;

    IEC_CU( cuEventDestroy( ev_t1 ) );
    IEC_CU( cuEventDestroy( ev_t0 ) );

    for( std::size_t i = 0; i < n; ++i ) {
        p[i].x[0] = hx[i];  p[i].x[1] = hy[i];  p[i].x[2] = hz[i];
        p[i].v[0] = hvx[i]; p[i].v[1] = hvy[i]; p[i].v[2] = hvz[i];
    }
}

double CGpuEnsembleTracker::last_kernel_time() const   { return( _impl->t_kernel ); }
double CGpuEnsembleTracker::last_transfer_time() const { return( _impl->t_transfer ); }

bool CGpuEnsembleTracker::compiled_in() { return( true ); }

bool CGpuEnsembleTracker::available( std::string *reason )
{
    CudaGlobal &g = cuda_global();
    if( !g.ok && reason )
        *reason = g.err;
    return( g.ok );
}

std::string CGpuEnsembleTracker::device_name()
{
    CudaGlobal &g = cuda_global();
    return( g.ok ? g.name : std::string() );
}

#else // IBSIMU_CYCL_HAVE_CUDA

// ---------------------------------------------------------------------------
// 未编译进 CUDA 支持时的桩实现
// ---------------------------------------------------------------------------

struct CGpuEnsembleTracker::Impl {
    std::string err;
};

static const char *NO_CUDA =
    "CUDA backend not built (rerun configure with CUDA headers/libs available)";

CGpuEnsembleTracker::CGpuEnsembleTracker( const CRingFieldMap3D &, double q, double m )
    : _relativistic(false)
{
    if( m <= 0.0 )
        throw( std::invalid_argument("CGpuEnsembleTracker: mass must be > 0") );
    _impl = new Impl();
    _impl->err = NO_CUDA;
}

CGpuEnsembleTracker::~CGpuEnsembleTracker()
{
    delete _impl;
    _impl = 0;
}

void CGpuEnsembleTracker::track( std::vector<EnsembleParticle> &, double, int )
{
    throw( std::runtime_error( std::string("CGpuEnsembleTracker: ") + NO_CUDA ) );
}

double CGpuEnsembleTracker::last_kernel_time() const   { return( 0.0 ); }
double CGpuEnsembleTracker::last_transfer_time() const { return( 0.0 ); }

bool CGpuEnsembleTracker::compiled_in() { return( false ); }

bool CGpuEnsembleTracker::available( std::string *reason )
{
    if( reason )
        *reason = NO_CUDA;
    return( false );
}

std::string CGpuEnsembleTracker::device_name() { return( std::string() ); }

#endif // IBSIMU_CYCL_HAVE_CUDA

} // namespace ibsimu_cycl
