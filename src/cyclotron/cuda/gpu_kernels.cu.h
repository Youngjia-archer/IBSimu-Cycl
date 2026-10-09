/** @file gpu_kernels.cu.h
 *  @brief IBSimu-Cycl: GPU (CUDA) 内核源码（以字符串形式内嵌，运行期编译）。
 *
 *  为什么内嵌字符串而不是 .cu 文件？
 *  本项目的构建系统是 automake + libtool。libtool 在编译未知后缀时无法推断
 *  编译配置，加上 `--tag=CXX` 后会强制注入 `-fPIC -DPIC`，而 nvcc 不接受裸
 *  `-fPIC`（会 fatal），`-Xcompiler` 形式的参数又会被 libtool 重排到命令首部。
 *  要让 nvcc 进入 libtool 流程需要额外的包装脚本，代价与风险都偏高。
 *
 *  因此这里把内核源码作为字符串内嵌，由 NVRTC 在**运行期**编译为 PTX 再用
 *  CUDA Driver API 加载。好处：
 *    - 构建期只需要 CUDA 头文件与 libnvrtc/libcuda，不需要 nvcc；
 *    - 不触自动 automake/libtool，CI（无 CUDA）不受影响；
 *    - 无 GPU 时 available() 返回 false，可优雅降级。
 *
 *  内核中的场求值与推进算法与 CPU 端**逐行对应**，便于交叉验证：
 *    - 场求值   <-> CRingFieldMap3D::operator()   (ringfield3d.cpp)
 *    - Boris 步 <-> CBorisPusher::step()          (boris.cpp)
 *
 *  当前 GPU 后端只处理**静态磁场**（E = 0）；时变/射频场仍在 CPU 端。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#ifndef IBSIMU_CYCL_GPU_KERNELS_CU_H
#define IBSIMU_CYCL_GPU_KERNELS_CU_H 1

namespace ibsimu_cycl {

/*! \brief CUDA 内核源码（由 NVRTC 在运行期编译）。 */
static const char *GPU_KERNEL_SOURCE = R"NVRTC_SRC(
// IBSimu-Cycl GPU kernels (NVRTC). 场图数据为按节点交错的 6 个 double：
//   [0]=b  [1]=dbr  [2]=dbth  [3]=tb  [4]=dtrb  [5]=dttb
// 与 CRingFieldMap3D::raw_data() 的布局一致。

#define IEC_PI2  6.283185307179586476925286766559
#define IEC_CL2  8.9875517873681764e16   // c^2

// ---- 场求值：与 CRingFieldMap3D::operator() 逐行对应 ----
__device__ __forceinline__ void iec_eval_b(
    const double *__restrict__ fld, int nr, int nt,
    double r0, double dr, double dtheta,
    double xx, double yy, double zz,
    double &bx, double &by, double &bz )
{
    const double rr = sqrt( xx*xx + yy*yy );

    if( rr == 0.0 ) {          // 轴上：Br/Btheta 取轴对称极限 0
        bx = 0.0; by = 0.0; bz = fld[0];
        return;
    }

    double th = atan2( yy, xx );
    if( th < 0.0 ) th += IEC_PI2;

    // ---- 单元与权重（只算一次，6 个分量共用）----
    int i0;
    double tr;
    const double fr = (rr - r0)/dr;
    if( fr <= 0.0 ) {
        i0 = 0; tr = 0.0;
    } else {
        i0 = (int)fr;
        if( i0 >= nr-1 ) { i0 = nr-1; tr = 0.0; }
        else             { tr = fr - (double)i0; }
    }
    const int i1 = ( i0 < nr-1 ) ? i0+1 : i0;

    const double ft = th/dtheta;
    int j = (int)ft;
    double tt = ft - (double)j;
    if( tt < 0.0 ) tt = 0.0;
    j %= nt;
    if( j < 0 ) j += nt;
    const int j1 = ( j+1 == nt ) ? 0 : j+1;

    const double w00 = (1.0-tr)*(1.0-tt);
    const double w01 = (1.0-tr)*tt;
    const double w10 = tr*(1.0-tt);
    const double w11 = tr*tt;

    const double *p00 = fld + 6*((long)i0*nt + j );
    const double *p01 = fld + 6*((long)i0*nt + j1);
    const double *p10 = fld + 6*((long)i1*nt + j );
    const double *p11 = fld + 6*((long)i1*nt + j1);

    const double f0 = w00*p00[0] + w01*p01[0] + w10*p10[0] + w11*p11[0];
    const double f1 = w00*p00[1] + w01*p01[1] + w10*p10[1] + w11*p11[1];
    const double f2 = w00*p00[2] + w01*p01[2] + w10*p10[2] + w11*p11[2];
    const double f3 = w00*p00[3] + w01*p01[3] + w10*p10[3] + w11*p11[3];
    const double f4 = w00*p00[4] + w01*p01[4] + w10*p10[4] + w11*p11[4];
    const double f5 = w00*p00[5] + w01*p01[5] + w10*p10[5] + w11*p11[5];

    const double z2 = zz*zz;
    const double z3 = z2*zz;

    const double bz_mid = f0 - 0.5*f3*z2;
    const double br_cyl = f1*zz - (1.0/6.0)*f4*z3;
    const double bt_cyl = (f2/rr)*zz - (1.0/(6.0*rr))*f5*z3;

    const double cs = xx/rr;
    const double sn = yy/rr;
    bx = br_cyl*cs - bt_cyl*sn;
    by = br_cyl*sn + bt_cyl*cs;
    bz = bz_mid;
}

// ---- Boris 推进：与 CBorisPusher::step() 对应（E = 0）----
__device__ __forceinline__ void iec_boris_step(
    double qm, double dt, int relativistic,
    double bx, double by, double bz,
    double &px, double &py, double &pz,
    double &vx, double &vy, double &vz )
{
    if( !relativistic ) {
        const double h  = 0.5*qm*dt;
        const double tx = bx*h, ty = by*h, tz = bz*h;
        const double t2 = tx*tx + ty*ty + tz*tz;

        const double upx = vx + (vy*tz - vz*ty);
        const double upy = vy + (vz*tx - vx*tz);
        const double upz = vz + (vx*ty - vy*tx);

        const double s  = 2.0/(1.0 + t2);
        const double sx = tx*s, sy = ty*s, sz = tz*s;

        vx = vx + (upy*sz - upz*sy);
        vy = vy + (upz*sx - upx*sz);
        vz = vz + (upx*sy - upy*sx);

        px += vx*dt; py += vy*dt; pz += vz*dt;
        return;
    }

    // 相对论：在 u = gamma*v 空间旋转
    const double v2 = vx*vx + vy*vy + vz*vz;
    const double g  = 1.0/sqrt( 1.0 - v2/IEC_CL2 );
    double ux = vx*g, uy = vy*g, uz = vz*g;

    const double u2 = ux*ux + uy*uy + uz*uz;
    const double gh = sqrt( 1.0 + u2/IEC_CL2 );

    const double hh = qm*dt/(2.0*gh);
    const double tx = bx*hh, ty = by*hh, tz = bz*hh;
    const double t2 = tx*tx + ty*ty + tz*tz;

    const double upx = ux + (uy*tz - uz*ty);
    const double upy = uy + (uz*tx - ux*tz);
    const double upz = uz + (ux*ty - uy*tx);

    const double s  = 2.0/(1.0 + t2);
    const double sx = tx*s, sy = ty*s, sz = tz*s;

    ux = ux + (upy*sz - upz*sy);
    uy = uy + (upz*sx - upx*sz);
    uz = uz + (upx*sy - upy*sx);

    const double u2n = ux*ux + uy*uy + uz*uz;
    const double gn  = sqrt( 1.0 + u2n/IEC_CL2 );
    const double ig  = 1.0/gn;

    vx = ux*ig; vy = uy*ig; vz = uz*ig;

    px += vx*dt; py += vy*dt; pz += vz*dt;
}

// ---- 主内核：一个线程推进一个粒子，粒子维并行 ----
extern "C" __global__ void iec_track(
    const double *__restrict__ fld, int nr, int nt,
    double r0, double dr, double dtheta,
    double qm, double dt, int relativistic, int nsteps,
    double *__restrict__ px, double *__restrict__ py, double *__restrict__ pz,
    double *__restrict__ vx, double *__restrict__ vy, double *__restrict__ vz,
    long n )
{
    const long i = (long)blockIdx.x*blockDim.x + threadIdx.x;
    if( i >= n ) return;

    double X  = px[i], Y  = py[i], Z  = pz[i];
    double VX = vx[i], VY = vy[i], VZ = vz[i];

    for( int s = 0; s < nsteps; ++s ) {
        double bx, by, bz;
        iec_eval_b( fld, nr, nt, r0, dr, dtheta, X, Y, Z, bx, by, bz );
        iec_boris_step( qm, dt, relativistic, bx, by, bz, X, Y, Z, VX, VY, VZ );
    }

    px[i] = X;  py[i] = Y;  pz[i] = Z;
    vx[i] = VX; vy[i] = VY; vz[i] = VZ;
}
)NVRTC_SRC";

} // namespace ibsimu_cycl

#endif // IBSIMU_CYCL_GPU_KERNELS_CU_H
