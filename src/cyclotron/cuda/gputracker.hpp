/** @file gputracker.hpp
 *  @brief IBSimu-Cycl: GPU (CUDA) 粒子系综跟踪器（R1-③）。
 *
 *  粒子维并行是回旋加速器三维跟踪的天然并行维，且 GPU 的吞吐优势远高于
 *  多核 CPU。本类把整个系综放到 GPU 上推进：
 *
 *    - 场数据（CRingFieldMap3D 的交错节点表）**一次性上传**显存，之后不再
 *      往返主机，避免每步传场；
 *    - 每个 GPU 线程推进一个粒子，跑完整的 nsteps 步（而不是每步启动一次
 *      kernel），因此没有"每步 fork/join"式的启动开销；
 *    - 粒子状态以 SoA（x/y/z/vx/vy/vz 六个数组）传输，保证合并访存。
 *
 *  代价是每次 track() 要做一次 H2D/D2H 传输，因此系综越大、步数越多，
 *  GPU 的相对优势越明显（见 tests/cycl_cuda_tracker.cpp）。
 *
 *  \note 当前只支持**静态磁场**（E = 0）；时变/射频场仍在 CPU 端处理。
 *  \note 本头文件不包含任何 CUDA 头，因此可用普通 C++ 编译器编译。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#ifndef IBSIMU_CYCL_GPUTRACKER_HPP
#define IBSIMU_CYCL_GPUTRACKER_HPP 1

#include <cstddef>
#include <string>
#include <vector>

#include "ensembletracker.hpp"
#include "ringfield3d.hpp"

namespace ibsimu_cycl {

/*! \brief GPU（CUDA）粒子系综跟踪器（静态磁场）。 */
class CGpuEnsembleTracker {

    struct Impl;
    Impl  *_impl;                 /*!< \brief 设备资源（Pimpl，隐藏 CUDA 依赖）。 */
    bool   _relativistic;

    CGpuEnsembleTracker( const CGpuEnsembleTracker & );
    CGpuEnsembleTracker &operator=( const CGpuEnsembleTracker & );

public:

    /*! \brief 构建时是否启用了 CUDA（configure 是否找到 CUDA 头文件与库）。 */
    static bool compiled_in();

    /*! \brief 运行期是否可用（编译进 CUDA + 有可用 GPU + 内核编译成功）。
     *
     *  \param reason 若非空且返回 false，则填入不可用原因
     */
    static bool available( std::string *reason = 0 );

    /*! \brief 设备名（不可用时返回空串）。 */
    static std::string device_name();

    /*! \brief 构造：绑定磁场图与粒子参数（此时不上传数据）。
     *
     *  \param B 由中平面 Bz 重建的三维回旋加速器磁场（需支持 raw_data()）
     *  \param q 粒子电荷 [C]
     *  \param m 粒子质量 [kg]
     */
    CGpuEnsembleTracker( const CRingFieldMap3D &B, double q, double m );

    ~CGpuEnsembleTracker();

    void set_relativistic( bool on ) { _relativistic = on; }
    bool is_relativistic() const { return( _relativistic ); }

    /*! \brief 在 GPU 上推进 \a nsteps 步，粒子状态原地更新。
     *
     *  \param p      粒子数组（输入/输出）
     *  \param dt     时间步 [s]
     *  \param nsteps 步数
     *
     *  \throws std::runtime_error GPU 不可用或 CUDA 调用失败时
     */
    void track( std::vector<EnsembleParticle> &p, double dt, int nsteps );

    /*! \brief 上一次 track() 的设备端执行时间 [s]（不含数据传输）。 */
    double last_kernel_time() const;

    /*! \brief 上一次 track() 的数据传输时间 [s]（H2D + D2H）。 */
    double last_transfer_time() const;
};

} // namespace ibsimu_cycl

#endif // IBSIMU_CYCL_GPUTRACKER_HPP
