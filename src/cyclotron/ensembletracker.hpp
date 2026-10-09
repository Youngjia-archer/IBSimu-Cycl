/** @file ensembletracker.hpp
 *  @brief IBSimu-Cycl: 粒子系综跟踪器（OpenMP 并行）。
 *
 *  回旋加速器的三维跟踪天然**数据并行**：每个粒子的推进只依赖只读的场
 *  （`VectorField::operator()` 为 const），粒子之间无相互作用（空间电荷除外）。
 *  因此可用 OpenMP 在粒子维并行。
 *
 *  时变场的时间是全局量：每步先在**串行区**设置时间，再进入并行区推进所有粒子，
 *  保证同一时间步内所有粒子看到相同的 RF 相位。
 *
 *  实现要点：只建立一个 OpenMP 并行区，步间用屏障同步（而不是每步一次
 *  `parallel for` 的 fork/join）。这样即使场求值很快（小场图/解析场），
 *  并行收益也不会被启动开销吃掉。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#ifndef IBSIMU_CYCL_ENSEMBLETRACKER_HPP
#define IBSIMU_CYCL_ENSEMBLETRACKER_HPP 1

#include <cstddef>
#include <vector>

#include "vec3d.hpp"
#include "vectorfield.hpp"
#include "timevaryingfield.hpp"

namespace ibsimu_cycl {

/*! \brief 系综中的一个粒子（位置 + 速度）。 */
struct EnsembleParticle {
    Vec3D x;
    Vec3D v;
};

/*! \brief 系综统计量。 */
struct EnsembleStats {
    double KE_min, KE_max, KE_mean;   /*!< \brief 动能 [J]。 */
    double r_min,  r_max;             /*!< \brief 半径 [m]。 */
    std::size_t n;                    /*!< \brief 粒子数。 */
};

/*! \brief 粒子系综跟踪器（OpenMP 并行）。 */
class CEnsembleTracker {

    const VectorField  *_B;   /*!< \brief 磁场（不拥有，可空）。 */
    CTimeVaryingField  *_E;   /*!< \brief 时变电场（不拥有，可空）。 */
    double              _q;   /*!< \brief 电荷 [C]。 */
    double              _m;   /*!< \brief 质量 [kg]。 */
    bool                _relativistic;

public:

    CEnsembleTracker( const VectorField *B, CTimeVaryingField *E,
                      double q, double m );

    void set_relativistic( bool on ) { _relativistic = on; }
    bool is_relativistic() const { return( _relativistic ); }

    /*! \brief 并行推进 \a nsteps 步，粒子状态原地更新。
     *
     *  \param p      粒子数组（输入/输出）
     *  \param dt     时间步 [s]
     *  \param nsteps 步数
     */
    void track( std::vector<EnsembleParticle> &p, double dt, int nsteps ) const;

    /*! \brief 同 \a track，但每步后记录系综统计量（用于诊断/可视化）。 */
    void track( std::vector<EnsembleParticle> &p, double dt, int nsteps,
                std::vector<EnsembleStats> &history ) const;

    /*! \brief 计算系综统计量。 */
    static EnsembleStats analyze( const std::vector<EnsembleParticle> &p,
                                  double m, bool relativistic );
};

} // namespace ibsimu_cycl

#endif // IBSIMU_CYCL_ENSEMBLETRACKER_HPP
