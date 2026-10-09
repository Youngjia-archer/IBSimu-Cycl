/** @file boris.hpp
 *  @brief IBSimu-Cycl: Boris 粒子推进器（支持时变 E/B、非相对论与相对论）。
 *
 *  Boris 算法（Birdsall & Langdon）是 PIC 中标准的、能量守恒的推进器。
 *  本类支持：
 *    - 时变电场（如 RF 场，见 CTimeVaryingField）；
 *    - 静态或时变磁场；
 *    - E、B 均可为 NULL（视为 0）；
 *    - **非相对论**（默认）与**相对论**两种模式。
 *
 *  时变场的时间由**调用者**在每步前设置（例如 `efield.set_time(t)`），
 *  这样推进器本身无需时间参数，与 IBSimu 的 `VectorField` 接口保持一致。
 *
 *  相对论模式使用标准 relativistic Boris（在 u = gamma*v 空间做磁场旋转），
 *  静磁场下 |u| 与 gamma 精确守恒。真实回旋加速器（如 PSI Ring）束流
 *  beta ~ 0.6，必须使用相对论模式。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#ifndef IBSIMU_CYCL_BORIS_HPP
#define IBSIMU_CYCL_BORIS_HPP 1

#include "vec3d.hpp"
#include "vectorfield.hpp"

namespace ibsimu_cycl {

/*! \brief Boris 推进器（非相对论 / 相对论）。 */
class CBorisPusher {

    double _q;             /*!< \brief 电荷 [C]。 */
    double _m;             /*!< \brief 质量 [kg]。 */
    bool   _relativistic;  /*!< \brief 是否使用相对论动力学。 */

public:

    /*! \brief 构造（默认非相对论）。
     *  \param q 电荷 [C]（含符号）
     *  \param m 质量 [kg]（> 0）
     */
    CBorisPusher( double q, double m );

    /*! \brief 切换相对论模式。 */
    void set_relativistic( bool on ) { _relativistic = on; }
    bool is_relativistic() const { return( _relativistic ); }

    double charge() const { return( _q ); }
    double mass() const { return( _m ); }
    double charge_to_mass() const { return( _q/_m ); }

    /*! \brief 洛伦兹因子 gamma（由速度计算，非相对论模式返回 1）。 */
    double gamma( const Vec3D &v ) const;

    /*! \brief 推进一步（原地更新 \a x, \a v）。
     *
     *  \param E  电场（时变者需在调用前 set_time(t)），可为 NULL
     *  \param B  磁场，可为 NULL
     *  \param x  位置 [m]（输入/输出）
     *  \param v  速度 [m/s]（输入/输出）
     *  \param dt 时间步 [s]
     */
    void step( const VectorField *E, const VectorField *B,
               Vec3D &x, Vec3D &v, double dt ) const;
};

} // namespace ibsimu_cycl

#endif // IBSIMU_CYCL_BORIS_HPP
