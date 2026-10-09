/** @file timevaryingfield.hpp
 *  @brief IBSimu-Cycl: 时变矢量场（RF 场）。
 *
 *  将空间模式场 A(x)（例如 RF 谐振腔的模式场图）包装为随时间振荡的场：
 *      F(x, t) = A(x) * ( offset + amplitude * cos( omega*t + phi0 ) )
 *
 *  IBSimu 的粒子步进器通过 VectorField::operator()(x) 取场，因此当前时间
 *  由 set_time() 显式设置（步进器每步调用一次）。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
 */

#ifndef IBSIMU_CYCL_TIMEVARYINGFIELD_HPP
#define IBSIMU_CYCL_TIMEVARYINGFIELD_HPP 1

#include "vectorfield.hpp"
#include "vec3d.hpp"

namespace ibsimu_cycl {

/*! \brief 时变矢量场。
 *
 *  空间模式由 \a spatial 提供（不接管所有权，可在多个时变场间共享）。
 *  时间行为由 \a omega（角频率, rad/s）、\a phi0（初相, rad）、
 *  \a amplitude（交流幅度）与 \a offset（直流偏置）决定。
 */
class CTimeVaryingField : public VectorField {

    const VectorField *_spatial;  /*!< \brief 空间模式场（不拥有）。 */
    double             _omega;    /*!< \brief 角频率 (rad/s)。 */
    double             _phi0;     /*!< \brief 初相位 (rad)。 */
    double             _amplitude;/*!< \brief 交流幅度因子。 */
    double             _offset;   /*!< \brief 直流偏置因子。 */
    double             _t;        /*!< \brief 当前时间 (s)。 */

    CTimeVaryingField( const CTimeVaryingField & );            /*!< \brief 禁止拷贝。 */
    CTimeVaryingField &operator=( const CTimeVaryingField & ); /*!< \brief 禁止赋值。 */

public:

    /*! \brief 构造函数。
     *
     *  \param spatial   空间模式场（不接管所有权，不得为 NULL）
     *  \param omega     角频率 (rad/s)
     *  \param phi0      初相位 (rad)
     *  \param amplitude 交流幅度因子（默认 1）
     *  \param offset    直流偏置因子（默认 0）
     */
    CTimeVaryingField( const VectorField *spatial, double omega,
                       double phi0 = 0.0, double amplitude = 1.0,
                       double offset = 0.0 );

    virtual ~CTimeVaryingField();

    /*! \brief 设置当前时间 (s)。 */
    void set_time( double t ) { _t = t; }

    /*! \brief 当前时间 (s)。 */
    double time() const { return( _t ); }

    /*! \brief 角频率 (rad/s)。 */
    double omega() const { return( _omega ); }

    /*! \brief 当前相位 omega*t + phi0 (rad)。 */
    double phase() const { return( _omega*_t + _phi0 ); }

    /*! \brief 当前调制因子 offset + amplitude*cos(phase)。 */
    double scale_factor() const;

    /*! \brief 场求值（返回直角坐标分量）。 */
    virtual const Vec3D operator()( const Vec3D &x ) const;
};

} // namespace ibsimu_cycl

#endif // IBSIMU_CYCL_TIMEVARYINGFIELD_HPP
