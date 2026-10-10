/** @file closedorbit.hpp
 *  @brief 闭合轨道求解：单扇区映射的不动点。
 *
 *  ## 为什么需要它
 *
 *  真实回旋加速器里束流走的是**闭合轨道**（单圈映射的不动点）。以「在 r=r_ref 处
 *  纯切向发射」作为初值起步得到的**不是**闭合轨道，而是围绕它做大幅 betatron
 *  振荡的轨迹——逐圈不重合，看起来像“轨道重叠且不稳定”。
 *
 *  实测（真实 PSI Ring 场图、r=3.3 m、γ=1.227）：
 *
 *  | 量 | 切向发射 | 闭合轨道 |
 *  | --- | --- | --- |
 *  | 逐圈重复性（同方位角半径差） | 0.24 m | ~1e-2 m |
 *  | 径向摆动 | 0.62 m | 0.13 m（即物理扇贝调制） |
 *  | 起始径向速度 | 0（切向） | -1.6e7 m/s（偏离切向约 5°） |
 *
 *  关键是：在强分离扇区场里（PSI Ring 的磁极只覆盖每个 45° 扇区中的约 19°），
 *  不动点在给定方位角上**本来就带非零径向速度**。
 *
 *  ## 方法
 *
 *  在方位角 \a th0 处把中平面轨道参数化为 (r, v_r)：|v| 由 γ 定，方位向速度
 *  $v_\theta=-\sqrt{|v|^2-v_r^2}$（该场配置下粒子顺时针运动）。利用场的 N 折对称，
 *  只需求**单扇区映射**的不动点（积分长度与噪声都降为 1/N，Jacobian 条件数也更好）：
 *  粗网格扫描定位 → 信赖域阻尼 Newton（数值 Jacobian）。
 *
 *  ## 精度限制
 *
 *  闭合残差的地板约 1e-3 m / 1e5 m/s，由**场图插值的间断**决定
 *  （`CRingFieldMap3D` 在 (r,θ) 上双线性插值，单元边界只有 C0 连续），
 *  不是算法问题。
 */

#ifndef IBSIMU_CYCL_CLOSEDORBIT_HPP
#define IBSIMU_CYCL_CLOSEDORBIT_HPP 1

#include "vec3d.hpp"
#include "vectorfield.hpp"


namespace ibsimu_cycl {


/*! \brief 闭合轨道在参考方位角处的状态与形状。 */
struct ClosedOrbit {
    double r;            /*!< 方位角 \a th0 处的半径 [m] */
    double vr;           /*!< 同名处的径向速度 [m/s] */
    double vmag;         /*!< 速率 [m/s]（由 γ 定） */
    double rmin;         /*!< 一整圈内的最小半径 [m] */
    double rmax;         /*!< 一整圈内的最大半径 [m] */
    double Bbar;         /*!< r_guess 处的方位平均 Bz [T]（诊断用） */
    int    iterations;   /*!< Newton 实际迭代次数 */
};


/*! \brief 求给定能量下的闭合轨道（成功返回；失败抛 std::runtime_error）。
 *
 *  \param B       三维磁场（中平面上的 Bz 决定轨道；单位 T）
 *  \param q       电荷 [C]
 *  \param m       质量 [kg]
 *  \param gamma   洛伦兹因子（决定速率 |v|）
 *  \param th0     参考方位角 [rad]（轨道在此处记录状态）
 *  \param r_guess 半径初值 [m]（一般取目标轨道半径）
 *  \param relativistic 是否用相对论推进器
 */
ClosedOrbit find_closed_orbit( const VectorField &B, double q, double m,
			       double gamma, double th0, double r_guess,
			       bool relativistic = true );


} // namespace ibsimu_cycl


#endif // IBSIMU_CYCL_CLOSEDORBIT_HPP
