/** @file cycl_info.hpp
 *  @brief IBSimu-Cycl: 模块版本与构建信息。
 *
 *  这是 IBSimu-Cycl 新增模块的入口头文件。P0 阶段仅提供版本/基线信息，
 *  用于验证"新增模块已正确接入上游 autotools 构建"。
 *
 *  IBSimu-Cycl 是 IBSimu 的派生作品，整体以 GPL-3.0-or-later 发布。
 *  上游 IBSimu 版权 (c) 2004-2013 Taneli Kalvas 及 The Regents of the
 *  University of California / LBNL，以 GPL-2.0-or-later 授权。
 */

#ifndef IBSIMU_CYCL_INFO_HPP
#define IBSIMU_CYCL_INFO_HPP

namespace ibsimu_cycl {

/// IBSimu-Cycl 模块版本号。
const char *version();

/// 本项目所基于的上游 IBSimu 基线提交。
const char *upstream_baseline();

} // namespace ibsimu_cycl

#endif // IBSIMU_CYCL_INFO_HPP
