/** @file cycl_info.cpp
 *  @brief IBSimu-Cycl: 模块版本与构建信息（实现）。
 *
 *  IBSimu-Cycl is a derivative work of IBSimu, licensed under
 *  GPL-3.0-or-later. See NOTICE and LICENSE in the project root.
 */

#include "cycl_info.hpp"

namespace ibsimu_cycl {

const char *version()
{
    return "0.0.1";
}

const char *upstream_baseline()
{
    // IBSimu upstream commit used as the starting baseline.
    return "09beedb (2026-07-29, branch master)";
}

} // namespace ibsimu_cycl
