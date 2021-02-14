//          __        ____ ____
// ___  ___/  |_ __ _/_   /_   |
// \  \/ /\   __\  |  \   ||   |
//  \   /  |  | |  |  /   ||   |
//   \_/   |__| |____/|___||___|
//
//  License: BSD License ; see LICENSE
//

#ifndef VTU11_ALIAS_HPP
#define VTU11_ALIAS_HPP

#include <string>
#include <map>
#include <utility>
#include <vector>

namespace vtu11
{

using StringStringMap = std::map<std::string, std::string>;

using DataSet = std::tuple<std::string, size_t, std::vector<double>>;

using VtkCellType = std::int8_t;

using HeaderType = size_t;

using Byte = unsigned char;

using VtkIndexType = std::int64_t;

} // namespace vtu11

#endif // VTU11_ALIAS_HPP
