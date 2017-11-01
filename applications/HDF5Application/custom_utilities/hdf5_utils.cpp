#include "custom_utilities/hdf5_utils.h"

namespace Kratos
{
namespace HDF5
{
namespace Detail
{

bool IsPath(std::string Path)
{
    return regex_match(Path, std::regex("(/\\w+)+"));
}

std::vector<std::string> Split(std::string Path, char Delimiter)
{
    std::vector<std::string> result;
    result.reserve(10);
    std::stringstream ss(Path);
    std::string sub_string;
    while (std::getline(ss, sub_string, Delimiter))
        if (sub_string.size() > 0)
            result.push_back(sub_string);

    return result;
}

} // namespace Detail.
} // namespace HDF5.
} // namespace Kratos.
