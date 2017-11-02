#include "custom_utilities/hdf5_utils.h"

#include "utilities/openmp_utils.h"

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

void GetRawPointers(ElementsContainerType const& rElementsIn,
                    std::vector<ElementType const*>& rElementsOut)
{
    const unsigned num_elems = rElementsIn.size();
    if (rElementsOut.size() != num_elems)
        rElementsOut.resize(num_elems);

    const int num_threads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector partition;
    OpenMPUtils::DivideInPartitions(num_elems, num_threads, partition);
#pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();
        ElementsContainerType::const_iterator it = rElementsIn.begin() + partition[thread_id];
        for (auto i = partition[thread_id]; i < partition[thread_id + 1]; ++i)
        {
            rElementsOut[i] = &(*it);
            ++it;
        }
    }
}

void GetRawPointers(ConditionsContainerType const& rConditionsIn,
                    std::vector<ConditionType const*>& rConditionsOut)
{
    const unsigned num_conds = rConditionsIn.size();
    if (rConditionsOut.size() != num_conds)
        rConditionsOut.resize(num_conds);

    const int num_threads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector partition;
    OpenMPUtils::DivideInPartitions(num_conds, num_threads, partition);
#pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();
        ConditionsContainerType::const_iterator it = rConditionsIn.begin() + partition[thread_id];
        for (auto i = partition[thread_id]; i < partition[thread_id + 1]; ++i)
        {
            rConditionsOut[i] = &(*it);
            ++it;
        }
    }
}

} // namespace Detail.
} // namespace HDF5.
} // namespace Kratos.
