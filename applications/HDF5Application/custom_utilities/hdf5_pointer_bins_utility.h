//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_HDF5_POINTER_BINS_UTILITY_H_INCLUDED)
#define KRATOS_HDF5_POINTER_BINS_UTILITY_H_INCLUDED

// System includes
#include <vector>
#include <typeinfo>

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/openmp_utils.h"

// Application includes

namespace Kratos
{
namespace HDF5
{
namespace Detail
{
///@addtogroup HDF5Application
///@{

// A class for sorting pointers into bins according to their dynamic types.
template <class DataType>
class PointerBinsUtility
{
public:
    ///@name Type Definitions
    ///@{
        typedef DataType* PointerType;
        typedef DataType const* ConstPointerType;
        typedef std::vector<ConstPointerType> BinType;
        typedef ConstPointerType KeyType;
    ///@}
    ///@name Life Cycle
    ///@{
    explicit PointerBinsUtility(const std::vector<KeyType>& rBinKeys) : mBinKeys{rBinKeys}
    {
        mBins.resize(mBinKeys.size());
    }
    ///@}
    ///@name Operations
    ///@{
    // Generate bins from a container of mixed dynamic types.
    template <class TOtherContainerType>
    void CreateBins(TOtherContainerType& rData);

    // Get a bin by its key.
    BinType& GetBin(KeyType key)
    {
        KRATOS_TRY;

        unsigned i;
        bool found = false;
        for (i = 0; i < mBinKeys.size(); ++i)
            if (typeid(mBinKeys[i]) == typeid(key))
            {
                found = true;
                break;
            }
        if (!found)
            KRATOS_ERROR << "No bin exists for key." << std::endl;

        return mBins[i];

        KRATOS_CATCH("");
    }

    void Clear()
    {
        for (auto& r_bin : mBins)
            r_bin.clear();
    }
    ///@}
private:
    ///@name Member Variables
    ///@{
    const std::vector<KeyType> mBinKeys;
    std::vector<BinType> mBins;
    ///@}
}; // class PointerBinsUtility

template <class DataType>
template <class TOtherContainerType>
void PointerBinsUtility<DataType>::CreateBins(TOtherContainerType& rData)
{
    KRATOS_TRY;

    for (auto& r_bin : mBins)
        r_bin.reserve(rData.size());

    const int num_threads = OpenMPUtils::GetNumThreads();
#ifdef _OPENMP
    std::vector<omp_lock_t> bin_locks(mBins.size());
    for (auto& lock : bin_locks)
        omp_init_lock(&lock);
#endif
    OpenMPUtils::PartitionVector partition;
    OpenMPUtils::DivideInPartitions(rData.size(), num_threads, partition);
#pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();
        typename TOtherContainerType::const_iterator it = rData.begin() + partition[thread_id];
        for (auto i = partition[thread_id]; i < partition[thread_id + 1]; ++i)
        {
            bool found = false;
            ConstPointerType p_item = &(*it);
            for (unsigned i = 0; i < mBinKeys.size(); ++i)
            {
                if (typeid(p_item) == typeid(mBinKeys[i]))
                {
#ifdef _OPENMP
                    omp_set_lock(&bin_locks[i]);
#endif
                    mBins[i].push_back(p_item);
#ifdef _OPENMP
                    omp_unset_lock(&bin_locks[i]);
#endif
                    found = true;
                    break;
                }
            }
            KRATOS_ERROR_IF(!found) << "Did not find bin for element #"
                                    << p_item->Id() << std::endl;
            ++it;
        }
    }

    KRATOS_CATCH("");
}

///@} addtogroup
} // namespace Detail.
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_POINTER_BINS_UTILITY_H_INCLUDED defined