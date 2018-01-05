//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
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
#include "hdf5_application_define.h"

namespace Kratos
{
namespace HDF5
{
namespace Internals
{
///@addtogroup HDF5Application
///@{
///@name Kratos Classes
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
            if (IsMatch(key, mBinKeys[i]))
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
    ///@name Private Operations
    ///@{
    template <class TOtherContainerType>
    void CreateSingleBin(TOtherContainerType& rData);

    template <class TOtherContainerType>
    void CreateMultipleBins(TOtherContainerType& rData);

    inline bool IsMatch(KeyType lhs, KeyType rhs);
    ///@}
}; // class PointerBinsUtility
///@} // Kratos Classes

template <class DataType>
template <class TOtherContainerType>
void PointerBinsUtility<DataType>::CreateBins(TOtherContainerType& rData)
{
    KRATOS_TRY;

    if (mBinKeys.size() == 1)
        CreateSingleBin(rData);
    else
        CreateMultipleBins(rData);

    KRATOS_CATCH("");
}

template <class DataType>
template <class TOtherContainerType>
void PointerBinsUtility<DataType>::CreateSingleBin(TOtherContainerType& rData)
{
    KRATOS_TRY;

    mBins[0].resize(rData.size());

    const int num_threads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector partition;
    OpenMPUtils::DivideInPartitions(rData.size(), num_threads, partition);
#pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();
        typename TOtherContainerType::const_iterator it = rData.begin() + partition[thread_id];
        for (auto i = partition[thread_id]; i < partition[thread_id + 1]; ++i)
        {
            ConstPointerType p_item = &(*it);
            if (IsMatch(p_item, mBinKeys[0]))
            mBins[0][i] = p_item;
            else
                KRATOS_ERROR << "Did not find bin for element #" << p_item->Id() << std::endl;
            ++it;
        }
    }

    KRATOS_CATCH("");
}

template <class DataType>
template <class TOtherContainerType>
void PointerBinsUtility<DataType>::CreateMultipleBins(TOtherContainerType& rData)
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
            for (unsigned j = 0; j < mBinKeys.size(); ++j)
            {
                if (IsMatch(p_item, mBinKeys[j]))
                {
#ifdef _OPENMP
                    omp_set_lock(&bin_locks[j]);
#endif
                    mBins[j].push_back(p_item);
#ifdef _OPENMP
                    omp_unset_lock(&bin_locks[j]);
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

template <class DataType>
bool PointerBinsUtility<DataType>::IsMatch(KeyType lhs, KeyType rhs)
{
    return (typeid(lhs) == typeid(rhs) &&
            lhs->GetGeometry().size() == rhs->GetGeometry().size() &&
            lhs->GetGeometry().WorkingSpaceDimension() == rhs->GetGeometry().WorkingSpaceDimension());
}

///@} addtogroup
} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_POINTER_BINS_UTILITY_H_INCLUDED defined