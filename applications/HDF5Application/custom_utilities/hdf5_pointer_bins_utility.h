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
    for (auto& r_bin : mBins)
        r_bin.reserve(rData.size());

    for (auto& r_item : rData)
    {
        ConstPointerType p_item = &r_item;
        for (unsigned i = 0; i < mBinKeys.size(); ++i)
            if (typeid(p_item) == typeid(mBinKeys[i]))
            {
                mBins[i].push_back(p_item);
                break;
            }
    }
}

///@} addtogroup
} // namespace Detail.
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_POINTER_BINS_UTILITY_H_INCLUDED defined