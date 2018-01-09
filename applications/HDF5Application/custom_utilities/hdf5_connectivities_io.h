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

#if !defined(KRATOS_HDF5_CONNECTIVITIES_IO_H_INCLUDED)
#define KRATOS_HDF5_CONNECTIVITIES_IO_H_INCLUDED

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "hdf5_application_define.h"
#include "custom_io/hdf5_file.h"
#include "custom_utilities/hdf5_connectivities_data.h"
#include "custom_utilities/hdf5_pointer_bins_utility.h"

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
struct FileItem
{
    const std::string Name;
    const std::string Path;
};

template <class TConnectivitiesType>
class ConnectivitiesInput
{
    struct InputItem
    {
        const std::string& Name;
        const std::string& Path;

        void ReadConnectivities(File& rFile,
                                NodesContainerType& rNodes,
                                PropertiesContainerType& rProperties,
                                unsigned StartIndex,
                                unsigned BlockSize,
                                PointerVectorSet<TConnectivitiesType, IndexedObject>& rConnectivities)
        {
            const TConnectivitiesType& r_creator =
                KratosComponents<TConnectivitiesType>::Get(Name);
            Internals::ConnectivitiesData connectivities;
            connectivities.ReadData(rFile, Path, StartIndex, BlockSize);
            connectivities.CreateEntities(r_creator, rNodes, rProperties, rConnectivities);
        }
    };

public:
    ///@name Type Definitions
    ///@{

    typedef typename std::vector<InputItem>::iterator iterator;

    ///@}
    ///@name Life Cycle
    ///@{

    ConnectivitiesInput(std::vector<FileItem>& rInputs)
    {
        for (const auto& r_item : rInputs)
        {
            mInputs.push_back({r_item.Name, r_item.Path});
        }
    }

    ///@}
    ///@name Operations
    ///@{

    inline iterator begin()
    {
        return mInputs.begin();
    }

    inline iterator end()
    {
        return mInputs.end();
    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    std::vector<InputItem> mInputs;

    ///@}
};

template <class TConnectivitiesType>
class ConnectivitiesOutput
{
    struct OutputItem
    {
        const std::string& Name;
        const std::string& Path;
        Internals::ConnectivitiesData Connectivities;

        void WriteConnectivities(File& rFile)
        {
            Connectivities.WriteData(rFile, Path);
        }
    };

public:
    ///@name Type Definitions
    ///@{

    typedef typename std::vector<OutputItem>::iterator iterator;

    ///@}
    ///@name Life Cycle
    ///@{

    ConnectivitiesOutput(std::vector<FileItem> const& rOutputItems, PointerVectorSet<TConnectivitiesType, IndexedObject> const& rConnectivities)
    {
        std::vector<const TConnectivitiesType*> bin_keys;
        for (const auto& r_item : rOutputItems)
        {
            mOutputs.push_back({r_item.Name, r_item.Path});
            bin_keys.push_back(&KratosComponents<TConnectivitiesType>::Get(r_item.Name));
        }

        Internals::PointerBinsUtility<TConnectivitiesType> bins(bin_keys);
        bins.CreateBins(rConnectivities);
        for (unsigned i = 0; i < mOutputs.size(); ++i)
            mOutputs[i].Connectivities.SetData(bins.GetBin(bin_keys[i]));
    }

    ///@}
    ///@name Operations
    ///@{

    inline iterator begin()
    {
        return mOutputs.begin();
    }

    inline iterator end()
    {
        return mOutputs.end();
    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    std::vector<OutputItem> mOutputs;

    ///@}
};

///@} // Kratos Classes
///@} addtogroup
} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_CONNECTIVITIES_IO_H_INCLUDED defined