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

#if !defined(KRATOS_HDF5_CONNECTIVITIES_DATA_H_INCLUDED)
#define KRATOS_HDF5_CONNECTIVITIES_DATA_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "hdf5_application_define.h"
#include "custom_io/hdf5_file.h"

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

/// A class representing connectivities in a mesh.
class ConnectivitiesData
{
public:
    ///@name Type Definitions
    ///@{
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ConnectivitiesData);

    ///@}
    ///@name Life Cycle
    ///@{
    ///@}
    ///@name Operations
    ///@{
    inline std::string const& Name() const
    {
        return mName;
    }

    inline Vector<int> const& GetIds() const
    {
        return mIds;
    }

    inline Vector<int> const& GetPropertiesIds() const
    {
        return mPropertiesIds;
    }

    inline Matrix<int> const& GetConnectivities() const
    {
        return mConnectivities;
    }

    inline unsigned size() const
    {
        return mIds.size();
    }

    void ReadData(File& rFile, const std::string& rPath, unsigned StartIndex, unsigned BlockSize);

    void WriteData(File& rFile, const std::string& rPath, WriteInfo& rInfo);

    void CreateEntities(NodesContainerType& rNodes,
                        PropertiesContainerType& rProperties,
                        ElementsContainerType& rElements);

    void CreateEntities(NodesContainerType& rNodes,
                        PropertiesContainerType& rProperties,
                        ConditionsContainerType& rConditions);

    // Fill data from elements of a single element type.
    void SetData(ElementsContainerType const& rElements);

    // Fill data from conditions of a single condition type.
    void SetData(ConditionsContainerType const& rConditions);

    void Clear();
    ///@}
private:
    ///@name Member Variables
    ///@{
    std::string mName;
    Vector<int> mIds;
    Vector<int> mPropertiesIds;
    Matrix<int> mConnectivities;
    ///@}
}; // class ConnectivitiesData

///@} // Kratos Classes
///@} addtogroup
} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_CONNECTIVITIES_DATA_H_INCLUDED defined