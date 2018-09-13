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

namespace Kratos
{
namespace HDF5
{

class File;
struct WriteInfo;

namespace Internals
{
///@addtogroup HDF5Application
///@{
///@name Kratos Classes
///@{

/// Represents connectivities information of a single element or condition type in a mesh.
/**
 * Acts as the intermediary between the HDF5 file and the Kratos elements and conditions.
 */
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

    inline unsigned size() const
    {
        return mIds.size();
    }

    /// Read data from a file.
    /**
     * Ensures valid element or condition data is read from the given path on
     * return. Previously stored data is replaced.
     */
    void ReadData(File& rFile, const std::string& rPath, unsigned StartIndex, unsigned BlockSize);

    /// Write data to a file.
    void WriteData(File& rFile, const std::string& rPath, WriteInfo& rInfo);

    // Create and append new elements to the container.
    void CreateEntities(NodesContainerType& rNodes,
                        PropertiesContainerType& rProperties,
                        ElementsContainerType& rElements) const;

    // Create and append new conditions to the container.
    void CreateEntities(NodesContainerType& rNodes,
                        PropertiesContainerType& rProperties,
                        ConditionsContainerType& rConditions) const;

    // Fill data from elements of a single element type.
    /**
     * Expects a uniform, non-empty container of a single element type.
     */
    void SetData(ElementsContainerType const& rElements);

    // Fill data from elements of a single element type.
    /**
     * Expects a registered element name and a uniform container of the
     * corresponding element type. The container may be empty.
     */
    void SetData(const std::string& rName, ElementsContainerType const& rElements);

    // Fill data from conditions of a single condition type.
    /**
     * Expects a uniform, non-empty container of a single condition type.
     */
    void SetData(ConditionsContainerType const& rConditions);

    // Fill data from conditions of a single condition type.
    /**
     * Expects a registered condition name and a uniform container of the
     * corresponding condition type. The container may be empty.
     */
    void SetData(const std::string& rName, ConditionsContainerType const& rConditions);

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