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

#if !defined(KRATOS_HDF5_CONNECTIVITIES_DATA_H_INCLUDED)
#define KRATOS_HDF5_CONNECTIVITIES_DATA_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "includes/model_part.h"

// Application includes
#include "custom_io/hdf5_file.h"
#include "custom_utilities/hdf5_utils.h"

namespace Kratos
{
namespace HDF5
{
///@addtogroup HDF5Application
///@{

class ConnectivitiesData
{
public:
    ///@name Type Definitions
    ///@{
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ConnectivitiesData);

    typedef ModelPart::NodesContainerType NodesContainerType;

    typedef ModelPart::ElementsContainerType ElementsContainerType;

    typedef ModelPart::ConditionsContainerType ConditionsContainerType;

    typedef ModelPart::PropertiesContainerType PropertiesContainerType;

    typedef ModelPart::NodeType NodeType;

    typedef ModelPart::ElementType ElementType;

    typedef ModelPart::ConditionType ConditionType;

    typedef ModelPart::PropertiesType PropertiesType;

    ///@}
    ///@name Life Cycle
    ///@{
    ///@}
    ///@name Operations
    ///@{
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

    void ReadData(File& rFile, std::string Path, unsigned StartIndex, unsigned BlockSize);

    void WriteData(File& rFile, std::string Path);

    void CreateElements(ElementType const& rElementType,
                        NodesContainerType& rNodes,
                        PropertiesContainerType& rProperties,
                        ElementsContainerType& rElements);

    void CreateConditions(ConditionType const& rConditionType,
                          NodesContainerType& rNodes,
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
    Vector<int> mIds;
    Vector<int> mPropertiesIds;
    Matrix<int> mConnectivities;
    ///@}
}; // class ConnectivitiesData

///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_CONNECTIVITIES_DATA_H_INCLUDED defined