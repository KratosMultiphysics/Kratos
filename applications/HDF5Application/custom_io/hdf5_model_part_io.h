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

#if !defined(KRATOS_HDF5_MODEL_PART_IO_H_INCLUDED)
#define KRATOS_HDF5_MODEL_PART_IO_H_INCLUDED

// System includes
#include <iostream>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/io.h"

// Application includes
#include "custom_io/hdf5_file.h"

namespace Kratos
{
///@addtogroup HDF5Application
///@{

class HDF5ModelPartIO : public IO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(HDF5ModelPartIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    HDF5ModelPartIO(Parameters& rParams, HDF5File::Pointer pFile);

    ///@}
    ///@name Operations
    ///@{

    std::size_t ReadNodesNumber() override;

    void WriteNodes(NodesContainerType const& rThisNodes) override;

    void ReadElements(NodesContainerType& rThisNodes,
                      PropertiesContainerType& rThisProperties,
                      ElementsContainerType& rThisElements) override;

    std::size_t ReadElementsConnectivities(ConnectivitiesContainerType& rElementsConnectivities) override;

    void WriteElements(ElementsContainerType const& rThisElements) override;

    void ReadConditions(NodesContainerType& rThisNodes,
                        PropertiesContainerType& rThisProperties,
                        ConditionsContainerType& rThisConditions) override;

    std::size_t ReadConditionsConnectivities(ConnectivitiesContainerType& rConditionsConnectivities) override;

    void ReadInitialValues(ModelPart& rThisModelPart) override;

    void ReadInitialValues(NodesContainerType& rThisNodes,
                           ElementsContainerType& rThisElements,
                           ConditionsContainerType& rThisConditions) override;

    void ReadModelPart(ModelPart& rThisModelPart) override;

    void WriteModelPart(ModelPart& rThisModelPart) override;

    ///@}

private:
    ///@name Member Variables
    ///@{
    HDF5File::Pointer mpFile;
    ///@}

    ///@name Private Operations
    ///@{

    ///@}
};

///@} addtogroup
} // namespace Kratos.

#endif // KRATOS_HDF5_MODEL_PART_IO_H_INCLUDED defined
