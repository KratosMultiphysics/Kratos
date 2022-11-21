//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{
///@addtogroup CoSimulationApplication
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class DataTransfer3D1DUtilities
 * @ingroup CoSimulationApplication
 * @brief This utility includes auxiliary methods to transfer from 3D domains to 1D domains and viceversa
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(CO_SIMULATION_APPLICATION) DataTransfer3D1DUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DataTransfer3D1DUtilities
    KRATOS_CLASS_POINTER_DEFINITION(DataTransfer3D1DUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DataTransfer3D1DUtilities() = delete;

    /// Assignment operator.
    DataTransfer3D1DUtilities& operator=(DataTransfer3D1DUtilities const& rOther) = delete;

    /// Copy constructor.
    DataTransfer3D1DUtilities(DataTransfer3D1DUtilities const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method copies from the 3D to the 1D
     * @param rModelPart3D The 3D model part
     * @param rModelPart1D The 1D model part
     */
    static void From3Dto1DDataTransfer(
        ModelPart& rModelPart3D,
        ModelPart& rModelPart1D
        );

    /**
     * @brief This method copies from the 1D to the 3D
     * @param rModelPart3D The 3D model part
     * @param rModelPart1D The 1D model part
     */
    static void From1Dto3DDataTransfer(
        ModelPart& rModelPart3D,
        ModelPart& rModelPart1D
        );

    ///@}
}; // Class DataTransfer3D1DUtilities

///@}

///@} addtogroup block

}  // namespace Kratos.