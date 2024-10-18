// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Armin Geiser
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

namespace Kratos {

///@name Kratos Globals
///@{

///@name Kratos Classes
///@{

/** 
 * @class DistributeLoadOnSurfaceProcess
 * @brief This process distributes a load on surface load conditions belonging to a modelpart
 * @details The load is distributed according to the surface area.
 * @author Armin Geiser
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) DistributeLoadOnSurfaceProcess 
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of DistributeLoadOnSurfaceProcess
    KRATOS_CLASS_POINTER_DEFINITION(DistributeLoadOnSurfaceProcess);

    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    DistributeLoadOnSurfaceProcess(ModelPart& rModelPart,
                                  Parameters Parameters);

    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitializeSolutionStep() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override {
        return "DistributeLoadOnSurfaceProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {
        rOStream << "DistributeLoadOnSurfaceProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart; /// The model part
    Parameters mParameters; /// The parameters

    ///@}
    ///@name Private Operations
    ///@{

    ///@}

}; // Class DistributeLoadOnSurfaceProcess

///@}

}  // namespace Kratos.