//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Armin Geiser
//

#if !defined(KRATOS_DISTRIBUTE_LOAD_ON_SURFACE_PROCESS_H_INCLUDED )
#define  KRATOS_DISTRIBUTE_LOAD_ON_SURFACE_PROCESS_H_INCLUDED

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

/// Process to create the animated Eigenvectors
/** This process distributes a load on surface load conditions belonging to a modelpart.
 *  The load is distributed according to the surface area.
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) DistributeLoadOnSurfaceProcess : public Process
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

    ModelPart& mrModelPart;
    Parameters mParameters;

    ///@}
    ///@name Private Operations
    ///@{

    ///@}

}; // Class DistributeLoadOnSurfaceProcess

///@}

}  // namespace Kratos.

#endif // KRATOS_DISTRIBUTE_LOAD_ON_SURFACE_PROCESS_H_INCLUDED  defined
