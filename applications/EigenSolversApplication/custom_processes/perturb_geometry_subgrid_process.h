// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Manuel Messmer
//

#if !defined(KRATOS_PERTURB_GEOMETRY_SUBGRID_PROCESS)
#define KRATOS_PERTURB_GEOMETRY_SUBGRID_PROCESS

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_processes/perturb_geometry_base_process.h"
#include "includes/model_part.h"
#include "custom_utilities/omp_node_search.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

/**
 * @class PerturbGeometrySubgridProcess
 *
 * @ingroup EigenSolversApplication
 *
 * @brief This class generates a random field based on a reduced correlation matrix
 * @details Random field is used to perturb initial geometry
 *
 * @author Manuel Messmer
 */
class KRATOS_API(EIGEN_SOLVERS_APPLICATION) PerturbGeometrySubgridProcess
    : public PerturbGeometryBaseProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PerturbGeometrySubgridProcess
    KRATOS_CLASS_POINTER_DEFINITION(PerturbGeometrySubgridProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PerturbGeometrySubgridProcess( ModelPart& rInitialModelPart, double MaximalDisplacement) :
        PerturbGeometryBaseProcess(rInitialModelPart, MaximalDisplacement){
    }

    /// Destructor.
    ~PerturbGeometrySubgridProcess() override
    = default;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    int CreateEigenvectors(ModelPart& rThisModelPart, double minDistance, double correlationLength, double truncationTolerance);

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "PerturbGeometrySubgridProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PerturbGeometrySubgridProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    PerturbGeometrySubgridProcess& operator=(PerturbGeometrySubgridProcess const& rOther) = delete;

    /// Copy constructor.
    PerturbGeometrySubgridProcess(PerturbGeometrySubgridProcess const& rOther) = delete;


    ///@}

}; // Class PerturbGeometrySubgridProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}
#endif /* KRATOS_PERTURB_GEOMETRY_SUBGRID_PROCESS defined */


