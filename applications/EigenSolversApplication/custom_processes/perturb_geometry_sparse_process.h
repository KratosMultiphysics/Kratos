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

#if !defined(KRATOS_PERTUB_GEOMETRY_SPARSE_PROCESS)
#define KRATOS_PERTUB_GEOMETRY_SPARSE_PROCESS

// System includes

// External includes

// Project includes
#include "custom_solvers/eigensystem_solver.h"
#include "custom_processes/perturb_geometry_base_process.h"
#include "includes/define.h"
#include "includes/model_part.h"

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
 * @class PerturbGeometrySparseProcess
 *
 * @ingroup EigenSolversApplication
 *
 * @brief This class generates a random field based on a sparse correlation matrix
 * @details Random field is used to perturb initial geometry
 *
 * @author Manuel Messmer
 */
class KRATOS_API(EIGEN_SOLVERS_APPLICATION) PerturbGeometrySparseProcess
    : public PerturbGeometryBaseProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PerturbGeometrySparseProcess
    KRATOS_CLASS_POINTER_DEFINITION(PerturbGeometrySparseProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PerturbGeometrySparseProcess( ModelPart& rInitialModelPart, double MaximalDisplacement) :
        PerturbGeometryBaseProcess(rInitialModelPart, MaximalDisplacement){
    }

    /// Destructor.
    ~PerturbGeometrySparseProcess() override
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

    int CreateEigenvectors(ModelPart& rThisModelPart, double correlationLength, double truncationTolerance, Parameters eigenvalueSettings );

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
        return "PerturbGeometrySparseProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PerturbGeometrySparseProcess";
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
    PerturbGeometrySparseProcess& operator=(PerturbGeometrySparseProcess const& rOther) = delete;

    /// Copy constructor.
    PerturbGeometrySparseProcess(PerturbGeometrySparseProcess const& rOther) = delete;


    ///@}

}; // Class PerturbGeometrySparseProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}
#endif /* KRATOS_PERTUBE_GEOMETRY_SPARSE_PROCESS defined */


