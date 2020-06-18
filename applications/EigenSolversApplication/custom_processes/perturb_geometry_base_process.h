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

#if !defined(KRATOS_PERTURB_GEOMETRY_PROCESS)
#define KRATOS_PERTURB_GEOMETRY_PROCESS

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "linear_solvers/direct_solver.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "custom_solvers/eigensystem_solver.h"
#include "custom_solvers/eigen_direct_solver.h"
#include "custom_utilities/omp_node_search.h"
#include "utilities/mortar_utilities.h"


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
 * @class PerturbGeometryBaseProcess
 *
 * @ingroup StructuralMechanicsApplication
 *
 * @brief Base class for geometry perturbation process
 * @details It takes into account all elements in the ModelPart
 *
 * @author Manuel Messmer
 */
class KRATOS_API(EIGEN_SOLVERS_APPLICATION) PerturbGeometryBaseProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PerturbGeometryBaseProcess
    KRATOS_CLASS_POINTER_DEFINITION(PerturbGeometryBaseProcess);

    typedef TUblasSparseSpace<double> TSparseSpaceType;
    typedef TUblasDenseSpace<double> TDenseSpaceType;

    typedef typename TSparseSpaceType::VectorPointerType SparseVectorPointerType;

    typedef typename TSparseSpaceType::MatrixPointerType SparseMatrixPointerType;

    typedef typename TDenseSpaceType::VectorType DenseVectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef typename TDenseSpaceType::MatrixPointerType DenseMatrixPointerType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PerturbGeometryBaseProcess( ModelPart& rInitialModelPart, double MaximalDisplacement) :
        mrInitialModelPart(rInitialModelPart),
        mMaximalDisplacement(MaximalDisplacement)
    {
        KRATOS_TRY
        MortarUtilities::ComputeNodesMeanNormalModelPart( mrInitialModelPart, false );
        mpPerturbationMatrix = TDenseSpaceType::CreateEmptyMatrixPointer();
        KRATOS_CATCH("")
    }

    /// Destructor.
    ~PerturbGeometryBaseProcess() override
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

    void AssembleEigenvectors(ModelPart& rThisModelPart, const std::vector<double>& variables );

    double CorrelationFunction( ModelPart::NodeIterator itNode1, ModelPart::NodeIterator itNode2, double CorrelationLenth);

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
        return "PerturbGeometryBaseProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PerturbGeometryBaseProcess";
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

    DenseMatrixPointerType mpPerturbationMatrix;

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

    ModelPart& mrInitialModelPart;

    double mMaximalDisplacement;

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
    PerturbGeometryBaseProcess& operator=(PerturbGeometryBaseProcess const& rOther) = delete;

    /// Copy constructor.
    PerturbGeometryBaseProcess(PerturbGeometryBaseProcess const& rOther) = delete;


    ///@}

}; // Class PerturbGeometryBaseProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}
#endif /* KRATOS_PERTURBE_GEOMETRY_PROCESS defined */


