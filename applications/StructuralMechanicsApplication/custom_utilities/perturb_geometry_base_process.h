// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Manuel Messmer
//

#if !defined(KRATOS_PERTURB_GEOMETRY_BASE_PROCESS)
#define KRATOS_PERTURB_GEOMETRY_BASE_PROCESS

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "processes/process.h"
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
 * @ingroup StructuralMechanicsApplication
 * @brief Base class for geometry perturbation process
 * @author Manuel Messmer
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) PerturbGeometryBaseProcess
    : public Process
{
public:

    ///@name Type Definitions
    ///@{

    typedef TUblasSparseSpace<double> TSparseSpaceType;
    typedef TUblasDenseSpace<double> TDenseSpaceType;

    typedef typename TDenseSpaceType::MatrixPointerType DenseMatrixPointerType;

    typedef typename TDenseSpaceType::VectorType DenseVectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    /// Pointer definition of PerturbGeometryBaseProcess
    KRATOS_CLASS_POINTER_DEFINITION(PerturbGeometryBaseProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    PerturbGeometryBaseProcess( ModelPart& rInitialModelPart, Parameters Settings) :
        mrInitialModelPart(rInitialModelPart),
        mCorrelationLength(Settings["correlation_length"].GetDouble()),
        mTruncationError(Settings["truncation_error"].GetDouble()),
        mEchoLevel(Settings["echo_level"].GetInt()),
        mMaximalDisplacement(Settings["max_displacement"].GetDouble())
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

    virtual int CreateRandomFieldVectors() = 0;

    /**
     * @brief Assemble random field and apply to initial geometry
     * @param rPerturbationMatrix Perturbation matrix. Stores eigenvectors of correlation matrix (colum-wise).
     * @param random_field Random field vector. Stores nodal deviations.
     */
    void ApplyRandomFieldVectorsToGeometry(ModelPart& rThisModelPart, const std::vector<double>& variables );

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

    ModelPart& mrInitialModelPart;

    double mCorrelationLength;

    double mTruncationError;

    int mEchoLevel;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Correlation function
     * @return Correlation value of two nodes
     */
    double CorrelationFunction( ModelPart::NodeIterator itNode1, ModelPart::NodeIterator itNode2, double CorrelationLenth);

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
#endif /* KRATOS_PERTURB_GEOMETRY_BASE_PROCESS defined */