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

#if !defined(KRATOS_PERTURB_GEOMETRY_BASE_UTILITY)
#define KRATOS_PERTURB_GEOMETRY_BASE_UTILITY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
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
 * @class PerturbGeometryBaseUtility
 * @ingroup StructuralMechanicsApplication
 * @brief Base class for geometry perturbation utility
 * @author Manuel Messmer
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) PerturbGeometryBaseUtility
{
public:

    ///@name Type Definitions
    ///@{

    typedef TUblasSparseSpace<double> TSparseSpaceType;
    typedef TUblasDenseSpace<double> TDenseSpaceType;

    typedef typename TDenseSpaceType::MatrixPointerType DenseMatrixPointerType;

    typedef typename TDenseSpaceType::VectorType DenseVectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    /// Pointer definition of PerturbGeometryBaseUtility
    KRATOS_CLASS_POINTER_DEFINITION(PerturbGeometryBaseUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    PerturbGeometryBaseUtility( ModelPart& rInitialModelPart, Parameters Settings) :
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
    virtual ~PerturbGeometryBaseUtility() {}

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
    virtual std::string Info() const
    {
        return "PerturbGeometryBaseUtility";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "PerturbGeometryBaseUtility";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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
    PerturbGeometryBaseUtility& operator=(PerturbGeometryBaseUtility const& rOther) = delete;

    /// Copy constructor.
    PerturbGeometryBaseUtility(PerturbGeometryBaseUtility const& rOther) = delete;


    ///@}

}; // Class PerturbGeometryBaseUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}
#endif /* KRATOS_PERTURB_GEOMETRY_BASE_UTILITY defined */