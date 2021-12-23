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

#if !defined(KRATOS_PERTUB_GEOMETRY_SPARSE_UTILITY)
#define KRATOS_PERTUB_GEOMETRY_SPARSE_UTILITY

// System includes

// External includes

// Project includes
#include "custom_utilities/perturb_geometry_base_utility.h"
#include "linear_solvers/linear_solver.h"

namespace Kratos {
///@addtogroup StructuralMechanicsApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class PerturbGeometrySparseUtility
 * @ingroup StructuralMechanicsApplication
 * @brief This class generates a random field based on a sparse correlation matrix.
 * @details Random field is used to perturb initial geometry.
 * @author Manuel Messmer
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) PerturbGeometrySparseUtility
    : public PerturbGeometryBaseUtility
{
public:

    ///@name Type Definitions
    ///@{

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType> LinearSolverType;

    typedef LinearSolverType::Pointer                       LinearSolverPointerType;

    typedef ModelPart::NodesContainerType::ContainerType    ResultNodesContainerType;

    typedef TSparseSpaceType::MatrixType                    SparseMatrixType;

    /// Pointer definition of PerturbGeometrySparseUtility
    KRATOS_CLASS_POINTER_DEFINITION(PerturbGeometrySparseUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    PerturbGeometrySparseUtility( ModelPart& rInitialModelPart, LinearSolverPointerType pEigenSolver, Parameters Settings) :
        PerturbGeometryBaseUtility(rInitialModelPart, Settings){
            mpEigenSolver = pEigenSolver;
    }

    /// Destructor.
    ~PerturbGeometrySparseUtility() override
    = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates Eigenvectors of a sparse correlation matrix
     * @details Generates sparse correlation matrix. Decomposes correlation matrix.
     * @param correlation_matrix Correlation matrix. Stores correlation value for all nodes.
     * @param rPerturbationMatrix Perturbation matrix. Stores eigenvectors of correlation matrix (colum-wise).
     */
    int CreateRandomFieldVectors() override;

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "PerturbGeometrySparseUtility";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PerturbGeometrySparseUtility";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    LinearSolverPointerType mpEigenSolver;

    /// Assignment operator.
    PerturbGeometrySparseUtility& operator=(PerturbGeometrySparseUtility const& rOther) = delete;

    /// Copy constructor.
    PerturbGeometrySparseUtility(PerturbGeometrySparseUtility const& rOther) = delete;

    ///@}

    }; // Class PerturbGeometrySparseUtility

///@}

///@} addtogroup block
}
#endif /* KRATOS_PERTUB_GEOMETRY_SPARSE_UTILITY defined */