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

#if !defined(KRATOS_PERTURB_GEOMETRY_SUBGRID_UTILITY)
#define KRATOS_PERTURB_GEOMETRY_SUBGRID_UTILITY

// System includes

// External includes

// Project includes
#include "custom_utilities/perturb_geometry_base_utility.h"
#include "linear_solvers/linear_solver.h"

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
 * @class PerturbGeometrySubgridUtility
 * @ingroup StructuralMechanicsApplication
 * @brief This class generates a random field based on a reduced correlation matrix
 * @details Random field is used to perturb initial geometry
 * @author Manuel Messmer
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) PerturbGeometrySubgridUtility
    : public PerturbGeometryBaseUtility
{
public:

    ///@name Type Definitions
    ///@{

    typedef LinearSolver<TDenseSpaceType, TDenseSpaceType>      LinearSolverType;

    typedef LinearSolverType::Pointer                           LinearSolverPointerType;

    typedef ModelPart::NodesContainerType::ContainerType        ResultNodesContainerType;

    /// Pointer definition of PerturbGeometrySubgridUtility
    KRATOS_CLASS_POINTER_DEFINITION(PerturbGeometrySubgridUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    PerturbGeometrySubgridUtility( ModelPart& rInitialModelPart, LinearSolverPointerType pEigenSolver, Parameters Settings) :
        PerturbGeometryBaseUtility(rInitialModelPart, Settings){
            mpEigenSolver = pEigenSolver;
            mMinDistanceSubgrid = Settings["min_distance_subgrid"].GetDouble();
    }

    /// Destructor.
    ~PerturbGeometrySubgridUtility() override
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

    /**
     * @brief Creates Eigenvectors of correlation matrix in a subgrid
     * @details Finds a subgrid (coarser mesh). Generates correlation matrix in subgrid. Decomposes correlation matrix.
     * @param correlation_matrix Correlation matrix. Stores correlation value for all nodes in the subgrid.
     * @param rPerturbationMatrix Perturbation matrix. Stores eigenvectors of correlation matrix (column-wise).
     */
    int CreateRandomFieldVectors() override;

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
        return "PerturbGeometrySubgridUtility";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PerturbGeometrySubgridUtility";
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

    LinearSolverPointerType mpEigenSolver;
    double mMinDistanceSubgrid;

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
    PerturbGeometrySubgridUtility& operator=(PerturbGeometrySubgridUtility const& rOther) = delete;

    /// Copy constructor.
    PerturbGeometrySubgridUtility(PerturbGeometrySubgridUtility const& rOther) = delete;

    ///@}

}; // Class PerturbGeometrySubgridUtility

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}
#endif /* KRATOS_PERTURB_GEOMETRY_SUBGRID_UTILITY defined */