// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:		 BSD License
//					 license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_CONTACT_UTILITIES)
#define KRATOS_CONTACT_UTILITIES

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "contact_structural_mechanics_application_variables.h"

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

///@}
///@name Kratos Classes
///@{

/**
 * @class ContactUtilities
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This class includes some utilities used for contact computations
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) ContactUtilities
{
public:
    ///@name Type Definitions
    ///@{

    // Some geometrical definitions
    typedef Node<3>                                              NodeType;
    typedef Point::CoordinatesArrayType              CoordinatesArrayType;

    /// Definition of geometries
    typedef Geometry<NodeType>                               GeometryType;

    /// The containers of the components of the model parts
    typedef ModelPart::NodesContainerType                  NodesArrayType;
    typedef ModelPart::ConditionsContainerType        ConditionsArrayType;

    /// Index type definition
    typedef std::size_t                                         IndexType;

    /// Size type definition
    typedef std::size_t                                          SizeType;

    /// Pointer definition of ExactMortarIntegrationUtility
    KRATOS_CLASS_POINTER_DEFINITION(ContactUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    /**
     * @brief This is the default constructor
     */
    ContactUtilities() = default;

    /// Destructor.
    virtual ~ContactUtilities() = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function computes the relative size of the mesh
     * @param rModelPart The modelpart to compute
     */
    static double CalculateRelativeSizeMesh(ModelPart& rModelPart);

    /**
     * @brief This method computes the maximal nodal H
     * @param rModelPart The modelpart to compute
     */
    static double CalculateMaxNodalH(ModelPart& rModelPart);

    /**
     * @brief This method computes the mean nodal H
     * @param rModelPart The modelpart to compute
     */
    static double CalculateMeanNodalH(ModelPart& rModelPart);

    /**
     * @brief This method computes the minimal nodal H
     * @param rModelPart The modelpart to compute
     */
    static double CalculateMinimalNodalH(ModelPart& rModelPart);

    /**
     * @brief This function scales the points according to a factor (to increase the bounding box)
     * @param rPointToScale The point to scale
     * @param rNormal The normal of the point
     * @param LengthSearch The factor considered to "grow" the node
     */
    template<class TPointType>
    static void ScaleNode(
        TPointType& rPointToScale,
        const array_1d<double, 3>& rNormal,
        const double LengthSearch
        )
    {
        noalias(rPointToScale.Coordinates()) = rPointToScale.Coordinates() + rNormal * LengthSearch;
    }

    /**
     * @brief Calculates the distance between nodes
     * @param rPointOrigin The first node
     * @param rPointDestiny The second node
     */
    static double DistancePoints(
        const GeometryType::CoordinatesArrayType& rPointOrigin,
        const GeometryType::CoordinatesArrayType& rPointDestiny
        );

    /**
     * @brief It calculates the center updated in u_n+1 or u_n+1/2
     * @param rModelPart The modelpart to update
     * @param DeltaTime The increment of time considered
     * @param HalfJump If the jumpt is just half dt
     */
    static void ComputeStepJump(
        ModelPart& rModelPart,
        const double DeltaTime,
        const bool HalfJump = true
        );

    /**
     * @brief It checks the activity of the current contact simulation
     * @param rModelPart The modelpart to check the activity
     * @param ThrowError If an error is thrown
     */
    static bool CheckActivity(
        ModelPart& rModelPart,
        const bool ThrowError = true
        );

    /**
     * @brief This method checks that the modelpart has a rotation DoF
     * @details Used in the convergence criteria
     * @param rModelPart The modelpart to clean up
     */
    static bool CheckModelPartHasRotationDoF(ModelPart& rModelPart);

    /**
     * @brief This method removes the model parts with computing conditions
     * @details So for example we can remove potential errors in remeshing processes
     * @param rModelPart The modelpart to clean up
     */
    static void CleanContactModelParts(ModelPart& rModelPart);

    /**
     * @brief It computes the explicit contributions of the conditions
     * @param rModelPart The modelpart to update
     */
    static void ComputeExplicitContributionConditions(ModelPart& rModelPart);

    /**
     * @brief It activates the conditions with active nodes
     * @param rModelPart The modelpart to check
     */
    static void ActivateConditionWithActiveNodes(ModelPart& rModelPart);

    /**
     * @brief It calculates the center updated in u_n+1/2
     * @param rThisGeometry The geometry to calculate
     * @return point: The center in u_n+1/2 (Newmark)
     */
    static array_1d<double, 3> GetHalfJumpCenter(GeometryType& rThisGeometry);

    /**
     * @brief It calculates the matrix containing the tangent vector of the r_gt (for frictional contact)
     * @param rGeometry The geometry to calculate
     * @param StepSlip The considered step slip
     * @return tangent_matrix The matrix containing the tangent vectors of the r_gt
     */
    template< std::size_t TDim, std::size_t TNumNodes>
    static BoundedMatrix<double, TNumNodes, TDim> ComputeTangentMatrixSlip(
        const GeometryType& rGeometry,
        const std::size_t StepSlip = 1
        )
    {
        /* DEFINITIONS */
        // Zero tolerance
        const double zero_tolerance = std::numeric_limits<double>::epsilon();
        // Tangent matrix
        BoundedMatrix<double, TNumNodes, TDim> tangent_matrix;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            const array_1d<double, 3>& r_gt = rGeometry[i_node].FastGetSolutionStepValue(WEIGHTED_SLIP, StepSlip);
            const double norm_slip = norm_2(r_gt);
            if (norm_slip > zero_tolerance) { // Non zero r_gt
                const array_1d<double, 3> tangent_slip = r_gt/norm_slip;
                for (std::size_t i_dof = 0; i_dof < TDim; ++i_dof)
                    tangent_matrix(i_node, i_dof) = tangent_slip[i_dof];
            } else { // We consider the tangent direction as auxiliar
                const array_1d<double, 3>& r_normal = rGeometry[i_node].FastGetSolutionStepValue(NORMAL);
                array_1d<double, 3> tangent_xi, tangent_eta;
                MathUtils<double>::OrthonormalBasis(r_normal, tangent_xi, tangent_eta);
                if (TDim == 3) {
                    for (std::size_t i_dof = 0; i_dof < 3; ++i_dof)
                        tangent_matrix(i_node, i_dof) = tangent_xi[i_dof];
                } else  {
                    if (std::abs(tangent_xi[2]) > std::numeric_limits<double>::epsilon()) {
                        for (std::size_t i_dof = 0; i_dof < 2; ++i_dof)
                            tangent_matrix(i_node, i_dof) = tangent_eta[i_dof];
                    } else {
                        for (std::size_t i_dof = 0; i_dof < 2; ++i_dof)
                            tangent_matrix(i_node, i_dof) = tangent_xi[i_dof];
                    }
                }
            }
        }

        return tangent_matrix;
    }

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
    /**
     * @brief It calculates the matrix of a variable of a geometry
     * @param rNodes The geometry to calculate
     * @param rVarName The name of the variable to calculate
     * @return var_matrix: The matrix containing the variables of the geometry
     */
    static Matrix GetVariableMatrix(
        const GeometryType& rNodes,
        const Variable<array_1d<double,3> >& rVarName
        );

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
};// class ContactUtilities

}
#endif /* KRATOS_CONTACT_UTILITIES defined */
