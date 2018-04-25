// KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_ALM_FRICTIONAL_MORTAR_CONTACT_CONDITION_H_INCLUDED )
#define  KRATOS_ALM_FRICTIONAL_MORTAR_CONTACT_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/ALM_mortar_contact_condition.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    typedef Point                                     PointType;
    typedef Node<3>                                    NodeType;
    typedef Geometry<NodeType>                     GeometryType;
    typedef Geometry<PointType>               GeometryPointType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod   IntegrationMethod;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{ *

/**
 * @class AugmentedLagrangianMethodFrictionalMortarContactCondition
 * @ingroup ContactStructuralMechanicsApplication
 * @brief AugmentedLagrangianMethodFrictionalMortarContactCondition
 * @details This is a contact condition which employes the mortar method with dual lagrange multiplier
 * The method has been taken from the Alexander Popps thesis:
 * Popp, Alexander: Mortar Methods for Computational Contact Mechanics and General Interface Problems, Technische Universität München, jul 2012
 * @author Vicente Mataix Ferrandiz
 */
template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) AugmentedLagrangianMethodFrictionalMortarContactCondition
    : public AugmentedLagrangianMethodMortarContactCondition<TDim, TNumNodes, FrictionalCase::FRICTIONAL, TNormalVariation>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of AugmentedLagrangianMethodFrictionalMortarContactCondition
    KRATOS_CLASS_POINTER_DEFINITION( AugmentedLagrangianMethodFrictionalMortarContactCondition );

    typedef AugmentedLagrangianMethodMortarContactCondition<TDim, TNumNodes, FrictionalCase::FRICTIONAL, TNormalVariation> BaseType;

    typedef Condition                                                                                             ConditionBaseType;

    typedef PairedCondition                                                                                 PairedConditionBaseType;

    typedef typename BaseType::MortarConditionMatrices                                                      MortarConditionMatrices;

    typedef typename BaseType::GeneralVariables                                                                    GeneralVariables;

    typedef typename BaseType::IntegrationUtility                                                                IntegrationUtility;

    typedef typename BaseType::DerivativesUtilitiesType                                                    DerivativesUtilitiesType;

    typedef typename BaseType::BelongType                                                                                BelongType;

    typedef typename BaseType::ConditionArrayListType                                                        ConditionArrayListType;

    typedef MortarOperator<TNumNodes>                                                                   MortarBaseConditionMatrices;

    typedef typename ConditionBaseType::VectorType                                                                       VectorType;

    typedef typename ConditionBaseType::MatrixType                                                                       MatrixType;

    typedef typename ConditionBaseType::IndexType                                                                         IndexType;

    typedef typename ConditionBaseType::GeometryType::Pointer                                                   GeometryPointerType;

    typedef typename ConditionBaseType::NodesArrayType                                                               NodesArrayType;

    typedef typename ConditionBaseType::PropertiesType                                                               PropertiesType;

    typedef typename ConditionBaseType::PropertiesType::Pointer                                               PropertiesPointerType;

    typedef typename ConditionBaseType::EquationIdVectorType                                                   EquationIdVectorType;

    typedef typename ConditionBaseType::DofsVectorType                                                               DofsVectorType;

    typedef Line2D2<Point>                                                                                                 LineType;

    typedef Triangle3D3<Point>                                                                                         TriangleType;

    typedef typename std::conditional<TDim == 2, LineType, TriangleType >::type                                   DecompositionType;

    typedef DerivativeDataFrictional<TDim, TNumNodes, TNormalVariation>                                          DerivativeDataType;

    static constexpr IndexType MatrixSize = TDim * (TNumNodes + TNumNodes + TNumNodes);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    AugmentedLagrangianMethodFrictionalMortarContactCondition()
        : BaseType()
    {
    }

    // Constructor 1
    AugmentedLagrangianMethodFrictionalMortarContactCondition(
        IndexType NewId,
        GeometryPointerType pGeometry
        ):BaseType(NewId, pGeometry)
    {
    }

    // Constructor 2
    AugmentedLagrangianMethodFrictionalMortarContactCondition(
        IndexType NewId,
        GeometryPointerType pGeometry,
        PropertiesPointerType pProperties
        ):BaseType( NewId, pGeometry, pProperties )
    {
    }

    // Constructor 3
    AugmentedLagrangianMethodFrictionalMortarContactCondition(
        IndexType NewId,
        GeometryPointerType pGeometry,
        PropertiesPointerType pProperties,
        GeometryType::Pointer pMasterGeometry
        ):BaseType( NewId, pGeometry, pProperties, pMasterGeometry )
    {
    }

    ///Copy constructor
    AugmentedLagrangianMethodFrictionalMortarContactCondition( AugmentedLagrangianMethodFrictionalMortarContactCondition const& rOther)
    {
    }

    /// Destructor.
    ~AugmentedLagrangianMethodFrictionalMortarContactCondition() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
    * @brief Called at the beginning of each solution step
    */
    void Initialize() override;

    /**
    * @brief Called at the ending of each solution step
    * @param rCurrentProcessInfo the current process info instance
    */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Creates a new element pointer from an arry of nodes
     * @param NewId the ID of the new element
     * @param rThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& rThisNodes,
        PropertiesPointerType pProperties
        ) const override;

    /**
     * @brief Creates a new element pointer from an existing geometry
     * @param NewId the ID of the new element
     * @param pGeom the  geometry taken to create the condition
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */

    Condition::Pointer Create(
        IndexType NewId,
        GeometryPointerType pGeom,
        PropertiesPointerType pProperties
        ) const override;

    /**
     * @brief Creates a new element pointer from an existing geometry
     * @param NewId the ID of the new element
     * @param pGeom the  geometry taken to create the condition
     * @param pProperties the properties assigned to the new element
     * @param pMasterGeom the paired geometry
     * @return a Pointer to the new element
     */
    Condition::Pointer Create(
        IndexType NewId,
        GeometryPointerType pGeom,
        PropertiesPointerType pProperties,
        GeometryPointerType pMasterGeom
        ) const override;

    /**
     * this is called during the assembling process in order
     * to calculate the condition contribution in explicit calculation.
     * NodalData is modified Inside the function, so the
     * The "AddEXplicit" FUNCTIONS THE ONLY FUNCTIONS IN WHICH A CONDITION
     * IS ALLOWED TO WRITE ON ITS NODES.
     * the caller is expected to ensure thread safety hence
     * SET/UNSETLOCK MUST BE PERFORMED IN THE STRATEGY BEFORE CALLING THIS FUNCTION
     * @param rCurrentProcessInfo the current process info instance
     */
    void AddExplicitContribution(ProcessInfo& rCurrentProcessInfo) override;

    /******************************************************************/
    /********** AUXILLIARY METHODS FOR GENERAL CALCULATIONS ***********/
    /******************************************************************/

    /**
     * @brief Sets on rResult the ID's of the element degrees of freedom
     * @param rResult The result vector with the ID's of the DOF
     * @param rCurrentProcessInfo the current process info instance
     */

    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets on ConditionalDofList the degrees of freedom of the considered element geometry
     * @param rConditionalDofList The list of DOFs
     * @param rCurrentProcessInfo The current process info instance
     */

    void GetDofList(
        DofsVectorType& rConditionalDofList,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo The current process information
     */
    int Check( const ProcessInfo& rCurrentProcessInfo ) override;

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

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    MortarBaseConditionMatrices mPreviousMortarOperators; /// These are the mortar operators from the previous converged step, necessary for a consistent definition of the slip

    // TODO: Define the "CL" or friction law to compute this. Or do it nodally

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /********************************************************************************/
    /**************** METHODS TO CALCULATE MORTAR CONDITION MATRICES ****************/
    /********************************************************************************/

    /**
     * @brief Calculates the local contibution of the LHS
     * @param rLocalLHS The local LHS to compute
     * @param rMortarConditionMatrices The mortar operators to be considered
     * @param rDerivativeData The class containing all the derivatives uses to compute the jacobian
     * @param rActiveInactive The integer that is used to identify which case is the currectly computed
     */

    void CalculateLocalLHS(
        Matrix& rLocalLHS,
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const IndexType rActiveInactive,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculates the local contibution of the RHS
     * @param rLocalRHS The local RHS to compute
     * @param rMortarConditionMatrices The mortar operators to be considered
     * @param rDerivativeData The class containing all the derivatives uses to compute the jacobian
     * @param rActiveInactive The integer that is used to identify which case is the currectly computed
     */

    void CalculateLocalRHS(
        Vector& rLocalRHS,
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const IndexType rActiveInactive,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /******************************************************************/
    /********** AUXILLIARY METHODS FOR GENERAL CALCULATIONS ***********/
    /******************************************************************/

    /**
     * @brief Returns a value depending of the active/inactive set
     * @param CurrentGeometry The geometry containing the nodes that are needed to be checked as active or inactive
     * @return The integer that can be used to identify the case to compute
     */

    IndexType GetActiveInactiveValue(GeometryType& CurrentGeometry) const override
    {
        IndexType value = 0;
        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            if (CurrentGeometry[i_node].Is(ACTIVE) == true) {
                if (CurrentGeometry[i_node].Is(SLIP) == true)
                    value += std::pow(3, i_node);
                else
                    value += 2 * std::pow(3, i_node);
            }
        }

        return value;
    }

    /**
     * @brief This method returns a vector containing the friction coefficients
     * @return The friction coefficient corresponding to each node
     */
    array_1d<double, TNumNodes> GetFrictionCoefficient()
    {
        // The friction coefficient
        array_1d<double, TNumNodes> friction_coeffient_vector;
        auto& geom = this->GetGeometry();

        for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node) {
            friction_coeffient_vector[i_node] = geom[i_node].GetValue(FRICTION_COEFFICIENT);
        }

        // TODO: Define the "CL" or friction law to compute this

        return friction_coeffient_vector;
    }

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
     * @brief It calculates the matrix containing the tangent vector of the slip (for frictional contact)
     * @param ThisNodes The geometry to calculate
     * @return tangent_matrix The matrix containing the tangent vectors of the slip
     */

    static inline BoundedMatrix<double, TNumNodes, TDim> ComputeTangentMatrixSlip(const GeometryType& ThisNodes) {
        /* DEFINITIONS */
        // Zero tolerance
        const double zero_tolerance = std::numeric_limits<double>::epsilon();
        // Tangent matrix
        BoundedMatrix<double, TNumNodes, TDim> tangent_matrix;

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            const array_1d<double, 3>& slip = ThisNodes[i_node].FastGetSolutionStepValue(WEIGHTED_SLIP);
            const double norm_slip = norm_2(slip);
            if (norm_slip > zero_tolerance) { // Non zero slip
                const array_1d<double, 3> tangent_slip = slip/norm_slip;
                for (std::size_t i_dof = 0; i_dof < TDim; ++i_dof)
                    tangent_matrix(i_node, i_dof) = tangent_slip[i_dof];
            } else { // We consider the tangent direction as auxiliar
                const array_1d<double, 3>& tangent_xi = ThisNodes[i_node].GetValue(TANGENT_XI);
                for (std::size_t i_dof = 0; i_dof < TDim; ++i_dof)
                    tangent_matrix(i_node, i_dof) = tangent_xi[i_dof];
            }
        }

        return tangent_matrix;
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
    }

    ///@}

}; // Class AugmentedLagrangianMethodFrictionalMortarContactCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.

#endif // KRATOS_ALM_FRICTIONAL_MORTAR_CONTACT_CONDITION_H_INCLUDED  defined
