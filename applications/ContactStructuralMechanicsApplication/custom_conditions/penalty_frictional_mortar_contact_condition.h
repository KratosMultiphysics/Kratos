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

#if !defined(KRATOS_PENALTY_FRICTIONAL_MORTAR_CONTACT_CONDITION_H_INCLUDED )
#define  KRATOS_PENALTY_FRICTIONAL_MORTAR_CONTACT_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "custom_conditions/mortar_contact_condition.h"

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
 * @class PenaltyMethodFrictionalMortarContactCondition
 * @ingroup ContactStructuralMechanicsApplication
 * @brief PenaltyMethodFrictionalMortarContactCondition
 * @details This is a frictional contact condition which employes the mortar method with penalty
 * The method has been taken from the Alexander Popps thesis:
 * Popp, Alexander: Mortar Methods for Computational Contact Mechanics and General Interface Problems, Technische Universität München, jul 2012
 * @author Vicente Mataix Ferrandiz
 * @tparam TDim The dimension of work
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TFrictional If we are solving a frictional or frictionless problem
 * @tparam TNormalVariation If we are consider normal variation
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster = TNumNodes>
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) PenaltyMethodFrictionalMortarContactCondition
    : public MortarContactCondition<TDim, TNumNodes, FrictionalCase::FRICTIONAL_PENALTY, TNormalVariation, TNumNodesMaster>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of PenaltyMethodFrictionalMortarContactCondition
    KRATOS_CLASS_POINTER_DEFINITION( PenaltyMethodFrictionalMortarContactCondition );

    typedef MortarContactCondition<TDim, TNumNodes, FrictionalCase::FRICTIONAL_PENALTY, TNormalVariation, TNumNodesMaster> BaseType;

    typedef Condition                                                                                             ConditionBaseType;

    typedef PairedCondition                                                                                 PairedConditionBaseType;

    typedef typename BaseType::MortarConditionMatrices                                                      MortarConditionMatrices;

    typedef typename BaseType::GeneralVariables                                                                    GeneralVariables;

    typedef typename BaseType::IntegrationUtility                                                                IntegrationUtility;

    typedef typename BaseType::DerivativesUtilitiesType                                                    DerivativesUtilitiesType;

    typedef typename BaseType::BelongType                                                                                BelongType;

    typedef typename BaseType::ConditionArrayListType                                                        ConditionArrayListType;

    typedef MortarOperator<TNumNodes, TNumNodesMaster>                                                  MortarBaseConditionMatrices;

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

    typedef DerivativeDataFrictional<TDim, TNumNodes, TNormalVariation, TNumNodesMaster>                         DerivativeDataType;

    static constexpr IndexType MatrixSize = TDim * (TNumNodes + TNumNodesMaster);

    static constexpr IndexType StepSlip = TNormalVariation ? 0 : 1;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    PenaltyMethodFrictionalMortarContactCondition()
        : BaseType()
    {
    }

    // Constructor 1
    PenaltyMethodFrictionalMortarContactCondition(
        IndexType NewId,
        GeometryPointerType pGeometry
        ):BaseType(NewId, pGeometry)
    {
    }

    // Constructor 2
    PenaltyMethodFrictionalMortarContactCondition(
        IndexType NewId,
        GeometryPointerType pGeometry,
        PropertiesPointerType pProperties
        ):BaseType( NewId, pGeometry, pProperties )
    {
    }

    // Constructor 3
    PenaltyMethodFrictionalMortarContactCondition(
        IndexType NewId,
        GeometryPointerType pGeometry,
        PropertiesPointerType pProperties,
        GeometryType::Pointer pMasterGeometry
        ):BaseType( NewId, pGeometry, pProperties, pMasterGeometry )
    {
    }

    ///Copy constructor
    PenaltyMethodFrictionalMortarContactCondition( PenaltyMethodFrictionalMortarContactCondition const& rOther)
    {
    }

    /// Destructor.
    ~PenaltyMethodFrictionalMortarContactCondition() override;

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
     * @brief Called at the begining of each solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

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
     * @brief This is called during the assembling process in order
     * to calculate the condition contribution in explicit calculation.
     * NodalData is modified Inside the function, so the
     * The "AddEXplicit" FUNCTIONS THE ONLY FUNCTIONS IN WHICH A CONDITION
     * IS ALLOWED TO WRITE ON ITS NODES.
     * the caller is expected to ensure thread safety hence
     * SET/UNSETLOCK MUST BE PERFORMED IN THE STRATEGY BEFORE CALLING THIS FUNCTION
     * @param rCurrentProcessInfo the current process info instance
     */
    void AddExplicitContribution(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This function is designed to make the element to assemble an rRHS vector identified by a variable rRHSVariable by assembling it to the nodes on the variable rDestinationVariable (double version)
     * @details The "AddEXplicit" FUNCTIONS THE ONLY FUNCTIONS IN WHICH AN ELEMENT IS ALLOWED TO WRITE ON ITS NODES.
     * The caller is expected to ensure thread safety hence SET/UNSETLOCK MUST BE PERFORMED IN THE STRATEGY BEFORE CALLING THIS FUNCTION
     * @param rRHSVector input variable containing the RHS vector to be assembled
     * @param rRHSVariable variable describing the type of the RHS vector to be assembled
     * @param rDestinationVariable variable in the database to which the rRHSVector will be assembled
     * @param rCurrentProcessInfo the current process info instance
     */
    void AddExplicitContribution(
        const VectorType& rRHSVector,
        const Variable<VectorType>& rRHSVariable,
        Variable<double >& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This function is designed to make the element to assemble an rRHS vector identified by a variable rRHSVariable by assembling it to the nodes on the variable (array_1d<double, 3>) version rDestinationVariable.
     * @details The "AddEXplicit" FUNCTIONS THE ONLY FUNCTIONS IN WHICH AN ELEMENT IS ALLOWED TO WRITE ON ITS NODES.
     * The caller is expected to ensure thread safety hence SET/UNSETLOCK MUST BE PERFORMED IN THE STRATEGY BEFORE CALLING THIS FUNCTION
     * @param rRHSVector input variable containing the RHS vector to be assembled
     * @param rRHSVariable variable describing the type of the RHS vector to be assembled
     * @param rDestinationVariable variable in the database to which the rRHSVector will be assembled
     * @param rCurrentProcessInfo the current process info instance
     */
    void AddExplicitContribution(const VectorType& rRHSVector,
        const Variable<VectorType>& rRHSVariable,
        Variable<array_1d<double, 3> >& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

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

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "PenaltyMethodFrictionalMortarContactCondition #" << this->Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PenaltyMethodFrictionalMortarContactCondition #" << this->Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        PrintInfo(rOStream);
        this->GetGeometry().PrintData(rOStream);
        this->GetPairedGeometry().PrintData(rOStream);
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

    bool mPreviousMortarOperatorsInitialized = false;     /// In order to know iw we need to initialize the previous operators

    MortarBaseConditionMatrices mPreviousMortarOperators; /// These are the mortar operators from the previous converged step, necessary for a consistent definition of the r_gt

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
     * @param rCurrentProcessInfo The current process information
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
     * @param rCurrentProcessInfo The current process information
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
     * @param rCurrentGeometry The geometry containing the nodes that are needed to be checked as active or inactive
     * @return The integer that can be used to identify the case to compute
     */
    IndexType GetActiveInactiveValue(GeometryType& rCurrentGeometry) const override
    {
        IndexType value = 0;
        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            if (rCurrentGeometry[i_node].Is(ACTIVE) == true) {
                if (rCurrentGeometry[i_node].Is(SLIP) == true)
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
        auto& r_geometry = this->GetGeometry();

        for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node) {
            friction_coeffient_vector[i_node] = r_geometry[i_node].GetValue(FRICTION_COEFFICIENT);
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
     * @brief It computes the previous mortar operators
     * @param rCurrentProcessInfo The current process information
     */
    void ComputePreviousMortarOperators( ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief It calculates the matrix containing the tangent vector of the r_gt (for frictional contact)
     * @param rGeometry The geometry to calculate
     * @return tangent_matrix The matrix containing the tangent vectors of the r_gt
     */
    static inline BoundedMatrix<double, TNumNodes, TDim> ComputeTangentMatrixSlip(const GeometryType& rGeometry) {
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

}; // Class PenaltyMethodFrictionalMortarContactCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.

#endif // KRATOS_PENALTY_FRICTIONAL_MORTAR_CONTACT_CONDITION_H_INCLUDED  defined
