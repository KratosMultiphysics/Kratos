// // KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_PENALTY_FRICTIONLESS_MORTAR_CONTACT_CONDITION_H_INCLUDED )
#define  KRATOS_PENALTY_FRICTIONLESS_MORTAR_CONTACT_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/mortar_contact_condition.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The definition of the size type
    typedef std::size_t SizeType;

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
 * @class PenaltyMethodFrictionlessMortarContactCondition
 * @ingroup ContactStructuralMechanicsApplication
 * @brief PenaltyMethodFrictionlessMortarContactCondition
 * @details This is a contact condition which employes the mortar method with penalty constraint
 * The method has been taken from the Alexander Popps thesis:
 * Popp, Alexander: Mortar Methods for Computational Contact Mechanics and General Interface Problems, Technische Universität München, jul 2012
 * @author Vicente Mataix Ferrandiz
 * @tparam TDim The dimension of work
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TNormalVariation If we are consider normal variation
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, const SizeType TNumNodesMaster = TNumNodes >
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) PenaltyMethodFrictionlessMortarContactCondition
    : public MortarContactCondition<TDim, TNumNodes, FrictionalCase::FRICTIONLESS_PENALTY, TNormalVariation, TNumNodesMaster>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of PenaltyMethodFrictionlessMortarContactCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( PenaltyMethodFrictionlessMortarContactCondition );

    typedef MortarContactCondition<TDim, TNumNodes, FrictionalCase::FRICTIONLESS_PENALTY, TNormalVariation, TNumNodesMaster> BaseType;

    typedef typename BaseType::MortarConditionMatrices                    MortarConditionMatrices;

    typedef typename BaseType::GeneralVariables                                  GeneralVariables;

    typedef typename BaseType::AeData                                                      AeData;

    typedef typename BaseType::IntegrationUtility                              IntegrationUtility;

    typedef typename BaseType::DerivativesUtilitiesType                  DerivativesUtilitiesType;

    typedef typename BaseType::BelongType                                              BelongType;

    typedef typename BaseType::ConditionArrayType                              ConditionArrayType;

    typedef typename BaseType::ConditionArrayListType                      ConditionArrayListType;

    typedef typename BaseType::DecompositionType                                DecompositionType;

    typedef typename BaseType::DerivativeDataType                              DerivativeDataType;

    typedef Condition                                                           ConditionBaseType;

    typedef PairedCondition                                               PairedConditionBaseType;

    typedef typename ConditionBaseType::VectorType                                     VectorType;

    typedef typename ConditionBaseType::MatrixType                                     MatrixType;

    typedef typename ConditionBaseType::IndexType                                       IndexType;

    typedef typename ConditionBaseType::GeometryType::Pointer                 GeometryPointerType;

    typedef typename ConditionBaseType::NodesArrayType                             NodesArrayType;

    typedef typename ConditionBaseType::PropertiesType                             PropertiesType;

    typedef typename ConditionBaseType::PropertiesType::Pointer             PropertiesPointerType;

    typedef typename ConditionBaseType::EquationIdVectorType                 EquationIdVectorType;

    typedef typename ConditionBaseType::DofsVectorType                             DofsVectorType;

    /// Point definition
    typedef Point                                                                       PointType;

    /// Node type definition
    typedef Node<3>                                                                      NodeType;

    /// Geoemtry type definition
    typedef Geometry<NodeType>                                                       GeometryType;

    // Type definition for integration methods
    typedef GeometryType::IntegrationPointsArrayType                        IntegrationPointsType;

    static constexpr IndexType MatrixSize = TDim * (TNumNodes + TNumNodesMaster);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    PenaltyMethodFrictionlessMortarContactCondition()
        : BaseType()
    {
    }

    // Constructor 1
    PenaltyMethodFrictionlessMortarContactCondition(
        IndexType NewId,
        GeometryPointerType pGeometry
        ):BaseType(NewId, pGeometry)
    {
    }

    // Constructor 2
    PenaltyMethodFrictionlessMortarContactCondition(
        IndexType NewId,
        GeometryPointerType pGeometry,
        PropertiesPointerType pProperties
        ):BaseType( NewId, pGeometry, pProperties )
    {
    }

    // Constructor 3
    PenaltyMethodFrictionlessMortarContactCondition(
        IndexType NewId,
        GeometryPointerType pGeometry,
        PropertiesPointerType pProperties,
        GeometryType::Pointer pMasterGeometry
        ):BaseType( NewId, pGeometry, pProperties, pMasterGeometry )
    {
    }

    ///Copy constructor
    PenaltyMethodFrictionlessMortarContactCondition( PenaltyMethodFrictionlessMortarContactCondition const& rOther)
    {
    }

    /// Destructor.
    ~PenaltyMethodFrictionlessMortarContactCondition() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

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
    void AddExplicitContribution(const ProcessInfo& rCurrentProcessInfo) override;

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
        const Variable<double >& rDestinationVariable,
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
        const Variable<array_1d<double, 3> >& rDestinationVariable,
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
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    /**
     * @brief Sets on ConditionalDofList the degrees of freedom of the considered element geometry
     * @param rConditionalDofList The list of DOFs
     * @param rCurrentProcessInfo The current process info instance
     */
    void GetDofList(
        DofsVectorType& rConditionalDofList,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo The current process information
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

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
        buffer << "PenaltyMethodFrictionlessMortarContactCondition #" << this->Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PenaltyMethodFrictionlessMortarContactCondition #" << this->Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        PrintInfo(rOStream);
        this->GetParentGeometry().PrintData(rOStream);
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
     * @param rCurrentGeometry The geometry containing the nodes that are needed to be checked as active or inactive
     * @return The integer that can be used to identify the case to compute
     */
    IndexType GetActiveInactiveValue(const GeometryType& rCurrentGeometry) const override
    {
        IndexType value = 0;
        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
            if (rCurrentGeometry[i_node].Is(ACTIVE) == true)
                value += 1 << i_node;

        return value;
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

}; // Class PenaltyMethodFrictionlessMortarContactCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.

#endif // KRATOS_PENALTY_FRICTIONLESS_MORTAR_CONTACT_CONDITION_H_INCLUDED  defined
