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

#if !defined(KRATOS_DALM_FRICTIONLESS_MORTAR_CONTACT_CONDITION_H_INCLUDED )
#define  KRATOS_DALM_FRICTIONLESS_MORTAR_CONTACT_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/ALM_frictionless_mortar_contact_condition.h"

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
///@{

/**
 * @class DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition
 * @ingroup ContactStructuralMechanicsApplication
 * @brief DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition
 * @details This is a contact condition which employes the mortar method with dual lagrange multiplier
 * The method has been taken from the Alexander Popps thesis:
 * Popp, Alexander: Mortar Methods for Computational Contact Mechanics and General Interface Problems, Technische Universität München, jul 2012
 * @author Vicente Mataix Ferrandiz
 */
template< std::size_t TDim, std::size_t TNumNodes, bool TNormalVariation >
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition
    : public AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim, TNumNodes, TNormalVariation>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition
    KRATOS_CLASS_POINTER_DEFINITION( DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition );

    typedef AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim, TNumNodes, TNormalVariation> BaseType;

    typedef typename BaseType::MortarConditionMatrices                                 MortarConditionMatrices;

    typedef Condition                                                                        ConditionBaseType;

    typedef PairedCondition                                                            PairedConditionBaseType;

    typedef typename ConditionBaseType::VectorType                                                  VectorType;

    typedef typename ConditionBaseType::MatrixType                                                  MatrixType;

    typedef typename ConditionBaseType::IndexType                                                    IndexType;

    typedef typename ConditionBaseType::GeometryType::Pointer                              GeometryPointerType;

    typedef typename ConditionBaseType::NodesArrayType                                          NodesArrayType;

    typedef typename ConditionBaseType::PropertiesType                                          PropertiesType;

    typedef typename ConditionBaseType::PropertiesType::Pointer                          PropertiesPointerType;

    typedef typename ConditionBaseType::EquationIdVectorType                              EquationIdVectorType;

    typedef typename ConditionBaseType::DofsVectorType                                          DofsVectorType;

    typedef typename std::vector<array_1d<PointType,TDim>>                              ConditionArrayListType;

    typedef Line2D2<Point>                                                                            LineType;

    typedef Triangle3D3<Point>                                                                    TriangleType;

    typedef typename std::conditional<TDim == 2, LineType, TriangleType >::type              DecompositionType;

    typedef DerivativeData<TDim, TNumNodes, TNormalVariation>                               DerivativeDataType;

    static constexpr std::size_t MatrixSize = TDim * (TNumNodes + TNumNodes) + 2 * TNumNodes;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition()
        : BaseType()
    {
    }

    // Constructor 1
    DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition(
        IndexType NewId,
        GeometryPointerType pGeometry
        ):BaseType(NewId, pGeometry)
    {
    }

    // Constructor 2
    DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition(
        IndexType NewId,
        GeometryPointerType pGeometry,
        PropertiesPointerType pProperties
        ):BaseType( NewId, pGeometry, pProperties )
    {
    }

    // Constructor 3
    DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition(
        IndexType NewId,
        GeometryPointerType pGeometry,
        PropertiesPointerType pProperties,
        GeometryType::Pointer pMasterGeometry
        ):BaseType( NewId, pGeometry, pProperties, pMasterGeometry )
    {
    }

    ///Copy constructor
    DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition( DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition const& rOther)
    {
    }

    /// Destructor.
    ~DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * Creates a new element pointer from an arry of nodes
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
     * Creates a new element pointer from an existing geometry
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
     * Creates a new element pointer from an existing geometry
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

    /******************************************************************/
    /********** AUXILLIARY METHODS FOR GENERAL CALCULATIONS ***********/
    /******************************************************************/

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     * @param rResult The result vector with the ID's of the DOF
     * @param rCurrentProcessInfo the current process info instance
     */

    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * Sets on ConditionalDofList the degrees of freedom of the considered element geometry
     * @param rConditionalDofList The list of DOFs
     * @param rCurrentProcessInfo The current process info instance
     */

    void GetDofList(
        DofsVectorType& rConditionalDofList,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
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
     * Calculates the local contibution of the LHS
     * @param rLocalLHS The local LHS to compute
     * @param rMortarConditionMatrices The mortar operators to be considered
     * @param rDerivativeData The class containing all the derivatives uses to compute the jacobian
     * @param rActiveInactive The integer that is used to identify which case is the currectly computed
     */

    void CalculateLocalLHS(
        Matrix& rLocalLHS,
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const IndexType rActiveInactive
        ) override;

    /**
     * Calculates the local contibution of the RHS
     * @param rLocalRHS The local RHS to compute
     * @param rMortarConditionMatrices The mortar operators to be considered
     * @param rDerivativeData The class containing all the derivatives uses to compute the jacobian
     * @param rActiveInactive The integer that is used to identify which case is the currectly computed
     */

    void CalculateLocalRHS(
        Vector& rLocalRHS,
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const IndexType rActiveInactive
        ) override;

    /**
     * Calculates the local DALM contibution of the LHS
     * @param rLocalLHS The local LHS to compute
     * @param rMortarConditionMatrices The mortar operators to be considered
     * @param rDerivativeData The class containing all the derivatives uses to compute the jacobian
     * @param rActiveInactive The integer that is used to identify which case is the currectly computed
     */

    void CalculateLocalLHSDALM(
        Matrix& rLocalLHS,
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const IndexType rActiveInactive
        );

    /**
     * Calculates the local DALM contibution of the RHS
     * @param rLocalRHS The local RHS to compute
     * @param rMortarConditionMatrices The mortar operators to be considered
     * @param rDerivativeData The class containing all the derivatives uses to compute the jacobian
     * @param rActiveInactive The integer that is used to identify which case is the currectly computed
     */

    void CalculateLocalRHSDALM(
        Vector& rLocalRHS,
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const IndexType rActiveInactive
        );

    /**
     * This method just resizes the LHS matrix
     * @param rLeftHandSideMatrix The LHS matrix
     */
    void ResizeLHS(MatrixType& rLeftHandSideMatrix) override;

    /**
     * This method just resizes the RHS vector
     * @param rRightHandSideVector The RHS vector
     */
    void ResizeRHS(VectorType& rRightHandSideVector) override;

    /**
     * This method just sets as zero the LHS matrix
     * @param rLeftHandSideMatrix The LHS matrix
     */
    void ZeroLHS(MatrixType& rLeftHandSideMatrix) override;

    /**
     *This method just sets as zero the RHS vector
     * @param rRightHandSideVector The RHS vector
     */
    void ZeroRHS(VectorType& rRightHandSideVector) override;

    /******************************************************************/
    /********** AUXILLIARY METHODS FOR GENERAL CALCULATIONS ***********/
    /******************************************************************/

    /**
     * Returns a value depending of the active/inactive set
     * @param CurrentGeometry The geometry containing the nodes that are needed to be checked as active or inactive
     * @return The integer that can be used to identify the case to compute
     */

    IndexType GetActiveInactiveValue(GeometryType& CurrentGeometry) const override
    {
        IndexType value = 0;
        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
            if (CurrentGeometry[i_node].Is(ACTIVE) == true)
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

}; // Class DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.

#endif // KRATOS_DALM_FRICTIONLESS_MORTAR_CONTACT_CONDITION_H_INCLUDED  defined
