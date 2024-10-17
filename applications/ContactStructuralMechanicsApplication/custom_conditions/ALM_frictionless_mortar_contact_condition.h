// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

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
    using SizeType = std::size_t;

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
 * @class AugmentedLagrangianMethodFrictionlessMortarContactCondition
 * @ingroup ContactStructuralMechanicsApplication
 * @brief AugmentedLagrangianMethodFrictionlessMortarContactCondition
 * @details This is a contact condition which employes the mortar method with dual lagrange multiplier
 * The method has been taken from the Alexander Popps thesis:
 * Popp, Alexander: Mortar Methods for Computational Contact Mechanics and General Interface Problems, Technische Universität München, jul 2012
 * @author Vicente Mataix Ferrandiz
 * @tparam TDim The dimension of work
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TNormalVariation If we are consider normal variation
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template< SizeType TDim, SizeType TNumNodes, bool TNormalVariation, const SizeType TNumNodesMaster = TNumNodes >
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) AugmentedLagrangianMethodFrictionlessMortarContactCondition
    : public MortarContactCondition<TDim, TNumNodes, FrictionalCase::FRICTIONLESS, TNormalVariation, TNumNodesMaster>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of AugmentedLagrangianMethodFrictionlessMortarContactCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( AugmentedLagrangianMethodFrictionlessMortarContactCondition );

    /// Base type for the mortar contact condition
    using BaseType = MortarContactCondition<TDim, TNumNodes, FrictionalCase::FRICTIONLESS, TNormalVariation, TNumNodesMaster>;

    /// Type for the matrices used in the mortar contact condition
    using MortarConditionMatrices = typename BaseType::MortarConditionMatrices;

    /// Base type for the condition
    using ConditionBaseType = Condition;

    /// Base type for the paired condition
    using PairedConditionBaseType = PairedCondition;

    /// Type for the vector used in the condition
    using VectorType = typename ConditionBaseType::VectorType;

    /// Type for the matrix used in the condition
    using MatrixType = typename ConditionBaseType::MatrixType;

    /// Type for the index used in the condition
    using IndexType = typename ConditionBaseType::IndexType;

    /// Pointer type for the geometry used in the condition
    using GeometryPointerType = typename ConditionBaseType::GeometryType::Pointer;

    /// Type for the array of nodes used in the condition
    using NodesArrayType = typename ConditionBaseType::NodesArrayType;

    /// Type for the properties used in the condition
    using PropertiesType = typename ConditionBaseType::PropertiesType;

    /// Pointer type for the properties used in the condition
    using PropertiesPointerType = typename PropertiesType::Pointer;

    /// Type for the vector of equation IDs used in the condition
    using EquationIdVectorType = typename ConditionBaseType::EquationIdVectorType;

    /// Type for the vector of DOFs used in the condition
    using DofsVectorType = typename ConditionBaseType::DofsVectorType;

    /// Definition of a point
    using PointType = Point;

    /// Definition of the geometry type
    using GeometryType = Geometry<Node>;

    /// Type for the integration points in the geometry
    using IntegrationPointsType = typename GeometryType::IntegrationPointsArrayType;

    /// Type for the array list of conditions with points
    using ConditionArrayListType = std::vector<array_1d<PointType, TDim>>;

    /// Line type definition
    using LineType = Line2D2<Point>;

    /// Triangle type definition
    using TriangleType = Triangle3D3<Point>;

    /// Type definition for the decomposition based on dimension
    using DecompositionType = typename std::conditional<TDim == 2, LineType, TriangleType>::type;

    /// Type definition for the derivative data based on dimension and number of nodes
    using DerivativeDataType = DerivativeData<TDim, TNumNodes, TNumNodesMaster>;

    /// Constant expression for matrix size
    static constexpr IndexType MatrixSize = TDim * (TNumNodes + TNumNodesMaster) + TNumNodes;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    AugmentedLagrangianMethodFrictionlessMortarContactCondition()
        : BaseType()
    {
    }

    // Constructor 1
    AugmentedLagrangianMethodFrictionlessMortarContactCondition(
        IndexType NewId,
        GeometryPointerType pGeometry
        ):BaseType(NewId, pGeometry)
    {
    }

    // Constructor 2
    AugmentedLagrangianMethodFrictionlessMortarContactCondition(
        IndexType NewId,
        GeometryPointerType pGeometry,
        PropertiesPointerType pProperties
        ):BaseType( NewId, pGeometry, pProperties )
    {
    }

    // Constructor 3
    AugmentedLagrangianMethodFrictionlessMortarContactCondition(
        IndexType NewId,
        GeometryPointerType pGeometry,
        PropertiesPointerType pProperties,
        GeometryType::Pointer pMasterGeometry
        ):BaseType( NewId, pGeometry, pProperties, pMasterGeometry )
    {
    }

    ///Copy constructor
    AugmentedLagrangianMethodFrictionlessMortarContactCondition( AugmentedLagrangianMethodFrictionlessMortarContactCondition const& rOther)
    {
    }

    /// Destructor.
    ~AugmentedLagrangianMethodFrictionlessMortarContactCondition() override;

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

    /********************************************************************************/
    /**************** METHODS TO CALCULATE MORTAR CONDITION MATRICES ****************/
    /********************************************************************************/

    /**
     * @brief Calculates the local contibution of the RHS
     * @param pCondition The condition pointer
     * @param rLocalRHS The local RHS to compute
     * @param rMortarConditionMatrices The mortar operators to be considered
     * @param rDerivativeData The class containing all the derivatives uses to compute the jacobian
     * @param rActiveInactive The integer that is used to identify which case is the currectly computed
     */
    static void StaticCalculateLocalRHS(
        PairedCondition* pCondition,
        Vector& rLocalRHS,
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const IndexType rActiveInactive,
        const ProcessInfo& rCurrentProcessInfo
        );

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
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
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
        buffer << "AugmentedLagrangianMethodFrictionlessMortarContactCondition #" << this->Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AugmentedLagrangianMethodFrictionlessMortarContactCondition #" << this->Id();
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
     * @param CurrentGeometry The geometry containing the nodes that are needed to be checked as active or inactive
     * @return The integer that can be used to identify the case to compute
     */
    IndexType GetActiveInactiveValue(const GeometryType& CurrentGeometry) const override
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

}; // Class AugmentedLagrangianMethodFrictionlessMortarContactCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.