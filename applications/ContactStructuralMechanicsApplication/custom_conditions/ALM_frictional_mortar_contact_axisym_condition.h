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
#include "custom_conditions/ALM_frictional_mortar_contact_condition.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    using PointType = Point;
    using GeometryType = Geometry<Node>;
    using GeometryPointType = Geometry<PointType>;
    ///Type definition for integration methods
    using IntegrationMethod = GeometryData::IntegrationMethod;

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
 * @class AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition
 * @ingroup ContactStructuralMechanicsApplication
 * @brief AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition
 * @todo Complete this
 * @author Vicente Mataix Ferrandiz
 */
template< std::size_t TNumNodes, bool TNormalVariation >
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition
    : public AugmentedLagrangianMethodFrictionalMortarContactCondition<2, TNumNodes, TNormalVariation>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition );

    /// The base type for the mortar contact condition
    using MortarBaseType = MortarContactCondition<2, TNumNodes, FrictionalCase::FRICTIONAL, TNormalVariation>;

    /// The base type for the augmented Lagrangian method with frictional mortar contact condition
    using BaseType = AugmentedLagrangianMethodFrictionalMortarContactCondition<2, TNumNodes, TNormalVariation>;

    /// Type for the matrices used in the mortar contact condition
    using MortarConditionMatrices = typename MortarBaseType::MortarConditionMatrices;

    /// Type for the general variables used in the mortar contact condition
    using GeneralVariables = typename MortarBaseType::GeneralVariables;

    /// Type for the augmented Lagrangian data used in the mortar contact condition
    using AeData = typename MortarBaseType::AeData;

    /// The base type for the condition
    using ConditionBaseType = Condition;

    /// The vector type used in the condition
    using VectorType = typename ConditionBaseType::VectorType;

    /// The matrix type used in the condition
    using MatrixType = typename ConditionBaseType::MatrixType;

    /// The index type used in the condition
    using IndexType = typename ConditionBaseType::IndexType;

    /// Pointer type for the geometry of the condition
    using GeometryPointerType = typename ConditionBaseType::GeometryType::Pointer;

    /// The array type for the nodes of the condition
    using NodesArrayType = typename ConditionBaseType::NodesArrayType;

    /// Pointer type for the properties of the condition
    using PropertiesPointerType = typename ConditionBaseType::PropertiesType::Pointer;

    /// Type for the vector of equation IDs of the condition
    using EquationIdVectorType = typename ConditionBaseType::EquationIdVectorType;

    /// Type for the vector of DOFs of the condition
    using DofsVectorType = typename ConditionBaseType::DofsVectorType;

    /// Type for the array of conditions with points
    using ConditionArrayListType = std::vector<array_1d<PointType, 2>>;

    /// The decomposition type for the condition
    using DecompositionType = Line2D2<Point>;

    /// Type for the derivative data used in frictional mortar contact condition
    using DerivativeDataType = DerivativeDataFrictional<2, TNumNodes>;

    /// Matrix size definition
    static constexpr IndexType MatrixSize = 2 * (TNumNodes + TNumNodes) + TNumNodes;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition(): BaseType()
    {
    }

    // Constructor 1
    AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition(
        IndexType NewId,
        GeometryPointerType pGeometry
        ):BaseType(NewId, pGeometry)
    {
    }

    // Constructor 2
    AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition(
        IndexType NewId,
        GeometryPointerType pGeometry,
        PropertiesPointerType pProperties
        ):BaseType( NewId, pGeometry, pProperties )
    {
    }

    // Constructor 3
    AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition(
        IndexType NewId,
        GeometryPointerType pGeometry,
        PropertiesPointerType pProperties,
        GeometryType::Pointer pMasterGeometry
        ):BaseType( NewId, pGeometry, pProperties, pMasterGeometry )
    {
    }

    ///Copy constructor
    AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition( AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition const& rOther)
    {
    }

    /// Destructor.
    ~AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a new element pointer from an arry of nodes
     * @param NewId The ID of the new element
     * @param rThisNodes The nodes of the new element
     * @param pProperties The properties assigned to the new element
     * @return a Pointer to the new element
     */
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& rThisNodes,
        PropertiesPointerType pProperties
        ) const override;

    /**
     * @brief Creates a new element pointer from an existing geometry
     * @param NewId The ID of the new element
     * @param pGeom The  geometry taken to create the condition
     * @param pProperties The properties assigned to the new element
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


    /******************************************************************/
    /********** AUXILLIARY METHODS FOR GENERAL CALCULATIONS ***********/
    /******************************************************************/

    /**
     * @brief This functions returns if the computation is axisymmetric or not
     * @return If axisymmetric or not
     */
    bool IsAxisymmetric() const override;

    /**
     * This functions computes the integration weight to consider
     * @param rVariables The kinematic variables
     */
    double GetAxisymmetricCoefficient(const GeneralVariables& rVariables) const override;

    /**
     * Calculates the radius of axisymmetry
     * @param rVariables Internal values
     * @return Radius The radius of axisymmetry
     */

    double CalculateRadius(const GeneralVariables& rVariables) const;

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
        buffer << "AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition #" << this->Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition #" << this->Id();
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

}; // Class AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.
