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
#include "custom_conditions/ALM_frictionless_mortar_contact_condition.h"

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
 * @class AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition
 * @ingroup ContactStructuralMechanicsApplication
 * @brief AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition
 * @todo Complete this
 * @author Vicente Mataix Ferrandiz
 */
template< std::size_t TNumNodes, bool TNormalVariation >
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition
    : public AugmentedLagrangianMethodFrictionlessMortarContactCondition<2, TNumNodes, TNormalVariation>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition );

    /// Base type for the mortar contact condition
    using MortarBaseType = MortarContactCondition<2, TNumNodes, FrictionalCase::FRICTIONLESS, TNormalVariation>;

    /// Base type for the condition
    using BaseType = AugmentedLagrangianMethodFrictionlessMortarContactCondition<2, TNumNodes, TNormalVariation>;

    /// Type for the matrices used in the mortar contact condition
    using MortarConditionMatrices = typename MortarBaseType::MortarConditionMatrices;

    /// Type for the general variables used in the mortar contact condition
    using GeneralVariables = typename MortarBaseType::GeneralVariables;

    /// Type for the AeData used in the mortar contact condition
    using AeData = typename MortarBaseType::AeData;

    /// Type for the base condition type
    using ConditionBaseType = Condition;

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

    /// Pointer type for the properties used in the condition
    using PropertiesPointerType = typename ConditionBaseType::PropertiesType::Pointer;

    /// Type for the vector of equation IDs used in the condition
    using EquationIdVectorType = typename ConditionBaseType::EquationIdVectorType;

    /// Type for the vector of DOFs used in the condition
    using DofsVectorType = typename ConditionBaseType::DofsVectorType;

    /// Type for the array list of conditions with points
    using ConditionArrayListType = std::vector<array_1d<PointType,2>>;

    /// Type for the line in 2D
    using DecompositionType = Line2D2<Point>;

    /// Type for the derivative data based on the dimension and number of nodes
    using DerivativeDataType = DerivativeData<2, TNumNodes>;

    /// Constant expression for matrix size
    static constexpr IndexType MatrixSize = 2 * (TNumNodes + TNumNodes) + TNumNodes;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition(): BaseType()
    {
    }

    // Constructor 1
    AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition(
        IndexType NewId,
        GeometryPointerType pGeometry
        ):BaseType(NewId, pGeometry)
    {
    }

    // Constructor 2
    AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition(
        IndexType NewId,
        GeometryPointerType pGeometry,
        PropertiesPointerType pProperties
        ):BaseType( NewId, pGeometry, pProperties )
    {
    }

    // Constructor 3
    AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition(
        IndexType NewId,
        GeometryPointerType pGeometry,
        PropertiesPointerType pProperties,
        GeometryType::Pointer pMasterGeometry
        ):BaseType( NewId, pGeometry, pProperties, pMasterGeometry )
    {
    }

    ///Copy constructor
    AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition( AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition const& rOther)
    {
    }

    /// Destructor.
    ~AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a new element pointer from an array of nodes
     * @param NewId The ID of the new element
     * @param ThisNodes tThe nodes of the new element
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

    /******************************************************************/
    /********** AUXILIARY METHODS FOR GENERAL CALCULATIONS ************/
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
     * @param rVariables: Internal values
     * @return Radius: The radius of axisymmetry
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
        buffer << "AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition #" << this->Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition #" << this->Id();
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

}; // Class AugmentedLagrangianMethodFrictionlessMortarContactAxisymCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.
