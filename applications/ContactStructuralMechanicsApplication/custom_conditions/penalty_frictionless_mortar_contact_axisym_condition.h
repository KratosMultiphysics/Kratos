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
#include "custom_conditions/penalty_frictionless_mortar_contact_condition.h"

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
 * @class PenaltyMethodFrictionlessMortarContactAxisymCondition
 * @ingroup ContactStructuralMechanicsApplication
 * @brief PenaltyMethodFrictionlessMortarContactAxisymCondition
 * @todo Complete this
 * @author Vicente Mataix Ferrandiz
 */
template< std::size_t TNumNodes, bool TNormalVariation >
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) PenaltyMethodFrictionlessMortarContactAxisymCondition
    : public PenaltyMethodFrictionlessMortarContactCondition<2, TNumNodes, TNormalVariation>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of PenaltyMethodFrictionlessMortarContactAxisymCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( PenaltyMethodFrictionlessMortarContactAxisymCondition );

    /// Base type for the mortar contact condition
    using MortarBaseType = MortarContactCondition<2, TNumNodes, FrictionalCase::FRICTIONLESS_PENALTY, TNormalVariation>;

    /// Base type for the penalty method frictionless mortar contact condition
    using BaseType = PenaltyMethodFrictionlessMortarContactCondition<2, TNumNodes, TNormalVariation>;

    /// Type definition for mortar condition matrices
    using MortarConditionMatrices = typename MortarBaseType::MortarConditionMatrices;

    /// Type definition for general variables
    using GeneralVariables = typename MortarBaseType::GeneralVariables;

    /// Type definition for AE data
    using AeData = typename MortarBaseType::AeData;

    /// Base type for the condition
    using ConditionBaseType = Condition;

    /// Vector type definition
    using VectorType = typename ConditionBaseType::VectorType;

    /// Matrix type definition
    using MatrixType = typename ConditionBaseType::MatrixType;

    /// Index type definition
    using IndexType = typename ConditionBaseType::IndexType;

    /// Geometry pointer definition
    using GeometryPointerType = typename ConditionBaseType::GeometryType::Pointer;

    /// Nodes array type definition
    using NodesArrayType = typename ConditionBaseType::NodesArrayType;

    /// Properties pointer definition
    using PropertiesPointerType = typename ConditionBaseType::PropertiesType::Pointer;

    /// Equation ID vector type definition
    using EquationIdVectorType = typename ConditionBaseType::EquationIdVectorType;

    /// Dofs vector type definition
    using DofsVectorType = typename ConditionBaseType::DofsVectorType;

    /// Type definition for an array list of conditions
    using ConditionArrayListType = typename std::vector<array_1d<PointType, 2>>;

    /// Line type definition
    using DecompositionType = Line2D2<PointType>;

    /// Type definition for derivative data
    using DerivativeDataType = DerivativeData<2, TNumNodes>;

    /// Constant expression for matrix size
    static constexpr IndexType MatrixSize = 2 * (TNumNodes + TNumNodes) + TNumNodes;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    PenaltyMethodFrictionlessMortarContactAxisymCondition(): BaseType()
    {
    }

    // Constructor 1
    PenaltyMethodFrictionlessMortarContactAxisymCondition(
        IndexType NewId,
        GeometryPointerType pGeometry
        ):BaseType(NewId, pGeometry)
    {
    }

    // Constructor 2
    PenaltyMethodFrictionlessMortarContactAxisymCondition(
        IndexType NewId,
        GeometryPointerType pGeometry,
        PropertiesPointerType pProperties
        ):BaseType( NewId, pGeometry, pProperties )
    {
    }

    // Constructor 3
    PenaltyMethodFrictionlessMortarContactAxisymCondition(
        IndexType NewId,
        GeometryPointerType pGeometry,
        PropertiesPointerType pProperties,
        GeometryType::Pointer pMasterGeometry
        ):BaseType( NewId, pGeometry, pProperties, pMasterGeometry )
    {
    }

    ///Copy constructor
    PenaltyMethodFrictionlessMortarContactAxisymCondition( PenaltyMethodFrictionlessMortarContactAxisymCondition const& rOther)
    {
    }

    /// Destructor.
    ~PenaltyMethodFrictionlessMortarContactAxisymCondition() override;

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
        buffer << "PenaltyMethodFrictionlessMortarContactAxisymCondition #" << this->Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PenaltyMethodFrictionlessMortarContactAxisymCondition #" << this->Id();
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

}; // Class PenaltyMethodFrictionlessMortarContactAxisymCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.