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

#if !defined(KRATOS_ALM_FRICTIONAL_MORTAR_CONTACT_AXISYM_CONDITION_H_INCLUDED )
#define  KRATOS_ALM_FRICTIONAL_MORTAR_CONTACT_AXISYM_CONDITION_H_INCLUDED

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
    KRATOS_CLASS_POINTER_DEFINITION( AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition );

    typedef AugmentedLagrangianMethodMortarContactCondition<2, TNumNodes, FrictionalCase::FRICTIONAL, TNormalVariation> MortarBaseType;

    typedef AugmentedLagrangianMethodFrictionalMortarContactCondition<2, TNumNodes, TNormalVariation>                         BaseType;

    typedef typename MortarBaseType::MortarConditionMatrices                                                   MortarConditionMatrices;

    typedef typename MortarBaseType::GeneralVariables                                                                 GeneralVariables;

    typedef typename MortarBaseType::AeData                                                                                     AeData;

    typedef Condition                                                                                                ConditionBaseType;

    typedef typename ConditionBaseType::VectorType                                                                          VectorType;

    typedef typename ConditionBaseType::MatrixType                                                                          MatrixType;

    typedef typename ConditionBaseType::IndexType                                                                            IndexType;

    typedef typename ConditionBaseType::GeometryType::Pointer                                                      GeometryPointerType;

    typedef typename ConditionBaseType::NodesArrayType                                                                  NodesArrayType;

    typedef typename ConditionBaseType::PropertiesType::Pointer                                                  PropertiesPointerType;

    typedef typename ConditionBaseType::EquationIdVectorType                                                      EquationIdVectorType;

    typedef typename ConditionBaseType::DofsVectorType                                                                  DofsVectorType;

    typedef typename std::vector<array_1d<PointType,2>>                                                         ConditionArrayListType;

    typedef Line2D2<Point>                                                                                           DecompositionType;

    typedef DerivativeDataFrictional<2, TNumNodes, TNormalVariation>                                                DerivativeDataType;

    static constexpr IndexType MatrixSize = 2 * (TNumNodes + TNumNodes) + TNumNodes;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition(): BaseType()
    {
    }

    // Constructor 1
    AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition(IndexType NewId, GeometryPointerType pGeometry):BaseType(NewId, pGeometry)
    {
    }

    // Constructor 2
    AugmentedLagrangianMethodFrictionalMortarContactAxisymCondition(IndexType NewId, GeometryPointerType pGeometry, PropertiesPointerType pProperties):BaseType( NewId, pGeometry, pProperties )
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
     * Creates a new element pointer from an arry of nodes
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
     * Creates a new element pointer from an existing geometry
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

    /******************************************************************/
    /********** AUXILLIARY METHODS FOR GENERAL CALCULATIONS ***********/
    /******************************************************************/

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

#endif // KRATOS_ALM_FRICTIONAL_MORTAR_CONTACT_AXISYM_CONDITION_H_INCLUDED  defined
