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

#if !defined(KRATOS_PENALTY_FRICTIONAL_MORTAR_CONTACT_AXISYM_CONDITION_H_INCLUDED )
#define  KRATOS_PENALTY_FRICTIONAL_MORTAR_CONTACT_AXISYM_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/penalty_frictional_mortar_contact_condition.h"

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
 * @class PenaltyMethodFrictionalMortarContactAxisymCondition
 * @ingroup ContactStructuralMechanicsApplication
 * @brief PenaltyMethodFrictionalMortarContactAxisymCondition
 * @todo Complete this
 * @author Vicente Mataix Ferrandiz
 */
template< std::size_t TNumNodes, bool TNormalVariation, std::size_t TNumNodesMaster = TNumNodes>
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) PenaltyMethodFrictionalMortarContactAxisymCondition
    : public PenaltyMethodFrictionalMortarContactCondition<2, TNumNodes, TNormalVariation, TNumNodesMaster>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of PenaltyMethodFrictionalMortarContactAxisymCondition
    KRATOS_CLASS_POINTER_DEFINITION( PenaltyMethodFrictionalMortarContactAxisymCondition );

    typedef MortarContactCondition<2, TNumNodes, FrictionalCase::FRICTIONAL_PENALTY, TNormalVariation, TNumNodesMaster> MortarBaseType;

    typedef PenaltyMethodFrictionalMortarContactCondition<2, TNumNodes, TNormalVariation, TNumNodesMaster>                    BaseType;

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

    typedef DerivativeDataFrictional<2, TNumNodes, TNormalVariation, TNumNodesMaster>                               DerivativeDataType;

    static constexpr IndexType MatrixSize = 2 * (TNumNodes + TNumNodesMaster);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    PenaltyMethodFrictionalMortarContactAxisymCondition(): BaseType()
    {
    }

    // Constructor 1
    PenaltyMethodFrictionalMortarContactAxisymCondition(IndexType NewId, GeometryPointerType pGeometry):BaseType(NewId, pGeometry)
    {
    }

    // Constructor 2
    PenaltyMethodFrictionalMortarContactAxisymCondition(IndexType NewId, GeometryPointerType pGeometry, PropertiesPointerType pProperties):BaseType( NewId, pGeometry, pProperties )
    {
    }

    ///Copy constructor
    PenaltyMethodFrictionalMortarContactAxisymCondition( PenaltyMethodFrictionalMortarContactAxisymCondition const& rOther)
    {
    }

    /// Destructor.
    ~PenaltyMethodFrictionalMortarContactAxisymCondition() override;

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
        buffer << "PenaltyMethodFrictionalMortarContactAxisymCondition #" << this->Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PenaltyMethodFrictionalMortarContactAxisymCondition #" << this->Id();
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

}; // Class PenaltyMethodFrictionalMortarContactAxisymCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.

#endif // KRATOS_PENALTY_FRICTIONAL_MORTAR_CONTACT_AXISYM_CONDITION_H_INCLUDED  defined
