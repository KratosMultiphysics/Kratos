//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//			 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_GEOMETRICAL_CONDITION_H_INCLUDED )
#define  KRATOS_GEOMETRICAL_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/condition.h"

namespace Kratos
{
///@name Kratos Globals
///@{
///@}
///@name Type Definitions
///@{
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
 * @class GeometricalCondition
 * @ingroup KratocCore
 * @brief This is pure geometric condition. The only pourpose for this definition is to create dummy conditions
 * @details Herits all method from base condition, and overrides the Clone() method
 * @author Vicente Mataix Ferrandiz
 */
class GeometricalCondition
    : public Condition
{
public:

    ///@name Type Definitions
    ///@{

    /// We define the base class Condition
    typedef Condition BaseType;
    
    /// Dfinition of the index type
    typedef BaseType::IndexType IndexType;

    /// Definition of the size type
    typedef BaseType::SizeType SizeType;
    
    /// Definition of the node type
    typedef BaseType::NodeType NodeType;

    /// Definition of the properties type
    typedef BaseType::PropertiesType PropertiesType;

    /// Definition of the geometry type with given NodeType
    typedef BaseType::GeometryType GeometryType;

    /// Definition of nodes container type, redefined from GeometryType
    typedef BaseType::NodesArrayType NodesArrayType;
    
    /// Counted pointer of GeometricalCondition
    KRATOS_CLASS_POINTER_DEFINITION( GeometricalCondition);
    
    ///@}

public:

    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor.
     * @param NewId The Id of the new created condition
     */
    GeometricalCondition(IndexType NewId = 0);

    /**
     * @brief Constructor using an array of nodes
     * @param NewId The Id of the new created condition
     * @param rThisNodes The array of nodes taht will define the geometry that will define the condition
     */
    GeometricalCondition(
        IndexType NewId, 
        const NodesArrayType& rThisNodes
        );

    /**
     * @brief Constructor using Geometry
     * @param NewId The Id of the new created condition
     * @param pGeometry The pointer to the geometry that will define the condition
     */
    GeometricalCondition(
        IndexType NewId, 
        GeometryType::Pointer pGeometry
        );

    /**
     * @brief Constructor using Properties
     * @param NewId The Id of the new created condition
     * @param pGeometry The pointer to the geometry that will define the condition
     * @param pProperties The pointer to the properties that will define the behaviour of the condition
     */
    GeometricalCondition(
        IndexType NewId, 
        GeometryType::Pointer pGeometry, 
        PropertiesType::Pointer pProperties
        );

    ///Copy constructor
    GeometricalCondition(GeometricalCondition const& rOther);

    /// Destructor.
    ~GeometricalCondition() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    GeometricalCondition& operator=(GeometricalCondition const& rOther);

    ///@}
    ///@name Operations
    ///@{
   
    /**
     * @brief Creates a new condition pointer
     * @param NewId the ID of the new condition
     * @param ThisNodes the nodes of the new condition
     * @param pProperties the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(
        IndexType NewId, 
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new condition pointer
     * @param NewId the ID of the new condition
     * @param pGeom the geometry to be employed
     * @param pProperties the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new condition pointer and clones the previous condition data
     * @param NewId the ID of the new condition
     * @param ThisNodes the nodes of the new condition
     * @param pProperties the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone (
        IndexType NewId, 
        NodesArrayType const& ThisNodes
        ) const override;

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
    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class GeometricalCondition

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_GEOMETRICAL_CONDITION_H_INCLUDED  defined
