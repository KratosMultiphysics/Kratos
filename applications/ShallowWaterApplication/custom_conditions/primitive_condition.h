//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#pragma once

// System includes


// External includes


// Project includes
#include "wave_condition.h"

namespace Kratos
{
///@addtogroup ShallowWaterApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @ingroup ShallowWaterApplication
 * @class PrimitiveCondition
 * @brief Implementation of a condition for shallow water waves problems
 * @author Miguel Maso Sotomayor
 */
template<std::size_t TNumNodes>
class PrimitiveCondition : public WaveCondition<TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;

    typedef Node NodeType;

    typedef Geometry<NodeType> GeometryType;

    typedef WaveCondition<TNumNodes> BaseType;

    typedef typename BaseType::NodesArrayType NodesArrayType;

    typedef typename BaseType::PropertiesType PropertiesType;

    typedef typename BaseType::EquationIdVectorType EquationIdVectorType;

    typedef typename BaseType::DofsVectorType DofsVectorType;

    typedef typename BaseType::ConditionData ConditionData;

    typedef typename BaseType::LocalVectorType LocalVectorType;

    typedef typename BaseType::LocalMatrixType LocalMatrixType;

    ///@}
    ///@name Pointer definition
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(PrimitiveCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    PrimitiveCondition() : BaseType(){}

    /**
     * @brief Constructor using an array of nodes
     */
    PrimitiveCondition(IndexType NewId, const NodesArrayType& ThisNodes) : BaseType(NewId, ThisNodes){}

    /**
     * @brief Constructor using Geometry
     */
    PrimitiveCondition(IndexType NewId, GeometryType::Pointer pGeometry) : BaseType(NewId, pGeometry){}

    /**
     * @brief Constructor using Geometry and Properties
     */
    PrimitiveCondition(IndexType NewId, GeometryType::Pointer pGeometry, typename PropertiesType::Pointer pProperties) : BaseType(NewId, pGeometry, pProperties){}

    /**
     * @brief Destructor
     */
    ~ PrimitiveCondition() override {};

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Create a new condition pointer
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<PrimitiveCondition<TNumNodes>>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    /**
     * @brief Create a new condition pointer
     * @param NewId: the ID of the new condition
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<PrimitiveCondition<TNumNodes>>(NewId, pGeom, pProperties);
    }

    /**
     * @brief Create a new condition pointer and clone the previous condition data
     * @param NewId the ID of the new condition
     * @param rThisNodes the nodes of the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override
    {
        Condition::Pointer p_new_elem = Create(NewId, this->GetGeometry().Create(ThisNodes), this->pGetProperties());
        p_new_elem->SetData(this->GetData());
        p_new_elem->Set(Flags(*this));
        return p_new_elem;
    }

    /**
     * @brief Turn back information as a string.
     */
    std::string Info() const override
    {
        return "PrimitiveCondition";
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    static constexpr IndexType mLocalSize = BaseType::mLocalSize;

    ///@}
    ///@name Protected Operations
    ///@{

    void CalculateGaussPointData(
        ConditionData& rData,
        const IndexType PointIndex,
        const array_1d<double,TNumNodes>& rN) override;

    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
    }

    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class PrimitiveCondition

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.
