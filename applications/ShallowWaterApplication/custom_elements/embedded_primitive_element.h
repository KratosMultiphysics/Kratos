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
//                   Ruben Zorrilla
//


#pragma once

// System includes


// External includes


// Project includes
#include "wave_element.h"

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
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

///@brief Implementation of a linear element for shallow water problems
template<std::size_t TNumNodes>
class EmbeddedPrimitiveElement : public WaveElement<TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;

    typedef Node<3> NodeType;

    typedef Geometry<NodeType> GeometryType;

    typedef WaveElement<TNumNodes> WaveElementType;

    typedef typename WaveElementType::NodesArrayType NodesArrayType;

    typedef typename WaveElementType::PropertiesType PropertiesType;

    typedef typename WaveElementType::ElementData ElementData;

    typedef typename WaveElementType::LocalMatrixType LocalMatrixType;

    typedef typename WaveElementType::LocalVectorType LocalVectorType;

    typedef typename WaveElementType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    ///@}
    ///@name Pointer definition
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(EmbeddedPrimitiveElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    EmbeddedPrimitiveElement() : WaveElementType(){}

    /**
     * @brief Constructor using an array of nodes
     */
    EmbeddedPrimitiveElement(IndexType NewId, const NodesArrayType& ThisNodes) : WaveElementType(NewId, ThisNodes){}

    /**
     * @brief Constructor using Geometry
     */
    EmbeddedPrimitiveElement(IndexType NewId, GeometryType::Pointer pGeometry) : WaveElementType(NewId, pGeometry){}

    /**
     * @brief Constructor using Geometry and Properties
     */
    EmbeddedPrimitiveElement(IndexType NewId, GeometryType::Pointer pGeometry, typename PropertiesType::Pointer pProperties) : WaveElementType(NewId, pGeometry, pProperties){}

    /**
     * @brief Destructor
     */
    ~ EmbeddedPrimitiveElement() override {};

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Create a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<EmbeddedPrimitiveElement<TNumNodes>>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    /**
     * @brief Create a new element pointer
     * @param NewId: the ID of the new element
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<EmbeddedPrimitiveElement<TNumNodes>>(NewId, pGeom, pProperties);
    }

    /**
     * @brief Create a new element pointer and clone the previous element data
     * @param NewId the ID of the new element
     * @param rThisNodes the nodes of the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override
    {
        Element::Pointer p_new_elem = Create(NewId, this->GetGeometry().Create(ThisNodes), this->pGetProperties());
        p_new_elem->SetData(this->GetData());
        p_new_elem->Set(Flags(*this));
        return p_new_elem;
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Turn back information as a string.
     */
    std::string Info() const override
    {
        return "EmbeddedPrimitiveElement";
    }

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    static constexpr IndexType mLocalSize = WaveElementType::mLocalSize;

    ///@}
    ///@name Protected Operations
    ///@{

    void UpdateGaussPointData(ElementData& rData, const array_1d<double,TNumNodes>& rN) override;

    double StabilizationParameter(const ElementData& rData) const override;

    void CalculateGeometryData(
        const GeometryType& rGeometry,
        Vector &rGaussWeights,
        Matrix &rNContainer,
        ShapeFunctionsGradientsType &rDN_DX) override;

    ///@}

private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class EmbeddedPrimitiveElement

///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.
