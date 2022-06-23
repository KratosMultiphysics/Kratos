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

    typedef WaveElement<TNumNodes> BaseType;

    using MatrixType = typename BaseType::MatrixType;

    using VectorType = typename BaseType::VectorType;

    typedef typename BaseType::NodesArrayType NodesArrayType;

    typedef typename BaseType::PropertiesType PropertiesType;

    typedef typename BaseType::ElementData ElementData;

    typedef typename BaseType::LocalMatrixType LocalMatrixType;

    typedef typename BaseType::LocalVectorType LocalVectorType;

    typedef typename BaseType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

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
    EmbeddedPrimitiveElement() : BaseType(){}

    /**
     * @brief Constructor using an array of nodes
     */
    EmbeddedPrimitiveElement(IndexType NewId, const NodesArrayType& ThisNodes) : BaseType(NewId, ThisNodes){}

    /**
     * @brief Constructor using Geometry
     */
    EmbeddedPrimitiveElement(IndexType NewId, GeometryType::Pointer pGeometry) : BaseType(NewId, pGeometry){}

    /**
     * @brief Constructor using Geometry and Properties
     */
    EmbeddedPrimitiveElement(IndexType NewId, GeometryType::Pointer pGeometry, typename PropertiesType::Pointer pProperties) : BaseType(NewId, pGeometry, pProperties){}

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

    /**
     * @brief Calculate the elemental contribution to the problem
     * @param rLeftHandSideMatrix Elemental left hand side matrix
     * @param rRightHandSideVector Elemental right hand side vector
     * @param rCurrentProcessInfo Reference to the ProcessInfo from the ModelPart containing the element
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

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

    static constexpr IndexType mLocalSize = BaseType::mLocalSize;

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

    int CalculateDistances(Vector& rDistances) const;

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
