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

#ifndef KRATOS_CONSERVATIVE_ELEMENT_RV_H_INCLUDED
#define KRATOS_CONSERVATIVE_ELEMENT_RV_H_INCLUDED

// System includes


// External includes


// Project includes
#include "conservative_element.h"

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
class ConservativeElementRV : public ConservativeElement<TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;

    typedef Node NodeType;

    typedef Geometry<NodeType> GeometryType;

    typedef ConservativeElement<TNumNodes> ConservativeElementType;

    typedef typename ConservativeElementType::WaveElementType WaveElementType;

    typedef typename ConservativeElementType::NodesArrayType NodesArrayType;

    typedef typename ConservativeElementType::PropertiesType PropertiesType;

    typedef typename ConservativeElementType::ElementData ElementData;

    typedef typename ConservativeElementType::LocalMatrixType LocalMatrixType;

    typedef typename ConservativeElementType::LocalVectorType LocalVectorType;

    typedef typename ConservativeElementType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    ///@}
    ///@name Pointer definition
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ConservativeElementRV);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    ConservativeElementRV() : ConservativeElementType(){}

    /**
     * @brief Constructor using an array of nodes
     */
    ConservativeElementRV(IndexType NewId, const NodesArrayType& ThisNodes) : ConservativeElementType(NewId, ThisNodes){}

    /**
     * @brief Constructor using Geometry
     */
    ConservativeElementRV(IndexType NewId, GeometryType::Pointer pGeometry) : ConservativeElementType(NewId, pGeometry){}

    /**
     * @brief Constructor using Geometry and Properties
     */
    ConservativeElementRV(IndexType NewId, GeometryType::Pointer pGeometry, typename PropertiesType::Pointer pProperties) : ConservativeElementType(NewId, pGeometry, pProperties){}

    /**
     * @brief Destructor
     */
    ~ ConservativeElementRV() override {};

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
        return Kratos::make_intrusive<ConservativeElementRV<TNumNodes>>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
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
        return Kratos::make_intrusive<ConservativeElementRV<TNumNodes>>(NewId, pGeom, pProperties);
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
    ///@name Input and output
    ///@{

    /**
     * @brief Turn back information as a string.
     */
    std::string Info() const override
    {
        return "ConservativeElementRV";
    }

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    static constexpr IndexType mLocalSize = ConservativeElementType::mLocalSize;

    ///@}
    ///@name Protected Operations
    ///@{

    void GetNodalData(ElementData& rData, const GeometryType& rGeometry, int Step = 0) override;

    void CalculateArtificialViscosity(
        BoundedMatrix<double,3,3>& rViscosity,
        BoundedMatrix<double,2,2>& rDiffusion,
        const ElementData& rData,
        const array_1d<double,TNumNodes>& rN,
        const BoundedMatrix<double,TNumNodes,2>& rDN_DX) override;

    void ShockCapturingParameters(
        double& rArtViscosity,
        double& rArtDiffusion,
        const ElementData& rData,
        const array_1d<double,TNumNodes>& rN,
        const BoundedMatrix<double,3,2>& rDN_DX);

    void AlgebraicResidual(
        array_1d<double,3>& rFlowResidual,
        double& rHeightresidual,
        BoundedMatrix<double,3,3>& rFlowGrad,
        array_1d<double,3>& rHeightGrad,
        const ElementData& rData,
        const array_1d<double,TNumNodes>& rN,
        const BoundedMatrix<double,3,2>& rDN_DX);

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

}; // Class ConservativeElementRV

///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CONSERVATIVE_ELEMENT_RV_H_INCLUDED  defined
