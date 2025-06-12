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

#ifndef KRATOS_CRANK_NICOLSON_WAVE_ELEMENT_H_INCLUDED
#define KRATOS_CRANK_NICOLSON_WAVE_ELEMENT_H_INCLUDED

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
class CrankNicolsonWaveElement : public WaveElement<TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;

    typedef Node NodeType;

    typedef Geometry<NodeType> GeometryType;

    typedef WaveElement<TNumNodes> WaveElementType;

    typedef typename WaveElementType::MatrixType MatrixType;

    typedef typename WaveElementType::VectorType VectorType;

    typedef typename WaveElementType::NodesArrayType NodesArrayType;

    typedef typename WaveElementType::PropertiesType PropertiesType;

    typedef typename WaveElementType::LocalMatrixType LocalMatrixType;

    typedef typename WaveElementType::LocalVectorType LocalVectorType;

    typedef typename WaveElementType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    ///@}
    ///@name Pointer definition
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(CrankNicolsonWaveElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    CrankNicolsonWaveElement() : WaveElementType(){}

    /**
     * @brief Constructor using an array of nodes
     */
    CrankNicolsonWaveElement(IndexType NewId, const NodesArrayType& ThisNodes) : WaveElementType(NewId, ThisNodes){}

    /**
     * @brief Constructor using Geometry
     */
    CrankNicolsonWaveElement(IndexType NewId, GeometryType::Pointer pGeometry) : WaveElementType(NewId, pGeometry){}

    /**
     * @brief Constructor using Geometry and Properties
     */
    CrankNicolsonWaveElement(IndexType NewId, GeometryType::Pointer pGeometry, typename PropertiesType::Pointer pProperties) : WaveElementType(NewId, pGeometry, pProperties){}

    /**
     * @brief Destructor
     */
    ~ CrankNicolsonWaveElement() override {};

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
        return Kratos::make_intrusive<CrankNicolsonWaveElement<TNumNodes>>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
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
        return Kratos::make_intrusive<CrankNicolsonWaveElement<TNumNodes>>(NewId, pGeom, pProperties);
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
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Calculate the elemental mass matrix
     * @param rMassMatrix the elemental mass matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{


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
        return "CrankNicolsonWaveElement";
    }

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

}; // Class CrankNicolsonWaveElement

///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CRANK_NICOLSON_WAVE_ELEMENT_H_INCLUDED  defined
