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

#ifndef KRATOS_BOUSSINESQ_ELEMENT_H_INCLUDED
#define KRATOS_BOUSSINESQ_ELEMENT_H_INCLUDED

// System includes


// External includes


// Project includes
#include "primitive_element.h"

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
class BoussinesqElement : public PrimitiveElement<TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;

    typedef Node NodeType;

    typedef Geometry<NodeType> GeometryType;

    typedef PrimitiveElement<TNumNodes> BaseType;

    typedef typename BaseType::VectorType VectorType;

    typedef typename BaseType::NodesArrayType NodesArrayType;

    typedef typename BaseType::PropertiesType PropertiesType;

    typedef typename BaseType::ElementData ElementData;

    typedef typename BaseType::LocalMatrixType LocalMatrixType;

    typedef typename BaseType::LocalVectorType LocalVectorType;

    typedef typename BaseType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    ///@}
    ///@name Pointer definition
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(BoussinesqElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    BoussinesqElement() : BaseType(){}

    /**
     * @brief Constructor using an array of nodes
     */
    BoussinesqElement(IndexType NewId, const NodesArrayType& ThisNodes) : BaseType(NewId, ThisNodes){}

    /**
     * @brief Constructor using Geometry
     */
    BoussinesqElement(IndexType NewId, GeometryType::Pointer pGeometry) : BaseType(NewId, pGeometry){}

    /**
     * @brief Constructor using Geometry and Properties
     */
    BoussinesqElement(IndexType NewId, GeometryType::Pointer pGeometry, typename PropertiesType::Pointer pProperties) : BaseType(NewId, pGeometry, pProperties){}

    /**
     * @brief Destructor
     */
    ~ BoussinesqElement() override {};

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
        return Kratos::make_intrusive<BoussinesqElement<TNumNodes>>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
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
        return Kratos::make_intrusive<BoussinesqElement<TNumNodes>>(NewId, pGeom, pProperties);
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
     * @brief Calculate the velocity laplacian projection
     * @param rCurrentProcessInfo Reference to the ProcessInfo from the ModelPart containing the elements
     */
    void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Calculate the rhs according to the Adams-Moulton scheme
     * @param rRightHandSideVector Elemental right hand side vector
     * @param rCurrentProcessInfo Reference to the ProcessInfo from the ModelPart containing the element
     * @see ResidualBasedAdamsMoultonScheme
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Add the explicit contribution according to the Adams-Bashforth scheme
     * @param rCurrentProcessInfo the current process info instance
     */
    void AddExplicitContribution(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief This method provides the specifications/requirements of the element
     * @return specifications The required specifications/requirements
     */
    const Parameters GetSpecifications() const override;

    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Turn back information as a string.
     */
    std::string Info() const override
    {
        return "BoussinesqElement";
    }

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    static constexpr IndexType mLocalSize = BaseType::mLocalSize;

    ///@}
    ///@name Protected Operations
    ///@{

    void AddRightHandSide(
        LocalVectorType& rRHS,
        ElementData& rData,
        const Matrix& rNContainer,
        const ShapeFunctionsGradientsType& rDN_DXContainer,
        const Vector& rWeights);

    void AddDispersionProjection(
        LocalVectorType& rDispersionH,
        LocalVectorType& rDispersionU,
        const ElementData& rData,
        const array_1d<double,TNumNodes>& rN,
        const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
        const double Weight = 1.0);

    void GetNodalData(ElementData& rData, const GeometryType& rGeometry, int Step = 0) override;

    void CalculateArtificialViscosity(
        BoundedMatrix<double,3,3>& rViscosity,
        BoundedMatrix<double,2,2>& rDiffusion,
        const ElementData& rData,
        const array_1d<double,TNumNodes>& rN,
        const BoundedMatrix<double,TNumNodes,2>& rDN_DX) override;

    void AlgebraicResidual(
        double& rMassResidual,
        array_1d<double,2>& rFreeSurfaceGradient,
        const ElementData& rData,
        const array_1d<double,TNumNodes>& rN,
        const BoundedMatrix<double,TNumNodes,2>& rDN_DX) const;

    void AddDispersiveTerms(
        LocalVectorType& rVector,
        const ElementData& rData,
        const array_1d<double,TNumNodes>& rN,
        const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
        const double Weight = 1.0) override;

    void AddMassTerms(
        LocalMatrixType& rMatrix,
        const ElementData& rData,
        const array_1d<double,TNumNodes>& rN,
        const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
        const double Weight = 1.0) override;

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

}; // Class BoussinesqElement

///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_BOUSSINESQ_ELEMENT_H_INCLUDED  defined
