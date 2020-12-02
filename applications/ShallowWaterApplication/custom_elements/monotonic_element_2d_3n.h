//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_MONOTONIC_ELEMENT_2D_3N_H_INCLUDED
#define KRATOS_MONOTONIC_ELEMENT_2D_3N_H_INCLUDED

// System includes

// External includes

// Project includes
#include "shallow_water_2d_3.h"

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

class MonotonicElement2D3N : public ShallowWater2D3
{
public:

    ///@name Type Definitions
    ///@{

    typedef ShallowWater2D3 BaseType;

    ///@}
    ///@name Pointer Definitions

    /// Pointer definition of MonotonicElement2D3N
    KRATOS_CLASS_POINTER_DEFINITION(MonotonicElement2D3N);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    MonotonicElement2D3N(IndexType NewId = 0)
    : BaseType(NewId)
    {}

    /**
     * Constructor using an array of nodes
     */
    MonotonicElement2D3N(IndexType NewId, const NodesArrayType& ThisNodes)
    : BaseType(NewId, ThisNodes)
    {}

    /**
     * Constructor using Geometry
     */
    MonotonicElement2D3N(IndexType NewId, GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
    {}

    /**
     * Constructor using Properties
     */
    MonotonicElement2D3N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
    {}

    /**
     * Destructor
     */
    ~MonotonicElement2D3N(){};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * ELEMENTS inherited from this class have to implement next
     * Create and Clone methods: MANDATORY
     */

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<MonotonicElement2D3N>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<MonotonicElement2D3N>(NewId, pGeom, pProperties);
    }

    /**
     * @brief It creates a new element pointer and clones the previous element data
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override
    {
        Element::Pointer p_new_elem = Kratos::make_intrusive<MonotonicElement2D3N>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
        p_new_elem->SetData(this->GetData());
        p_new_elem->Set(Flags(*this));
        return p_new_elem;
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental consistent mass matrix
     * @param rMassMatrix the elemental mass matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * In the monotonic element this is called during the assembling
     * process in order to calculate the elemental diffusion matrix
     * to ensure monotonicity
     * @param rDampingMatrix the elemental damping matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

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
        return "Monotonic shallow water element";
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

    void AddShockCapturingTerm(
        MatrixType& rLHS,
        const ElementData& rData,
        const BoundedMatrix<double,3,2>& rDN_DX) override {}

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
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class MonotonicElement2D3N

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_MONOTONIC_ELEMENT_2D_3N_H_INCLUDED  defined
