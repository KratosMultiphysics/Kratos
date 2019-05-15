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

#ifndef KRATOS_CONSERVED_VAR_ELEM_H_INCLUDED
#define KRATOS_CONSERVED_VAR_ELEM_H_INCLUDED

// System includes


// External includes


// Project includes
#include "primitive_var_element.hpp"

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

/// Implementation of a linear element for shallow water interpolating coserved variables
template< unsigned int TNumNodes >
class ConservedVarElement : public PrimitiveVarElement<TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of ConservedVarElement
    KRATOS_CLASS_POINTER_DEFINITION( ConservedVarElement );

    typedef PrimitiveVarElement<TNumNodes>                             BaseType;

    typedef typename BaseType::IndexType                              IndexType;

    typedef typename BaseType::GeometryType                        GeometryType;

    typedef typename BaseType::GeometryType::Pointer        GeometryPointerType;

    typedef typename BaseType::PropertiesType                    PropertiesType;

    typedef typename BaseType::PropertiesType::Pointer    PropertiesPointerType;

    typedef typename BaseType::NodesArrayType                    NodesArrayType;

    typedef typename BaseType::ElementVariables                ElementVariables;

    typedef typename BaseType::EquationIdVectorType        EquationIdVectorType;

    typedef typename BaseType::DofsVectorType                    DofsVectorType;

    typedef typename BaseType::VectorType                            VectorType;

    typedef typename BaseType::MatrixType                            MatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ConservedVarElement() :
        BaseType()
    {}

    /// Constructor using a Geometry instance
    ConservedVarElement(IndexType NewId, GeometryPointerType pGeometry) :
        BaseType(NewId, pGeometry)
    {}

    /// Constructor using geometry and properties
    ConservedVarElement(IndexType NewId, GeometryPointerType pGeometry, PropertiesPointerType pProperties) :
        BaseType(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~ ConservedVarElement() override {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& rThisNodes, PropertiesPointerType pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive< ConservedVarElement <TNumNodes> >(NewId, this->GetGeometry().Create(rThisNodes), pProperties);
        KRATOS_CATCH("")
    }

    Element::Pointer Create(IndexType NewId, GeometryPointerType pGeom, PropertiesPointerType pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive< ConservedVarElement <TNumNodes> >(NewId, pGeom, pProperties);
        KRATOS_CATCH("")
    }

    /**
     * It clones the selected element variables, creating a new one
     * @param NewId the ID of the new element
     * @param rThisNodes the nodes of the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override
    {
        KRATOS_TRY
        Element::Pointer p_new_elem = Create(NewId, this->GetGeometry().Create(ThisNodes), this->pGetProperties());
        p_new_elem->SetData(this->GetData());
        p_new_elem->Set(Flags(*this));
        return p_new_elem;
        KRATOS_CATCH("")
    }

    /// Check that all required data containers are properly initialized and registered in Kratos
    /**
     * @return 0 if no errors are detected.
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    /// Fill given vector with the linear system row index for the element's degrees of freedom
    /**
     * @param rResult
     * @param rCurrentProcessInfo
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    /// Fill given array with containing the element's degrees of freedom
    /**
     * @param rElementalDofList
     * @param rCurrentProcessInfo
     */
    void GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo) override;

    /// Evaluate the elemental contribution to the problem for turbulent viscosity.
    /**
     * @param rLeftHandSideMatrix Elemental left hand side matrix
     * @param rRightHandSideVector Elemental right hand side vector
     * @param rCurrentProcessInfo Reference to the ProcessInfo from the ModelPart containg the element
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

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

    void GetNodalValues(ElementVariables& rVariables) override;

    void GetElementValues(const BoundedMatrix<double,TNumNodes, 2>& rDN_DX, ElementVariables& rVariables) override;

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
    ///@name Serialization
    ///@{

    friend class Serializer;

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


    ///@}


}; // Class ConservedVarElement

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CONSERVED_VAR_ELEM_H_INCLUDED  defined
