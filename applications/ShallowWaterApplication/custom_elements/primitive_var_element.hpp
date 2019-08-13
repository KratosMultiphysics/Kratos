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

#ifndef KRATOS_PRIMITIVE_VAR_ELEM_H_INCLUDED
#define KRATOS_PRIMITIVE_VAR_ELEM_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"

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

/// Implementation of a linear element for shallow water problems
template< unsigned int TNumNodes >
class PrimitiveVarElement : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of PrimitiveVarElement
    KRATOS_CLASS_POINTER_DEFINITION( PrimitiveVarElement );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PrimitiveVarElement() :
        Element()
    {}

    /// Constructor using a Geometry instance
    PrimitiveVarElement(IndexType NewId, GeometryType::Pointer pGeometry) :
        Element(NewId, pGeometry)
    {}

    /// Constructor using geometry and properties
    PrimitiveVarElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
        Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~ PrimitiveVarElement() {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive< PrimitiveVarElement <TNumNodes> >(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
        KRATOS_CATCH("")
    }

    virtual Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive< PrimitiveVarElement <TNumNodes> >(NewId, pGeom, pProperties);
        KRATOS_CATCH("")
    }

    /**
     * It clones the selected element variables, creating a new one
     * @param NewId the ID of the new element
     * @param rThisNodes the nodes of the new element
     * @return a Pointer to the new element
     */
    virtual Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override
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
    virtual int Check(const ProcessInfo& rCurrentProcessInfo) override;

    /// Fill given vector with the linear system row index for the element's degrees of freedom
    /**
     * @param rResult
     * @param rCurrentProcessInfo
     */
    virtual void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    /// Fill given array with containing the element's degrees of freedom
    /**
     * @param rElementalDofList
     * @param rCurrentProcessInfo
     */
    virtual void GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo) override;

    /// Evaluate the elemental contribution to the problem for turbulent viscosity.
    /**
     * @param rLeftHandSideMatrix Elemental left hand side matrix
     * @param rRightHandSideVector Elemental right hand side vector
     * @param rCurrentProcessInfo Reference to the ProcessInfo from the ModelPart containg the element
     */
    virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    virtual void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    virtual void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

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

    struct ElementVariables
    {
        double dt_inv;
        double lumping_factor;
        double dyn_tau;
        double gravity;
        double manning2;
        double height_units;

        double height;
        array_1d<double,2> velocity;
        array_1d<double,2> momentum;
        array_1d<double,2> height_grad;
        BoundedMatrix<double,2,2> velocity_grad;
        double velocity_div;

        array_1d<double, TNumNodes*3> depth;
        array_1d<double, TNumNodes*3> rain;
        array_1d<double, TNumNodes*3> unknown;
        array_1d<double, TNumNodes*3> proj_unk;
    };

    virtual void InitializeElement(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateGeometry(BoundedMatrix<double, TNumNodes, 2>& rDN_DX, double& rArea);

    virtual double ComputeElemSize(const BoundedMatrix<double, TNumNodes, 2>& rDN_DX);

    virtual void GetNodalValues(ElementVariables& rVariables);

    virtual void GetElementValues(const BoundedMatrix<double,TNumNodes, 2>& rDN_DX, ElementVariables& rVariables);

    virtual void ComputeStabilizationParameters(const ElementVariables& rVariables,
                                        const double& rElemSize,
                                        double& rTauU,
                                        double& rTauH,
                                        double& rKdc);

    virtual void ComputeAuxMatrices(
            const BoundedMatrix<double,TNumNodes, TNumNodes>& rNcontainer,
            const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
            const ElementVariables& rVariables,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rMassMatrixScalar,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rMassMatrixVector,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rScalarGrad,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rVectorDiv,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rScalarDiff,
            BoundedMatrix<double,TNumNodes*3,TNumNodes*3>& rVectorDiff );

    virtual void CalculateLumpedMassMatrix(BoundedMatrix<double, TNumNodes*3, TNumNodes*3>& rMassMatrix);

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

}; // Class PrimitiveVarElement

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_PRIMITIVE_VAR_ELEM_H_INCLUDED  defined
