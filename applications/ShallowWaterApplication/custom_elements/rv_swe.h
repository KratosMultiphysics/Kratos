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

#ifndef KRATOS_RV_SWE_H_INCLUDED
#define KRATOS_RV_SWE_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "custom_utilities/element_framework.h"

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
template< size_t TNumNodes, ElementFramework TFramework >
class RV_SWE : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of RV_SWE
    KRATOS_CLASS_POINTER_DEFINITION( RV_SWE );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RV_SWE() : Element(){}

    /// Constructor using a Geometry instance
    RV_SWE(IndexType NewId, GeometryType::Pointer pGeometry) : Element(NewId, pGeometry){}

    /// Constructor using geometry and properties
    RV_SWE(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : Element(NewId, pGeometry, pProperties){}

    /// Destructor.
    virtual ~ RV_SWE(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive< RV_SWE <TNumNodes, TFramework> >(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive< RV_SWE <TNumNodes, TFramework> >(NewId, pGeom, pProperties);
    }

    /**
     * It clones the selected element variables, creating a new one
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

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

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
        static constexpr size_t LocalSize = TNumNodes*3;
        const double epsilon = 1e-4;

        double dt_inv;
        double lumping_factor;
        double dyn_tau;
        double gravity;
        double manning2;
        double porosity;
        double height_units;

        double height;
        array_1d<double, 2> velocity;
        array_1d<double, 2> projected_velocity;
        array_1d<double, 2> height_grad;
        BoundedMatrix<double, 2, 2> velocity_grad;
        double velocity_div;
        int sign;

        array_1d<double, LocalSize> depth;
        array_1d<double, LocalSize> rain;
        array_1d<double, LocalSize> unknown;
        array_1d<double, LocalSize> prev_unk;
        array_1d<double, LocalSize> proj_unk;

        BoundedMatrix<double, LocalSize, LocalSize> MassMatrixScalar;
        BoundedMatrix<double, LocalSize, LocalSize> MassMatrixVector;
        BoundedMatrix<double, LocalSize, LocalSize> ScalarGrad;
        BoundedMatrix<double, LocalSize, LocalSize> VectorDiv;
        BoundedMatrix<double, LocalSize, LocalSize> ScalarDiff;
        BoundedMatrix<double, LocalSize, LocalSize> VectorDiff;
        BoundedMatrix<double, LocalSize, LocalSize> Convection;
        BoundedMatrix<double, LocalSize, LocalSize> ScalarConvectionStabilization;
        BoundedMatrix<double, LocalSize, LocalSize> VectorConvectionStabilization;
        BoundedMatrix<double, LocalSize, LocalSize> FrictionStabilization;
    };

    void CheckVariableKey();

    void CheckVariableInNodalData(Node<3>& rNode);

    virtual void InitializeElementVariables(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateGeometry(BoundedMatrix<double, TNumNodes, 2>& rDN_DX, double& rArea);

    virtual void GetNodalValues(ElementVariables& rVariables);

    virtual void CalculateElementValues(const BoundedMatrix<double,TNumNodes, 2>& rDN_DX, ElementVariables& rVariables);

    virtual void ComputeStabilizationParameters(
        const ElementVariables& rVariables,
        double& rTauU,
        double& rTauH,
        double& rKappaU,
        double& rKappaH);

    virtual void BuildMassMatrices(
        array_1d<double,TNumNodes>& rN,
        ElementVariables& rVariables);

    virtual void BuildGradientMatrices(
        array_1d<double,TNumNodes>& rN,
        BoundedMatrix<double, TNumNodes, 2>& rDN_DX,
        ElementVariables& rVariables);

    virtual void BuildDiffusivityMatrices(
        BoundedMatrix<double, TNumNodes, 2>& rDN_DX,
        ElementVariables& rVariables);

    virtual void BuildConvectionMatrices(
        array_1d<double, TNumNodes>& rN,
        BoundedMatrix<double, TNumNodes, 2>& rDN_DX,
        ElementVariables& rVariables);

    virtual void AddInertiaTerms(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ElementVariables& rVariables);
    
    virtual void AddConvectiveTerms(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ElementVariables& rVariables);
    
    virtual void AddWaveTerms(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ElementVariables& rVariables);
    
    virtual void AddFrictionTerms(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ElementVariables& rVariables);
    
    virtual void AddStabilizationTerms(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ElementVariables& rVariables);
    
    virtual void AddSourceTerms(
        VectorType& rRightHandSideVector,
        ElementVariables& rVariables);

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

}; // Class RV_SWE

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_RV_SWE_H_INCLUDED  defined
