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

#ifndef KRATOS_SWE_H_INCLUDED
#define KRATOS_SWE_H_INCLUDED

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
class SWE : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SWE);

    typedef std::size_t IndexType;

    typedef Node<3> NodeType;

    typedef array_1d<double, 3 * TNumNodes> LocalVectorType;

    typedef BoundedMatrix<double, 3*TNumNodes, 3*TNumNodes> LocalMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SWE() : Element(){}

    /// Constructor using a Geometry instance
    SWE(IndexType NewId, GeometryType::Pointer pGeometry) : Element(NewId, pGeometry){}

    /// Constructor using geometry and properties
    SWE(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : Element(NewId, pGeometry, pProperties){}

    /// Destructor.
    virtual ~ SWE(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive< SWE <TNumNodes, TFramework> >(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive< SWE <TNumNodes, TFramework> >(NewId, pGeom, pProperties);
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
        // Constants
        static constexpr size_t LocalSize = TNumNodes * 3;

        // Element values
        double epsilon;
        double dt_inv;
        double lumping_factor;
        double dyn_tau;
        double gravity;
        double manning2;
        double porosity;
        double height_units;
        double permeability;
        double discharge_penalty;

        // Element variables
        array_1d<double, 2> projected_momentum; // It is used to compute friction terms
        double momentum_div; // It is used to compute shock capturing
        double velocity_div; // It is used to compute the convective term
        double height;
        double wave_vel_2;
        array_1d<double, 2> surface_grad; // Is used to compute shock capturing
        array_1d<double, 2> velocity; // It is used to compute the convective stabilization

        // Unknowns and nodal values
        LocalVectorType rain;
        LocalVectorType unknown;
        LocalVectorType prev_unk;

        // Shape functions and derivatives
        BoundedMatrix<double, 2, LocalSize> N_q;
        array_1d<double, LocalSize> N_h;
        array_1d<double, LocalSize> DN_DX_q;
        BoundedMatrix<double, 2, LocalSize> DN_DX_h;
        BoundedMatrix<double, 2, LocalSize> Grad_q1;
        BoundedMatrix<double, 2, LocalSize> Grad_q2;
    };

    void InitializeElementVariables(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);

    void CalculateGeometry(BoundedMatrix<double, TNumNodes, 2>& rDN_DX, double& rArea);

    void GetNodalValues(ElementVariables& rVariables);

    void CalculateElementValues(const BoundedMatrix<double, TNumNodes, 2>& rDN_DX, ElementVariables& rVariables);

    void ComputeStabilizationParameters(
        const ElementVariables& rVariables,
        double& rTauU,
        double& rTauF);

    void ComputeConvectionStabilizationParameters(
        const ElementVariables& rVariables,
        double& rTau);

    void BuildAuxiliaryMatrices(
        const array_1d<double, TNumNodes>& rN,
        const BoundedMatrix<double, TNumNodes, 2>& rDN_DX,
        ElementVariables& rVariables);

    void AddInertiaTerms(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ElementVariables& rVariables);

    void AddConvectiveTerms(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ElementVariables& rVariables);

    void AddWaveTerms(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ElementVariables& rVariables);

    void AddFrictionTerms(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ElementVariables& rVariables);

    void AddStabilizationTerms(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ElementVariables& rVariables);

    void AddSourceTerms(
        VectorType& rRightHandSideVector,
        ElementVariables& rVariables);

    void CalculateLumpedMassMatrix(LocalMatrixType& rMassMatrix);

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

}; // Class SWE

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SWE_H_INCLUDED  defined
