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

#ifndef KRATOS_CV_SWE_H_INCLUDED
#define KRATOS_CV_SWE_H_INCLUDED

// System includes


// External includes


// Project includes
#include "rv_swe.h"

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
class CV_SWE : public RV_SWE<TNumNodes, TFramework>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of CV_SWE
    KRATOS_CLASS_POINTER_DEFINITION( CV_SWE );

    typedef RV_SWE<TNumNodes, TFramework>               BaseType;
	typedef Properties                            PropertiesType;
    typedef Geometry<Node<3>>                       GeometryType;
    typedef Geometry<Node<3>>::PointsArrayType    NodesArrayType;
    typedef Vector                                    VectorType;
    typedef Matrix                                    MatrixType;

    // typedef typename BaseType::GeometryType                        GeometryType;

    // typedef typename BaseType::GeometryType::Pointer        GeometryPointerType;

    // typedef typename BaseType::PropertiesType                    PropertiesType;

    // typedef typename BaseType::PropertiesType::Pointer    PropertiesPointerType;

    // typedef typename BaseType::NodesArrayType                    NodesArrayType;

    // typedef typename BaseType::ElementVariables                ElementVariables;

    typedef typename BaseType::EquationIdVectorType        EquationIdVectorType;

    typedef typename BaseType::DofsVectorType                    DofsVectorType;

    // typedef typename BaseType::VectorType                            VectorType;

    // typedef typename BaseType::MatrixType                            MatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CV_SWE() : BaseType(){}

    /// Constructor using a Geometry instance
    CV_SWE(IndexType NewId, GeometryType::Pointer pGeometry) : BaseType(NewId, pGeometry){}

    /// Constructor using geometry and properties
    CV_SWE(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : BaseType(NewId, pGeometry, pProperties){}

    /// Destructor.
    virtual ~ CV_SWE(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive< CV_SWE <TNumNodes, TFramework> >(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive< CV_SWE <TNumNodes, TFramework> >(NewId, pGeom, pProperties);
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

    struct ConservativeElementVariables : public BaseType::ElementVariables
    {
        array_1d<double, 2> momentum;
        array_1d<double, 2> projected_momentum;
        BoundedMatrix<double, 2, 2> momentum_grad;
        double momentum_div;
        double wave_vel_2;
        array_1d<double, BaseType::ElementVariables::LocalSize> nodal_velocity;
    };

    void GetNodalValues(ConservativeElementVariables& rVariables);

    void CalculateElementValues(const BoundedMatrix<double,TNumNodes, 2>& rDN_DX, ConservativeElementVariables& rVariables);

    void ComputeStabilizationParameters(
        const ConservativeElementVariables& rVariables,
        double& rTauU,
        double& rTauH,
        double& rKappaU,
        double& rKappaH);

    void BuildConvectionMatrices(
        array_1d<double, TNumNodes>& rN,
        BoundedMatrix<double, TNumNodes, 2>& rDN_DX,
        ConservativeElementVariables& rVariables);

    void AddWaveTerms(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ConservativeElementVariables& rVariables);

    void AddFrictionTerms(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ConservativeElementVariables& rVariables);

    void AddStabilizationTerms(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ConservativeElementVariables& rVariables);

    void AddSourceTerms(
        VectorType& rRightHandSideVector,
        ConservativeElementVariables& rVariables);

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

}; // Class CV_SWE

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CV_SWE_H_INCLUDED  defined
