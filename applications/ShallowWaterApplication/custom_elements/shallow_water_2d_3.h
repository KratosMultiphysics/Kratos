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

#ifndef KRATOS_SHALLOW_WATER_2D_3_H_INCLUDED
#define KRATOS_SHALLOW_WATER_2D_3_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"

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

class ShallowWater2D3 : public Element
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of ShallowWater2D3
    KRATOS_CLASS_POINTER_DEFINITION(ShallowWater2D3);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    ShallowWater2D3(IndexType NewId = 0)
    : Element(NewId)
    {}

    /**
     * Constructor using an array of nodes
     */
    ShallowWater2D3(IndexType NewId, const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes)
    {}

    /**
     * Constructor using Geometry
     */
    ShallowWater2D3(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {}

    /**
     * Constructor using Properties
     */
    ShallowWater2D3(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {}

    /**
     * Destructor
     */
    ~ShallowWater2D3(){};

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
        return Kratos::make_intrusive<ShallowWater2D3>(NewId, GetGeometry().Create(ThisNodes), pProperties);
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
        return Kratos::make_intrusive<ShallowWater2D3>(NewId, pGeom, pProperties);
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
        Element::Pointer p_new_elem = Kratos::make_intrusive<ShallowWater2D3>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
        p_new_elem->SetData(this->GetData());
        p_new_elem->Set(Flags(*this));
        return p_new_elem;
    }

    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult: the elemental equation ID vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const override;

    /**
     * determines the elemental list of DOFs
     * @param ElementalDofList: the list of DOFs
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& CurrentProcessInfo) const override;

    /**
     * Getting method to obtain the variable which defines the degrees of freedom
     */
    void GetValuesVector(Vector& rValues, int Step = 0) const override;

    /**
     * Getting method to obtain the time derivative of variable which defines the degrees of freedom
     */
    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) const override;

    /**
     * Getting method to obtain the second time derivative of variable which defines the degrees of freedom
     */
    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override;

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side matrix only
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental mass matrix
     * @param rMassMatrix the elemental mass matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * In the flux corrected scheme this is called during the assembling
     * process in order to calculate the elemental diffusion matrix
     * to ensure monotonicity.
     * This method should not be called by the stabilized scheme.
     * @param rDampingMatrix the elemental damping matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * This method provides the place to perform checks on the completeness of the input
     * and the compatibility with the problem options as well as the constitutive laws selected
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * Access for variables on Integration points.
     * This gives access to variables stored in the constitutive law on each integration point.
     * Specializations of element must specify the actual interface to the integration points!
     * Note, that these functions expect a std::vector of values for the specified variable type that
     * contains a value for each integration point!
     * @param rVariable: the specified variable
     * @param rValues: where to store the values for the specified variable type at each integration point
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
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
        return "Shallow water element";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info() << Id();
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:

    ///@name Protected type definitions
    ///@{

    struct ElementData
    {
        double dt_inv;
        double stab_factor;
        double shock_stab_factor;
        double rel_dry_height;
        double gravity;

        double height;
        array_1d<double,3> flow_rate;
        array_1d<double,3> velocity;
        double manning2;

        array_1d<double,3> topography;
        array_1d<double,3> rain;
        array_1d<double,9> unknown;
        array_1d<double,9> mesh_acc;

        void InitializeData(const ProcessInfo& rCurrentProcessInfo);
        void GetNodalData(const GeometryType& rGeometry, const BoundedMatrix<double,3,2>& rDN_DX);
    };

    ///@}
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

    virtual void AddGradientTerms(
        MatrixType& rLHS,
        VectorType& rRHS,
        const ElementData& rData,
        const array_1d<double,3>& rN,
        const BoundedMatrix<double,3,2>& rDN_DX);

    virtual void AddSourceTerms(
        MatrixType& rLHS,
        VectorType& rRHS,
        const ElementData& rData,
        const array_1d<double,3>& rN,
        const BoundedMatrix<double,3,2>& rDN_DX);

    virtual void AddShockCapturingTerm(
        MatrixType& rLHS,
        const ElementData& rData,
        const BoundedMatrix<double,3,2>& rDN_DX);

    virtual void AddDesingularizationTerm(
        MatrixType& rLHS,
        const ElementData& rData);

    void ComputeMassMatrix(
        BoundedMatrix<double,9,9>& rMatrix,
        const ElementData& rData,
        const array_1d<double,3>& rN,
        const BoundedMatrix<double,3,2>& rDN_DX);

    void ComputeMassMatrix(
        BoundedMatrix<double,9,9>& rFlowMatrix,
        BoundedMatrix<double,9,9>& rHeightMatrix,
        const ElementData& rData,
        const array_1d<double,3>& rN,
        const BoundedMatrix<double,3,2>& rDN_DX);

    void ComputeGradientMatrix(
        BoundedMatrix<double,9,9>& rMatrix,
        const ElementData& rData,
        const array_1d<double,3>& rN,
        const BoundedMatrix<double,3,2>& rDN_DX);

    void ComputeDiffusionMatrix(
        BoundedMatrix<double,9,9>& rMatrix,
        const ElementData& rData,
        const BoundedMatrix<double,3,2>& rDN_DX,
        const BoundedMatrix<double,2,2>& rK1,
        const BoundedMatrix<double,2,2>& rK2,
        const BoundedMatrix<double,2,2>& rKh);

    void ComputeGradientVector(
        array_1d<double,9>& rVector,
        const ElementData& rData,
        const array_1d<double,3>& rN,
        const BoundedMatrix<double,3,2>& rDN_DX);

    void ComputeCrossWindDiffusivityTensors(
        BoundedMatrix<double,2,2>& rK1,
        BoundedMatrix<double,2,2>& rK2,
        BoundedMatrix<double,2,2>& rKh,
        const ElementData& rData,
        const BoundedMatrix<double,3,2>& rDN_DX);

    void AlgebraicResidual(
        array_1d<double,3>& rFlowResidual,
        double& rHeightresidual,
        BoundedMatrix<double,3,3> rFlowGrad,
        array_1d<double,3> rHeightGrad,
        const ElementData& rData,
        const BoundedMatrix<double,3,2>& rDN_DX);

    void StreamLineTensor(
        BoundedMatrix<double,2,2>& rTensor,
        const array_1d<double,3>& rVector);

    void CrossWindTensor(
        BoundedMatrix<double,2,2>& rTensor,
        const array_1d<double,3>& rVeector);

    double StabilizationParameter(const ElementData& rData);

    double WetFraction(double Height, double Epsilon);

    array_1d<double,3> CharacteristicLength(const ElementData& rData);

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

}; // Class ShallowWater2D3

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_SHALLOW_WATER_2D_3_H_INCLUDED  defined
