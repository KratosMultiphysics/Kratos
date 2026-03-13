//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicoló Antonelli
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "iga_application_variables.h"


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

class StokesElement : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of StokesElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(StokesElement);

    typedef Element BaseType;

    /// Type for shape function values container
    typedef Kratos::Vector ShapeFunctionsType;

    /// Type for shape function derivatives container
    typedef Kratos::Matrix ShapeDerivativesType;

    /// Type for geometry calls
    typedef Kratos::Geometry< Node > GeometryType;
    // static constexpr std::size_t NumNodes = TDim + 1;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    StokesElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry);

    StokesElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~StokesElement();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override;


    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    StokesElement() : Element()
    {
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{
    
    IntegrationMethod GetIntegrationMethod() const override;

    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void GetFirstDerivativesVector(Vector &rValues, int Step) const override;

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

    template< class TVariableType >
    void EvaluateInPoint(TVariableType& rResult,
                         const Kratos::Variable<TVariableType>& Var,
                         const ShapeFunctionsType& rShapeFunc,
                         GeometryType& rGeom)
    {
        rResult = rShapeFunc[0] * rGeom[0].FastGetSolutionStepValue(Var);

        for(SizeType i = 1; i < rShapeFunc.size(); i++)
        {
            rResult += rShapeFunc[i] * rGeom[i].FastGetSolutionStepValue(Var);
        }
    }

    virtual void CalculateTau(double MuEffective);

    double ElementSize();

    void AddMomentumTerms(MatrixType &rLHS,
            VectorType &rRHS,
            const array_1d<double,3>& rBodyForce,
            const double TauTwo,
            const ShapeFunctionsType &rN,
            const ShapeDerivativesType &rDN_DX,
            const double Weigth,
            const Matrix& rD,
            Vector& rStressVector);

    void AddContinuityTerms(MatrixType &rLHS,
            VectorType &rRHS,
            const array_1d<double,3>& rBodyForce,
            const double TauOne,
            const ShapeFunctionsType &rN,
            const ShapeDerivativesType &rDN_DX,
            const double Weight);

    void AddSecondOrderStabilizationTerms(MatrixType &rLHS,
            VectorType &rRHS,
            const array_1d<double,3>& rBodyForce,
            const double TauOne,
            const ShapeFunctionsType &rN,
            const ShapeDerivativesType &rDN_DX,
            const double Weight,
            const Matrix& rD);
    
    Vector CalculateStressAtIntegrationPoint(const ProcessInfo& rCurrentProcessInfo);
    
    void InitializeMaterial();

    ConstitutiveLaw::Pointer mpConstitutiveLaw; /// The pointer containing the constitutive laws


    ///@}
    ///@name Protected member Variables
    ///@{

    /**
     * Internal variables used in the kinematic calculations
     */
    struct ConstitutiveVariables
    {
        ConstitutiveLaw::StrainVectorType StrainVector;
        ConstitutiveLaw::StressVectorType StressVector;
        ConstitutiveLaw::VoigtSizeMatrixType D;

        /**
         * The default constructor
         * @param StrainSize The size of the strain vector in Voigt notation
         */
        ConstitutiveVariables(const SizeType StrainSize)
        {
            if (StrainVector.size() != StrainSize)
                StrainVector.resize(StrainSize);

            if (StressVector.size() != StrainSize)
                StressVector.resize(StrainSize);

            if (D.size1() != StrainSize || D.size2() != StrainSize)
                D.resize(StrainSize, StrainSize);

            noalias(StrainVector) = ZeroVector(StrainSize);
            noalias(StressVector) = ZeroVector(StrainSize);
            noalias(D)            = ZeroMatrix(StrainSize, StrainSize);
        }
    };


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{

    // Protected default constructor necessary for serialization

    ///@}

private:
    ///@name Static Member Variables

    // Dimension
    IndexType mDim;

    // Basis function order
    IndexType mBasisFunctionsOrder;

    /// Order of integration.
    Element::IntegrationMethod mIntegrationMethod;

    /// Shape function derivatives at each integration point.
    std::vector< ShapeDerivativesType > mDN_DX;

    /// Integration weights at each integration point (as fraction of the element's area).
    std::vector< double > mGaussWeight;

    void ApplyConstitutiveLaw(
        const Matrix& rB, 
        ConstitutiveLaw::Parameters& rValues,
        ConstitutiveVariables& rConstitutiveVariables) const;

    void CalculateB(
        Matrix& rB,
        const ShapeDerivativesType& r_DN_DX) const;
    
    void CalculateBDerivativeDx(
        Matrix& BDerivativeDx,
        const ShapeDerivativesType& r_DDN_DDX) const;
    
    void CalculateBDerivativeDy(
        Matrix& BDerivativeDy,
        const ShapeDerivativesType& r_DDN_DDX) const;

    void CalculateBDerivativeDx3D(
        Matrix& BDerivativeDx,
        const ShapeDerivativesType& r_DDN_DDX) const;

    void CalculateBDerivativeDy3D(
        Matrix& BDerivativeDy,
        const ShapeDerivativesType& r_DDN_DDX) const;

    void CalculateBDerivativeDz3D(
        Matrix& BDerivativeDz,
        const ShapeDerivativesType& r_DDN_DDX) const;

    void GetSolutionCoefficientVector(
        Vector& rValues) const;

    ///@}
    ///@name Member Variables
    ///@{

    double mTauOne = 0; ///< Stabilization parameter for momentum equation
    double mTauTwo = 0; ///< Stabilization parameter for continuity equation

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
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    // std::vector<std::size_t> GetSurrogateFacesIds();

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

}; // Class StokesElement

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}  // namespace Kratos.
