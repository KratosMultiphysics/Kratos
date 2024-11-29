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

#if !defined(KRATOS_STOKES_ELEMENT_H_INCLUDED )
#define  KRATOS_STOKES_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"


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

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

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
    
    void CalculateOnIntegrationPoints(
        const Variable<Matrix>& rVariable,  
        std::vector<Matrix>& D_constitutive_matrix, 
        const ProcessInfo& rCurrentProcessInfo) override;   

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

    virtual void CalculateTau(double& TauOne,
                              double& TauTwo,
                              const double DynViscosity,
                              const double h_element_size,
                              double mu_effective);

    double ElementSize();

    void AddMomentumTerms(MatrixType &rLHS,
            VectorType &rRHS,
            const double Density,
            const double Viscosity,
            const array_1d<double,3>& BodyForce,
            const double TauTwo,
            const ShapeFunctionsType &N,
            const ShapeDerivativesType &DN_DX,
            const double Weigth,
            const Matrix& r_D,
            Vector& r_stress_vector);

    void AddContinuityTerms(MatrixType &rLHS,
            VectorType &rRHS,
            const double Density,
            const array_1d<double,3>& BodyForce,
            const double TauOne,
            const ShapeFunctionsType &N,
            const ShapeDerivativesType &DN_DX,
            const double Weight);

    void AddSecondOrderStabilizationTerms(MatrixType &rLHS,
            VectorType &rRHS,
            const double Density,
            const array_1d<double,3>& BodyForce,
            const double TauOne,
            const ShapeFunctionsType &N,
            const ShapeDerivativesType &DN_DX,
            const double Weight,
            const Matrix& r_D,
            Vector& r_stress_vector);
    
    Vector CalculateStressAtIntegrationPoint(const ProcessInfo& rCurrentProcessInfo);

    Matrix CalculateConstitutiveMatrixAtIntegrationPoint(const ProcessInfo& rCurrentProcessInfo);
    
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
    int mDim;

    // Basis function order
    IndexType mBasisFunctionsOrder;

    /// Order of integration.
    Element::IntegrationMethod mIntegrationMethod;

    /// Shape function derivatives at each integration point.
    std::vector< ShapeDerivativesType > mDN_DX;

    /// Integration weights at each integration point (as fraction of the element's area).
    std::vector< double > mGaussWeight;

    void CalculateB(
        Matrix& rB,
        const ShapeDerivativesType& r_DN_DX) const;

    void CalculateI(
        const ShapeFunctionsType &N,
        Matrix& rI) const;
    
    void CalculateB_second_order(
        Matrix& B_second_order,
        const ShapeDerivativesType& r_DDN_DDX) const;
    
    void CalculateB_derivative_x(
        Matrix& B_derivative_x,
        const ShapeDerivativesType& r_DDN_DDX) const;
    
    void CalculateB_derivative_y(
        Matrix& B_derivative_y,
        const ShapeDerivativesType& r_DDN_DDX) const;

    void GetValuesVector(
        Vector& rValues) const;

    ///@}
    ///@name Member Variables
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

    /// Assignment operator.
    //StokesElement& operator=(const StokesElement& rOther);

    /// Copy constructor.
    //StokesElement(const StokesElement& rOther);

    ///@}
    Parameters ReadParamatersFile(
        const std::string& rDataFileName) const;

}; // Class StokesElement

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    StokesElement& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const StokesElement& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);
      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_LAPLACIAN_SHIFTED_BOUNDARY_ELEMENT_H_INCLUDED  defined