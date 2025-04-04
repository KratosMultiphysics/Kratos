//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela, Riccardo Rossi
//

#if !defined(KRATOS_SPALART_ALLMARAS_H_INCLUDED )
#define  KRATOS_SPALART_ALLMARAS_H_INCLUDED

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/element.h"
#include "geometries/geometry_data.h"


namespace Kratos
{
///@addtogroup FluidDynamicsApplication
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

/// Short class definition.
/** Detail class definition.
  */
class SpalartAllmaras : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SpalartAllmaras
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SpalartAllmaras);

    /// Type for shape function values container
    typedef Kratos::Vector ShapeFunctionsType;

    /// Type for shape function derivatives container
    typedef Kratos::Matrix ShapeDerivativesType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SpalartAllmaras(IndexType NewId = 0) :
        Element(NewId),
        mIntegrationMethod(GeometryData::IntegrationMethod::GI_GAUSS_1)
    {}

    /// Constructor using a Geometry instance
    SpalartAllmaras(IndexType NewId, GeometryType::Pointer pGeometry) :
        Element(NewId, pGeometry),
        mIntegrationMethod(pGeometry->GetDefaultIntegrationMethod())
    {}

    /// Constructor using geometry and properties
    SpalartAllmaras(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
        Element(NewId, pGeometry, pProperties),
        mIntegrationMethod(pGeometry->GetDefaultIntegrationMethod())
    {}

    /// Additional constructor using a specific quadrature
    SpalartAllmaras(IndexType NewId, GeometryType::Pointer pGeometry, const Element::IntegrationMethod& ThisIntegrationMethod) :
        Element(NewId, pGeometry),
        mIntegrationMethod(ThisIntegrationMethod)
    {}

    /// Additional constructor using a specific quadrature
    SpalartAllmaras(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, const Element::IntegrationMethod& ThisIntegrationMethod) :
        Element(NewId,pGeometry,pProperties),
        mIntegrationMethod(ThisIntegrationMethod)
    {}

    /// Destructor.
    ~SpalartAllmaras() override
    {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Create a new SpalartAllmaras element and return a pointer to it
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties)  const override;

    /// Create a new SpalartAllmaras element and return a pointer to it
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties)  const override;

    /// Check that all required data containers are properly initialized and registered in Kratos
    /** @return 0 if no errors are detected.
      */
    int Check(const ProcessInfo &rCurrentProcessInfo) const override;

    /// Calculate Shape function derivatives for the element
    void Initialize(const ProcessInfo &rCurrentProcessInfo) override;

    /// Calculates the projection term for stabilization
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

//    /// Compute projection of convective term for stabilization
//    virtual void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo);

    /// Evaluate the elemental contribution to the problem for turbulent viscosity.
    /**
     * @param rLeftHandSideMatrix Elemental left hand side matrix
     * @param rRightHandSideVector Elemental right hand side vector
     * @param rCurrentProcessInfo Reference to the ProcessInfo from the ModelPart containing the element
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;
//    {
//        KRATOS_THROW_ERROR(std::logic_error, "SplartAllmaras::CalculateRightHandSide method not implemented", "");
//    }

    /// Fill given array with containing the element's degrees of freedom
    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const override;

    /// Fill given vector with the linear system row index for the element's degrees of freedom
    void EquationIdVector(Element::EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    /// Fill given vector with nodal values of the problem variable (TURBULENT_VISCOSITY)
    void GetValuesVector(Vector& rValues, int Step = 0) const override;

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
        std::stringstream buffer;
        buffer << "SpalartAllmaras" << this->GetGeometry().WorkingSpaceDimension() << "D #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SpalartAllmaras" << this->GetGeometry().WorkingSpaceDimension() << "D #" << Id() << std::endl;
        rOStream << "Number of Nodes: " << this->GetGeometry().PointsNumber() << std::endl;
        rOStream << "Integration method: " << static_cast<int>(this->mIntegrationMethod);
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        this->PrintInfo(rOStream);
        rOStream << "Geometry Data: " << std::endl;
        this->GetGeometry().PrintData(rOStream);
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

    void InitializeElementData();

    void AddMassTerm(MatrixType &rMassMatrix,
                     const ShapeFunctionsType& N,
                     const double Weight);

    void AddConvection(MatrixType& rLHS,
                       const ShapeFunctionsType& N,
                       const Vector& UGradN,
                       const double Weight);

    void AddModelTerms(MatrixType& rLHS,
                       const double MolecularViscosity,
                       const double LastEddyViscosity,
                       const array_1d<double,3>& rLastEddyViscosityGradient,
                       const double Distance,
                       const array_1d<double,3>& rVelocity,
                       const ShapeFunctionsType& N,
                       const ShapeDerivativesType& DN_DX,
                       const double Weight);

    void AddStabilization(MatrixType& rLHS,
                          VectorType& rRHS,
                          const double Tau,
                          const ShapeFunctionsType& N,
                          const ShapeDerivativesType& DN_DX,
                          const Vector& UGradN,
                          const double Weight);

    void EvaluateConvection(Vector& rResult,
                            const array_1d<double,3>& rConvVel,
                            const ShapeDerivativesType& DN_DX);

    template< class TVariableType >
    void EvaluateInPoint(TVariableType& rResult,
                         const Kratos::Variable<TVariableType>& Var,
                         const ShapeFunctionsType& rShapeFunc)
    {
        rResult = rShapeFunc[0] * this->GetGeometry()[0].FastGetSolutionStepValue(Var);

        for(SizeType i = 1; i < rShapeFunc.size(); i++)
        {
            rResult += rShapeFunc[i] * this->GetGeometry()[i].FastGetSolutionStepValue(Var);
        }
    }

    double MeasureOfVorticity(const ShapeDerivativesType& DN_DX);

    void VelocityGradientNorms(double& rNormS,
                               double& rNormOmega,
                               const ShapeDerivativesType& DN_DX);

    double CalculateTau(double ElementSize, const ProcessInfo& rCurrentProcessInfo);

    double ElementSize();


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

    /// Order of integration for the element.
    Element::IntegrationMethod mIntegrationMethod;

    /// Shape function derivatives at each integration point
    std::vector< ShapeDerivativesType > mDN_DX;

    /// Determinant of the Jacobian.
    /** Used to calculate the integration weight at each integration point.
      * Note that we store it once as we are assuming straight-edged elements,
      * otherwise, the Jacobian is not constant and should be evaluated at each Gauss Point.
      */
    double mDetJ;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_TRY;
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
        //enum IntegrationMethod {GI_GAUSS_1, GI_GAUSS_2, GI_GAUSS_3, GI_GAUSS_4, GI_GAUSS_5,NumberOfIntegrationMethods};
        unsigned int IntMethod = 0;
        switch(mIntegrationMethod)
        {
        case GeometryData::IntegrationMethod::GI_GAUSS_1:
            IntMethod = 1;
            break;
        case GeometryData::IntegrationMethod::GI_GAUSS_2:
            IntMethod = 2;
            break;
        case GeometryData::IntegrationMethod::GI_GAUSS_3:
            IntMethod = 3;
            break;
        case GeometryData::IntegrationMethod::GI_GAUSS_4:
            IntMethod = 4;
            break;
        case GeometryData::IntegrationMethod::GI_GAUSS_5:
            IntMethod = 5;
            break;
        default:
            KRATOS_ERROR << "Unknown integration method encountered on serializer save for SpalartAllmaras element: " << static_cast<int>(mIntegrationMethod) << std::endl;
            break;
        }
        rSerializer.save("IntMethod",IntMethod);
        rSerializer.save("mDN_DX",mDN_DX);
        rSerializer.save("mDetJ",mDetJ);
        KRATOS_CATCH("");
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_TRY;
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer,Element);
        unsigned int IntMethod = 0;
        rSerializer.load("IntMethod",IntMethod);
        switch(mIntegrationMethod)
        {
        case GeometryData::IntegrationMethod::GI_GAUSS_1:
            mIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
            break;
        case GeometryData::IntegrationMethod::GI_GAUSS_2:
            mIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
            break;
        case GeometryData::IntegrationMethod::GI_GAUSS_3:
            mIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_3;
            break;
        case GeometryData::IntegrationMethod::GI_GAUSS_4:
            mIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_4;
            break;
        case GeometryData::IntegrationMethod::GI_GAUSS_5:
            mIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
            break;
        default:
            KRATOS_ERROR << "Unknown integration method encountered on serializer load for SpalartAllmaras element: " << static_cast<int>(IntMethod) << std::endl;
            break;
        }
        rSerializer.load("mDN_DX",mDN_DX);
        rSerializer.load("mDetJ",mDetJ);
        KRATOS_CATCH("");
    }

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

    /// Assignment operator.
    SpalartAllmaras& operator=(SpalartAllmaras const& rOther);

    /// Copy constructor.
    SpalartAllmaras(SpalartAllmaras const& rOther);


    ///@}

}; // Class SpalartAllmaras

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  SpalartAllmaras& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SpalartAllmaras& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SPALART_ALLMARAS_H_INCLUDED  defined
