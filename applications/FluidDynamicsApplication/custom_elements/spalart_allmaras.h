/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

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
    KRATOS_CLASS_POINTER_DEFINITION(SpalartAllmaras);

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
        mIntegrationMethod(GeometryData::GI_GAUSS_1)
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
    virtual ~SpalartAllmaras()
    {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Create a new SpalartAllmaras element and return a pointer to it
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties)  const;

    /// Check that all required data containers are properly initialized and registered in Kratos
    /** @return 0 if no errors are detected.
      */
    virtual int Check(const ProcessInfo &rCurrentProcessInfo);

    /// Calculate Shape function derivatives for the element
    virtual void Initialize();

    /// Calculates the projection term for stabilization
    virtual void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);

//    /// Compute projection of convective term for stabilization
//    virtual void InitializeNonLinearIteration(ProcessInfo &CurrentProcessInfo);

    /// Evaluate the elemental contribution to the problem for turbulent viscosity.
    /**
     * @param rLeftHandSideMatrix Elemental left hand side matrix
     * @param rRightHandSideVector Elemental right hand side vector
     * @param rCurrentProcessInfo Reference to the ProcessInfo from the ModelPart containg the element
     */
    virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR(std::logic_error, "SplartAllmaras::CalculateRightHandSide method not implemented", "");
    }

    /// Fill given array with containing the element's degrees of freedom
    virtual void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo);

    /// Fill given vector with the linear system row index for the element's degrees of freedom
    virtual void EquationIdVector(Element::EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    /// Fill given vector with nodal values of the problem variable (TURBULENT_VISCOSITY)
    virtual void GetValuesVector(Vector& rValues, int Step = 0);

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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "SpalartAllmaras" << this->GetGeometry().WorkingSpaceDimension() << "D #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "SpalartAllmaras" << this->GetGeometry().WorkingSpaceDimension() << "D #" << Id() << std::endl;
        rOStream << "Number of Nodes: " << this->GetGeometry().PointsNumber() << std::endl;
        rOStream << "Integration method: " << this->mIntegrationMethod;
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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

    void AddMassTerm(MatrixType &rMassMatrix,
                     const ShapeFunctionsType N,
                     const double Weight);

    void AddConvection(MatrixType& rLHS,
                      const ShapeFunctionsType& N,
                      const Vector& UGradN,
                      const double Weight);

    void AddModelTerms(MatrixType& rLHS,
                       const double MolecularViscosity,
                       const double LastEddyViscosity,
                       const array_1d<double,3> rLastEddyViscosityGradient,
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
                         const Kratos::Variable<TVariableType> Var,
                         const ShapeFunctionsType& rShapeFunc)
    {
        rResult = rShapeFunc[0] * this->GetGeometry()[0].FastGetSolutionStepValue(Var);

        for(SizeType i = 1; i < rShapeFunc.size(); i++)
        {
            rResult += rShapeFunc[i] * this->GetGeometry()[i].FastGetSolutionStepValue(Var);
        }
    }

    double MeasureOfVorticity(const ShapeDerivativesType& DN_DX);


    double CalculateTau(const ProcessInfo& rCurrentProcessInfo);


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

    /// Total elemental area or volume
    double mElementSize;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_TRY;
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
        enum IntegrationMethod {GI_GAUSS_1, GI_GAUSS_2, GI_GAUSS_3, GI_GAUSS_4, GI_GAUSS_5,
                                NumberOfIntegrationMethods
                           };
        unsigned int IntMethod = 0;
        switch(mIntegrationMethod)
        {
        case GeometryData::GI_GAUSS_1:
            IntMethod = 1;
            break;
        case GeometryData::GI_GAUSS_2:
            IntMethod = 2;
            break;
        case GeometryData::GI_GAUSS_3:
            IntMethod = 3;
            break;
        case GeometryData::GI_GAUSS_4:
            IntMethod = 4;
            break;
        case GeometryData::GI_GAUSS_5:
            IntMethod = 5;
            break;
        default:
            KRATOS_ERROR(std::invalid_argument,"Unknown integration method encountered on serializer save for SpalartAllmaras element: ",mIntegrationMethod);
            break;
        }
        rSerializer.save("IntMethod",IntMethod);
        rSerializer.save("mDN_DX",mDN_DX);
        rSerializer.save("mElementSize",mElementSize);
        KRATOS_CATCH("");
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_TRY;
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer,Element);
        unsigned int IntMethod = 0;
        rSerializer.load("IntMethod",IntMethod);
        switch(mIntegrationMethod)
        {
        case 1:
            mIntegrationMethod = GeometryData::GI_GAUSS_1;
            break;
        case 2:
            mIntegrationMethod = GeometryData::GI_GAUSS_2;
            break;
        case 3:
            mIntegrationMethod = GeometryData::GI_GAUSS_3;
            break;
        case 4:
            mIntegrationMethod = GeometryData::GI_GAUSS_4;
            break;
        case 5:
            mIntegrationMethod = GeometryData::GI_GAUSS_5;
            break;
        default:
            KRATOS_ERROR(std::invalid_argument,"Unknown integration method encountered on serializer load for SpalartAllmaras element: ",IntMethod);
            break;
        }
        rSerializer.load("mDN_DX",mDN_DX);
        rSerializer.load("mElementSize",mElementSize);
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


