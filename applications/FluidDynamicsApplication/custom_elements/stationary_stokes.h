#ifndef KRATOS_STATIONARY_STOKES_H
#define KRATOS_STATIONARY_STOKES_H

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

#ifndef KRATOS_FLUID_TAYLOR_HOOD_H
#define KRATOS_FLUID_TAYLOR_HOOD_H
// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/element.h"
#include "geometries/geometry_data.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"


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

/// GLS-stabilized element for the solution of the stationary Stokes problem.
template< unsigned int TDim >
class StationaryStokes : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of StationaryStokes
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(StationaryStokes);

    /// Type for shape function values container
    typedef Kratos::Vector ShapeFunctionsType;

    /// Type for shape function derivatives container
    typedef Kratos::Matrix ShapeDerivativesType;

    /// Type for geometry calls
    typedef Kratos::Geometry< Node<3> > GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    StationaryStokes(IndexType NewId = 0) :
        Element(NewId)
    {}

    /// Constructor using a Geometry instance
    StationaryStokes(IndexType NewId, GeometryType::Pointer pGeometry) :
        Element(NewId, pGeometry)
    {}

    /// Constructor using geometry and properties
    StationaryStokes(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
        Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    ~StationaryStokes() override
    {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Create a new StationaryStokes element and return a pointer to it
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties)  const override;


    /// Create a new element of this type.
	/**
	 @param NewId Index of the new element
     @param pGeom A pointer to the geometry of the new element
	 @param pProperties Pointer to the element's properties
	 */
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    /// Check that all required data containers are properly initialized and registered in Kratos
    /** @return 0 if no errors are detected.
      */
    int Check(const ProcessInfo &rCurrentProcessInfo) override;

    /// Calculate Shape function derivatives and Jacobian at each integration point
    void Initialize() override;

    /// Evaluate the elemental contribution to the problem.
    /**
     * @param rLeftHandSideMatrix Elemental left hand side matrix
     * @param rRightHandSideVector Elemental right hand side vector
     * @param rCurrentProcessInfo Reference to the ProcessInfo from the ModelPart containg the element
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override
    {
        MatrixType TmpLHS;
        this->CalculateLocalSystem(TmpLHS,rRightHandSideVector,rCurrentProcessInfo);
    }

    /// Fill given array with containing the element's degrees of freedom
    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo) override;

    /// Fill given vector with the linear system row index for the element's degrees of freedom
    void EquationIdVector(Element::EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void GetFirstDerivativesVector(Vector &rValues, int Step = 0) override;

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
        buffer << "StationaryStokes" << this->GetGeometry().WorkingSpaceDimension() << "D #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "StationaryStokes" << this->GetGeometry().WorkingSpaceDimension() << "D #" << Id() << std::endl;
        rOStream << "Number of Nodes: " << this->GetGeometry().PointsNumber() << std::endl;
        rOStream << "Integration method: " << this->mIntegrationMethod;
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

    void AddMomentumTerms(MatrixType &rLHS,
                          VectorType &rRHS,
                          const double Density,
                          const double Viscosity,
                          const array_1d<double,3>& BodyForce,
                          const double TauTwo,
                          const ShapeFunctionsType &N,
                          const ShapeDerivativesType &DN_DX,
                          const double Weigth);

    void AddContinuityTerms(MatrixType &rLHS,
                            VectorType &rRHS,
                            const double Density,
                            const array_1d<double,3>& BodyForce,
                            const double TauOne,
                            const ShapeFunctionsType &N,
                            const ShapeDerivativesType &DN_DX,
                            const double Weight);

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
                              const double DynViscosity);

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

    /// Order of integration.
    Element::IntegrationMethod mIntegrationMethod;

    /// Shape function derivatives at each integration point.
    std::vector< ShapeDerivativesType > mDN_DX;

    /// Integration weights at each integration point (as fraction of the element's area).
    std::vector< double > mGaussWeight;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_TRY;
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );

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
            KRATOS_THROW_ERROR(std::invalid_argument,"Unknown integration method encountered on serializer save for StationaryStokes element: ",mIntegrationMethod);
            break;
        }
        rSerializer.save("IntMethod",IntMethod);
        rSerializer.save("mDN_DX",mDN_DX);
        rSerializer.save("mGaussWeight",mGaussWeight);
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
            KRATOS_THROW_ERROR(std::invalid_argument,"Unknown integration method encountered on serializer load for StationaryStokes element: ",IntMethod);
            break;
        }
        rSerializer.load("mDN_DX",mDN_DX);
        rSerializer.load("mGaussWeight",mGaussWeight);
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
    StationaryStokes& operator=(StationaryStokes const& rOther);

    /// Copy constructor.
    StationaryStokes(StationaryStokes const& rOther);


    ///@}

}; // Class StationaryStokes

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim >
inline std::istream& operator >> (std::istream& rIStream,
                                  StationaryStokes<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim >
inline std::ostream& operator << (std::ostream& rOStream,
                                  const StationaryStokes<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FLUID_TAYLOR_HOOD_H


#endif // KRATOS_STATIONARY_STOKES_H
