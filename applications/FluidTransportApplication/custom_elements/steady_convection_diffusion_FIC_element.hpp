//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Albert Puigferrat Perez
//                   Ignasi de Pouplana
//

#if !defined(KRATOS_STEADY_CONVECTION_DIFFUSION_FIC_ELEMENT_H_INCLUDED )
#define  KRATOS_STEADY_CONVECTION_DIFFUSION_FIC_ELEMENT_H_INCLUDED

// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"
#include "includes/convection_diffusion_settings.h"

// Application includes
// TODO: Create base element
#include "fluid_transport_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class SteadyConvectionDiffusionFICElement : public Element
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( SteadyConvectionDiffusionFICElement );

    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;


///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    SteadyConvectionDiffusionFICElement(IndexType NewId = 0) : Element( NewId ) {}

    /// Constructor using an array of nodes
    SteadyConvectionDiffusionFICElement(IndexType NewId, const NodesArrayType& ThisNodes) : Element(NewId, ThisNodes) {}

    /// Constructor using Geometry
    SteadyConvectionDiffusionFICElement(IndexType NewId, GeometryType::Pointer pGeometry) : Element(NewId, pGeometry) {}

    /// Constructor using Properties
    SteadyConvectionDiffusionFICElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : Element( NewId, pGeometry, pProperties ) {}

    /// Destructor
    virtual ~SteadyConvectionDiffusionFICElement() {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief GetIntegrationMethod Return the integration order to be used.
     * @return Gauss Order
     */
    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

    void GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo ) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo ) override;

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    struct ElementVariables
    {
        ///Properties variables
        double rho_dot_c;
        double Peclet;
        double AlphaV;
        double lv;

        array_1d<double,TDim> HVector;
        BoundedMatrix<double,TDim,TDim> DifMatrixK;

        ///ProcessInfo variables

        ///Nodal variables
        array_1d<double,TNumNodes> NodalPhi;
        array_1d<double,TNumNodes> NodalQSource;
        array_1d<array_1d<double,3>, TNumNodes> NodalVel;

        ///Variables computed at each GP
        double IntegrationCoefficient;
        double QSource;
        array_1d<double,TNumNodes> N;
        array_1d<double,TDim> VelInter;
        BoundedMatrix<double,TNumNodes,TDim> GradNT;

        //Auxiliary
        BoundedMatrix<double,TNumNodes,TDim> AdvMatrixAux;
        BoundedMatrix<double,TNumNodes,TDim> DifMatrixAux;
        BoundedMatrix<double,TDim,TDim> FICMatrixAuxOne;
        BoundedMatrix<double,TDim,TNumNodes> FICMatrixAuxTwo;
        BoundedMatrix<double,TNumNodes,TNumNodes> AdvMatrixAuxTwo;
        BoundedMatrix<double,TNumNodes,TNumNodes> DifMatrixAuxTwo;
        BoundedMatrix<double,TNumNodes,TNumNodes> FICMatrixAuxThree;


    };

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo );

    void InitializeElementVariables(ElementVariables& rVariables, const GeometryType& Geom, const PropertiesType& Prop, const ProcessInfo& CurrentProcessInfo);

    void CalculateHVector(ElementVariables& rVariables, const PropertiesType& Prop, const ProcessInfo& CurrentProcessInfo);


    double ProjectedElementSize(const Geometry<Node<3> >& rGeometry, const array_1d<double,3>& rVelocity);

    double AverageElementSize(const Geometry<Node<3> >& rGeometry);

    void InterpolateVariableWithComponents(array_1d<double,TDim>& rVector,const Matrix& Ncontainer,
                                        const array_1d<array_1d<double,TDim>, TNumNodes>& VariableWithComponents,const unsigned int& GPoint);


    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddAdvectionMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddDiffusiveMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddFICMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);



    void CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddRHSAdvection(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddRHSDiffusive(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddRHSFIC(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddSourceForce(VectorType& rRightHandSideVector, ElementVariables& rVariables);


    void CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const double& detJ, const double& weight);


    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo );


///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    }

    /// Assignment operator.
    SteadyConvectionDiffusionFICElement & operator=(SteadyConvectionDiffusionFICElement const& rOther);

    /// Copy constructor.
    SteadyConvectionDiffusionFICElement(SteadyConvectionDiffusionFICElement const& rOther);

}; // Class SteadyConvectionDiffusionFICElement

} // namespace Kratos

#endif // KRATOS_STEADY_CONVECTION_DIFFUSION_FIC_ELEMENT_H_INCLUDED  defined
