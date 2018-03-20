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
//#include "custom_elements/U_Pw_element.hpp"  <-- Modificar per fer un element basic
//#include "custom_utilities/element_utilities.hpp"
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
    //using UPwElement<TDim,TNumNodes>::mThisIntegrationMethod;

    
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
    
    //void CalculateRightHandSide(VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo ) override;
    
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);
            
    void GetValuesVector(Vector& rValues, int Step = 0);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;
    
    // void GetValueOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3>>& rValues, const ProcessInfo& rCurrentProcessInfo) override;
    
    // void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;
    
    // void CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3>>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;
    
    // void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:
        
    struct ElementVariables
    {
        ///Properties variables
        double alpha;
        array_1d<double,TDim> HVector;
        boost::numeric::ublas::bounded_matrix<double,TDim,TDim> AlphaMatrix;

        ///ProcessInfo variables
    
        ///Nodal variables
        array_1d<double,TNumNodes> NodalPhi;
        array_1d<double,TNumNodes> NodalQSource;
        array_1d<array_1d<double,TDim>, TNumNodes> NodalVel;
        
        ///Variables computed at each GP
        double IntegrationCoefficient;
        double QSource;
        array_1d<double,TNumNodes> N;
        array_1d<double,TDim> VelInter;
        boost::numeric::ublas::bounded_matrix<double,TNumNodes,TDim> GradNT;

        //Auxiliary
        boost::numeric::ublas::bounded_matrix<double,TNumNodes,TDim> AdvMatrixAux;
        boost::numeric::ublas::bounded_matrix<double,TNumNodes,TDim> DifMatrixAux;
        boost::numeric::ublas::bounded_matrix<double,TDim,TDim> FICMatrixAuxOne;
        boost::numeric::ublas::bounded_matrix<double,TDim,TNumNodes> FICMatrixAuxTwo;
        boost::numeric::ublas::bounded_matrix<double,TNumNodes,TNumNodes> AdvMatrixAuxTwo;
        boost::numeric::ublas::bounded_matrix<double,TNumNodes,TNumNodes> DifMatrixAuxTwo;
        boost::numeric::ublas::bounded_matrix<double,TNumNodes,TNumNodes> FICMatrixAuxThree;

        ///Constitutive Law parameters
        //Vector StrainVector;
        //Vector StressVector;
        //Matrix ConstitutiveMatrix;
        //Vector Np;
        //Matrix GradNpT;
        //Matrix F;
        //double detF;

    };
    
    /// Member Variables
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------       
    
    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo );

    void InitializeElementVariables(ElementVariables& rVariables, const GeometryType& Geom, const PropertiesType& Prop, const ProcessInfo& CurrentProcessInfo);

    void CalculateHVector(ElementVariables& rVariables);

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
    
    void CalculateAndAddFICForce(VectorType& rRightHandSideVector, ElementVariables& rVariables);
    

    void CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const double& detJ, const double& weight);

    
    //void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo );


///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    /// Serialization
    
    friend class Serializer;
    
    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    }

    virtual void load(Serializer& rSerializer)
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
