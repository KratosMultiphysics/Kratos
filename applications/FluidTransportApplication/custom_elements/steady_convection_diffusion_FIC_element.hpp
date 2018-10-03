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
#include <cmath>
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"
#include "includes/convection_diffusion_settings.h"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/element_size_calculator.h"

// Application includes
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

    virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo ) override;

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
        double AlphaVBar;
        double AlphaV;
        double AlphaR;
        double lv;
        double lsc;
        double OmegaV;
        double SigmaV;
        double LambdaV;
        double XiV;
        double Residual;
        double Beta;
        double NormGradPhi;
        double absorption;
        double DifSC;
        double AuxDiffusion;
        double CosinusNormals;
        double CosinusGradPhi;
        double LowTolerance;
        double HighTolerance;

        //Transient variables
        double TransientAbsorption;
        double TransientResidual;

        int IterationNumber;

        array_1d<double,TDim> HVector;
        array_1d<double,TDim> HvVector;
        array_1d<double,TDim> HrVector;
        array_1d<double,TDim> HscVector;
        array_1d<double,TDim> GradPhi;
        array_1d<double,TDim> FICVectorAuxOne;

        BoundedMatrix<double,TDim,TDim> DifMatrix;
        BoundedMatrix<double,TDim,TDim> DifMatrixK;
        BoundedMatrix<double,TDim,TDim> DifMatrixV;
        BoundedMatrix<double,TDim,TDim> DifMatrixS;
        BoundedMatrix<double,TDim,TDim> DifMatrixR;
        BoundedMatrix<double,TDim,TDim> DifMatrixSC;

        BoundedMatrix<double,TDim,TDim> IdentityMatrix;

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
        array_1d<double,TDim> VelInterHat;
        BoundedMatrix<double,TNumNodes,TDim> GradNT;

        //Auxiliary
        BoundedMatrix<double,TNumNodes,TDim> AdvMatrixAux;
        BoundedMatrix<double,TNumNodes,TDim> AbpMatrixAux;
        BoundedMatrix<double,TNumNodes,TDim> DifMatrixAux;
        BoundedMatrix<double,TNumNodes,TDim> MatrixAux;
        BoundedMatrix<double,TDim,TNumNodes> FICMatrixAuxOne;
        BoundedMatrix<double,TNumNodes,TNumNodes> AdvMatrixAuxTwo;
        BoundedMatrix<double,TNumNodes,TNumNodes> DifMatrixAuxTwo;
        BoundedMatrix<double,TNumNodes,TNumNodes> FICMatrixAuxTwo;

    };

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo );

    virtual void InitializeElementVariables(ElementVariables& rVariables, const GeometryType& Geom, const PropertiesType& Prop, const ProcessInfo& CurrentProcessInfo);

    void CalculateNormalsAngle(ElementVariables& rVariables);

    void CalculateBoundaryLv(ElementVariables& rVariables);

    virtual void CalculateDiffusivityVariables(ElementVariables& rVariables, const PropertiesType& Prop, const ProcessInfo& CurrentProcessInfo);

    virtual void CalculateHVector(ElementVariables& rVariables, const PropertiesType& Prop, const ProcessInfo& CurrentProcessInfo);

    void CalculatePeclet(ElementVariables& rVariables, const Geometry<Node<3> >& rGeom, const double& NormVel, const ProcessInfo& CurrentProcessInfo,
                                                    const PropertiesType& Prop);

    void CalculateFICBeta(ElementVariables& rVariables);

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3>>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3>>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;


    virtual void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddAdvectionMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddDiffusiveMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddAbsorptionMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddFICMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);



    void CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddRHSAdvection(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddRHSDiffusive(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddRHSAbsorption(VectorType& rRightHandSideVector, ElementVariables& rVariables);

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
