// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//

// Application includes
#define _USE_MATH_DEFINES
#include "custom_elements/steady_state_Pw_piping_element.hpp"
#include <math.h>
namespace Kratos
{


template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer SteadyStatePwPipingElement<TDim,TNumNodes>::
    Create( IndexType NewId,
            NodesArrayType const& ThisNodes,
            PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new SteadyStatePwPipingElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer SteadyStatePwPipingElement<TDim,TNumNodes>::
    Create(IndexType NewId,
           GeometryType::Pointer pGeom,
           PropertiesType::Pointer pProperties) const
{
    return Element::Pointer( new SteadyStatePwPipingElement( NewId, pGeom, pProperties ) );
}
template< unsigned int TDim, unsigned int TNumNodes >
int SteadyStatePwPipingElement<TDim,TNumNodes>::
    Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY
    int ierr = SteadyStatePwInterfaceElement<TDim, TNumNodes>::Check(rCurrentProcessInfo);
    if (ierr != 0) return ierr;
    // todo check piping parameters
    return ierr;
    KRATOS_CATCH( "" );

}

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyStatePwPipingElement<TDim, TNumNodes>::
Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
        // KRATOS_INFO("0-UPwSmallStrainInterfaceElement::Initialize()") << std::endl;

        SteadyStatePwInterfaceElement<TDim, TNumNodes>::Initialize(rCurrentProcessInfo);
    //const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    this->CalculateLength(Geom);

    this->Set(ACTIVE, false);

    KRATOS_CATCH("");
}


template< >
void SteadyStatePwPipingElement<2, 4>::CalculateLength(const GeometryType& Geom)
{
    KRATOS_TRY;
    // KRATOS_INFO("0-UPwSmallStrainInterfaceElement<2,4>:::CalculateInitialGap()") << std::endl;


    array_1d<double, 3> Vx;
    noalias(Vx) = Geom.GetPoint(1) - Geom.GetPoint(0);
    this->Length = norm_2(Vx);


    // KRATOS_INFO("1-UPwSmallStrainInterfaceElement<2,4>:::CalculateInitialGap()") << std::endl;
    KRATOS_CATCH("")
}

template< >
void SteadyStatePwPipingElement<3, 6>::CalculateLength(const GeometryType& Geom)
{
    KRATOS_ERROR << " Length of SteadyStatePwPipingElement3D6N element is not implemented" << std::endl;
}
template< >
void SteadyStatePwPipingElement<3, 8>::CalculateLength(const GeometryType& Geom)
{
    KRATOS_ERROR << " Length of SteadyStatePwPipingElement3D8N element is not implemented" << std::endl;
}


//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void SteadyStatePwPipingElement<TDim,TNumNodes>::
    CalculateAll(MatrixType& rLeftHandSideMatrix,
                 VectorType& rRightHandSideVector,
                 const ProcessInfo& CurrentProcessInfo,
                 const bool CalculateStiffnessMatrixFlag,
                 const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY
    // KRATOS_INFO("0-SteadyStatePwInterfaceElement:::CalculateAll()") << std::endl;

    //Previous definitions
    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = Geom.IntegrationPoints( mThisIntegrationMethod );
    const unsigned int NumGPoints = IntegrationPoints.size();

    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer = Geom.ShapeFunctionsLocalGradients( mThisIntegrationMethod );
    GeometryType::JacobiansType JContainer(NumGPoints);
    Geom.Jacobian( JContainer, mThisIntegrationMethod );
    Vector detJContainer(NumGPoints);
    Geom.DeterminantOfJacobian(detJContainer,mThisIntegrationMethod);

    //Element variables
    InterfaceElementVariables Variables;
    
    this->InitializeElementVariables(Variables,
                                     Geom,
                                     Prop,
                                     CurrentProcessInfo);

   // this->CalculateLength(Geom);
    // VG: TODO
    // Perhaps a new parameter to get join width and not minimum joint width
    Variables.JointWidth = Prop[MINIMUM_JOINT_WIDTH];
	
    double eq_pipe_height =  this->CalculateEquilibriumPipeHeight(Geom, Prop);

    //Auxiliary variables
    array_1d<double,TDim> RelDispVector;
    SFGradAuxVariables SFGradAuxVars;

    // create general parametes of retention law
    RetentionLaw::Parameters RetentionParameters(Geom, this->GetProperties(), CurrentProcessInfo);

    //Loop over integration points
    for ( unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        //Compute Np, StrainVector, JointWidth, GradNpT
        noalias(Variables.Np) = row(NContainer,GPoint);

        this->template 
            CalculateShapeFunctionsGradients< Matrix >(Variables.GradNpT,
                                                       SFGradAuxVars,
                                                       JContainer[GPoint],
                                                       Variables.RotationMatrix,
                                                       DN_DeContainer[GPoint],
                                                       NContainer,
                                                       Variables.JointWidth,
                                                       GPoint);

        //Compute BodyAcceleration and Permeability Matrix
        GeoElementUtilities::
            InterpolateVariableWithComponents<TDim, TNumNodes>(Variables.BodyAcceleration,
                                                                NContainer,
                                                                Variables.VolumeAcceleration,
                                                                GPoint );

        InterfaceElementUtilities::FillPermeabilityMatrix( Variables.LocalPermeabilityMatrix,
                                                           Variables.JointWidth,
                                                           Prop[TRANSVERSAL_PERMEABILITY] );

        CalculateRetentionResponse( Variables,
                                    RetentionParameters,
                                    GPoint );

        //Compute weighting coefficient for integration
        Variables.IntegrationCoefficient = 
            this->CalculateIntegrationCoefficient(IntegrationPoints,
                                                  GPoint,
                                                  detJContainer[GPoint]);

        //Contributions to the left hand side
        if (CalculateStiffnessMatrixFlag) this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        //Contributions to the right hand side
        if (CalculateResidualVectorFlag)  this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);

    }

    // KRATOS_INFO("1-SteadyStatePwInterfaceElement:::CalculateAll()") << std::endl;
    KRATOS_CATCH( "" )
}

template< >
double SteadyStatePwPipingElement<2, 4>::CalculateWaterPressureGradient(const GeometryType& Geom)
{
    return abs(Geom[1].FastGetSolutionStepValue(WATER_PRESSURE) - Geom[0].FastGetSolutionStepValue(WATER_PRESSURE)) / this->Length;
}
template< >
double SteadyStatePwPipingElement<3, 6>::CalculateWaterPressureGradient(const GeometryType& Geom)
{
    KRATOS_ERROR << " pressure gradient calculation of SteadyStatePwPipingElement3D6N element is not implemented" << std::endl;
    return 0;
}
template< >
double SteadyStatePwPipingElement<3, 8>::CalculateWaterPressureGradient(const GeometryType& rVariables)
{
    KRATOS_ERROR << " pressure gradient calculation of SteadyStatePwPipingElement3D8N element is not implemented" << std::endl;
    return 0;
}


/// <summary>
/// Calculates the equilibrium pipe height of a piping element according to Sellmeijers rule
/// </summary>
/// <param name="rVariables"></param>
/// <param name="Geom"></param>
/// <returns></returns>
template< unsigned int TDim, unsigned int TNumNodes >
double SteadyStatePwPipingElement<TDim,TNumNodes>::
    CalculateEquilibriumPipeHeight(const GeometryType& Geom, const PropertiesType& Prop)
{

    // todo add modelFactor input and calculate slope of pipe
    const double modelFactor = 1;
    const double pipeSlope = 0;
    const double d70 = Prop[PIPE_D_70];
    const double eta = Prop[PIPE_ETA];
    const double theta = Prop[PIPE_THETA];
    const double SolidDensity = Prop[DENSITY_SOLID];
    const double FluidDensity = Prop[DENSITY_WATER];
 
    // calculate pressure gradient over element
    double dpdx = CalculateWaterPressureGradient(Geom);
    
    // return infinite when dpdx is 0
    if (dpdx < DBL_EPSILON)
    { 
        return 1e10;
    }

    // gravity is taken from first node
    array_1d<double, 3> gravity_array= Geom[0].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    const double gravity = norm_2(gravity_array);
    return modelFactor * M_PI / 3.0 * d70 * (SolidDensity - FluidDensity) * gravity * eta  * sin((theta  + pipeSlope) * M_PI / 180.0) / cos(theta * M_PI / 180.0) / dpdx;

}
    
template class SteadyStatePwPipingElement<2,4>;
template class SteadyStatePwPipingElement<3,6>;
template class SteadyStatePwPipingElement<3,8>;

} // Namespace Kratos
