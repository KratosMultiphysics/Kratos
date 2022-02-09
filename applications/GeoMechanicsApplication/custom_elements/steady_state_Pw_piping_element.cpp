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
    int ierr = SteadyStatePwInterfaceElement::Check(rCurrentProcessInfo);
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
    const GeometryType& Geom = this->GetGeometry();
    this->CalculateLength(Geom);
    //this->gravity = norm_2(rVariables.InterfaceVariables.VolumeAcceleration);
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
    PipingElementVariables PipingVariables;
    
    this->InitializeElementVariables(PipingVariables,
                                     Geom,
                                     Prop,
                                     CurrentProcessInfo);
    InterfaceElementVariables Variables = PipingVariables.InterfaceVariables;
   // this->CalculateLength(Geom);
    // VG: TODO
    // Perhaps a new parameter to get join width and not minimum joint width
    Variables.JointWidth = Prop[MINIMUM_JOINT_WIDTH];

    double pipe_height =  this->CalculateEquilibriumPipeHeight(PipingVariables, Geom);

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
double SteadyStatePwPipingElement<2, 4>::CalculateWaterPressureGradient(InterfaceElementVariables& rVariables)
{
    return abs(rVariables.PressureVector[1] - rVariables.PressureVector[0]) / this->Length;
}
template< >
double SteadyStatePwPipingElement<3, 6>::CalculateWaterPressureGradient(InterfaceElementVariables& rVariables)
{
    return 0;
}
template< >
double SteadyStatePwPipingElement<3, 8>::CalculateWaterPressureGradient(InterfaceElementVariables& rVariables)
{
    return 0;
}


template< unsigned int TDim, unsigned int TNumNodes >
double SteadyStatePwPipingElement<TDim,TNumNodes>::
    CalculateEquilibriumPipeHeight(PipingElementVariables& rVariables, const GeometryType Geom)
{

    // todo add modelFactor input and calculate slope of pipe
    const double modelFactor = 1;
    const double pipeSlope = 0;
    const double dhdx = 1;
 

    double dpdx = CalculateWaterPressureGradient(rVariables.InterfaceVariables);
    
    if (dpdx < DBL_EPSILON)
    { 
        return 1e10;
    }
    //noalias(NodeVolumeAcceleration) = rGeom[node].FastGetSolutionStepValue(VOLUME_ACCELERATION, 0);
    //const double g = norm_2(NodeVolumeAcceleration);

    array_1d<double, 3> gravity_array = rVariables.InterfaceVariables.VolumeAcceleration;
    const double gravity = norm_2(gravity_array);
    //double length = geom.Length();
    return modelFactor * M_PI / 3.0 * rVariables.d70  * (rVariables.InterfaceVariables.SolidDensity - rVariables.InterfaceVariables.FluidDensity) * gravity * rVariables.eta  * sin((rVariables.theta  + pipeSlope) * M_PI / 180.0) / cos(rVariables.theta * M_PI / 180.0) / dpdx;

}


template< unsigned int TDim, unsigned int TNumNodes >
void SteadyStatePwPipingElement<TDim,TNumNodes>::
    InitializeElementVariables( PipingElementVariables& rVariables,
                                const GeometryType& Geom,
                                const PropertiesType& Prop,
                                const ProcessInfo& CurrentProcessInfo )
{
    SteadyStatePwInterfaceElement::InitializeElementVariables(rVariables.InterfaceVariables, Geom, Prop, CurrentProcessInfo);

    rVariables.d70                = Prop[PIPE_D_70];
    rVariables.eta                = Prop[PIPE_ETA];
    rVariables.theta              = Prop[PIPE_THETA];
}
    


template class SteadyStatePwPipingElement<2,4>;
template class SteadyStatePwPipingElement<3,6>;
template class SteadyStatePwPipingElement<3,8>;

} // Namespace Kratos
