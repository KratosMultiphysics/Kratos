// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//

// Application includes
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
    int ierr = SteadyStatePwInterfaceElement<TDim, TNumNodes>::Check(rCurrentProcessInfo);
    if (ierr != 0) return ierr;

    KRATOS_TRY
    const PropertiesType& rProp = this->GetProperties();
    // Verify properties
    if (rProp.Has(PIPE_ETA) == false ||
        rProp[PIPE_ETA] < 0.0)
        KRATOS_ERROR << "PIPE_ETA has Key zero, is not defined or has an invalid value at element " << this->Id() << std::endl;

    if (rProp.Has(PIPE_THETA) == false ||
        rProp[PIPE_THETA] < 0.0)
        KRATOS_ERROR << "PIPE_THETA has Key zero, is not defined or has an invalid value at element " << this->Id() << std::endl;

    if (rProp.Has(PIPE_D_70) == false ||
        rProp[PIPE_D_70] < 0.0)
        KRATOS_ERROR << "PIPE_D_70 has Key zero, is not defined or has an invalid value at element " << this->Id() << std::endl;

    if (rProp.Has(PIPE_START_ELEMENT) == false)
        KRATOS_ERROR << "PIPE_START_ELEMENT has Key zero, is not defined or has an invalid value at element " << this->Id() << std::endl;

    if (rProp.Has(PIPE_MODIFIED_D) == false)
        KRATOS_ERROR << "PIPE_MODIFIED_D has Key zero, is not defined or has an invalid value at element " << this->Id() << std::endl;

    if (rProp.Has(PIPE_MODEL_FACTOR) == false)
        KRATOS_ERROR << "PIPE_MODEL_FACTOR has Key zero, is not defined or has an invalid value at element " << this->Id() << std::endl;
    KRATOS_CATCH( "" );

    return ierr;
}

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyStatePwPipingElement<TDim, TNumNodes>::
Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    SteadyStatePwInterfaceElement<TDim, TNumNodes>::Initialize(rCurrentProcessInfo);

    this->CalculateLength(this->GetGeometry());
    this->CalculateSlope(this->GetGeometry());

    // initialse pipe parameters if not initalised, (important for staged analysis. 
    if (!this->pipe_initialised)
    { 
        this->pipe_initialised = true;
        this->SetValue(PIPE_EROSION, false);

        // initialise pipe height with a small value
        double smallPipeHeight = 1e-10;
        this->SetValue(PIPE_HEIGHT, smallPipeHeight);
        this->SetValue(PREV_PIPE_HEIGHT, smallPipeHeight);
        this->SetValue(DIFF_PIPE_HEIGHT, 0);

        this->SetValue(PIPE_ACTIVE, false);
    }

    KRATOS_CATCH("");
}

template< unsigned int TDim, unsigned int TNumNodes >
array_1d<double, TNumNodes> SteadyStatePwPipingElement<TDim, TNumNodes>::HydraulicHeadsFromWaterPressures() {

    const double NumericalLimit = std::numeric_limits<double>::epsilon();
    const PropertiesType& rProp = this->GetProperties();
    const GeometryType& rGeom = this->GetGeometry();

    //Defining necessary variables
    array_1d<double, TNumNodes> NodalHydraulicHead;
    for (unsigned int node = 0; node < TNumNodes; ++node) {
        array_1d<double, 3> NodeVolumeAcceleration;
        noalias(NodeVolumeAcceleration) = rGeom[node].FastGetSolutionStepValue(VOLUME_ACCELERATION, 0);
        const double g = norm_2(NodeVolumeAcceleration);
        if (g > NumericalLimit) {
            const double FluidWeight = g * rProp[DENSITY_WATER];

            array_1d<double, 3> NodeCoordinates;
            noalias(NodeCoordinates) = rGeom[node].Coordinates();
            array_1d<double, 3> NodeVolumeAccelerationUnitVector;
            noalias(NodeVolumeAccelerationUnitVector) = NodeVolumeAcceleration / g;

            const double WaterPressure = rGeom[node].FastGetSolutionStepValue(WATER_PRESSURE);
            NodalHydraulicHead[node] = -inner_prod(NodeCoordinates, NodeVolumeAccelerationUnitVector)
                - PORE_PRESSURE_SIGN_FACTOR * WaterPressure / FluidWeight;
        }
        else {
            NodalHydraulicHead[node] = 0.0;
        }
    }
    return NodalHydraulicHead;
}

template< >
void SteadyStatePwPipingElement<2, 4>::CalculateLength(const GeometryType& Geom)
{
    KRATOS_TRY

	double dx = Geom.GetPoint(1)[0] - Geom.GetPoint(0)[0];
    double dy = Geom.GetPoint(1)[1] - Geom.GetPoint(0)[1];

	this->SetValue(PIPE_ELEMENT_LENGTH, sqrt(pow(dx, 2) + pow(dy, 2)));
    
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

template< >
void SteadyStatePwPipingElement<2, 4>::CalculateSlope(const GeometryType& Geom)
{
    KRATOS_TRY

	double dx = Geom.GetPoint(1)[0] - Geom.GetPoint(0)[0];
    double dy = Geom.GetPoint(1)[1] - Geom.GetPoint(0)[1];
    this->SetValue(PIPE_ELEMENT_SLOPE, atan(dy / dx) / (M_PI / 180.0));
    KRATOS_INFO("PipeAngle") << this->GetValue(PIPE_ELEMENT_SLOPE) << std::endl;

    KRATOS_CATCH("")
}

template< >
void SteadyStatePwPipingElement<3, 6>::CalculateSlope(const GeometryType& Geom)
{
    KRATOS_ERROR << " Slope of SteadyStatePwPipingElement3D6N element is not implemented" << std::endl;
}
template< >
void SteadyStatePwPipingElement<3, 8>::CalculateSlope(const GeometryType& Geom)
{
    KRATOS_ERROR << " Slope of SteadyStatePwPipingElement3D8N element is not implemented" << std::endl;
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyStatePwPipingElement<TDim, TNumNodes>::
CalculateOnIntegrationPoints(const Variable<bool>& rVariable,
    std::vector<bool>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // KRATOS_INFO("0-TransientPwInterfaceElement:::CalculateOnIntegrationPoints<double>()") << std::endl;

    if (rVariable == PIPE_ACTIVE)
    {
        const GeometryType& Geom = this->GetGeometry();
        const unsigned int OutputGPoints = Geom.IntegrationPointsNumber(this->GetIntegrationMethod());
        if (rValues.size() != OutputGPoints)
            rValues.resize(OutputGPoints);

        bool pipe_active = this->GetValue(rVariable);
        for (unsigned int GPoint = 0; GPoint < OutputGPoints; ++GPoint) {
            rValues[GPoint] = pipe_active;
        }
    }
    
    KRATOS_CATCH("")
}

template< unsigned int TDim, unsigned int TNumNodes >
void SteadyStatePwPipingElement<TDim, TNumNodes>::
CalculateOnIntegrationPoints(const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    const GeometryType& Geom = this->GetGeometry();
    const unsigned int OutputGPoints = Geom.IntegrationPointsNumber(this->GetIntegrationMethod());

    if (rVariable == PIPE_HEIGHT)
    {
        if (rValues.size() != OutputGPoints)
            rValues.resize(OutputGPoints);

        double pipe_height = this->GetValue(rVariable);
        for (unsigned int GPoint = 0; GPoint < OutputGPoints; ++GPoint) {
            rValues[GPoint] = pipe_height;
        }
    }
    else
    {
        TransientPwInterfaceElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }
    KRATOS_CATCH("")
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

    // Set joint width as pipe height
    Variables.JointWidth = this->GetValue(PIPE_HEIGHT);
	
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

    KRATOS_CATCH( "" )
}

template< >
double SteadyStatePwPipingElement<2, 4>::CalculateWaterPressureGradient(const GeometryType& Geom)
{
    double length = this->GetValue(PIPE_ELEMENT_LENGTH);
	return abs((Geom[3].FastGetSolutionStepValue(WATER_PRESSURE) + Geom[0].FastGetSolutionStepValue(WATER_PRESSURE))/2 
        - (Geom[2].GetSolutionStepValue(WATER_PRESSURE)+ Geom[1].GetSolutionStepValue(WATER_PRESSURE))/2) / length;

    
}
template< >
double SteadyStatePwPipingElement<3, 6>::CalculateWaterPressureGradient(const GeometryType& Geom)
{
    KRATOS_ERROR << " pressure gradient calculation of SteadyStatePwPipingElement3D6N element is not implemented" << std::endl;
}
template< >
double SteadyStatePwPipingElement<3, 8>::CalculateWaterPressureGradient(const GeometryType& Geom)
{
    KRATOS_ERROR << " pressure gradient calculation of SteadyStatePwPipingElement3D8N element is not implemented" << std::endl;
}


template< >
double SteadyStatePwPipingElement<2, 4>::CalculateHeadGradient(const GeometryType& Geom)
{
    double length = this->GetValue(PIPE_ELEMENT_LENGTH);
    array_1d<double, 4> NodalHydraulicHead = HydraulicHeadsFromWaterPressures();
    return abs((NodalHydraulicHead[3] + NodalHydraulicHead[0]) / 2
        - (NodalHydraulicHead[2] + NodalHydraulicHead[1]) / 2) / length;
}
template< >
double SteadyStatePwPipingElement<3, 6>::CalculateHeadGradient(const GeometryType& Geom)
{
    KRATOS_ERROR << " hydraulic head gradient calculation of SteadyStatePwPipingElement3D6N element is not implemented" << std::endl;
}
template< >
double SteadyStatePwPipingElement<3, 8>::CalculateHeadGradient(const GeometryType& Geom)
{
    KRATOS_ERROR << " hydraulic head gradient calculation of SteadyStatePwPipingElement3D8N element is not implemented" << std::endl;
}
/// <summary>
///  Calculate the particle diameter for the particles in the pipe. The particle diameter equals d70, 
/// when the unmodified sellmeijer piping rule is used. 
/// </summary>
/// <param name="Prop"></param>
/// <param name="Geom"></param>
/// <returns></returns>
template< unsigned int TDim, unsigned int TNumNodes >
double SteadyStatePwPipingElement<TDim, TNumNodes>::CalculateParticleDiameter(const PropertiesType& Prop)
{
    double diameter;

    if (Prop[PIPE_MODIFIED_D])
        diameter = 2.08e-4 * pow((Prop[PIPE_D_70] / 2.08e-4), 0.4);
    else
        diameter = Prop[PIPE_D_70];
    return diameter;
}

/// <summary>
/// Calculates the equilibrium pipe height of a piping element according to Sellmeijers rule
/// </summary>
/// <param name="Prop"></param>
/// <param name="Geom"></param>
/// <returns></returns>
template< unsigned int TDim, unsigned int TNumNodes >
double SteadyStatePwPipingElement<TDim,TNumNodes>:: CalculateEquilibriumPipeHeight(const PropertiesType& Prop, const GeometryType& Geom)
{
    const double modelFactor = Prop[PIPE_MODEL_FACTOR];
    const double eta = Prop[PIPE_ETA];
    const double theta = Prop[PIPE_THETA];
    const double SolidDensity = Prop[DENSITY_SOLID];
    const double FluidDensity = Prop[DENSITY_WATER];
    const double pipeSlope = this->GetValue(PIPE_ELEMENT_SLOPE);

    // calculate pressure gradient over element
    double dhdl = CalculateHeadGradient(Geom);
    // calculate particle diameter
    double particle_d = CalculateParticleDiameter(Prop);

    // return infinite when dpdxy is 0
	if (dhdl < std::numeric_limits<double>::epsilon())
    { 
        return 1e10;
    }

    // gravity is taken from first node
    array_1d<double, 3> gravity_array= Geom[0].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    const double gravity = norm_2(gravity_array);

    // KRATOS_INFO("modelFactor") << modelFactor << std::endl;
    // KRATOS_INFO("particle_d") << particle_d << std::endl;
    // KRATOS_INFO("solidDensity") << SolidDensity << std::endl;
    // KRATOS_INFO("fluidDensity") << FluidDensity << std::endl;
    // KRATOS_INFO("eta") << eta << std::endl;
    // KRATOS_INFO("theta") << theta << std::endl;
    double pHeight = this->GetValue(PIPE_HEIGHT);
    double equilibriumHeight = modelFactor * M_PI / 3.0 * particle_d * (SolidDensity / FluidDensity - 1) * eta * sin((theta + pipeSlope) * M_PI / 180.0) / cos((theta * M_PI / 180.0)) / dhdl;
	KRATOS_INFO("Equilibrium") << pipeSlope << " : " << dhdl << " : " << particle_d << " : " << equilibriumHeight << " : " << pHeight << std::endl;
    return equilibriumHeight;
	
}

template class SteadyStatePwPipingElement<2,4>;
template class SteadyStatePwPipingElement<3,6>;
template class SteadyStatePwPipingElement<3,8>;

} // Namespace Kratos
  