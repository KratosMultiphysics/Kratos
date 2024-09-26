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
#include "custom_elements/steady_state_Pw_piping_element.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "utilities/math_utils.h"
#include <cmath>

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer SteadyStatePwPipingElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                                     NodesArrayType const& ThisNodes,
                                                                     PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SteadyStatePwPipingElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer SteadyStatePwPipingElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                     GeometryType::Pointer pGeom,
                                                                     PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SteadyStatePwPipingElement(NewId, pGeom, pProperties));
}

template <unsigned int TDim, unsigned int TNumNodes>
int SteadyStatePwPipingElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    int ierr = TransientPwLineElement<TDim, TNumNodes>::Check(rCurrentProcessInfo);
    if (ierr != 0) return ierr;

    KRATOS_TRY
    const PropertiesType& rProp = this->GetProperties();
    // Verify properties
    if (rProp.Has(PIPE_ETA) == false || rProp[PIPE_ETA] < 0.0)
        KRATOS_ERROR << "PIPE_ETA has Key zero, is not defined or has an "
                        "invalid value at element "
                     << this->Id() << std::endl;

    if (rProp.Has(PIPE_THETA) == false || rProp[PIPE_THETA] < 0.0)
        KRATOS_ERROR << "PIPE_THETA has Key zero, is not defined or has an "
                        "invalid value at element "
                     << this->Id() << std::endl;

    if (rProp.Has(PIPE_D_70) == false || rProp[PIPE_D_70] < 0.0)
        KRATOS_ERROR << "PIPE_D_70 has Key zero, is not defined or has an "
                        "invalid value at element "
                     << this->Id() << std::endl;

    if (rProp.Has(PIPE_START_ELEMENT) == false)
        KRATOS_ERROR << "PIPE_START_ELEMENT has Key zero, is not defined or "
                        "has an invalid value at element "
                     << this->Id() << std::endl;

    if (rProp.Has(PIPE_MODIFIED_D) == false)
        KRATOS_ERROR << "PIPE_MODIFIED_D has Key zero, is not defined or has "
                        "an invalid value at element "
                     << this->Id() << std::endl;

    if (rProp.Has(PIPE_MODEL_FACTOR) == false)
        KRATOS_ERROR << "PIPE_MODEL_FACTOR has Key zero, is not defined or has "
                        "an invalid value at element "
                     << this->Id() << std::endl;
    KRATOS_CATCH("");

    return ierr;
}

template <unsigned int TDim, unsigned int TNumNodes>
void SteadyStatePwPipingElement<TDim, TNumNodes>::SetValueAtElement( Variable<double>& rVariable,
                                                                     double& rValue)
{
    // Copies properties
    Properties::Pointer p_new_prop = Kratos::make_shared<Properties>(this->GetProperties());
    // Adds new properties to the element
    p_new_prop->SetValue(rVariable, rValue);
    this->SetProperties(p_new_prop);
}

template <unsigned int TDim, unsigned int TNumNodes>
BoundedMatrix<double, TNumNodes, TNumNodes> SteadyStatePwPipingElement<TDim, TNumNodes>::CalculatePermeabilityMatrix(
    const GeometryType::ShapeFunctionsGradientsType& rShapeFunctionGradients,
    const Vector&                                    rIntegrationCoefficients,
    const ProcessInfo&                               rCurrentProcessInfo) const
{

    const auto& r_properties = this->GetProperties();
    RetentionLaw::Parameters RetentionParameters(r_properties, rCurrentProcessInfo);
    BoundedMatrix<double, 1, 1> constitutive_matrix;
    GeoElementUtilities::FillPermeabilityMatrix(constitutive_matrix, r_properties);

    auto result = BoundedMatrix<double, TNumNodes, TNumNodes>{ZeroMatrix{TNumNodes, TNumNodes}};
    for (unsigned int integration_point_index = 0;
         integration_point_index < this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
         ++integration_point_index) {
        const double RelativePermeability = mRetentionLawVector[integration_point_index]->CalculateRelativePermeability(RetentionParameters);
        double dynamic_viscosity_inverse = 1.0 / r_properties[DYNAMIC_VISCOSITY];
        result += GeoTransportEquationUtilities::CalculatePermeabilityMatrix<TDim, TNumNodes>(
            rShapeFunctionGradients[integration_point_index], dynamic_viscosity_inverse, constitutive_matrix,
            RelativePermeability * r_properties[PIPE_HEIGHT], rIntegrationCoefficients[integration_point_index]);
    }
    return result;
}

template <unsigned int TDim, unsigned int TNumNodes>
void SteadyStatePwPipingElement<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    TransientPwLineElement<TDim, TNumNodes>::Initialize(rCurrentProcessInfo);

    this->CalculateLength(this->GetGeometry());

    // initialse pipe parameters if not initalised, (important for staged analysis.)
    if (!this->pipe_initialised) {
        this->pipe_initialised = true;
        this->SetValue(PIPE_EROSION, false);

        // initialise pipe height with a small value
        double smallPipeHeight = 1e-10;
        this->SetValue(PIPE_HEIGHT, smallPipeHeight);
        this->SetValue(PREV_PIPE_HEIGHT, smallPipeHeight);
        this->SetValue(DIFF_PIPE_HEIGHT, 0);

        const PropertiesType& rProp = this->GetProperties();
        this->SetValue(PERMEABILITY_XX, rProp[PERMEABILITY_XX]);

        this->SetValue(PIPE_ACTIVE, false);
        //this->Set(ACTIVE, false);
    }

    //KRATOS_INFO("SteadyStatePwPipingElement") << "Initialised Piping Element" << std::endl;
    //KRATOS_INFO("SteadyStatePwPipingElement") << "Pipe Id: " << this->Id() << std::endl;
    //KRATOS_INFO("SteadyStatePwPipingElement") << "Pipe Length: " << this->GetValue(PIPE_ELEMENT_LENGTH) << std::endl;
    //KRATOS_INFO("SteadyStatePwPipingElement") << "Pipe Height: " << this->GetValue(PIPE_HEIGHT) << std::endl;
    //KRATOS_INFO("SteadyStatePwPipingElement") << "Pipe Active: " << this->Is(ACTIVE) << std::endl;
    //KRATOS_INFO("SteadyStatePwPipingElement") << "Pipe Permeability: " << this->GetValue(PERMEABILITY_XX) << std::endl;


    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void SteadyStatePwPipingElement<TDim, TNumNodes>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    TransientPwLineElement<TDim, TNumNodes>::InitializeSolutionStep(rCurrentProcessInfo);

    // todo This doesnt seem right
    //auto onElement = this->GetProperties().GetValue(PERMEABILITY_XX);
    //KRATOS_INFO("SteadyStatePwPipingElement")<< "Resetting Permeability: " << onElement << " to " << this->GetValue(PERMEABILITY_XX) << std::endl;
    //this->SetValueAtElement(PERMEABILITY_XX, onElement);


}
template <unsigned int TDim, unsigned int TNumNodes>
void SteadyStatePwPipingElement<TDim, TNumNodes>::CalculateLength(const GeometryType& Geom)
{
    // currently length is only calculated in x direction
    KRATOS_TRY
    this->SetValue(PIPE_ELEMENT_LENGTH, std::abs(Geom.GetPoint(1)[0] - Geom.GetPoint(0)[0]));
    KRATOS_INFO("Element length") << this->GetValue(PIPE_ELEMENT_LENGTH) << std::endl;
    KRATOS_CATCH("")
}




//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void SteadyStatePwPipingElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<bool>& rVariable, std::vector<bool>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if (rVariable == PIPE_ACTIVE) {
        const GeometryType& Geom = this->GetGeometry();
        const unsigned int OutputGPoints = Geom.IntegrationPointsNumber(this->GetIntegrationMethod());
        if (rValues.size() != OutputGPoints) rValues.resize(OutputGPoints);

        bool pipe_active = this->GetValue(rVariable);
        for (unsigned int GPoint = 0; GPoint < OutputGPoints; ++GPoint) {
            if (pipe_active == false) {
                rValues[GPoint] = 0;
            } else {
                rValues[GPoint] = 1;
            }
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void SteadyStatePwPipingElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if (rVariable == PIPE_HEIGHT) {
        const GeometryType& Geom = this->GetGeometry();
        const unsigned int OutputGPoints = Geom.IntegrationPointsNumber(this->GetIntegrationMethod());
        if (rValues.size() != OutputGPoints) rValues.resize(OutputGPoints);

        double pipe_height = this->GetValue(rVariable);
        for (unsigned int GPoint = 0; GPoint < OutputGPoints; ++GPoint) {
            rValues[GPoint] = pipe_height;
        }
    } else if (rVariable == PERMEABILITY_XX) {
        const GeometryType& Geom = this->GetGeometry();
        const unsigned int OutputGPoints = Geom.IntegrationPointsNumber(this->GetIntegrationMethod());
        if (rValues.size() != OutputGPoints) rValues.resize(OutputGPoints);

        double permeability = this->GetValue(rVariable);
        for (unsigned int GPoint = 0; GPoint < OutputGPoints; ++GPoint) {
            rValues[GPoint] = permeability;
        }
    } else {
        TransientPwLineElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    }
    KRATOS_CATCH("")
}


template <unsigned int TDim, unsigned int TNumNodes>
double SteadyStatePwPipingElement<TDim, TNumNodes>::CalculateHeadGradient(const PropertiesType& Prop,
                                                               const GeometryType&   Geom,
                                                               double                dx)
{
    const auto nodalHead = GeoElementUtilities::CalculateNodalHydraulicHeadFromWaterPressures(Geom, Prop);
    auto headGradient =(std::abs(nodalHead[0] - nodalHead[1])) / dx;
    KRATOS_INFO("SteadyStatePwPipingElement") << "Head Gradient: " << headGradient << std::endl;
    KRATOS_INFO("SteadyStatePwPipingElement") << "dx: " << dx << std::endl;
    KRATOS_INFO("SteadyStatePwPipingElement") << "Head 0: " << nodalHead[0] << std::endl;
    KRATOS_INFO("SteadyStatePwPipingElement") << "Head 1: " << nodalHead[1] << std::endl;
    return headGradient;
}

/// <summary>
///  Calculate the particle diameter for the particles in the pipe. The particle diameter equals d70,
/// when the unmodified sellmeijer piping rule is used.
/// </summary>
/// <param name="Prop"></param>
/// <param name="Geom"></param>
/// <returns></returns>
template <unsigned int TDim, unsigned int TNumNodes>
double SteadyStatePwPipingElement<TDim, TNumNodes>::CalculateParticleDiameter(const PropertiesType& Prop)
{
    double diameter;

    if (Prop[PIPE_MODIFIED_D]) diameter = 2.08e-4 * pow((Prop[PIPE_D_70] / 2.08e-4), 0.4);
    else diameter = Prop[PIPE_D_70];
    return diameter;
}

/// <summary>
/// Calculates the equilibrium pipe height of a piping element according to Sellmeijers rule
/// </summary>
/// <param name="Prop"></param>
/// <param name="Geom"></param>
/// <returns></returns>
template <unsigned int TDim, unsigned int TNumNodes>
double SteadyStatePwPipingElement<TDim, TNumNodes>::CalculateEquilibriumPipeHeight(const PropertiesType& Prop,
                                                                                   const GeometryType& Geom,
                                                                                   double pipe_length)
{
    const double modelFactor  = Prop[PIPE_MODEL_FACTOR];
    const double eta          = Prop[PIPE_ETA];
    const double theta        = Prop[PIPE_THETA];
    const double SolidDensity = Prop[DENSITY_SOLID];
    const double FluidDensity = Prop[DENSITY_WATER];

    // calculate head gradient over element
    double dhdx = CalculateHeadGradient(Prop, Geom, pipe_length);
    KRATOS_INFO("SteadyStatePwPipingElement") << "Head Gradient: " << dhdx << std::endl;

    // calculate particle diameter
    double particle_d = CalculateParticleDiameter(Prop);
    KRATOS_INFO("SteadyStatePwPipingElement") << "Particle Diameter: " << particle_d << std::endl;

    // todo calculate slope of pipe, currently pipe is assumed to be horizontal
    const double pipeSlope = 0;

    // return infinite when dhdx is 0
    if (dhdx < std::numeric_limits<double>::epsilon()) {
        return 1e10;
    }

    double EquilibriumPipeHeight = modelFactor * Globals::Pi / 3.0 * particle_d * (SolidDensity / FluidDensity - 1) * eta *
                                   std::sin(MathUtils<>::DegreesToRadians(theta + pipeSlope)) /
                                   std::cos(MathUtils<>::DegreesToRadians(theta)) / dhdx;



    KRATOS_INFO("SteadyStatePwPipingElement") << "Equilibrium Pipe Height: " << EquilibriumPipeHeight << std::endl;

    return EquilibriumPipeHeight;
}

template <unsigned int TDim, unsigned int TNumNodes>
bool SteadyStatePwPipingElement<TDim, TNumNodes>::InEquilibrium(const PropertiesType& Prop, const GeometryType& Geom)
{
    // Calculation if Element in Equilibrium
    // double pipeEquilibriumPipeHeight = CalculateEquilibriumPipeHeight(Prop, Geom);

    // Logic if in equilibrium

    const bool inEquilibrium = false;
    return inEquilibrium;
}

template class SteadyStatePwPipingElement<2, 2>;
template class SteadyStatePwPipingElement<2, 3>;
template class SteadyStatePwPipingElement<2, 4>;
template class SteadyStatePwPipingElement<2, 5>;
template class SteadyStatePwPipingElement<3, 2>;
template class SteadyStatePwPipingElement<3, 3>;
} // Namespace Kratos
