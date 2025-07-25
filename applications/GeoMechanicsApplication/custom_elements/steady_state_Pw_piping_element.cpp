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
#include "custom_utilities/transport_equation_utilities.hpp"
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
        NewId, this->GetGeometry().Create(ThisNodes), pProperties,
        this->GetStressStatePolicy().Clone(), this->CloneIntegrationCoefficientModifier()));
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer SteadyStatePwPipingElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                     GeometryType::Pointer pGeom,
                                                                     PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SteadyStatePwPipingElement(NewId, pGeom, pProperties,
                                                           this->GetStressStatePolicy().Clone(),
                                                           this->CloneIntegrationCoefficientModifier()));
}

template <unsigned int TDim, unsigned int TNumNodes>
int SteadyStatePwPipingElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    int ierr = SteadyStatePwInterfaceElement<TDim, TNumNodes>::Check(rCurrentProcessInfo);
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
    KRATOS_CATCH("")

    return ierr;
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string SteadyStatePwPipingElement<TDim, TNumNodes>::Info() const
{
    return "SteadyStatePwPipingElement";
}

template <unsigned int TDim, unsigned int TNumNodes>
void SteadyStatePwPipingElement<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    SteadyStatePwInterfaceElement<TDim, TNumNodes>::Initialize(rCurrentProcessInfo);

    this->CalculateLength(this->GetGeometry());

    // initialse pipe parameters if not initalised, (important for staged analysis.
    if (!this->pipe_initialised) {
        this->pipe_initialised = true;
        this->SetValue(PIPE_EROSION, false);

        // initialise pipe height with a small value
        double smallPipeHeight = 1e-10;
        this->SetValue(PIPE_HEIGHT, smallPipeHeight);
        this->SetValue(PREV_PIPE_HEIGHT, smallPipeHeight);
        this->SetValue(DIFF_PIPE_HEIGHT, 0);

        this->SetValue(PIPE_ACTIVE, false);
    }

    KRATOS_CATCH("")
}

template <>
void SteadyStatePwPipingElement<2, 4>::CalculateLength(const GeometryType& Geom)
{
    // currently length is only calculated in x direction
    KRATOS_TRY
    this->SetValue(PIPE_ELEMENT_LENGTH, std::abs(Geom.GetPoint(1)[0] - Geom.GetPoint(0)[0]));
    KRATOS_CATCH("")
}

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
            rValues[GPoint] = pipe_active;
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
    } else {
        TransientPwInterfaceElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
            rVariable, rValues, rCurrentProcessInfo);
    }
    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void SteadyStatePwPipingElement<TDim, TNumNodes>::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                                               VectorType& rRightHandSideVector,
                                                               const ProcessInfo& CurrentProcessInfo,
                                                               bool CalculateStiffnessMatrixFlag,
                                                               bool CalculateResidualVectorFlag)
{
    KRATOS_TRY
    // Previous definitions
    const PropertiesType&                           Prop = this->GetProperties();
    const GeometryType&                             Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        Geom.IntegrationPoints(mThisIntegrationMethod);
    const unsigned int NumGPoints = IntegrationPoints.size();

    // Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues(mThisIntegrationMethod);
    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer =
        Geom.ShapeFunctionsLocalGradients(mThisIntegrationMethod);
    GeometryType::JacobiansType JContainer(NumGPoints);
    Geom.Jacobian(JContainer, mThisIntegrationMethod);
    Vector detJContainer(NumGPoints);
    Geom.DeterminantOfJacobian(detJContainer, mThisIntegrationMethod);

    // Element variables
    InterfaceElementVariables Variables;

    this->InitializeElementVariables(Variables, Geom, Prop, CurrentProcessInfo);

    // Set joint width as pipe height
    Variables.JointWidth = this->GetValue(PIPE_HEIGHT);

    // Auxiliary variables
    array_1d<double, TDim> RelDispVector;
    SFGradAuxVariables     SFGradAuxVars;

    RetentionLaw::Parameters RetentionParameters(this->GetProperties());

    const auto integration_coefficients =
        this->CalculateIntegrationCoefficients(IntegrationPoints, detJContainer);

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        // Compute Np, StrainVector, JointWidth, GradNpT
        noalias(Variables.Np) = row(NContainer, GPoint);

        this->template CalculateShapeFunctionsGradients<Matrix>(
            Variables.GradNpT, SFGradAuxVars, JContainer[GPoint], Variables.RotationMatrix,
            DN_DeContainer[GPoint], NContainer, Variables.JointWidth, GPoint);

        // Compute BodyAcceleration and Permeability Matrix
        GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
            Variables.BodyAcceleration, NContainer, Variables.VolumeAcceleration, GPoint);

        InterfaceElementUtilities::FillPermeabilityMatrix(
            Variables.LocalPermeabilityMatrix, Variables.JointWidth, Prop[TRANSVERSAL_PERMEABILITY]);

        this->CalculateRetentionResponse(Variables, RetentionParameters, GPoint);

        Variables.IntegrationCoefficient = integration_coefficients[GPoint];

        // Contributions to the left hand side
        if (CalculateStiffnessMatrixFlag) this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        // Contributions to the right hand side
        if (CalculateResidualVectorFlag)
            this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);
    }

    KRATOS_CATCH("")
}

template <>
double SteadyStatePwPipingElement<2, 4>::CalculateHeadGradient(const PropertiesType& Prop,
                                                               const GeometryType&   Geom,
                                                               double                dx)
{
    const auto nodalHead = GeoElementUtilities::CalculateNodalHydraulicHeadFromWaterPressures(Geom, Prop);
    return std::abs((nodalHead[3] + nodalHead[0]) / 2 - (nodalHead[2] + nodalHead[1]) / 2) / dx;
}

/// <summary>
/// Calculates the equilibrium pipe height of a piping element according to Sellmeijers rule
/// </summary>
/// <param name="rProperties"></param>
/// <param name="rGeometry"></param>
/// <returns></returns>
template <unsigned int TDim, unsigned int TNumNodes>
double SteadyStatePwPipingElement<TDim, TNumNodes>::CalculateEquilibriumPipeHeight(
    const PropertiesType& rProperties, const GeometryType& rGeometry, double PipeLength)
{
    // calculate head gradient over element
    const auto dhdx = CalculateHeadGradient(rProperties, rGeometry, PipeLength);

    // return infinite when dhdx is 0
    if (dhdx < std::numeric_limits<double>::epsilon()) return 1e10;

    // todo calculate slope of pipe, currently pipe is assumed to be horizontal
    constexpr auto pipe_slope = 0.0;

    return rProperties[PIPE_MODEL_FACTOR] * Globals::Pi / 3.0 *
           GeoTransportEquationUtilities::CalculateParticleDiameter(rProperties) *
           (rProperties[DENSITY_SOLID] / rProperties[DENSITY_WATER] - 1) * rProperties[PIPE_ETA] *
           std::sin(MathUtils<>::DegreesToRadians(rProperties[PIPE_THETA] + pipe_slope)) /
           std::cos(MathUtils<>::DegreesToRadians(rProperties[PIPE_THETA])) / dhdx;
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

template class SteadyStatePwPipingElement<2, 4>;
} // Namespace Kratos
