#include <cmath>
#include <tuple>
#include <numeric>

#include "utilities/parallel_utilities.h"
#include "processes/calculate_nodal_area_process.h"

#include "fluid_dynamics_application_variables.h"
#include "shock_capturing_entropy_viscosity_process.h"

#define ENTROPY YOUNG_MODULUS

namespace Kratos {


ShockCapturingEntropyViscosityProcess::ShockCapturingEntropyViscosityProcess(
    ModelPart& rModelPart,
    Parameters rParameters)
    : Process()
    , mrModelPart(rModelPart)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mComputeAreasEveryStep = rParameters["calculate_nodal_area_at_each_step"].GetBool();
    mTunableConstant = rParameters["tunable_constant"].GetDouble();
    mTunableConstantMax = rParameters["tunable_constant_max"].GetDouble();
    mArtificialBulkViscosityPrandtl = rParameters["artificial_bulk_viscosity_Prandtl"].GetDouble();
    mArtificialConductivityPrandtl = rParameters["artificial_conductivity_Prandtl"].GetDouble();

    KRATOS_CATCH("")
}



const Parameters ShockCapturingEntropyViscosityProcess::GetDefaultParameters() const
{
    return Parameters(R"(
    {
        "model_part_name" : "",
        "calculate_nodal_area_at_each_step" : false,
        "tunable_constant"       : 0.0,
        "tunable_constant_max"   : 0.0,
        "artificial_bulk_viscosity_Prandtl"   : 0.1,
        "artificial_conductivity_Prandtl"     : 0.1
    }
    )");
}


double ShockCapturingEntropyViscosityProcess::TotalDerivativeUtil::
Divergence(const Matrix& rShapeFunGradients, const Matrix& rNodalValues)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(rShapeFunGradients.size1() != rNodalValues.size1())
        << "Matrices rShapeFunGradients and rNodalValues must be the same size (Size1: " 
        << rShapeFunGradients.size1() << " != " << rNodalValues.size1() <<")" << std::endl;
    KRATOS_DEBUG_ERROR_IF(rShapeFunGradients.size2() != rNodalValues.size2())
        << "Matrices rShapeFunGradients and rNodalValues must be the same size (Size2: " 
        << rShapeFunGradients.size2() << " != " << rNodalValues.size2() <<")" << std::endl; 

    double divergence = 0.0;
    for(std::size_t i=0; i<rShapeFunGradients.size1(); ++i)
    {
        for(std::size_t j=0; j<rShapeFunGradients.size2(); ++j)
        {
            divergence += rShapeFunGradients(i,j) * rNodalValues(i,j);
        }
    }
    return divergence;

    KRATOS_CATCH("")
}

void ShockCapturingEntropyViscosityProcess::ExecuteBeforeSolutionLoop()
{
    UpdateNodalAreaProcess();
}


void ShockCapturingEntropyViscosityProcess::ExecuteFinalizeSolutionStep()
{

    if (mComputeAreasEveryStep)
    {
        UpdateNodalAreaProcess();
    }

    ComputeNodalEntropies();
    ComputeArtificialMagnitudes();
}


void ShockCapturingEntropyViscosityProcess::UpdateNodalAreaProcess()
{
    CalculateNodalAreaProcess<false> nodal_area_process(mrModelPart);
    nodal_area_process.Execute();
}


void ShockCapturingEntropyViscosityProcess::ComputeNodalEntropies()
{
    const double heat_capacity_ratio = mrModelPart.GetProcessInfo().GetValue(HEAT_CAPACITY_RATIO);
    
    block_for_each(mrModelPart.Nodes(), [heat_capacity_ratio](NodeType& r_node)
    {
        const double density = r_node.GetValue(DENSITY);
        const double pressure = r_node.GetValue(PRESSURE);

        const auto entropy = ComputeEntropy(density, pressure, heat_capacity_ratio);
        r_node.FastGetSolutionStepValue(ENTROPY) = entropy;

        r_node.SetValue(ARTIFICIAL_DYNAMIC_VISCOSITY, 0.0);
        r_node.SetValue(ARTIFICIAL_CONDUCTIVITY, 0.0);
        r_node.SetValue(ARTIFICIAL_BULK_VISCOSITY, 0.0);
    });
}

void ShockCapturingEntropyViscosityProcess::ComputeArtificialMagnitudes()
{
    KRATOS_TRY

    const double delta_time = mrModelPart.GetProcessInfo().GetValue(DELTA_TIME);
    const double heat_capacity_ratio = mrModelPart.GetProcessInfo().GetValue(HEAT_CAPACITY_RATIO);

    block_for_each(mrModelPart.Elements(), [&](Element& r_element)
    {
        const auto inf_norm = ComputeElementalInfNormData(r_element, delta_time, heat_capacity_ratio);
        const double h2 = ComputeHSquared(r_element);

        const double mu_e = mTunableConstant * h2 * inf_norm.Density * inf_norm.EntropyResidual;
        const double mu_max = mTunableConstantMax * std::sqrt(h2) * inf_norm.Density * inf_norm.TotalVelocity;

        const double mu_h = std::min(mu_e, mu_max);
        const double mu_bulk = mArtificialBulkViscosityPrandtl * mu_h / inf_norm.Density;
        const double kappa = mArtificialConductivityPrandtl * mu_h / (heat_capacity_ratio - 1.0);

        DistributeVariablesToNodes(r_element, mu_h, mu_bulk, kappa);
    });

    KRATOS_CATCH("")
}

void ShockCapturingEntropyViscosityProcess::DistributeVariablesToNodes(
    Element& rElement,
    const double ArtificialDynamicViscosity,
    const double ArtificialBulkViscosity,
    const double ArtificialConductivity) const
{
    auto& r_geometry = rElement.GetGeometry();
    const double element_volume = r_geometry.LocalSpaceDimension() == 3 ? r_geometry.Volume() : r_geometry.Area();
    
    for(unsigned int i=0; i<r_geometry.size(); ++i)
    {
        auto& r_node = r_geometry[i];
        const double weight = element_volume / r_node.GetValue(NODAL_AREA);
        r_node.SetLock();
        r_node.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY) += weight * ArtificialDynamicViscosity;
        r_node.GetValue(ARTIFICIAL_BULK_VISCOSITY) += weight * ArtificialBulkViscosity;
        r_node.GetValue(ARTIFICIAL_CONDUCTIVITY) += weight * ArtificialConductivity;
        r_node.UnSetLock();
    }
}


double ShockCapturingEntropyViscosityProcess::ComputeEntropy(
    const double Density,
    const double Pressure,
    const double Gamma)
{
    return Density / (Gamma - 1.0) * std::log(Pressure / std::pow(Density, Gamma));
}


double ShockCapturingEntropyViscosityProcess::ComputeHSquared(const Element& rElement)
{
    double h_squared = 0.0;
    const auto& r_geometry = rElement.GetGeometry();

    // H is the shortest edge of the element
    for(unsigned int i=0; i<r_geometry.size(); ++i)
    {
        const unsigned j = (i + 1) % r_geometry.size();
        const auto edge = r_geometry[j] - r_geometry[i];
        h_squared = std::min(h_squared, inner_prod(edge, edge));
    }

    return h_squared;
}


ShockCapturingEntropyViscosityProcess::InfNormData ShockCapturingEntropyViscosityProcess::ComputeElementalInfNormData(
    const Element& rElement,
    const double DeltaTime,
    const double HeatCapacityRatio)
{
    TotalDerivativeUtil entropy_td;
    TotalDerivativeUtil density_td;
    Vector total_velocities;

    std::tie(entropy_td, density_td, total_velocities) = BuildTotalDerivativeUtils(rElement, DeltaTime, HeatCapacityRatio);
    return ComputeInfNorms(rElement.GetGeometry(), entropy_td, density_td, total_velocities);
}

std::tuple<ShockCapturingEntropyViscosityProcess::TotalDerivativeUtil, ShockCapturingEntropyViscosityProcess::TotalDerivativeUtil, Vector>
ShockCapturingEntropyViscosityProcess::BuildTotalDerivativeUtils(const Element& rElement, const double DeltaTime, const double HeatCapacityRatio)
{
    KRATOS_TRY

    const auto& r_geometry = rElement.GetGeometry();

    // Loading nodal values
    TotalDerivativeUtil entropy_total_derivative(r_geometry.size(), r_geometry.LocalSpaceDimension());
    TotalDerivativeUtil density_total_derivative(r_geometry.size(), r_geometry.LocalSpaceDimension());
    Vector total_velocities(r_geometry.size(), 0.0);

    for(unsigned int i=0; i<r_geometry.size(); ++i)
    {
        const auto velocity = r_geometry[i].FastGetSolutionStepValue(VELOCITY);
        const auto temperature = r_geometry[i].FastGetSolutionStepValue(TEMPERATURE);
        
        total_velocities[i] = norm_2(velocity) + std::sqrt(HeatCapacityRatio * temperature);

        const auto& r_node = r_geometry[i];
        entropy_total_derivative.LoadNodalValues(ENTROPY, r_node, i, velocity, DeltaTime);
        density_total_derivative.LoadNodalValues(DENSITY, r_node, i, velocity, DeltaTime);

        
        entropy_total_derivative.Value[i] /= density_total_derivative.Value[i];
    }

    // Computing inf norm over gauss points
    return std::tie(entropy_total_derivative, density_total_derivative, total_velocities);

    KRATOS_CATCH("")
}


ShockCapturingEntropyViscosityProcess::InfNormData ShockCapturingEntropyViscosityProcess::ComputeInfNorms(
    const Geometry<NodeType>& rGeometry,
    const TotalDerivativeUtil& rEntropyTotalDerivative,
    const TotalDerivativeUtil& rDensityTotalDerivative,
    const Vector& rTotalVelocities)
{
    KRATOS_TRY

    constexpr auto integration_method = GeometryData::IntegrationMethod::GI_GAUSS_1;
    const auto n_gauss_points = rGeometry.IntegrationPointsNumber(integration_method);
    
    const auto& r_shape_functions = rGeometry.ShapeFunctionsValues(integration_method);
    GeometryData::ShapeFunctionsGradientsType r_shape_functions_gradients;
    rGeometry.ShapeFunctionsIntegrationPointsGradients(r_shape_functions_gradients, integration_method);

    double max_residual = 0.0;  // max(Dh1, Dh2) in the paper
    double max_density = 0.0;
    double max_total_velocity = 0.0;

    for(unsigned int g = 0; g<n_gauss_points; ++g)
    {
        const auto N = row(r_shape_functions, g);
        const auto& G = r_shape_functions_gradients[g];

        const double entropy_residual = rEntropyTotalDerivative.ComputeAtGaussPoint(N, G);
        const double density_residual = rDensityTotalDerivative.ComputeAtGaussPoint(N, G);

        const double specific_entropy_residual = inner_prod(rEntropyTotalDerivative.Value, N) * density_residual;

        const double max_gp_residual = std::max(entropy_residual, specific_entropy_residual);
        const double density = inner_prod(rDensityTotalDerivative.Value, N);
        const double total_velocity = inner_prod(rTotalVelocities, N);
        
        max_residual =  std::max(max_residual, max_gp_residual);
        max_density = std::max(max_density, density);
        max_total_velocity = std::max(max_total_velocity, total_velocity);
    }

    return {max_residual, max_density, max_total_velocity};

    KRATOS_CATCH("")
}


} // namespace Kratos 