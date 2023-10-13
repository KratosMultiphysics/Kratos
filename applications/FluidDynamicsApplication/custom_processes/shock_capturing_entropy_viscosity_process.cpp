//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Eduard Gómez
//

#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "processes/calculate_nodal_area_process.h"

#include "fluid_dynamics_application_variables.h"
#include "shock_capturing_entropy_viscosity_process.h"


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
    mEntropyConstant = rParameters["entropy_constant"].GetDouble();     // In the article: c_e
    mEnergyConstant = rParameters["energy_constant"].GetDouble();       // In the article: c_max

    mArtificialMassDiffusivityPrandtl = rParameters["artificial_mass_viscosity_Prandtl"].GetDouble();   // In the article: p rho
    mArtificialConductivityPrandtl = rParameters["artificial_conductivity_Prandtl"].GetDouble();        // In the article: p

    KRATOS_CATCH("")
}


int ShockCapturingEntropyViscosityProcess::Check()
{
    KRATOS_TRY

    const int err_code = Process::Check();

    KRATOS_ERROR_IF_NOT(mrModelPart.GetProcessInfo().Has(DOMAIN_SIZE))
        << "Missing variable DOMAIN_SIZE in model part process info." << std::endl;
    const auto domain_size = mrModelPart.GetProcessInfo().GetValue(DOMAIN_SIZE);
    KRATOS_ERROR_IF(domain_size !=2 && domain_size !=3)
        << "ShockCapturingEntropyViscosityProcess is only implemented for 2D and 3D domains." << std::endl;


    block_for_each(mrModelPart.Nodes(), [](const NodeType& r_node)
    {
        KRATOS_ERROR_IF_NOT(r_node.SolutionStepsDataHas(NUMERICAL_ENTROPY))
            << "Missing NUMERICAL_ENTROPY variable from node #" << r_node.Id() << std::endl;

        KRATOS_ERROR_IF_NOT(r_node.SolutionStepsDataHas(DENSITY))
            << "Missing DENSITY variable from node #" << r_node.Id() << std::endl;

        KRATOS_ERROR_IF_NOT(r_node.SolutionStepsDataHas(PRESSURE))
            << "Missing PRESSURE variable from node #" << r_node.Id() << std::endl;

        KRATOS_ERROR_IF_NOT(r_node.SolutionStepsDataHas(TEMPERATURE))
            << "Missing TEMPERATURE variable from node #" << r_node.Id() << std::endl;

        KRATOS_ERROR_IF_NOT(r_node.SolutionStepsDataHas(VELOCITY))
            << "Missing VELOCITY variable from node #" << r_node.Id() << std::endl;;
    });

    if(mrModelPart.GetCommunicator().LocalMesh().NumberOfElements() != 0)
    {
        KRATOS_ERROR_IF_NOT(mrModelPart.ElementsBegin()->GetProperties().Has(HEAT_CAPACITY_RATIO))
            << "Variable HEAT_CAPACITY_RATIO missing from elemental properties." << std::endl;

#ifdef KRATOS_DEBUG
        const double first_gamma = mrModelPart.ElementsBegin()->GetProperties().GetValue(HEAT_CAPACITY_RATIO);
        static constexpr double tolerance = 1e-8;

        block_for_each(mrModelPart.Elements(), [&](Element& r_element)
        {
            const double this_gamma = r_element.GetProperties().GetValue(HEAT_CAPACITY_RATIO);
            KRATOS_ERROR_IF((first_gamma - this_gamma) > tolerance)
                << "HEAT_CAPACITY_RATIO is not constant in the domain" << std::endl;
        });
#endif
    }

    return err_code;

    KRATOS_CATCH("")
}


const Parameters ShockCapturingEntropyViscosityProcess::GetDefaultParameters() const
{
    return Parameters(R"(
    {
        "model_part_name" : "",
        "calculate_nodal_area_at_each_step" : false,
        "entropy_constant"  : 1.0,
        "energy_constant"   : 0.25,
        "artificial_mass_viscosity_Prandtl"   : 0.1,
        "artificial_conductivity_Prandtl"     : 0.1
    }
    )");
}


double ShockCapturingEntropyViscosityProcess::TotalDerivativeUtil::Divergence(
    const Matrix& rShapeFunGradients,
    const Matrix& rNodalValues)
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


void ShockCapturingEntropyViscosityProcess::ExecuteInitializeSolutionStep()
{
    if (mComputeAreasEveryStep || !mIsInitialized)
    {
        UpdateNodalAreaProcess();
    }

    if(!mIsInitialized)
    {
        ComputeNodalEntropies(1);
        /* ^ Necessary in order to compute derivative in first step.
         *   Stored in buffer index 1 to prevent it from being overwritten at the
         *      end of this same time-step.
         *   Ideally would be computed at ExecuteInitialize but initial condition
         *      processes don't run until ExecuteInitializeSolutionStep
         */
        mIsInitialized = true;
    }

    ComputeNodalEntropies();
    ComputeArtificialMagnitudes();
}


void ShockCapturingEntropyViscosityProcess::UpdateNodalAreaProcess()
{
    CalculateNodalAreaProcess<false> nodal_area_process(mrModelPart);
    nodal_area_process.Execute();
}


void ShockCapturingEntropyViscosityProcess::ComputeNodalEntropies(const unsigned int WriteBufferIndex)
{
    if(mrModelPart.GetCommunicator().LocalMesh().NumberOfElements() == 0) return; // empty mpdelpart

    const double heat_capacity_ratio = mrModelPart.ElementsBegin()->GetProperties().GetValue(HEAT_CAPACITY_RATIO);

    block_for_each(mrModelPart.Nodes(), [&](NodeType& r_node)
    {
        const double density = r_node.FastGetSolutionStepValue(DENSITY);
        const double pressure = r_node.FastGetSolutionStepValue(PRESSURE);

        const auto entropy = ComputeEntropy(density, pressure, heat_capacity_ratio);

        r_node.FastGetSolutionStepValue(NUMERICAL_ENTROPY, WriteBufferIndex) = entropy;

        r_node.SetValue(ARTIFICIAL_MASS_DIFFUSIVITY, 0.0);
        r_node.SetValue(ARTIFICIAL_DYNAMIC_VISCOSITY, 0.0);
        r_node.SetValue(ARTIFICIAL_CONDUCTIVITY, 0.0);
        r_node.SetValue(ARTIFICIAL_BULK_VISCOSITY, 0.0);
    });
}


/**
 * @brief 2D specialization. Assumes nodes are contiguous.
 */
template<>
double ShockCapturingEntropyViscosityProcess::MinimumEdgeLengthSquared<2>(
    const Element& rElement)
{
    const auto& r_geometry = rElement.GetGeometry();

    const auto edge_length_squared = 
        [&](const std::size_t j, const std::size_t k) -> double
    {
        const array_1d<double, 3> edge = r_geometry[k] - r_geometry[j];
        return inner_prod(edge, edge);
    };

    double h_squared = edge_length_squared(0, r_geometry.size()-1);

    for(unsigned int i=1; i < r_geometry.size(); ++i)
    {
        const double length_2 = edge_length_squared(i-1, i);
        h_squared = std::min(h_squared, length_2);
    }

    return h_squared;
}


/**
 * @brief 3D specialization. Works properly only for tetrahedra
 */
template<>
double ShockCapturingEntropyViscosityProcess::MinimumEdgeLengthSquared<3>(
    const Element& rElement)
{
    const auto& r_geometry = rElement.GetGeometry();

    KRATOS_WARNING_IF("ShockCapturingEntropyViscosityProcess", r_geometry.size() != 4)
        << "This process is only properly implemented in 2D (any geometry) and 3D tetrahedra" << std::endl;
    // Should also work for reasonably non-deformed hexahedra (since diagonals are counted as well)

    double h_squared = std::numeric_limits<double>::max();
    for(unsigned int i=0; i<r_geometry.size(); ++i)
    {
        for(unsigned int j=0; j<i; ++j)
        {
            const array_1d<double, 3> edge = r_geometry[j] - r_geometry[i];
            h_squared = std::min(h_squared, inner_prod(edge, edge));
        }
    }

    return h_squared;
}



void ShockCapturingEntropyViscosityProcess::ComputeArtificialMagnitudes()
{
    KRATOS_TRY

    const double delta_time = mrModelPart.GetProcessInfo().GetValue(DELTA_TIME);

    if(mrModelPart.GetCommunicator().LocalMesh().NumberOfElements() == 0) return; // Empty modelpart

    const double heat_capacity_ratio = mrModelPart.ElementsBegin()->GetProperties().GetValue(HEAT_CAPACITY_RATIO);
    const double specific_heat_c_v = mrModelPart.ElementsBegin()->GetProperties().GetValue(SPECIFIC_HEAT);

    const unsigned int ndim = mrModelPart.ElementsBegin()->GetGeometry().LocalSpaceDimension();
    
    const auto geometry_size = [&ndim]() -> std::function<double(Geometry<Node>*)>
    {
        if(ndim == 2) return [](const Geometry<Node> * const p_geom) { return p_geom->Area(); };
        if(ndim == 3) return [](const Geometry<Node> * const p_geom) { return p_geom->Volume(); };
        KRATOS_ERROR << "Invalid number of dimensions (" << ndim <<"). Only 2D and 3D are supported" << std::endl;
    }(); // The simpler "const auto var = condition ? lambda1 : lambda2;" does not compile with MSVC

    const auto minimum_edge_squared = [&]() -> std::function<double(const Element&)>
    {
        if(ndim == 2) return [](const Element& r_element) { return MinimumEdgeLengthSquared<2>(r_element); };
        if(ndim == 3) return [](const Element& r_element) { return MinimumEdgeLengthSquared<3>(r_element); };
        KRATOS_ERROR << "Invalid number of dimensions (" << ndim <<"). Only 2D and 3D are supported" << std::endl;
    }();

    block_for_each(mrModelPart.Elements(), [&](Element& r_element)
    {
        const auto inf_norm = ComputeElementalInfNormData(r_element, delta_time, heat_capacity_ratio, specific_heat_c_v);
        const double h2 = minimum_edge_squared(r_element);

        r_element.SetValue(SHOCK_SENSOR, inf_norm.EntropyResidual);

        const double mu_e = mEntropyConstant * h2 * inf_norm.Density * inf_norm.EntropyResidual;
        const double mu_max = mEnergyConstant * std::sqrt(h2) * inf_norm.Density * inf_norm.TotalVelocity;
        const double mu_h = std::min(mu_e, mu_max);

        const double mu_rho = mArtificialMassDiffusivityPrandtl * mu_h / inf_norm.Density;
        const double kappa  = mArtificialConductivityPrandtl * mu_h / (heat_capacity_ratio - 1.0);

        DistributeVariablesToNodes(r_element, mu_h, mu_rho, kappa, geometry_size);
    });

    KRATOS_CATCH("")
}


void ShockCapturingEntropyViscosityProcess::DistributeVariablesToNodes(
    Element& rElement,
    const double ArtificialDynamicViscosity,
    const double ArtificialMassDiffusivity,
    const double ArtificialConductivity,
    const std::function<double(Geometry<Node>*)>& rGeometrySize) const
{
    auto& r_geometry = rElement.GetGeometry();
    const double element_volume = rGeometrySize(&r_geometry);

    for(unsigned int i=0; i<r_geometry.size(); ++i)
    {
        auto& r_node = r_geometry[i];
        const double weight = element_volume / r_node.GetValue(NODAL_AREA);
        r_node.SetLock();
        r_node.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY) += weight * ArtificialDynamicViscosity;
        r_node.GetValue(ARTIFICIAL_MASS_DIFFUSIVITY) += weight * ArtificialMassDiffusivity;
        r_node.GetValue(ARTIFICIAL_CONDUCTIVITY) += weight * ArtificialConductivity;
        r_node.UnSetLock();
    }
}


ShockCapturingEntropyViscosityProcess::InfNormData ShockCapturingEntropyViscosityProcess::ComputeElementalInfNormData(
    const Element& rElement,
    const double DeltaTime,
    const double HeatCapacityRatio,
    const double SpecificHeatCV)
{
    TotalDerivativeUtil entropy_td;
    TotalDerivativeUtil density_td;
    Vector total_velocities;

    std::tie(entropy_td, density_td, total_velocities) = BuildTotalDerivativeUtils(rElement, DeltaTime, HeatCapacityRatio, SpecificHeatCV);
    return ComputeInfNorms(rElement.GetGeometry(), entropy_td, density_td, total_velocities);
}


std::tuple<ShockCapturingEntropyViscosityProcess::TotalDerivativeUtil, ShockCapturingEntropyViscosityProcess::TotalDerivativeUtil, Vector>
ShockCapturingEntropyViscosityProcess::BuildTotalDerivativeUtils(
    const Element& rElement,
    const double DeltaTime,
    const double HeatCapacityRatio,
    const double SpecificHeatCV)
{
    KRATOS_TRY

    const auto& r_geometry = rElement.GetGeometry();

    // Loading nodal values
    TotalDerivativeUtil entropy_total_derivative(r_geometry.LocalSpaceDimension(), r_geometry.size());
    TotalDerivativeUtil density_total_derivative(r_geometry.LocalSpaceDimension(), r_geometry.size());
    Vector total_velocities(r_geometry.size(), 0.0);

    for(unsigned int i=0; i<r_geometry.size(); ++i)
    {
        const auto& r_node = r_geometry[i];

        const auto& velocity = r_node.FastGetSolutionStepValue(VELOCITY);
        const double temperature = r_node.FastGetSolutionStepValue(TEMPERATURE);

        total_velocities[i] = norm_2(velocity) + std::sqrt(HeatCapacityRatio * (HeatCapacityRatio - 1) * SpecificHeatCV * temperature); 

        entropy_total_derivative.LoadNodalValues(NUMERICAL_ENTROPY, r_node, i, velocity, DeltaTime);
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

        const double max_local_residual = std::max(
            std::fabs(entropy_residual),
            std::fabs(specific_entropy_residual));

        const double density = inner_prod(rDensityTotalDerivative.Value, N);
        const double total_velocity = inner_prod(rTotalVelocities, N);

        max_residual =  std::max(max_residual, max_local_residual);
        max_density = std::max(max_density, density);
        max_total_velocity = std::max(max_total_velocity, total_velocity);
    }

    return {max_residual, max_density, max_total_velocity};

    KRATOS_CATCH("")
}


} // namespace Kratos
