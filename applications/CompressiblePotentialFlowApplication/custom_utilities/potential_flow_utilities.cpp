//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez

#include "custom_utilities/potential_flow_utilities.h"
#include "compressible_potential_flow_application_variables.h"
#include "fluid_dynamics_application_variables.h"
#include "includes/model_part.h"

namespace Kratos {
namespace PotentialFlowUtilities {
template <int Dim, int NumNodes>
array_1d<double, NumNodes> GetWakeDistances(const Element& rElement)
{
    const auto distances = rElement.GetValue(WAKE_ELEMENTAL_DISTANCES);
    KRATOS_ERROR_IF(distances.size() < NumNodes)
        << "Wake element with Id #" <<  rElement.Id() << " has no distances " << distances << std::endl;
    return rElement.GetValue(WAKE_ELEMENTAL_DISTANCES);
}

template <int Dim, int NumNodes>
void GetEquationIdVectorNormalElement(const Element& rElement, EquationIdVectorType& rElementalIdList)
{
    const auto r_geometry = rElement.GetGeometry();
    for (unsigned int i = 0; i < NumNodes; i++){
        rElementalIdList[i] = r_geometry[i].GetDof(VELOCITY_POTENTIAL).EquationId();
    }
}

template <int Dim, int NumNodes>
void GetEquationIdVectorKuttaElement(const Element& rElement, EquationIdVectorType& rElementalIdList)
{
    const auto& r_geometry = rElement.GetGeometry();
    // Kutta elements have only negative part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (!r_geometry[i].GetValue(TRAILING_EDGE))
            rElementalIdList[i] = r_geometry[i].GetDof(VELOCITY_POTENTIAL).EquationId();
        else
            rElementalIdList[i] = r_geometry[i].GetDof(AUXILIARY_VELOCITY_POTENTIAL).EquationId();
    }
}

template <int Dim, int NumNodes>
void GetDofListNormalElement(const Element& rElement, DofsVectorType& rElementalDofList)
{
    const auto r_geometry = rElement.GetGeometry();
    for (unsigned int i = 0; i < NumNodes; i++){
        rElementalDofList[i] = r_geometry[i].pGetDof(VELOCITY_POTENTIAL);
    }
}

template <int Dim, int NumNodes>
BoundedVector<double, NumNodes> GetPotentialOnNormalElement(const Element& rElement)
{
    const int kutta = rElement.GetValue(KUTTA);
    array_1d<double, NumNodes> potentials;

    const auto r_geometry = rElement.GetGeometry();

    if (kutta == 0) {
        for (unsigned int i = 0; i < NumNodes; i++) {
            potentials[i] = r_geometry[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
        }
    }
    else {
        for (unsigned int i = 0; i < NumNodes; i++) {
            if (!r_geometry[i].GetValue(TRAILING_EDGE)) {
                potentials[i] = r_geometry[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
            }
            else {
                potentials[i] = r_geometry[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL);
            }
        }
    }

    return potentials;
}

template <int Dim, int NumNodes>
BoundedVector<double, 2 * NumNodes> GetPotentialOnWakeElement(
    const Element& rElement, const array_1d<double, NumNodes>& rDistances)
{
    const auto upper_potentials =
        GetPotentialOnUpperWakeElement<Dim, NumNodes>(rElement, rDistances);

    const auto lower_potentials =
        GetPotentialOnLowerWakeElement<Dim, NumNodes>(rElement, rDistances);

    BoundedVector<double, 2 * NumNodes> split_element_values;
    for (unsigned int i = 0; i < NumNodes; i++) {
        split_element_values[i] = upper_potentials[i];
        split_element_values[NumNodes + i] = lower_potentials[i];
    }

    return split_element_values;
}

template <int Dim, int NumNodes>
BoundedVector<double, NumNodes> GetPotentialOnUpperWakeElement(
    const Element& rElement, const array_1d<double, NumNodes>& rDistances)
{
    array_1d<double, NumNodes> upper_potentials;
    const auto r_geometry = rElement.GetGeometry();

    for (unsigned int i = 0; i < NumNodes; i++){
        if (rDistances[i] > 0.0){
            upper_potentials[i] = r_geometry[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
        }
        else{
            upper_potentials[i] = r_geometry[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL);
        }
    }

    return upper_potentials;
}

template <int Dim, int NumNodes>
BoundedVector<double, NumNodes> GetPotentialOnLowerWakeElement(
    const Element& rElement, const array_1d<double, NumNodes>& rDistances)
{
    array_1d<double, NumNodes> lower_potentials;
    const auto r_geometry = rElement.GetGeometry();

    for (unsigned int i = 0; i < NumNodes; i++){
        if (rDistances[i] < 0.0){
            lower_potentials[i] = r_geometry[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
        }
        else{
            lower_potentials[i] = r_geometry[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL);
        }
    }

    return lower_potentials;
}

template <int Dim, int NumNodes>
array_1d<double, Dim> ComputeVelocity(const Element& rElement)
{
    const int wake = rElement.GetValue(WAKE);

    if (wake == 0)
        return ComputeVelocityNormalElement<Dim,NumNodes>(rElement);
    else
        return ComputeVelocityUpperWakeElement<Dim,NumNodes>(rElement);
}

template <int Dim, int NumNodes>
array_1d<double, Dim> ComputeVelocityNormalElement(const Element& rElement)
{
    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(rElement.GetGeometry(), data.DN_DX, data.N, data.vol);

    data.potentials = GetPotentialOnNormalElement<Dim,NumNodes>(rElement);

    return prod(trans(data.DN_DX), data.potentials);
}

template <int Dim, int NumNodes>
array_1d<double, Dim> ComputeVelocityUpperWakeElement(const Element& rElement)
{
    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(rElement.GetGeometry(), data.DN_DX, data.N, data.vol);

    const auto& r_distances = GetWakeDistances<Dim,NumNodes>(rElement);

    data.potentials = GetPotentialOnUpperWakeElement<Dim,NumNodes>(rElement, r_distances);

    return prod(trans(data.DN_DX), data.potentials);
}

template <int Dim, int NumNodes>
array_1d<double, Dim> ComputeVelocityLowerWakeElement(const Element& rElement)
{
    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(rElement.GetGeometry(), data.DN_DX, data.N, data.vol);

    const auto& r_distances = GetWakeDistances<Dim,NumNodes>(rElement);

    data.potentials = GetPotentialOnLowerWakeElement<Dim,NumNodes>(rElement, r_distances);

    return prod(trans(data.DN_DX), data.potentials);
}

template <int Dim, int NumNodes>
double ComputeIncompressiblePressureCoefficient(const Element& rElement, const ProcessInfo& rCurrentProcessInfo)
{
    const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    const double free_stream_velocity_norm = inner_prod(free_stream_velocity, free_stream_velocity);

    KRATOS_ERROR_IF(free_stream_velocity_norm < std::numeric_limits<double>::epsilon())
        << "Error on element -> " << rElement.Id() << "\n"
        << "free_stream_velocity_norm must be larger than zero." << std::endl;

    array_1d<double, Dim> v = ComputeVelocity<Dim,NumNodes>(rElement);

    double pressure_coefficient = (free_stream_velocity_norm - inner_prod(v, v)) /
               free_stream_velocity_norm; // 0.5*(norm_2(free_stream_velocity) - norm_2(v));
    return pressure_coefficient;
}

template <int Dim, int NumNodes>
double ComputePerturbationIncompressiblePressureCoefficient(const Element& rElement, const ProcessInfo& rCurrentProcessInfo)
{
    const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    const double free_stream_velocity_norm = inner_prod(free_stream_velocity, free_stream_velocity);

    KRATOS_ERROR_IF(free_stream_velocity_norm < std::numeric_limits<double>::epsilon())
        << "Error on element -> " << rElement.Id() << "\n"
        << "free_stream_velocity_norm must be larger than zero." << std::endl;

    array_1d<double, Dim> velocity = ComputeVelocity<Dim,NumNodes>(rElement);
    for (unsigned int i = 0; i < Dim; i++){
        velocity[i] += free_stream_velocity[i];
    }

    double pressure_coefficient = (free_stream_velocity_norm - inner_prod(velocity, velocity)) /
               free_stream_velocity_norm; // 0.5*(norm_2(free_stream_velocity) - norm_2(v));
    return pressure_coefficient;
}


template <int Dim, int NumNodes>
double ComputeCompressiblePressureCoefficient(const Element& rElement, const ProcessInfo& rCurrentProcessInfo)
{
    // Reading free stream conditions
    const array_1d<double, 3>& vinfinity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    const double M_inf = rCurrentProcessInfo[FREE_STREAM_MACH];
    const double heat_capacity_ratio = rCurrentProcessInfo[HEAT_CAPACITY_RATIO];

    // Computing local velocity
    array_1d<double, Dim> v = ComputeVelocity<Dim, NumNodes>(rElement);

    // Computing squares
    const double v_inf_2 = inner_prod(vinfinity, vinfinity);
    const double M_inf_2 = M_inf * M_inf;
    double v_2 = inner_prod(v, v);

    KRATOS_ERROR_IF(v_inf_2 < std::numeric_limits<double>::epsilon())
        << "Error on element -> " << rElement.Id() << "\n"
        << "v_inf_2 must be larger than zero." << std::endl;

    const double base = 1 + (heat_capacity_ratio - 1) * M_inf_2 * (1 - v_2 / v_inf_2) / 2;

    return 2 * (pow(base, heat_capacity_ratio / (heat_capacity_ratio - 1)) - 1) /
           (heat_capacity_ratio * M_inf_2);
}

template <int Dim, int NumNodes>
double ComputePerturbationCompressiblePressureCoefficient(const Element& rElement, const ProcessInfo& rCurrentProcessInfo)
{
    // Reading free stream conditions
    const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    const double M_inf = rCurrentProcessInfo[FREE_STREAM_MACH];
    const double heat_capacity_ratio = rCurrentProcessInfo[HEAT_CAPACITY_RATIO];

    // Computing local velocity
    array_1d<double, Dim> velocity = ComputeVelocity<Dim,NumNodes>(rElement);
    for (unsigned int i = 0; i < Dim; i++){
        velocity[i] += free_stream_velocity[i];
    }

    // Computing squares
    const double v_inf_2 = inner_prod(free_stream_velocity, free_stream_velocity);
    const double M_inf_2 = M_inf * M_inf;
    double v_2 = inner_prod(velocity, velocity);

    KRATOS_ERROR_IF(v_inf_2 < std::numeric_limits<double>::epsilon())
        << "Error on element -> " << rElement.Id() << "\n"
        << "v_inf_2 must be larger than zero." << std::endl;

    const double base = 1 + (heat_capacity_ratio - 1) * M_inf_2 * (1 - v_2 / v_inf_2) / 2;

    return 2 * (pow(base, heat_capacity_ratio / (heat_capacity_ratio - 1)) - 1) /
           (heat_capacity_ratio * M_inf_2);
}

template <int Dim, int NumNodes>
double ComputeLocalSpeedOfSound(const Element& rElement, const ProcessInfo& rCurrentProcessInfo)
{
    // Implemented according to Equation 8.7 of Drela, M. (2014) Flight Vehicle
    // Aerodynamics, The MIT Press, London
    // Reading free stream conditions
    const array_1d<double, 3>& v_inf = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    const double M_inf = rCurrentProcessInfo[FREE_STREAM_MACH];
    const double heat_capacity_ratio = rCurrentProcessInfo[HEAT_CAPACITY_RATIO];
    const double a_inf = rCurrentProcessInfo[SOUND_VELOCITY];

    // Computing local velocity
    array_1d<double, Dim> v = ComputeVelocity<Dim, NumNodes>(rElement);

    // Computing squares
    const double v_inf_2 = inner_prod(v_inf, v_inf);
    const double M_inf_2 = M_inf * M_inf;
    const double v_2 = inner_prod(v, v);

    KRATOS_ERROR_IF(v_inf_2 < std::numeric_limits<double>::epsilon())
        << "Error on element -> " << rElement.Id() << "\n"
        << "v_inf_2 must be larger than zero." << std::endl;

    return a_inf * sqrt(1 + (heat_capacity_ratio - 1) * M_inf_2 * (1 - v_2 / v_inf_2) / 2);
}

template <int Dim, int NumNodes>
double ComputePerturbationLocalSpeedOfSound(const Element& rElement, const ProcessInfo& rCurrentProcessInfo)
{
    // Implemented according to Equation 8.7 of Drela, M. (2014) Flight Vehicle
    // Aerodynamics, The MIT Press, London
    // Reading free stream conditions
    const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    const double M_inf = rCurrentProcessInfo[FREE_STREAM_MACH];
    const double heat_capacity_ratio = rCurrentProcessInfo[HEAT_CAPACITY_RATIO];
    const double a_inf = rCurrentProcessInfo[SOUND_VELOCITY];

    // Computing local velocity
    array_1d<double, Dim> velocity = ComputeVelocity<Dim,NumNodes>(rElement);
    for (unsigned int i = 0; i < Dim; i++){
        velocity[i] += free_stream_velocity[i];
    }

    // Computing squares
    const double v_inf_2 = inner_prod(free_stream_velocity, free_stream_velocity);
    const double M_inf_2 = M_inf * M_inf;
    const double v_2 = inner_prod(velocity, velocity);

    KRATOS_ERROR_IF(v_inf_2 < std::numeric_limits<double>::epsilon())
        << "Error on element -> " << rElement.Id() << "\n"
        << "v_inf_2 must be larger than zero." << std::endl;

    return a_inf * sqrt(1 + (heat_capacity_ratio - 1) * M_inf_2 * (1 - v_2 / v_inf_2) / 2);
}

template <int Dim, int NumNodes>
double ComputeLocalMachNumber(const Element& rElement, const ProcessInfo& rCurrentProcessInfo)
{
    // Implemented according to Equation 8.8 of Drela, M. (2014) Flight Vehicle
    // Aerodynamics, The MIT Press, London

    array_1d<double, Dim> velocity = ComputeVelocity<Dim, NumNodes>(rElement);
    const double velocity_module = sqrt(inner_prod(velocity, velocity));
    const double local_speed_of_sound = ComputeLocalSpeedOfSound<Dim, NumNodes>(rElement, rCurrentProcessInfo);

    return velocity_module / local_speed_of_sound;
}

template <int Dim, int NumNodes>
double ComputePerturbationLocalMachNumber(const Element& rElement, const ProcessInfo& rCurrentProcessInfo)
{
    // Implemented according to Equation 8.8 of Drela, M. (2014) Flight Vehicle
    // Aerodynamics, The MIT Press, London

    const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    array_1d<double, Dim> velocity = ComputeVelocity<Dim,NumNodes>(rElement);
    for (unsigned int i = 0; i < Dim; i++){
        velocity[i] += free_stream_velocity[i];
    }
    const double velocity_module = sqrt(inner_prod(velocity, velocity));
    const double local_speed_of_sound = ComputePerturbationLocalSpeedOfSound<Dim, NumNodes>(rElement, rCurrentProcessInfo);

    return velocity_module / local_speed_of_sound;
}

template <int Dim, int NumNodes>
double ComputePerturbationDensity(const Element& rElement, const ProcessInfo& rCurrentProcessInfo)
{
    // Reading free stream conditions
    const double rho_inf = rCurrentProcessInfo[FREE_STREAM_DENSITY];
    const double M_inf = rCurrentProcessInfo[FREE_STREAM_MACH];
    const double heat_capacity_ratio = rCurrentProcessInfo[HEAT_CAPACITY_RATIO];

    // Computing local mach number
    double local_mach_number = ComputePerturbationLocalMachNumber<Dim, NumNodes>(rElement, rCurrentProcessInfo);

    // Computing squares
    const double M_inf_2 = M_inf * M_inf;
    const double M_2 = local_mach_number * local_mach_number;

    // Computing density according to Equation 8.9 of Drela, M. (2014) Flight Vehicle
    // Aerodynamics, The MIT Press, London
    const double numerator = 1 + (heat_capacity_ratio - 1) * M_inf_2 / 2;
    const double denominator = 1 + (heat_capacity_ratio - 1) * M_2 / 2;
    const double base = numerator / denominator;

    if (base > 0.0)
    {
        return rho_inf * pow(base, 1 / (heat_capacity_ratio - 1));
    }
    else
    {
        KRATOS_WARNING("ComputeDensity") << "Using density correction" << std::endl;
        return rho_inf * 0.00001;
    }
}

template <int Dim, int NumNodes>
double ComputeDensityDerivative(const double& rDensity, const ProcessInfo& rCurrentProcessInfo)
{
    // Reading free stream conditions
    const double rho_inf = rCurrentProcessInfo[FREE_STREAM_DENSITY];
    const double heat_capacity_ratio = rCurrentProcessInfo[HEAT_CAPACITY_RATIO];
    const double a_inf = rCurrentProcessInfo[SOUND_VELOCITY];

    return -pow(rho_inf, heat_capacity_ratio - 1) *
           pow(rDensity, 2 - heat_capacity_ratio) / (2 * a_inf * a_inf);
}

template <int Dim, int NumNodes>
double ComputeUpwindDensity(const Element& rElement, const Element& rUpstreamElement, const ProcessInfo& rCurrentProcessInfo)
{
    const double upwind_factor = ComputeSwitchingOperator<Dim, NumNodes>(rElement, rUpstreamElement, rCurrentProcessInfo);
    const double density = ComputePerturbationDensity<Dim, NumNodes>(rElement, rCurrentProcessInfo);
    const double upstream_density = ComputePerturbationDensity<Dim, NumNodes>(rUpstreamElement, rCurrentProcessInfo);

    return density - upwind_factor * (density - upstream_density);
}

template <int Dim, int NumNodes>
double ComputeSwitchingOperator(const Element& rElement, const Element& rUpstreamElement, const ProcessInfo& rCurrentProcessInfo)
{
    const double upwind_factor = ComputeUpwindFactor<Dim, NumNodes>(rElement, rCurrentProcessInfo);
    const double upstream_upwind_factor = ComputeUpwindFactor<Dim, NumNodes>(rUpstreamElement, rCurrentProcessInfo);

    // Subsonic flow (local_mach_number < mach_number_limit)
    if(upwind_factor < 0.0){
        return 0.0;
    }
    // Supersonic flow and accelerating (local_mach_number > upstream_mach_number)
    else if( upwind_factor > upstream_upwind_factor){
        return upwind_factor;
    }
    // Supersonic flow and decelerating (local_mach_number < upstream_mach_number)
    else{
        return upstream_upwind_factor;
    }
}

template <int Dim, int NumNodes>
double ComputeUpwindFactor(const Element& rElement, const ProcessInfo& rCurrentProcessInfo)
{
    const double mach_number_limit = rCurrentProcessInfo[MACH_LIMIT];
    const double M_c_2 = mach_number_limit * mach_number_limit;

    const double local_mach_number = ComputePerturbationLocalMachNumber<Dim, NumNodes>(rElement, rCurrentProcessInfo);
    const double M_2 = local_mach_number * local_mach_number;

    return 1 - M_c_2 / M_2;
}

template <int Dim, int NumNodes>
double ComputeDerivativeUpwindFactorWRTMachNumberSquared(const Element& rElement, const ProcessInfo& rCurrentProcessInfo)
{
    const double mach_number_limit = rCurrentProcessInfo[MACH_LIMIT];
    const double M_c_2 = mach_number_limit * mach_number_limit;

    const double local_mach_number = ComputePerturbationLocalMachNumber<Dim, NumNodes>(rElement, rCurrentProcessInfo);
    const double M_4 = local_mach_number * local_mach_number * local_mach_number * local_mach_number;

    return M_c_2 / M_4;
}

template <int Dim, int NumNodes>
double ComputeDerivativeMachNumberSquaredWRTVelocitySquared(const Element& rElement, const ProcessInfo& rCurrentProcessInfo)
{
    const double free_stream_mach_number = rCurrentProcessInfo[FREE_STREAM_MACH];
    const double hcr = rCurrentProcessInfo[HEAT_CAPACITY_RATIO];
    const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];

    array_1d<double, Dim> velocity = ComputeVelocity<Dim,NumNodes>(rElement);
    for (unsigned int i = 0; i < Dim; i++){
        velocity[i] += free_stream_velocity[i];
    }
    // Computing squares
    const double v_2 = inner_prod(velocity, velocity);
    const double v_inf_2 = inner_prod(free_stream_velocity, free_stream_velocity);
    const double M_inf_2 = free_stream_mach_number * free_stream_mach_number;

    const double numerator = 1 + (hcr - 1) / 2.0 * M_inf_2;
    const double denominator = pow(1 + (hcr - 1) / 2.0 * M_inf_2 * (1 - v_2 / v_inf_2), 2);

    return M_inf_2 * numerator / (v_inf_2 * denominator);
}

template <int Dim, int NumNodes>
BoundedVector<double, NumNodes> ComputeDrhoDphiSupersonicAccelerating(const Element& rElement, const Element& rUpstreamElement, const ProcessInfo& rCurrentProcessInfo)
{
    const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(rElement.GetGeometry(), data.DN_DX, data.N, data.vol);
    array_1d<double, Dim> velocity = ComputeVelocity<Dim,NumNodes>(rElement);
    for (unsigned int i = 0; i < Dim; i++){
        velocity[i] += free_stream_velocity[i];
    }
    const BoundedVector<double, NumNodes> DNV = prod(data.DN_DX, velocity);

    const double upwind_factor = ComputeUpwindFactor<Dim, NumNodes>(rElement, rCurrentProcessInfo);
    const double density = ComputePerturbationDensity<Dim, NumNodes>(rElement, rCurrentProcessInfo);
    const double upstream_density = ComputePerturbationDensity<Dim, NumNodes>(rUpstreamElement, rCurrentProcessInfo);
    const double Drho_Dv2 = ComputeDensityDerivative<Dim, NumNodes>(density, rCurrentProcessInfo);
    const double Dmu_DM2 = ComputeDerivativeUpwindFactorWRTMachNumberSquared<Dim, NumNodes>(rElement, rCurrentProcessInfo);
    const double DM2_Dv2 = ComputeDerivativeMachNumberSquaredWRTVelocitySquared<Dim, NumNodes>(rElement, rCurrentProcessInfo);

    const double factor = 2 * (Drho_Dv2 * (1 - upwind_factor) - Dmu_DM2 * DM2_Dv2 * (density - upstream_density));
    const BoundedVector<double, NumNodes> Drho_DPhi = factor * DNV;

    return Drho_DPhi;
}

template <int Dim, int NumNodes>
BoundedVector<double, NumNodes> ComputeDrhoDphiUpSupersonicAccelerating(const Element& rElement, const Element& rUpstreamElement, const ProcessInfo& rCurrentProcessInfo)
{
    const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(rUpstreamElement.GetGeometry(), data.DN_DX, data.N, data.vol);
    array_1d<double, Dim> velocity = ComputeVelocity<Dim,NumNodes>(rUpstreamElement);
    for (unsigned int i = 0; i < Dim; i++){
        velocity[i] += free_stream_velocity[i];
    }
    const BoundedVector<double, NumNodes> DNV = prod(data.DN_DX, velocity);

    const double upwind_factor = ComputeUpwindFactor<Dim, NumNodes>(rElement, rCurrentProcessInfo);
    const double density_up = ComputePerturbationDensity<Dim, NumNodes>(rUpstreamElement, rCurrentProcessInfo);
    const double Drho_Dv2_up = ComputeDensityDerivative<Dim, NumNodes>(density_up, rCurrentProcessInfo);

    const BoundedVector<double, NumNodes> Drho_DPhiUp = 2 * Drho_Dv2_up * upwind_factor * DNV;

    return Drho_DPhiUp;
}

template <int Dim, int NumNodes>
bool CheckIfElementIsCutByDistance(const BoundedVector<double, NumNodes>& rNodalDistances)
{
    // Initialize counters
    unsigned int number_of_nodes_with_positive_distance = 0;
    unsigned int number_of_nodes_with_negative_distance = 0;

    // Count how many element nodes are above and below the wake
    for (unsigned int i = 0; i < rNodalDistances.size(); i++) {
        if (rNodalDistances(i) < 0.0) {
            number_of_nodes_with_negative_distance += 1;
        }
        else {
            number_of_nodes_with_positive_distance += 1;
        }
    }

    // Elements with nodes above and below the wake are wake elements
    return number_of_nodes_with_negative_distance > 0 &&
           number_of_nodes_with_positive_distance > 0;
}

bool CheckIfElementIsTrailingEdge(const Element& rElement)
{
    const auto& r_geometry = rElement.GetGeometry();
    bool is_trailing_edge = false;

    for(unsigned int i_node = 0; i_node<r_geometry.size(); i_node++){

        if (r_geometry[i_node].GetValue(TRAILING_EDGE)) {
            is_trailing_edge = true;
        }
    }

    return is_trailing_edge;
}

template <int Dim>
void CheckIfWakeConditionsAreFulfilled(const ModelPart& rWakeModelPart, const double& rTolerance, const int& rEchoLevel)
{
    unsigned int number_of_unfulfilled_wake_conditions = 0;
    for (auto& r_element : rWakeModelPart.Elements()){
        const bool wake_condition_is_fulfilled =
            CheckWakeCondition<Dim, Dim + 1>(r_element, rTolerance, rEchoLevel);
        if (!wake_condition_is_fulfilled)
        {
            number_of_unfulfilled_wake_conditions += 1;
        }
    }
    KRATOS_WARNING_IF("CheckIfWakeConditionsAreFulfilled", number_of_unfulfilled_wake_conditions > 0)
        << "THE WAKE CONDITION IS NOT FULFILLED IN " << number_of_unfulfilled_wake_conditions
        << " ELEMENTS WITH AN ABSOLUTE TOLERANCE OF " << rTolerance << std::endl;
}

template <int Dim, int NumNodes>
bool CheckWakeCondition(const Element& rElement, const double& rTolerance, const int& rEchoLevel)
{
    const auto upper_velocity = ComputeVelocityUpperWakeElement<Dim,NumNodes>(rElement);
    const auto lower_velocity = ComputeVelocityLowerWakeElement<Dim,NumNodes>(rElement);

    bool wake_condition_is_fulfilled = true;
    for (unsigned int i = 0; i < upper_velocity.size(); i++){
        if(std::abs(upper_velocity[i] - lower_velocity[i]) > rTolerance){
            wake_condition_is_fulfilled = false;
            break;
        }
    }

    KRATOS_WARNING_IF("CheckWakeCondition", !wake_condition_is_fulfilled && rEchoLevel > 0)
        << "WAKE CONDITION NOT FULFILLED IN ELEMENT # " << rElement.Id() << std::endl;
    KRATOS_WARNING_IF("CheckWakeCondition", !wake_condition_is_fulfilled && rEchoLevel > 1)
        << "WAKE CONDITION NOT FULFILLED IN ELEMENT # " << rElement.Id()
        << " upper_velocity  = " << upper_velocity
        << " lower_velocity  = " << lower_velocity << std::endl;

    return wake_condition_is_fulfilled;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template instantiation

// 2D
template array_1d<double, 3> GetWakeDistances<2, 3>(const Element& rElement);
template void GetEquationIdVectorNormalElement<2, 3>(const Element& rElement, EquationIdVectorType& rElementalDofList);
template void GetEquationIdVectorKuttaElement<2, 3>(const Element& rElement, EquationIdVectorType& rElementalDofList);
template void GetDofListNormalElement<2, 3>(const Element& rElement, DofsVectorType& rElementalDofList);
template BoundedVector<double, 3> GetPotentialOnNormalElement<2, 3>(const Element& rElement);
template BoundedVector<double, 2 * 3> GetPotentialOnWakeElement<2, 3>(
    const Element& rElement, const array_1d<double, 3>& rDistances);
template BoundedVector<double, 3> GetPotentialOnUpperWakeElement<2, 3>(
    const Element& rElement, const array_1d<double, 3>& rDistances);
template BoundedVector<double, 3> GetPotentialOnLowerWakeElement<2, 3>(
    const Element& rElement, const array_1d<double, 3>& rDistances);
template array_1d<double, 2> ComputeVelocityNormalElement<2, 3>(const Element& rElement);
template array_1d<double, 2> ComputeVelocityUpperWakeElement<2, 3>(const Element& rElement);
template array_1d<double, 2> ComputeVelocityLowerWakeElement<2, 3>(const Element& rElement);
template array_1d<double, 2> ComputeVelocity<2, 3>(const Element& rElement);
template double ComputeIncompressiblePressureCoefficient<2, 3>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputePerturbationIncompressiblePressureCoefficient<2, 3>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeCompressiblePressureCoefficient<2, 3>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputePerturbationCompressiblePressureCoefficient<2, 3>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeLocalSpeedOfSound<2, 3>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputePerturbationLocalSpeedOfSound<2, 3>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeLocalMachNumber<2, 3>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputePerturbationLocalMachNumber<2, 3>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputePerturbationDensity<2, 3>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeDensityDerivative<2, 3>(const double& rDensity, const ProcessInfo& rCurrentProcessInfo);
template double ComputeUpwindDensity<2, 3>(const Element& rElement, const Element& rUpstreamElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeSwitchingOperator<2, 3>(const Element& rElement, const Element& rUpstreamElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeUpwindFactor<2, 3>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeDerivativeUpwindFactorWRTMachNumberSquared<2, 3>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeDerivativeMachNumberSquaredWRTVelocitySquared<2, 3>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template BoundedVector<double, 3> ComputeDrhoDphiSupersonicAccelerating<2, 3>(const Element& rElement, const Element& rUpstreamElement, const ProcessInfo& rCurrentProcessInfo);
template BoundedVector<double, 3> ComputeDrhoDphiUpSupersonicAccelerating<2, 3>(const Element& rElement, const Element& rUpstreamElement, const ProcessInfo& rCurrentProcessInfo);
template bool CheckIfElementIsCutByDistance<2, 3>(const BoundedVector<double, 3>& rNodalDistances);
template void KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) CheckIfWakeConditionsAreFulfilled<2>(const ModelPart&, const double& rTolerance, const int& rEchoLevel);
template bool CheckWakeCondition<2, 3>(const Element& rElement, const double& rTolerance, const int& rEchoLevel);

// 3D
template array_1d<double, 4> GetWakeDistances<3, 4>(const Element& rElement);
template void GetEquationIdVectorNormalElement<3, 4>(const Element& rElement, EquationIdVectorType& rElementalDofList);
template void GetEquationIdVectorKuttaElement<3, 4>(const Element& rElement, EquationIdVectorType& rElementalDofList);
template void GetDofListNormalElement<3, 4>(const Element& rElement, DofsVectorType& rElementalDofList);
template BoundedVector<double, 4> GetPotentialOnNormalElement<3, 4>(const Element& rElement);
template BoundedVector<double, 2 * 4> GetPotentialOnWakeElement<3, 4>(
    const Element& rElement, const array_1d<double, 4>& rDistances);
template BoundedVector<double, 4> GetPotentialOnUpperWakeElement<3, 4>(
    const Element& rElement, const array_1d<double, 4>& rDistances);
template BoundedVector<double, 4> GetPotentialOnLowerWakeElement<3, 4>(
    const Element& rElement, const array_1d<double, 4>& rDistances);
template array_1d<double, 3> ComputeVelocityNormalElement<3, 4>(const Element& rElement);
template array_1d<double, 3> ComputeVelocityUpperWakeElement<3, 4>(const Element& rElement);
template array_1d<double, 3> ComputeVelocityLowerWakeElement<3, 4>(const Element& rElement);
template array_1d<double, 3> ComputeVelocity<3, 4>(const Element& rElement);
template double ComputeIncompressiblePressureCoefficient<3, 4>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputePerturbationIncompressiblePressureCoefficient<3, 4>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeCompressiblePressureCoefficient<3, 4>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputePerturbationCompressiblePressureCoefficient<3, 4>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeLocalSpeedOfSound<3, 4>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputePerturbationLocalSpeedOfSound<3, 4>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeLocalMachNumber<3, 4>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputePerturbationLocalMachNumber<3, 4>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputePerturbationDensity<3, 4>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeDensityDerivative<3, 4>(const double& rDensity, const ProcessInfo& rCurrentProcessInfo);
template double ComputeUpwindDensity<3, 4>(const Element& rElement, const Element& rUpstreamElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeSwitchingOperator<3, 4>(const Element& rElement, const Element& rUpstreamElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeUpwindFactor<3, 4>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeDerivativeUpwindFactorWRTMachNumberSquared<3, 4>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeDerivativeMachNumberSquaredWRTVelocitySquared<3, 4>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template BoundedVector<double, 4> ComputeDrhoDphiSupersonicAccelerating<3, 4>(const Element& rElement, const Element& rUpstreamElement, const ProcessInfo& rCurrentProcessInfo);
template BoundedVector<double, 4> ComputeDrhoDphiUpSupersonicAccelerating<3, 4>(const Element& rElement, const Element& rUpstreamElement, const ProcessInfo& rCurrentProcessInfo);
template bool CheckIfElementIsCutByDistance<3, 4>(const BoundedVector<double, 4>& rNodalDistances);
template void  KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) CheckIfWakeConditionsAreFulfilled<3>(const ModelPart&, const double& rTolerance, const int& rEchoLevel);
template bool CheckWakeCondition<3, 4>(const Element& rElement, const double& rTolerance, const int& rEchoLevel);
} // namespace PotentialFlow
} // namespace Kratos
