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
    return rElement.GetValue(WAKE_ELEMENTAL_DISTANCES);
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
array_1d<double, Dim> ComputePerturbedVelocity(
    const Element& rElement, 
    const ProcessInfo& rCurrentProcessInfo)
{
    const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    array_1d<double, Dim> velocity = ComputeVelocity<Dim,NumNodes>(rElement);
    for (unsigned int i = 0; i < Dim; i++)
    {
        velocity[i] += free_stream_velocity[i];
    }

    return velocity;
}

template <int Dim, int NumNodes>
double ComputeMaximumVelocitySquared(const ProcessInfo& rCurrentProcessInfo)
{
    // Following Fully Simulataneous Coupling of the Full Potential Equation
    //           and the Integral Boundary Layer Equations in Three Dimensions
    //           by Brian Nishida (1996), Section A.2 and Section 2.5

    // maximum local squared mach number (user defined, 3.0 used as default)
    const double max_local_mach_squared = rCurrentProcessInfo[MACH_SQUARED_LIMIT];

    // read free stream values
    const double heat_capacity_ratio = rCurrentProcessInfo[HEAT_CAPACITY_RATIO];
    const double free_stream_mach = rCurrentProcessInfo[FREE_STREAM_MACH];
    const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];

    KRATOS_ERROR_IF(free_stream_mach < std::numeric_limits<double>::epsilon())
        << "ComputeMaximumVelocitySquared: free_stream_mach must be larger than zero." << std::endl;

    // make squares of values
    const double free_stream_mach_squared = std::pow(free_stream_mach, 2.0);
    const double free_stream_velocity_squared = inner_prod(free_stream_velocity, free_stream_velocity);

    // calculate velocity
    const double numerator = (2.0 + (heat_capacity_ratio - 1.0) * free_stream_mach_squared );
    const double denominator = (2.0 + (heat_capacity_ratio - 1.0) * max_local_mach_squared );
    const double factor = free_stream_velocity_squared * max_local_mach_squared / free_stream_mach_squared;

    KRATOS_ERROR_IF(denominator < std::numeric_limits<double>::epsilon())
        << "ComputeMaximumVelocitySquared: denominatior must be larger than zero." << std::endl;

    return factor * numerator / denominator;
}

// This function returns the square of the magnitude of the velocity,
// clamping it if it is over the maximum allowed
template <int Dim, int NumNodes>
double ComputeClampedVelocitySquared(
    const array_1d<double, Dim>& rVelocity,
    const ProcessInfo& rCurrentProcessInfo)
{
    // compute max velocity allowed by limit Mach number
    const double max_velocity_squared = ComputeMaximumVelocitySquared<Dim, NumNodes>(rCurrentProcessInfo);
    double local_velocity_squared = inner_prod(rVelocity, rVelocity);

    // check if local velocity should be changed
    if (local_velocity_squared > max_velocity_squared)
    {
        KRATOS_WARNING("Clamped local velocity") <<
        "SQUARE OF LOCAL VELOCITY ABOVE ALLOWED SQUARE OF VELOCITY"
        << " local_velocity_squared  = " << local_velocity_squared
        << " max_velocity_squared  = " << max_velocity_squared << std::endl;

        local_velocity_squared = max_velocity_squared;
    }

    return local_velocity_squared;
}

template <int Dim, int NumNodes>
double ComputeVelocityMagnitude(
    const double localMachNumberSquared,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Following Fully Simulataneous Coupling of the Full Potential Equation
    //           and the Integral Boundary Layer Equations in Three Dimensions
    //           by Brian Nishida (1996), Section A.2 and Section 2.5

    // read free stream values
    const double heat_capacity_ratio = rCurrentProcessInfo[HEAT_CAPACITY_RATIO];
    const double free_stream_mach = rCurrentProcessInfo[FREE_STREAM_MACH];
    const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];

    KRATOS_ERROR_IF(free_stream_mach < std::numeric_limits<double>::epsilon())
        << "ComputeVelocityMagnitude: free_stream_mach must be larger than zero." << std::endl;

    // make squares of values
    const double free_stream_mach_squared = std::pow(free_stream_mach, 2);
    const double free_stream_velocity_squared = inner_prod(free_stream_velocity, free_stream_velocity);

    // calculate velocity
    const double numerator = (2.0 + (heat_capacity_ratio - 1.0) * free_stream_mach_squared );
    const double denominator = (2.0 + (heat_capacity_ratio - 1.0) * localMachNumberSquared );
    const double factor = free_stream_velocity_squared * localMachNumberSquared / free_stream_mach_squared;

    KRATOS_ERROR_IF(denominator < std::numeric_limits<double>::epsilon())
        << "ComputeVelocityMagnitude: denominator must be larger than zero." << std::endl;

    return factor * numerator / denominator;
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

    const double base = 1.0 + (heat_capacity_ratio - 1.0) * M_inf_2 * (1.0 - v_2 / v_inf_2) / 2.0;

    return 2.0 * (std::pow(base, heat_capacity_ratio / (heat_capacity_ratio - 1.0)) - 1.0) /
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

    const double base = 1.0 + (heat_capacity_ratio - 1.0) * M_inf_2 * (1.0 - v_2 / v_inf_2) / 2.0;

    return 2.0 * (std::pow(base, heat_capacity_ratio / (heat_capacity_ratio - 1.0)) - 1.0) /
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

    return a_inf * std::sqrt(1.0 + (heat_capacity_ratio - 1.0) * M_inf_2 * (1.0 - v_2 / v_inf_2) / 2.0);
}

template <int Dim, int NumNodes>
double ComputeLocalSpeedofSoundSquared(
    const array_1d<double, Dim>& rVelocity,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Implemented according to Equation 8.7 of Drela, M. (2014) Flight Vehicle
    // Aerodynamics, The MIT Press, London

    // read free stream values
    const double free_stream_speed_sound = rCurrentProcessInfo[SOUND_VELOCITY];

    // make squares of value
    const double free_stream_speed_sound_squared = std::pow(free_stream_speed_sound,2.0);

    // computes square of velocity including clamping according to MACH_SQUARED_LIMIT
    const double local_velocity_squared = ComputeClampedVelocitySquared<Dim, NumNodes>(rVelocity, rCurrentProcessInfo);

    // computes square bracket term with clamped velocity squared
    const double speed_of_sound_factor = ComputeSquaredSpeedofSoundFactor<Dim, NumNodes>(local_velocity_squared, rCurrentProcessInfo);

    return free_stream_speed_sound_squared * speed_of_sound_factor;
}

template <int Dim, int NumNodes>
double ComputeSquaredSpeedofSoundFactor(
    const double localVelocitySquared,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Implemented according to Equation 8.7 of Drela, M. (2014) Flight Vehicle
    // Aerodynamics, The MIT Press, London

    // read free stream values
    const double heat_capacity_ratio = rCurrentProcessInfo[HEAT_CAPACITY_RATIO];
    const double free_stream_mach = rCurrentProcessInfo[FREE_STREAM_MACH];
    const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];

    // make squares of values
    const double free_stream_mach_squared = std::pow(free_stream_mach, 2.0);
    const double free_stream_velocity_squared = inner_prod(free_stream_velocity, free_stream_velocity);

    return 1.0 + 0.5*(heat_capacity_ratio - 1.0)*
            free_stream_mach_squared*(1.0 - (localVelocitySquared/free_stream_velocity_squared));
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

    return a_inf * std::sqrt(1 + (heat_capacity_ratio - 1) * M_inf_2 * (1 - v_2 / v_inf_2) / 2);
}

template <int Dim, int NumNodes>
double ComputeLocalMachNumber(const Element& rElement, const ProcessInfo& rCurrentProcessInfo)
{
    // Implemented according to Equation 8.8 of Drela, M. (2014) Flight Vehicle
    // Aerodynamics, The MIT Press, London

    array_1d<double, Dim> velocity = ComputeVelocity<Dim, NumNodes>(rElement);
    const double velocity_module = std::sqrt(inner_prod(velocity, velocity));
    const double local_speed_of_sound = ComputeLocalSpeedOfSound<Dim, NumNodes>(rElement, rCurrentProcessInfo);

    return velocity_module / local_speed_of_sound;
}

template <int Dim, int NumNodes>
double ComputeLocalMachNumberSquared(
    const array_1d<double, Dim>& rVelocity,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Implemented according to Equation 8.8 of
    // Drela, M. (2014) Flight VehicleAerodynamics, The MIT Press, London

    const double local_speed_of_sound_squared = ComputeLocalSpeedofSoundSquared<Dim, NumNodes>(rVelocity, rCurrentProcessInfo);

    KRATOS_ERROR_IF(local_speed_of_sound_squared < std::numeric_limits<double>::epsilon())
        << "ComputeLocalMachNumberSquared: local speed of sound squared squared is less than zero." << std::endl;

    // computes square of velocity including clamping according to MACH_SQUARED_LIMIT
    const double local_velocity_squared = ComputeClampedVelocitySquared<Dim, NumNodes>(rVelocity, rCurrentProcessInfo);

    return local_velocity_squared / local_speed_of_sound_squared;
}

template <int Dim, int NumNodes>
double ComputeDerivativeLocalMachSquaredWRTVelocitySquared(
    const array_1d<double, Dim>& rVelocity,
    const double localMachNumberSquared,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Following Fully Simulataneous Coupling of the Full Potential Equation
    //           and the Integral Boundary Layer Equations in Three Dimensions
    //           by Brian Nishida (1996), Section A.2.3

    // read free stream values
    const double heat_capacity_ratio = rCurrentProcessInfo[HEAT_CAPACITY_RATIO];
    const double free_stream_mach = rCurrentProcessInfo[FREE_STREAM_MACH];
    const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];

    // make squares of values
    const double free_stream_mach_squared = std::pow(free_stream_mach, 2.0);
    const double free_stream_velocity_squared = inner_prod(free_stream_velocity, free_stream_velocity);

    KRATOS_ERROR_IF(free_stream_velocity_squared < std::numeric_limits<double>::epsilon())
        << "ComputeDerivativeLocalMachSquaredWRTVelocitySquared: free stream velocity squared squared is less than zero." << std::endl;

    // computes square of velocity including clamping according to MACH_SQUARED_LIMIT
    const double local_velocity_squared = ComputeClampedVelocitySquared<Dim, NumNodes>(rVelocity, rCurrentProcessInfo);

    KRATOS_ERROR_IF(local_velocity_squared < std::numeric_limits<double>::epsilon())
        << "ComputeDerivativeLocalMachSquaredWRTVelocitySquared: local velocity squared squared is less than zero." << std::endl;

    // square bracket term
    const double speed_of_sound_factor = ComputeSquaredSpeedofSoundFactor<Dim, NumNodes>(local_velocity_squared, rCurrentProcessInfo);

    KRATOS_ERROR_IF(speed_of_sound_factor < std::numeric_limits<double>::epsilon())
        << "ComputeDerivativeLocalMachSquaredWRTVelocitySquared: speed of sound factor must be larger than zero." << std::endl;

    const double second_term_factor = 0.5 * (heat_capacity_ratio - 1.0) / free_stream_velocity_squared * free_stream_mach_squared;

    return localMachNumberSquared * ((1.0 / local_velocity_squared) + second_term_factor / speed_of_sound_factor);
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
    const double velocity_module = std::sqrt(inner_prod(velocity, velocity));
    const double local_speed_of_sound = ComputePerturbationLocalSpeedOfSound<Dim, NumNodes>(rElement, rCurrentProcessInfo);

    return velocity_module / local_speed_of_sound;
}

template <int Dim, int NumNodes>
double ComputeUpwindFactor(
        double localMachNumberSquared, 
        const ProcessInfo& rCurrentProcessInfo)
{
    // Following Fully Simulataneous Coupling of the Full Potential Equation
    //           and the Integral Boundary Layer Equations in Three Dimensions
    //           by Brian Nishida (1996), Equation 2.13

    // read free stream values
    // default CRITICAL_MACH - 0.99 and default UPWIND_FACTOR_CONSTANT - 1.0
    const double critical_mach = rCurrentProcessInfo[CRITICAL_MACH];
    const double upwind_factor_constant = rCurrentProcessInfo[UPWIND_FACTOR_CONSTANT];

    if(localMachNumberSquared < 1e-3){
        localMachNumberSquared = 1e-3;
        KRATOS_WARNING("ComputeUpwindFactor") << "localMachNumberSquared is smaller than 1-3 and is being clamped to 1e-3"  <<  std::endl;
    }

    return upwind_factor_constant * (1.0 - std::pow(critical_mach, 2.0) / localMachNumberSquared);
}

template <int Dim, int NumNodes>
double SelectMaxUpwindFactor(
        const array_1d<double, Dim>& rCurrentVelocity, 
        const array_1d<double, Dim>& rUpwindVelocity, 
        const ProcessInfo& rCurrentProcessInfo)
{
    // Following Fully Simulataneous Coupling of the Full Potential Equation
    //           and the Integral Boundary Layer Equations in Three Dimensions
    //           by Brian Nishida (1996), Equation 2.13
    const double current_element_mach_squared = ComputeLocalMachNumberSquared<Dim,NumNodes>(rCurrentVelocity, rCurrentProcessInfo);
    const double upwind_element_mach_squared = ComputeLocalMachNumberSquared<Dim,NumNodes>(rUpwindVelocity, rCurrentProcessInfo);

    array_1d<double, 3> upwind_factor_options(3, 0.0);

    upwind_factor_options[1] = ComputeUpwindFactor<Dim, NumNodes>(current_element_mach_squared, rCurrentProcessInfo);
    upwind_factor_options[2] = ComputeUpwindFactor<Dim, NumNodes>(upwind_element_mach_squared, rCurrentProcessInfo);

    const auto case_option = ComputeUpwindFactorCase<Dim, NumNodes>(upwind_factor_options);
    return upwind_factor_options[case_option];
}

template <int Dim, int NumNodes>
size_t ComputeUpwindFactorCase(array_1d<double, 3>& rUpwindFactorOptions)
{
    // Following Fully Simulataneous Coupling of the Full Potential Equation
    //           and the Integral Boundary Layer Equations in Three Dimensions
    //           by Brian Nishida (1996), Equation 2.13

    // a subsonic current element should always return case 0
    if (rUpwindFactorOptions[1] < 0.0)
    {
        rUpwindFactorOptions[1] = 0.0;
        rUpwindFactorOptions[2] = 0.0;
    }
    const auto max_upwind_factor_opt = std::max_element(rUpwindFactorOptions.begin(), rUpwindFactorOptions.end());
    
    // Case 0: Subsonic flow
    // Case 1: Supersonic and accelerating flow (M^2 > M^2_up)
    // Case 2: Supersonic and decelerating flow (M^2 < M^2_up)
    // If Case 1 = Case 2, returns Case 1
    return std::distance(rUpwindFactorOptions.begin(), max_upwind_factor_opt);
}

template <int Dim, int NumNodes>
double ComputeUpwindFactorDerivativeWRTMachSquared(
    const double localMachNumberSquared,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Following Fully Simulataneous Coupling of the Full Potential Equation
    //           and the Integral Boundary Layer Equations in Three Dimensions
    //           by Brian Nishida (1996), section A.2

    // read free stream values
    // default CRITICAL_MACH - 0.99 and default UPWIND_FACTOR_CONSTANT - 1.0
    const double critical_mach = rCurrentProcessInfo[CRITICAL_MACH];
    const double upwind_factor_constant = rCurrentProcessInfo[UPWIND_FACTOR_CONSTANT];

    KRATOS_ERROR_IF(std::pow(localMachNumberSquared, 2.0) < std::numeric_limits<double>::epsilon())
        << "ComputeUpwindFactorDerivativeWRTMachSquared: local mach number squared is less than zero." << std::endl;

    return upwind_factor_constant * std::pow(critical_mach, 2.0) / std::pow(localMachNumberSquared, 2.0);
}

template <int Dim, int NumNodes>
double ComputeUpwindFactorDerivativeWRTVelocitySquared(
    const array_1d<double, Dim>& rVelocity,
    const ProcessInfo& rCurrentProcessInfo)
{

    // Following Fully Simulataneous Coupling of the Full Potential Equation
    //           and the Integral Boundary Layer Equations in Three Dimensions
    //           by Brian Nishida (1996), section A.2
    const double mach_squared = ComputeLocalMachNumberSquared<Dim,NumNodes>(rVelocity, rCurrentProcessInfo);
    const double upwind_factor_derivative = ComputeUpwindFactorDerivativeWRTMachSquared<Dim, NumNodes>(mach_squared,rCurrentProcessInfo);
    const double mach_number_derivative = ComputeDerivativeLocalMachSquaredWRTVelocitySquared<Dim, NumNodes>(rVelocity, mach_squared, rCurrentProcessInfo);
    return upwind_factor_derivative * mach_number_derivative;
}

template <int Dim, int NumNodes>
double ComputeDensity(
    const double localMachNumberSquared, 
    const ProcessInfo& rCurrentProcessInfo)
{
    // Implemented according to Equation 8.9 of Drela, M. (2014) Flight Vehicle
    // Aerodynamics, The MIT Press, London

    // reading free stream values
    const double free_stream_density = rCurrentProcessInfo[FREE_STREAM_DENSITY];
    const double free_stream_mach = rCurrentProcessInfo[FREE_STREAM_MACH];
    const double heat_capacity_ratio = rCurrentProcessInfo[HEAT_CAPACITY_RATIO];

    // density calculation
    const double numerator = 1.0 + (0.5 * (heat_capacity_ratio - 1.0)) * std::pow(free_stream_mach, 2.0);
    const double denominator = 1.0 + (0.5 * (heat_capacity_ratio - 1.0)) * localMachNumberSquared;

    KRATOS_ERROR_IF(denominator < std::numeric_limits<double>::epsilon())
        << "ComputeDensity: denominatior must be larger than zero." << std::endl;

    KRATOS_ERROR_IF((heat_capacity_ratio - 1.0) < std::numeric_limits<double>::epsilon())
        << "ComputeDensity: heat capacity ratio is smaller than 1." << std::endl;

    return free_stream_density * std::pow((numerator/denominator), 1.0/(heat_capacity_ratio - 1.0));
}

template <int Dim, int NumNodes>
double ComputeUpwindedDensity(
    const array_1d<double, Dim>& rCurrentVelocity, 
    const array_1d<double, Dim>& rUpwindVelocity, 
    const ProcessInfo& rCurrentProcessInfo)
{
    // Following Fully Simulataneous Coupling of the Full Potential Equation
    //           and the Integral Boundary Layer Equations in Three Dimensions
    //           by Brian Nishida (1996), Equation 2.12

    const double upwind_factor = SelectMaxUpwindFactor<Dim,NumNodes>(rCurrentVelocity, rUpwindVelocity, rCurrentProcessInfo);

    const double current_element_mach_squared = ComputeLocalMachNumberSquared<Dim,NumNodes>(rCurrentVelocity, rCurrentProcessInfo);
    const double upwind_element_mach_squared = ComputeLocalMachNumberSquared<Dim,NumNodes>(rUpwindVelocity, rCurrentProcessInfo);

    const double current_element_density = ComputeDensity<Dim,NumNodes>(current_element_mach_squared, rCurrentProcessInfo);
    const double upwind_element_density = ComputeDensity<Dim,NumNodes>(upwind_element_mach_squared, rCurrentProcessInfo);

    return current_element_density - upwind_factor * (current_element_density - upwind_element_density);
}

template <int Dim, int NumNodes>
double ComputeDensityDerivativeWRTVelocitySquared(
    const double localMachNumberSquared, 
    const ProcessInfo& rCurrentProcessInfo)
{
    // Following Fully Simulataneous Coupling of the Full Potential Equation
    //           and the Integral Boundary Layer Equations in Three Dimensions
    //           by Brian Nishida (1996), Section A.2.5

    // read free stream values
    const double free_stream_density = rCurrentProcessInfo[FREE_STREAM_DENSITY];
    const double free_stream_mach = rCurrentProcessInfo[FREE_STREAM_MACH];
    const double heat_capacity_ratio = rCurrentProcessInfo[HEAT_CAPACITY_RATIO];
    const BoundedVector<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];

    KRATOS_ERROR_IF(free_stream_mach < std::numeric_limits<double>::epsilon())
        << "ComputeDensityDerivativeWRTVelocitySquared: free stream mach number must be larger than zero." << std::endl;

    const double local_velocity_squared = ComputeVelocityMagnitude<Dim, NumNodes>(localMachNumberSquared, rCurrentProcessInfo);

    // ratio of speed of sound to free stream speed of sound
    const double speed_of_sound_ratio = ComputeSquaredSpeedofSoundFactor<Dim, NumNodes>(local_velocity_squared, rCurrentProcessInfo);

    const double free_stream_values_const = -0.5 * free_stream_density * std::pow(free_stream_mach,2.0) / inner_prod(free_stream_velocity, free_stream_velocity);
    const double speed_of_sound_power = (2.0 - heat_capacity_ratio) / (heat_capacity_ratio - 1.0);

    KRATOS_ERROR_IF((heat_capacity_ratio - 1.0) < std::numeric_limits<double>::epsilon())
        << "ComputeDensityDerivativeWRTVelocitySquared: heat capacity ratio must not be 1.0." << std::endl;

    return free_stream_values_const * std::pow(speed_of_sound_ratio, speed_of_sound_power);
}

template <int Dim, int NumNodes>
double ComputeUpwindedDensityDerivativeWRTVelocitySquaredSupersonicAccelerating(
    const array_1d<double, Dim>& rCurrentVelocity, 
    const double currentMachNumberSquared,
    const double upwindMachNumberSquared, 
    const ProcessInfo& rCurrentProcessInfo)
{
    // Following Fully Simulataneous Coupling of the Full Potential Equation
    //           and the Integral Boundary Layer Equations in Three Dimensions
    //           by Brian Nishida (1996), Section A.2.6
    const double upwind_factor = ComputeUpwindFactor<Dim, NumNodes>(currentMachNumberSquared, rCurrentProcessInfo);

    const double upwind_factor_derivative = ComputeUpwindFactorDerivativeWRTVelocitySquared<Dim, NumNodes>(rCurrentVelocity, rCurrentProcessInfo);

    const double Drho_Dq2 = ComputeDensityDerivativeWRTVelocitySquared<Dim, NumNodes>(currentMachNumberSquared, rCurrentProcessInfo);

    const double current_density = ComputeDensity<Dim, NumNodes>(currentMachNumberSquared, rCurrentProcessInfo);
    const double upwind_density = ComputeDensity<Dim, NumNodes>(upwindMachNumberSquared, rCurrentProcessInfo);

    return Drho_Dq2 * (1.0 - upwind_factor) - upwind_factor_derivative * (current_density - upwind_density);
}

template <int Dim, int NumNodes>
double ComputeUpwindedDensityDerivativeWRTVelocitySquaredSupersonicDeaccelerating(
    const double currentMachNumberSquared,
    const double upwindMachSquared,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Following Fully Simulataneous Coupling of the Full Potential Equation
    //           and the Integral Boundary Layer Equations in Three Dimensions
    //           by Brian Nishida (1996), Section A.2.6    
    // const double current_mach_sq = ComputeLocalMachNumberSquared<Dim, NumNodes>(rCurrentVelocity, rCurrentProcessInfo);
    const double Drho_Dq2 = ComputeDensityDerivativeWRTVelocitySquared<Dim, NumNodes>(currentMachNumberSquared, rCurrentProcessInfo);

    // const double upwind_mach_sq = ComputeLocalMachNumberSquared<Dim, NumNodes>(rUpwindVelocity, rCurrentProcessInfo);
    const double upwind_factor = ComputeUpwindFactor<Dim, NumNodes>(upwindMachSquared, rCurrentProcessInfo);

    return Drho_Dq2 * (1.0 - upwind_factor);
}

template <int Dim, int NumNodes>
double ComputeUpwindedDensityDerivativeWRTUpwindVelocitySquaredSupersonicAccelerating(
    const double currentMachNumberSquared,
    const double upwindMachNumberSquared, 
    const ProcessInfo& rCurrentProcessInfo)
{
    // Following Fully Simulataneous Coupling of the Full Potential Equation
    //           and the Integral Boundary Layer Equations in Three Dimensions
    //           by Brian Nishida (1996), Section A.2.6
    const double Drho_Dq2 = ComputeDensityDerivativeWRTVelocitySquared<Dim, NumNodes>(upwindMachNumberSquared, rCurrentProcessInfo);

    const double upwind_factor = ComputeUpwindFactor<Dim, NumNodes>(currentMachNumberSquared, rCurrentProcessInfo);

    return upwind_factor * Drho_Dq2;
}

template <int Dim, int NumNodes>
double ComputeUpwindedDensityDerivativeWRTUpwindVelocitySquaredSupersonicDeaccelerating(
    const array_1d<double, Dim>& rUpwindVelocity,
    const double currentMachNumberSquared,
    const double upwindMachNumberSquared, 
    const ProcessInfo& rCurrentProcessInfo)
{
    // Following Fully Simulataneous Coupling of the Full Potential Equation
    //           and the Integral Boundary Layer Equations in Three Dimensions
    //           by Brian Nishida (1996), Section A.2.6
    const double upwind_factor = ComputeUpwindFactor<Dim, NumNodes>(upwindMachNumberSquared, rCurrentProcessInfo);

    const double upwind_factor_derivative = ComputeUpwindFactorDerivativeWRTVelocitySquared<Dim, NumNodes>(rUpwindVelocity, rCurrentProcessInfo);

    const double Drho_Dq2 = ComputeDensityDerivativeWRTVelocitySquared<Dim, NumNodes>(upwindMachNumberSquared, rCurrentProcessInfo);

    const double current_density = ComputeDensity<Dim, NumNodes>(currentMachNumberSquared, rCurrentProcessInfo);
    const double upwind_density = ComputeDensity<Dim, NumNodes>(upwindMachNumberSquared, rCurrentProcessInfo);

    return upwind_factor * Drho_Dq2 - upwind_factor_derivative * (current_density - upwind_density);
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

template <int Dim,int NumNodes>
void GetSortedIds(std::vector<size_t>& Ids,
    const GeometryType& rGeom)
{
    Ids.resize(rGeom.PointsNumber());
    for (unsigned int i = 0; i < Ids.size(); i++)
    {
        Ids[i] = rGeom[i].Id();
    }
    std::sort(Ids.begin(), Ids.end());
}

template <int Dim, int NumNodes>
void GetNodeNeighborElementCandidates(GlobalPointersVector<Element>& ElementCandidates, const GeometryType& rGeom)
{
    for (int i = 0; i < Dim; i++)
    {
        const GlobalPointersVector<Element>& rNodeElementCandidates =
            rGeom[i].GetValue(NEIGHBOUR_ELEMENTS);
        for (unsigned int j = 0; j < rNodeElementCandidates.size(); j++)
        {
            ElementCandidates.push_back(rNodeElementCandidates(j));
        }
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// Template instantiation

// 2D
template array_1d<double, 3> GetWakeDistances<2, 3>(const Element& rElement);
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
template array_1d<double, 2> ComputePerturbedVelocity<2,3>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeMaximumVelocitySquared<2, 3>(const ProcessInfo& rCurrentProcessInfo);
template double ComputeClampedVelocitySquared<2, 3>(const array_1d<double, 2>& rVelocity, const ProcessInfo& rCurrentProcessInfo);
template double ComputeVelocityMagnitude<2, 3>(const double localMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);
template double ComputeIncompressiblePressureCoefficient<2, 3>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputePerturbationIncompressiblePressureCoefficient<2, 3>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeCompressiblePressureCoefficient<2, 3>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputePerturbationCompressiblePressureCoefficient<2, 3>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeLocalSpeedOfSound<2, 3>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeLocalSpeedofSoundSquared<2, 3>(const array_1d<double, 2>& rVelocity,const ProcessInfo& rCurrentProcessInfo);
template double ComputeSquaredSpeedofSoundFactor<2, 3>(const double localVelocitySquared, const ProcessInfo& rCurrentProcessInfo);
template double ComputePerturbationLocalSpeedOfSound<2, 3>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeLocalMachNumber<2, 3>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeLocalMachNumberSquared<2, 3>(const array_1d<double, 2>& rVelocity, const ProcessInfo& rCurrentProcessInfo);
template double ComputeDerivativeLocalMachSquaredWRTVelocitySquared<2, 3>(const array_1d<double, 2>& rVelocity, const double localMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);
template double ComputePerturbationLocalMachNumber<2, 3>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeUpwindFactor<2,3>(double localMachNumberSquared,const ProcessInfo& rCurrentProcessInfo);
template double SelectMaxUpwindFactor<2, 3>(const array_1d<double, 2>& rCurrentVelocity, const array_1d<double, 2>& rUpwindVelocity, const ProcessInfo& rCurrentProcessInfo);
template size_t ComputeUpwindFactorCase<2, 3>(array_1d<double, 3>& rUpwindFactorOptions);
template double ComputeUpwindFactorDerivativeWRTMachSquared<2,3>(const double localMachNumberSquared,const ProcessInfo& rCurrentProcessInfo);
template double ComputeUpwindFactorDerivativeWRTVelocitySquared<2,3>(const array_1d<double, 2>& rVelocity,const ProcessInfo& rCurrentProcessInfo);
template double ComputeDensity<2, 3>(const double localMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);
template double ComputeUpwindedDensity<2,3>(const array_1d<double, 2>& rCurrentVelocity, const array_1d<double, 2>& rUpwindVelocity, const ProcessInfo& rCurrentProcessInfo);
template double ComputeDensityDerivativeWRTVelocitySquared<2,3>(const double localVelocitySquared, const ProcessInfo& rCurrentProcessInfo);
template double ComputeUpwindedDensityDerivativeWRTVelocitySquaredSupersonicAccelerating<2,3>(const array_1d<double, 2>& rCurrentVelocity, const double currentMachNumberSquared, const double upwindMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);
template double ComputeUpwindedDensityDerivativeWRTVelocitySquaredSupersonicDeaccelerating<2,3>(const double currentMachNumberSquared, const double upwindMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);
template double ComputeUpwindedDensityDerivativeWRTUpwindVelocitySquaredSupersonicAccelerating<2,3>(const double currentMachNumberSquared, const double upwindMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);
template double ComputeUpwindedDensityDerivativeWRTUpwindVelocitySquaredSupersonicDeaccelerating<2,3>(const array_1d<double, 2>& rUpwindVelocity, const double currentMachNumberSquared, const double upwindMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);
template bool CheckIfElementIsCutByDistance<2, 3>(const BoundedVector<double, 3>& rNodalDistances);
template void KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) CheckIfWakeConditionsAreFulfilled<2>(const ModelPart&, const double& rTolerance, const int& rEchoLevel);
template bool CheckWakeCondition<2, 3>(const Element& rElement, const double& rTolerance, const int& rEchoLevel);
template void GetSortedIds<2, 3>(std::vector<size_t>& Ids, const GeometryType& rGeom);
template void GetNodeNeighborElementCandidates<2, 3>(GlobalPointersVector<Element>& ElementCandidates, const GeometryType& rGeom);
// 3D
template array_1d<double, 4> GetWakeDistances<3, 4>(const Element& rElement);
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
template array_1d<double, 3> ComputePerturbedVelocity<3,4>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeMaximumVelocitySquared<3, 4>(const ProcessInfo& rCurrentProcessInfo);
template double ComputeClampedVelocitySquared<3, 4>(const array_1d<double, 3>& rVelocity, const ProcessInfo& rCurrentProcessInfo);
template double ComputeVelocityMagnitude<3, 4>(const double localMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);
template double ComputeIncompressiblePressureCoefficient<3, 4>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputePerturbationIncompressiblePressureCoefficient<3, 4>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeCompressiblePressureCoefficient<3, 4>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputePerturbationCompressiblePressureCoefficient<3, 4>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeLocalSpeedOfSound<3, 4>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeLocalSpeedofSoundSquared<3, 4>(const array_1d<double, 3>& rVelocity,const ProcessInfo& rCurrentProcessInfo);
template double ComputeSquaredSpeedofSoundFactor<3, 4>(const double localVelocitySquared, const ProcessInfo& rCurrentProcessInfo);
template double ComputePerturbationLocalSpeedOfSound<3, 4>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeLocalMachNumber<3, 4>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeLocalMachNumberSquared<3, 4>(const array_1d<double, 3>& rVelocity, const ProcessInfo& rCurrentProcessInfo);
template double ComputeDerivativeLocalMachSquaredWRTVelocitySquared<3, 4>(const array_1d<double, 3>& rVelocity, const double localMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);
template double ComputePerturbationLocalMachNumber<3, 4>(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);
template double ComputeUpwindFactor<3, 4>(double localMachNumberSquared,const ProcessInfo& rCurrentProcessInfo);
template double SelectMaxUpwindFactor<3, 4>(const array_1d<double, 3>& rCurrentVelocity, const array_1d<double, 3>& rUpwindVelocity, const ProcessInfo& rCurrentProcessInfo);
template size_t ComputeUpwindFactorCase<3, 4>(array_1d<double, 3>& rUpwindFactorOptions);
template double ComputeUpwindFactorDerivativeWRTMachSquared<3,4>(const double localMachNumberSquared,const ProcessInfo& rCurrentProcessInfo);
template double ComputeUpwindFactorDerivativeWRTVelocitySquared<3,4>(const array_1d<double, 3>& rVelocity,const ProcessInfo& rCurrentProcessInfo);
template double ComputeDensity<3, 4>(const double localMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);
template double ComputeUpwindedDensity<3, 4>(const array_1d<double, 3>& rCurrentVelocity, const array_1d<double, 3>& rUpwindVelocity, const ProcessInfo& rCurrentProcessInfo);
template double ComputeDensityDerivativeWRTVelocitySquared<3,4>(const double localVelocitySquared, const ProcessInfo& rCurrentProcessInfo);
template double ComputeUpwindedDensityDerivativeWRTVelocitySquaredSupersonicAccelerating<3,4>(const array_1d<double, 3>& rCurrentVelocity, const double currentMachNumberSquared, const double upwindMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);
template double ComputeUpwindedDensityDerivativeWRTVelocitySquaredSupersonicDeaccelerating<3,4>(const double currentMachNumberSquared, const double upwindMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);
template double ComputeUpwindedDensityDerivativeWRTUpwindVelocitySquaredSupersonicAccelerating<3,4>(const double currentMachNumberSquared, const double upwindMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);
template double ComputeUpwindedDensityDerivativeWRTUpwindVelocitySquaredSupersonicDeaccelerating<3,4>(const array_1d<double, 3>& rUpwindVelocity, const double currentMachNumberSquared, const double upwindMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);
template bool CheckIfElementIsCutByDistance<3, 4>(const BoundedVector<double, 4>& rNodalDistances);
template void  KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) CheckIfWakeConditionsAreFulfilled<3>(const ModelPart&, const double& rTolerance, const int& rEchoLevel);
template bool CheckWakeCondition<3, 4>(const Element& rElement, const double& rTolerance, const int& rEchoLevel);
template void GetSortedIds<3, 4>(std::vector<size_t>& Ids, const GeometryType& rGeom);
template void GetNodeNeighborElementCandidates<3, 4>(GlobalPointersVector<Element>& ElementCandidates, const GeometryType& rGeom);
} // namespace PotentialFlow
} // namespace Kratos
