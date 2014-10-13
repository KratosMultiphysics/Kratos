//
//   Project Name:        Kratos
//   Last Modified by:    $Author: G.Casas $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "spheric_swimming_particle.h"
#include "../applications/DEM_application/custom_utilities/GeometryFunctions.h"
#include "../applications/DEM_application/custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application.h"

namespace Kratos
{

SphericSwimmingParticle::SphericSwimmingParticle(): SphericParticle(){}

SphericSwimmingParticle::SphericSwimmingParticle(IndexType NewId, GeometryType::Pointer pGeometry): SphericParticle(NewId, pGeometry){}

SphericSwimmingParticle::SphericSwimmingParticle( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
: SphericParticle(NewId, pGeometry, pProperties){}

SphericSwimmingParticle::SphericSwimmingParticle(IndexType NewId, NodesArrayType const& ThisNodes)
: SphericParticle(NewId, ThisNodes){}

Element::Pointer SphericSwimmingParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
 return Element::Pointer(new SphericSwimmingParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

/// Destructor.
SphericSwimmingParticle::~SphericSwimmingParticle(){}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

void SphericSwimmingParticle::ComputeAdditionalForces(array_1d<double, 3>& additionally_applied_force,
                                                    array_1d<double, 3>& additionally_applied_moment,
                                                    ProcessInfo& rCurrentProcessInfo,
                                                    const array_1d<double,3>& gravity)
{
    KRATOS_TRY

    //const array_1d<double, 3>& gravity = rCurrentProcessInfo[GRAVITY];
    const double& fluid_density        = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_DENSITY_PROJECTED);
    const double mass                  = mSqrtOfRealMass * mSqrtOfRealMass;

    array_1d<double, 3> buoyancy;
    array_1d<double, 3> drag_force;
    array_1d<double, 3> lift_force;
    array_1d<double, 3> virtual_mass_force;

    // The decomposition of forces that is considered here follows Jackson (The Dynamics of Fluidized Particles, 2000);
    // so that the role of f_n1 therein is played by additionally_applied_force here

    ComputeBuoyancy(buoyancy, fluid_density, gravity, rCurrentProcessInfo);
    ComputeDragForce(drag_force, fluid_density, rCurrentProcessInfo);
    ComputeVirtualMassForce(virtual_mass_force, fluid_density, rCurrentProcessInfo);
    ComputeLiftForce(lift_force, fluid_density, rCurrentProcessInfo);

    additionally_applied_force[0] = drag_force[0] + virtual_mass_force[0] + lift_force[0];
    additionally_applied_force[1] = drag_force[1] + virtual_mass_force[1] + lift_force[1];
    additionally_applied_force[2] = drag_force[2] + virtual_mass_force[2] + lift_force[2];

    UpdateNodalValues(additionally_applied_force, buoyancy, drag_force, virtual_mass_force, lift_force);

    additionally_applied_force[0] += buoyancy[0] + mass * gravity[0];
    additionally_applied_force[1] += buoyancy[1] + mass * gravity[1];
    additionally_applied_force[2] += buoyancy[2] + mass * gravity[2];

    KRATOS_CATCH( "" )
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
// Here nodal values are modified to record DEM forces that we want to print. In Kratos this is an exception since nodal values are meant to be modified only outside the element. Here it was not possible.

void SphericSwimmingParticle::UpdateNodalValues(const array_1d<double, 3>& hydrodynamic_force,
                                                const array_1d<double, 3>& buoyancy,
                                                const array_1d<double, 3>& drag_force,
                                                const array_1d<double, 3>& virtual_mass_force,
                                                const array_1d<double, 3>& lift_force)
{
    GetGeometry()(0)->FastGetSolutionStepValue(HYDRODYNAMIC_FORCE)     = hydrodynamic_force;
    GetGeometry()(0)->FastGetSolutionStepValue(BUOYANCY)               = buoyancy;

    if (mHasDragForceNodalVar){
        GetGeometry()(0)->FastGetSolutionStepValue(DRAG_FORCE)         = drag_force;
    }

    if (mHasVirtualMassForceNodalVar){
        GetGeometry()(0)->FastGetSolutionStepValue(VIRTUAL_MASS_FORCE) = virtual_mass_force;
    }

    if (mHasLiftForceNodalVar){
        GetGeometry()(0)->FastGetSolutionStepValue(LIFT_FORCE)         = lift_force;
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

void SphericSwimmingParticle::ComputeBuoyancy(array_1d<double, 3>& buoyancy, const double& fluid_density, const array_1d<double, 3>& gravity, ProcessInfo& rCurrentProcessInfo)
{

    if (mBuoyancyForceType == 0 || GetGeometry()[0].IsNot(INSIDE) || GetGeometry()[0].Is(BLOCKED)){ // case of identically null buoyancy
        noalias(buoyancy) = ZeroVector(3);
        return;
    }

    else {
        const double volume = 1.333333333333333 * KRATOS_M_PI * mRadius * mRadius * mRadius;

        if (mDragForceType == 2){ // Weatherford
            noalias(buoyancy) =  - gravity * fluid_density * volume;
        }

        else {
            const array_1d<double, 3>& pressure_grad = GetGeometry()(0)->FastGetSolutionStepValue(PRESSURE_GRAD_PROJECTED);
            noalias(buoyancy) = - volume * pressure_grad;
        }
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

void SphericSwimmingParticle::ComputeDragForce(array_1d<double, 3>& drag_force, const double& fluid_density, ProcessInfo& rCurrentProcessInfo)
{
    if (mDragForceType == 0 || GetGeometry()[0].IsNot(INSIDE) || GetGeometry()[0].Is(BLOCKED)){ // case of identically null drag force
        noalias(drag_force) = ZeroVector(3);
        return;
    }

    else {
        const array_1d<double, 3> fluid_vel     = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
        const array_1d<double, 3>& particle_vel = GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3> slip_vel;

        if (mFluidModelType == 0){ // fluid velocity is modified as a post-process
            const double fluid_fraction         = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_FRACTION_PROJECTED);
            noalias(slip_vel)                   = fluid_vel / fluid_fraction - particle_vel;
        }

        else {
            noalias(slip_vel)                   = fluid_vel - particle_vel;
        }

        const double norm_of_slip_vel           = MathUtils<double>::Norm3(slip_vel);
        double drag_coeff;

        // calculating the 'dimensional' drag coefficient, i.e., the factor by which the slip velocity must be multiplied to yield the drag force

        if (mDragForceType == 1){
            drag_coeff = ComputeStokesDragCoefficient(fluid_density, rCurrentProcessInfo);
        }

        else if (mDragForceType == 2){ // formulations of Haider (1989) and Chien (1994)
            drag_coeff = ComputeWeatherfordDragCoefficient(norm_of_slip_vel, fluid_density, rCurrentProcessInfo);
        }

        else if (mDragForceType == 3){ // formulation of Ganser (1993)
            drag_coeff = ComputeGanserDragCoefficient(norm_of_slip_vel, fluid_density, rCurrentProcessInfo);
        }

        else if (mDragForceType == 4){ // formulation of Ishii and Zuber (1979)
            drag_coeff = ComputeIshiiDragCoefficient(norm_of_slip_vel, fluid_density, rCurrentProcessInfo);
        }

        else if (mDragForceType == 5){ // Newton regime (Re ~ 1000 - 250000), formulation of Haider and Levenspiel (1989)
            drag_coeff = ComputeNewtonRegimeDragCoefficient(norm_of_slip_vel, fluid_density);
        }

        else {
            std::cout << "The integer value designating the drag coefficient calculation model" << std::endl;
            std::cout << " (mDragForceType = " << mDragForceType << "), is not supported" << std::endl << std::flush;
            return;
        }

        noalias(drag_force) = drag_coeff * slip_vel;
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

void SphericSwimmingParticle::ComputeVirtualMassForce(array_1d<double, 3>& virtual_mass_force, const double& fluid_density, ProcessInfo& rCurrentProcessInfo)
{
    if (mVirtualMassForceType == 0 || GetGeometry()[0].IsNot(INSIDE) || GetGeometry()[0].Is(BLOCKED)){ // case of identically null virtual mass force
        noalias(virtual_mass_force) = ZeroVector(3);
        return;
    }

    else {
        const double delta_t_inv                = 1 / rCurrentProcessInfo[DELTA_TIME];
        array_1d<double, 3> fluid_acc           = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_ACCEL_PROJECTED);
        const array_1d<double, 3>& particle_acc = delta_t_inv * (GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY) - GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY, 1));
        const double fluid_fraction             = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_FRACTION_PROJECTED);
        array_1d<double, 3> slip_acc;

    if (mFluidModelType == 0){ // fluid velocity is modified as a post-process
        noalias(slip_acc) = fluid_acc / fluid_fraction - particle_acc;
    }

    else {
        noalias(slip_acc) = fluid_acc - particle_acc;
    }

    double virtual_mass_coeff = 0.5; // inviscid case

    if (mVirtualMassForceType == 2) { // Zuber (1964) (moderate values of solid fraction)
        virtual_mass_coeff = 0.5 + 1.5 * (1 - fluid_fraction);
    }

    noalias(virtual_mass_force) = virtual_mass_coeff * fluid_density * slip_acc;
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

void SphericSwimmingParticle::ComputeLiftForce(array_1d<double, 3>& lift_force, const double& fluid_density, ProcessInfo& rCurrentProcessInfo)
{
    if (mLiftForceType == 0 || GetGeometry()[0].IsNot(INSIDE) || GetGeometry()[0].Is(BLOCKED)){ // case of identically null lift force
        noalias(lift_force) = ZeroVector(3);

        return;
    }

    else {
        const array_1d<double, 3>& fluid_vel    = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
        const array_1d<double, 3>& particle_vel = GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3> minus_slip_vel;

        if (mFluidModelType == 0){ // fluid velocity is modified as a post-process
            const double fluid_fraction = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_FRACTION_PROJECTED);
            noalias(minus_slip_vel)     = particle_vel - fluid_vel / fluid_fraction;
        }

        else {
            noalias(minus_slip_vel) = particle_vel - fluid_vel;
        }

        const double& shear_rate                      = GetGeometry()(0)->FastGetSolutionStepValue(SHEAR_RATE_PROJECTED);
        const array_1d<double, 3> vorticity           = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_VORTICITY_PROJECTED);
        const double vorticity_norm                   = MathUtils<double>::Norm3(vorticity);
        const double norm_of_slip_vel                 = MathUtils<double>::Norm3(minus_slip_vel);
        const array_1d<double, 3> vort_cross_slip_vel = MathUtils<double>::CrossProduct(vorticity, minus_slip_vel);
        double lift_coeff;

        if (mLiftForceType == 1){ // El Samni, E.A. (1949), see paper by R. K. Clark (1994)
            lift_coeff = ComputeSaffmanLiftCoefficient(norm_of_slip_vel, fluid_density, shear_rate, vorticity_norm, rCurrentProcessInfo);
        }

        else {
            std::cout << "The integer value designating the lift coefficient calculation model" << std::endl;
            std::cout << " (mLiftForceType = " << mLiftForceType << "), is not supported" << std::endl << std::flush;
            return;
        }

        noalias(lift_force) = lift_coeff * vort_cross_slip_vel; // the direction is given by the vorticity x (- slip_vel) (Jackson, 2000), which is normalized here
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

void SphericSwimmingParticle::ComputeParticleReynoldsNumber(double norm_of_slip_vel, double kinematic_viscosity, double& reynolds)
{
    reynolds = 2 * mRadius * norm_of_slip_vel / kinematic_viscosity;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

void SphericSwimmingParticle::AdditionalCalculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == REYNOLDS_NUMBER){

        if (GetGeometry()[0].IsNot(INSIDE)){
            Output = 0.0;
        }

        else {
            const array_1d<double, 3>& fluid_vel    = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
            const array_1d<double, 3>& particle_vel = GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3> slip_vel;

        if (mFluidModelType == 0){ // fluid velocity is modified as a post-process
            const double fluid_fraction         = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_FRACTION_PROJECTED);
            noalias(slip_vel)                   = fluid_vel / fluid_fraction - particle_vel;
        }

        else {
            noalias(slip_vel)                   = fluid_vel - particle_vel;
        }

        const double norm_of_slip_vel           = MathUtils<double>::Norm3(slip_vel);
        const double kinematic_viscosity        = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_VISCOSITY_PROJECTED);
        ComputeParticleReynoldsNumber(norm_of_slip_vel, kinematic_viscosity, Output);
        }
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

double SphericSwimmingParticle::ComputeStokesDragCoefficient(const double fluid_density, ProcessInfo& rCurrentProcessInfo)
{
    const double kinematic_viscosity = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_VISCOSITY_PROJECTED);
    const double fluid_fraction      = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_FRACTION_PROJECTED);

    double  drag_coeff  = 6.0 * KRATOS_M_PI * kinematic_viscosity * fluid_density * fluid_fraction * mRadius;

    return drag_coeff;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

double SphericSwimmingParticle::ComputeWeatherfordDragCoefficient(const double& norm_of_slip_vel, const double fluid_density, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //const double& particle_density             = GetGeometry()(0)->GetSolutionStepValue(PARTICLE_DENSITY);
    const double particle_density = GetDensity();
    const double kinematic_viscosity           = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_VISCOSITY_PROJECTED);
    const double sphericity                    = GetGeometry()(0)->GetSolutionStepValue(PARTICLE_SPHERICITY);
    //const array_1d<double, 3>& buoyancy       = GetGeometry()(0)->FastGetSolutionStepValue(BUOYANCY);//S
    const array_1d<double, 3>& gravity         = rCurrentProcessInfo[GRAVITY];
    const int manually_imposed_drag_law_option = rCurrentProcessInfo[MANUALLY_IMPOSED_DRAG_LAW_OPTION];
    const int drag_modifier_type               = rCurrentProcessInfo[DRAG_MODIFIER_TYPE];
    const double gel_strength                  = GetGeometry()(0)->FastGetSolutionStepValue(GEL_STRENGTH);
    const double power_law_n                   = GetGeometry()(0)->FastGetSolutionStepValue(POWER_LAW_N);
    const double power_law_K                   = GetGeometry()(0)->FastGetSolutionStepValue(POWER_LAW_K);
    const double yield_stress                  = GetGeometry()(0)->FastGetSolutionStepValue(YIELD_STRESS);

    int non_newtonian_option = 1;

    if (fabs(power_law_n - 1.0) < 0.00001  ||  fabs(yield_stress) < 0.00001) {
        non_newtonian_option = 0;
    }

    const double initial_drag_force            = rCurrentProcessInfo[INIT_DRAG_FORCE];
    const double drag_law_slope                = rCurrentProcessInfo[DRAG_LAW_SLOPE];
    const double power_law_tol                 = rCurrentProcessInfo[POWER_LAW_TOLERANCE];

    const double area                          = KRATOS_M_PI * mRadius * mRadius;
    const array_1d<double, 3> weight           = mSqrtOfRealMass * mSqrtOfRealMass * gravity;
    const array_1d<double, 3> buoyancy         = fluid_density / particle_density * weight; // hydrostatic case!! (only for Weatherford)

    double shahs_term_vel                      = 0.0;
    double beta                                = 0.0;
    double F0                                  = 0.0;
    double regularization_v                    = 0.02 * mRadius;
    double reynolds;
    double drag_coeff;

    if (!non_newtonian_option){ // Newtonian
        ComputeParticleReynoldsNumber(norm_of_slip_vel, kinematic_viscosity, reynolds);

        if (!non_newtonian_option && reynolds < 0.01){
            reynolds = 0.01;
        }

        CalculateNewtonianDragCoefficient(non_newtonian_option, reynolds, sphericity, drag_coeff, drag_modifier_type);
        drag_coeff = 0.5 * fluid_density * area * drag_coeff * norm_of_slip_vel;
    }

    else {
        shahs_term_vel = CalculateShahsTerm(power_law_n, power_law_K, power_law_tol, fluid_density, particle_density, sphericity, drag_modifier_type);

        if (!manually_imposed_drag_law_option){
            F0 = 4.0 * gel_strength * area; //initial value
            beta = (MathUtils<double>::Norm3(weight) - MathUtils<double>::Norm3(buoyancy) - F0) / shahs_term_vel; //S
        }

        else {
            F0 = initial_drag_force; //initial value
            beta = drag_law_slope; //slope
        }

        if (norm_of_slip_vel >= regularization_v){
            drag_coeff = (F0 + beta * norm_of_slip_vel) / norm_of_slip_vel;
        }

        else {
            drag_coeff = (F0 + beta * regularization_v) / regularization_v;
        }
    }

    return drag_coeff;

    KRATOS_CATCH("")
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

void SphericSwimmingParticle::CalculateNewtonianDragCoefficient(int non_newtonian_option,
                                                                const double reynolds,
                                                                double sphericity,
                                                                double& drag_coeff,
                                                                int drag_modifier_type)
{
    if (reynolds < 1){
        drag_coeff = 24.0; // Reynolds;
    }

    else {

        if (reynolds > 1000){
            drag_coeff = 0.44;
        }

        else{
            drag_coeff = 24.0 / reynolds * (1.0 + 0.15 * pow(reynolds, 0.687));
        }
    }

    if (!non_newtonian_option){ // Newtonian

        if (drag_coeff > 2.0){
            drag_coeff = 2.0; // watch out!
        }
    }

    if (sphericity < 0.9999){
        drag_coeff = CalculateDragCoeffFromSphericity(reynolds, sphericity, drag_modifier_type);
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

double SphericSwimmingParticle::CalculateDragCoeffFromSphericity(const double reynolds,
                                                                 const double sphericity,
                                                                 const int drag_modifier_type)
{
    double cdrag = 0.0;

    if (drag_modifier_type == 1){ // visual-Red Book
        double interpolator   = (1 - sphericity) / (1 - 0.806);
        double cdrag_modifier = 1 + 0.97 * interpolator + 0.715 * interpolator * log10(reynolds);

        if (reynolds < 1){
            cdrag_modifier += 0.3 * interpolator * pow(- 1.0 * log10(reynolds), 1.6);
        }

        cdrag = cdrag_modifier * cdrag;
    }

    if (drag_modifier_type == 2){ // Hayder
        cdrag = 24 / reynolds * (1 + exp(2.3288 - 6.4581 * sphericity + 2.4486 * sphericity * sphericity) * pow(reynolds, 0.0964 + 0.5565 * sphericity)) + 73.69 * reynolds * exp(- 5.0748 * sphericity) / (reynolds + 5.378 * exp(6.2122 * sphericity));
    }

    if (drag_modifier_type == 3){ // Chien
        cdrag = 30 / reynolds + 67.289 * exp(- 5.03 * sphericity);
    }

    return cdrag;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

double SphericSwimmingParticle::CalculateShahsTerm(double power_law_N,
                                                   double power_law_K,
                                                   double power_law_tol,
                                                   const double& fluid_density,
                                                   const double& particle_density,
                                                   double sphericity,
                                                   int drag_modifier_type)
{
    if (fabs(power_law_N) < power_law_tol || fabs(power_law_K) < power_law_tol){
        std::cout << "WARNING: Shah's method is being used with Power Law data being zero!!" << std::endl << std::flush;
    }

    double shah_A_i = 1 / (6.9148 * power_law_N * power_law_N - 24.838 * power_law_N + 22.642);
    double shah_B_i = 1 / (-0.5067 * power_law_N * power_law_N + 1.3234 * power_law_N - 0.1744);

    double dimensionless_shah = sqrt(pow(13.08, 2 - power_law_N) * pow(2 * mRadius, power_law_N + 2) * pow(fluid_density, power_law_N) * pow(particle_density - fluid_density, 2 - power_law_N) / (pow(2, 2 * (power_law_N - 1)) * power_law_K * power_law_K));
    double reynolds = pow(dimensionless_shah * shah_A_i, shah_B_i);
    double fi_i = CalculateDragCoeffFromSphericity(reynolds, 1.0, drag_modifier_type) / CalculateDragCoeffFromSphericity(reynolds, sphericity, drag_modifier_type);
    dimensionless_shah = sqrt(pow(fi_i, 2 - power_law_N)) * dimensionless_shah;
    reynolds = pow(dimensionless_shah * shah_A_i, shah_B_i);

    double terminal_vel =  pow(pow(2, power_law_N - 1) * power_law_K * reynolds / (pow(2 * mRadius, power_law_N) * fluid_density), 1 / (2 - power_law_N)) ;

    return terminal_vel;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

double SphericSwimmingParticle::ComputeGanserDragCoefficient(const double& norm_of_slip_vel, const double fluid_density, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const double kinematic_viscosity         = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_VISCOSITY_PROJECTED);
    const double sphericity                  = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_SPHERICITY);
    const int isometric_shape                = 1; // TEMPORARY!! yes (1) or no (0); shold be given as data
    const double surface_area                = 4 * KRATOS_M_PI * mRadius * mRadius; // TEMPORARY!! corresponding to a sphere; should be generalized b taking it as a parameter
    const double surface_area_circular_diam  = sqrt(4.0 * surface_area / KRATOS_M_PI);

    double equiv_reynolds;
    double k_1;
    double k_2;
    double drag_coeff;

    ComputeGanserParameters(isometric_shape, sphericity, surface_area_circular_diam, k_1, k_2);
    ComputeParticleReynoldsNumber(norm_of_slip_vel, kinematic_viscosity, equiv_reynolds);
    equiv_reynolds *= k_1 * k_2;

    // calculating adimensional drag coefficient

    drag_coeff =  k_2 * (24 * (1 + 0.1118 * pow((equiv_reynolds), 0.6567)) / (equiv_reynolds) + 0.4305 / (1 + 3305 / equiv_reynolds));

    // and then the dimensional drag coefficient

    drag_coeff *= 0.5 * fluid_density * surface_area * norm_of_slip_vel;

    return drag_coeff;

    KRATOS_CATCH("")
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

double SphericSwimmingParticle::ComputeIshiiDragCoefficient(const double& norm_of_slip_vel, const double fluid_density, ProcessInfo& rCurrentProcessInfo)
{
    const double kinematic_viscosity = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_VISCOSITY_PROJECTED);
    double coeff = 0.45;
    double reynolds;
    ComputeParticleReynoldsNumber(norm_of_slip_vel, kinematic_viscosity, reynolds);

    if (reynolds <= 1000){
        coeff = (24 + 2.4 * pow(reynolds, 0.75)) / reynolds;
    }

    double drag_coeff = 0.5 * coeff * KRATOS_M_PI * mRadius * mRadius;

    return drag_coeff;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

double SphericSwimmingParticle::ComputeNewtonRegimeDragCoefficient(const double& norm_of_slip_vel, const double fluid_density)
{
    const double sphericity = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_SPHERICITY);
    const double radius     = GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
    double drag_coeff       = 0.25 * KRATOS_M_PI * radius * radius * fluid_density * norm_of_slip_vel;

    if (sphericity >= 1.0){
        drag_coeff *= 44.0;

        return drag_coeff;
    }

    else {
        const double kinematic_viscosity = GetGeometry()(0)->FastGetSolutionStepValue(FLUID_VISCOSITY_PROJECTED);
        double A = exp(2.3288 - 6.4581 * sphericity + 2.4486 * sphericity * sphericity);
        double B = 0.0964 + 0.5565 * sphericity;
        double C = exp(4.905  - 13.8944 * sphericity + 18.4222 * sphericity * sphericity - 10.2599 * sphericity * sphericity * sphericity);
        double D = exp(1.4681 + 12.2584 * sphericity - 20.7322 * sphericity * sphericity + 15.8855 * sphericity * sphericity * sphericity);
        double particle_reynolds;
        ComputeParticleReynoldsNumber(norm_of_slip_vel, kinematic_viscosity, particle_reynolds);
        drag_coeff *= (24.0 * (1.0 + A * pow(particle_reynolds, B))) / particle_reynolds + C * particle_reynolds / (particle_reynolds + D);

        return drag_coeff;
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

 void SphericSwimmingParticle::ComputeGanserParameters(const int isometric_shape, const double sphericity, const double dn, double& k_1, double& k_2)
 {
     if (isometric_shape){
         k_1 = 3 / (1 + 2 / sqrt(sphericity));
     }

     else {
         k_1 = 3 / (0.5 * dn / mRadius + 2 / sqrt(sphericity));
     }

     k_2 = pow(10.0, 1.8148 * pow(- log10(sphericity), 0.5743));
 }

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

 double SphericSwimmingParticle::ComputeSaffmanLiftCoefficient(const double& norm_of_slip_vel,
                                                               const double fluid_density,
                                                               const double norm_of_shear_rate,
                                                               const double vorticity_norm,
                                                               ProcessInfo& rCurrentProcessInfo)
 {
     if (vorticity_norm > 0.000000000001 && norm_of_slip_vel > 0.000000000001){
         const double yield_stress   = 0.0; // we are considering a Bingham type fluid
         const double power_law_K    = GetGeometry()(0)->FastGetSolutionStepValue(POWER_LAW_K);
         const double power_law_n    = GetGeometry()(0)->FastGetSolutionStepValue(POWER_LAW_N);
         const double shear_rate_p   = norm_of_slip_vel / mRadius * (4.5 / power_law_n - 3.5); // graphic model by Unhlherr et al. (fit by Wallis, G.B. and Dobson, J.E., 1973)
         double equivalent_viscosity = yield_stress / shear_rate_p + power_law_K * pow(shear_rate_p, power_law_n - 1);
         const double coeff          = std::max(0.09 * norm_of_slip_vel, 5.82 * sqrt(0.5 * norm_of_shear_rate * equivalent_viscosity / fluid_density));
         const double lift_coeff     = 0.5 * KRATOS_M_PI * mRadius * mRadius * fluid_density * coeff * norm_of_slip_vel / vorticity_norm;
         return(lift_coeff);
     }

     else {
         return 0.0;
     }
 }

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

 void SphericSwimmingParticle::CustomInitialize()
 {
     mHasDragForceNodalVar        = GetGeometry()(0)->SolutionStepsDataHas(DRAG_FORCE);
     mHasVirtualMassForceNodalVar = GetGeometry()(0)->SolutionStepsDataHas(VIRTUAL_MASS_FORCE);
     mHasLiftForceNodalVar        = GetGeometry()(0)->SolutionStepsDataHas(LIFT_FORCE);
 }

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

 void SphericSwimmingParticle::AdditionalMemberDeclarationFirstStep(const ProcessInfo& r_process_info)
 {
     mBuoyancyForceType    = r_process_info[BUOYANCY_FORCE_TYPE];
     mDragForceType        = r_process_info[DRAG_FORCE_TYPE];
     mVirtualMassForceType = r_process_info[VIRTUAL_MASS_FORCE_TYPE];
     mLiftForceType        = r_process_info[LIFT_FORCE_TYPE];
     mFluidModelType       = r_process_info[FLUID_MODEL_TYPE];
 }

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

}  // namespace Kratos.

