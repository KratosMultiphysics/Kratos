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

#define SWIMMING_COPY_SECOND_TO_FIRST_3(a, b)            a[0]  = b[0]; a[1]  = b[1]; a[2]  = b[2];
#define SWIMMING_ADD_SECOND_TO_FIRST(a, b)               a[0] += b[0]; a[1] += b[1]; a[2] += b[2];
#define SWIMMING_SET_COMPONENTS_TO_ZERO_3(a)             a[0]  = 0.0;  a[1]  = 0.0;  a[2]  = 0.0;
#define SWIMMING_SET_COMPONENTS_TO_ZERO_3x3(a)           a[0][0] = 0.0; a[0][1] = 0.0; a[0][2] = 0.0; a[1][0] = 0.0; a[1][1] = 0.0; a[1][2] = 0.0; a[2][0] = 0.0; a[2][1] = 0.0; a[2][2] = 0.0;
#define SWIMMING_MULTIPLY_BY_SCALAR_3(a, b)              a[0] = b * a[0]; a[1] = b * a[1]; a[2] = b * a[2];
#define SWIMMING_MODULUS_3(a)                            sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])
#define SWIMMING_INNER_PRODUCT_3(a, b)                       (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
#define SWIMMING_SET_TO_CROSS_OF_FIRST_TWO_3(a, b, c)    c[0] = a[1] * b[2] - a[2] * b[1]; c[1] = a[2] * b[0] - a[0] * b[2]; c[2] = a[0] * b[1] - a[1] * b[0];
#define SWIMMING_POW_2(a)                                a * a
#define SWIMMING_POW_3(a)                                a * a * a
#define SWIMMING_POW_4(a)                                a * a * a * a
#define SWIMMING_POW_5(a)                                a * a * a * a * a

namespace Kratos
{

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

template< class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeAdditionalForces(array_1d<double, 3>& additionally_applied_force,
                                                      array_1d<double, 3>& additionally_applied_moment,
                                                      ProcessInfo& r_current_process_info,
                                                      const array_1d<double,3>& gravity)
{
    KRATOS_TRY

    if(!mCouplingType) {
        TBaseElement::ComputeAdditionalForces(additionally_applied_force, additionally_applied_moment, r_current_process_info, gravity);
        return;
    }

    NodeType& node = GetGeometry()[0];
    mFluidDensity                           = node.FastGetSolutionStepValue(FLUID_DENSITY_PROJECTED);
    mKinematicViscosity                     = node.FastGetSolutionStepValue(FLUID_VISCOSITY_PROJECTED);
    mFluidFraction                          = node.FastGetSolutionStepValue(FLUID_FRACTION_PROJECTED);
    const array_1d<double, 3>& fluid_vel    = node.FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
    const array_1d<double, 3>& particle_vel = node.FastGetSolutionStepValue(VELOCITY);
    mLastTimeStep = 0.0;

    if (mFluidModelType == 0){ // fluid velocity is modified as a post-process
        noalias(mSlipVel) = fluid_vel / mFluidFraction - particle_vel;
    }

    else {
        noalias(mSlipVel) = fluid_vel - particle_vel;
    }

    mNormOfSlipVel = SWIMMING_MODULUS_3(mSlipVel);

    array_1d<double, 3> buoyancy;
    array_1d<double, 3> drag_force;
    array_1d<double, 3> virtual_mass_force;
    array_1d<double, 3> saffman_lift_force;
    array_1d<double, 3> magnus_lift_force;
    array_1d<double, 3> brownian_motion_force;

    // The decomposition of forces that is considered here follows Jackson (The Dynamics of Fluidized Particles, 2000);
    // so that the role of f_n1 therein is played by additionally_applied_force here

    ComputeBuoyancy(buoyancy, gravity, r_current_process_info);
    ComputeDragForce(drag_force, r_current_process_info);
    ComputeVirtualMassForce(virtual_mass_force, r_current_process_info);
    ComputeSaffmanLiftForce(saffman_lift_force, r_current_process_info);
    ComputeMagnusLiftForce(magnus_lift_force, r_current_process_info);
    ComputeHydrodynamicTorque(additionally_applied_moment, r_current_process_info);
    ComputeBrownianMotionForce(brownian_motion_force, r_current_process_info);

    noalias(additionally_applied_force) += drag_force + virtual_mass_force + saffman_lift_force + magnus_lift_force + brownian_motion_force;

    UpdateNodalValues(additionally_applied_force, additionally_applied_moment, buoyancy, drag_force, virtual_mass_force, saffman_lift_force, magnus_lift_force, r_current_process_info);
    
    //Now add the contribution of base class function (gravity or other forces added in upper levels):
    TBaseElement::ComputeAdditionalForces(additionally_applied_force, additionally_applied_moment, r_current_process_info, gravity);
	additionally_applied_force += buoyancy;

    KRATOS_CATCH( "" )
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
// Here nodal values are modified to record DEM forces that we want to print. In Kratos this is an exception since nodal values are meant to be modified only outside the element. Here it was not possible.
template< class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::UpdateNodalValues(const array_1d<double, 3>& hydrodynamic_force,
                                                const array_1d<double, 3>& hydrodynamic_moment,
                                                const array_1d<double, 3>& buoyancy,
                                                const array_1d<double, 3>& drag_force,
                                                const array_1d<double, 3>& virtual_mass_force,
                                                const array_1d<double, 3>& saffman_lift_force,
                                                const array_1d<double, 3>& magnus_lift_force,
                                                ProcessInfo& r_current_process_info)
{
    noalias(GetGeometry()[0].FastGetSolutionStepValue(HYDRODYNAMIC_FORCE))      = hydrodynamic_force;
    noalias(GetGeometry()[0].FastGetSolutionStepValue(BUOYANCY))                = buoyancy;

    if (mHasHydroMomentNodalVar){
        noalias(GetGeometry()[0].FastGetSolutionStepValue(HYDRODYNAMIC_MOMENT)) = hydrodynamic_moment;
    }

    if (mHasDragForceNodalVar){
        noalias(GetGeometry()[0].FastGetSolutionStepValue(DRAG_FORCE))          = drag_force;
    }

    if (mHasVirtualMassForceNodalVar){
        noalias(GetGeometry()[0].FastGetSolutionStepValue(VIRTUAL_MASS_FORCE))  = virtual_mass_force;
    }

    if (mHasLiftForceNodalVar){
        noalias(GetGeometry()[0].FastGetSolutionStepValue(LIFT_FORCE))          = saffman_lift_force + magnus_lift_force;
    }

    if (mHasDragCoefficientVar){
        double drag_coefficient = ComputeDragCoefficient(r_current_process_info);
        GetGeometry()[0].FastGetSolutionStepValue(DRAG_COEFFICIENT)    = drag_coefficient;
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

template< class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeBuoyancy(array_1d<double, 3>& buoyancy, const array_1d<double, 3>& gravity, ProcessInfo& r_current_process_info)
{
    if (mBuoyancyForceType == 0 || GetGeometry()[0].IsNot(INSIDE) || GetGeometry()[0].Is(BLOCKED)){ // case of identically null buoyancy
        noalias(buoyancy) = ZeroVector(3);
        return;
    }

    else {
        const double volume = TBaseElement::CalculateVolume();

        if (mDragForceType == 2){ // Weatherford
            noalias(buoyancy) =  - gravity *  mFluidDensity * volume;
        }

        else {
            const array_1d<double, 3>& pressure_grad = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_GRAD_PROJECTED);
            noalias(buoyancy) = - volume * pressure_grad;
        }
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

template< class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeDragForce(array_1d<double, 3>& drag_force, ProcessInfo& r_current_process_info)
{
    if (mDragForceType == 0 || GetGeometry()[0].IsNot(INSIDE) || GetGeometry()[0].Is(BLOCKED)){ // case of identically null drag force
        noalias(drag_force) = ZeroVector(3);
        return;
    }

    else { // calculating the 'dimensional' drag coefficient, i.e., the factor by which the slip velocity must be multiplied to yield the drag force
        ProcessInfo const& const_current_process_info = r_current_process_info;
        double drag_coeff = ComputeDragCoefficient(const_current_process_info);

        ApplyDragPorosityModification(drag_coeff);

        noalias(drag_force) = drag_coeff * mSlipVel;
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeVirtualMassForce(array_1d<double, 3>& virtual_mass_force, ProcessInfo& r_current_process_info)
{
    if (mVirtualMassForceType == 0 || GetGeometry()[0].IsNot(INSIDE) || GetGeometry()[0].Is(BLOCKED)){ // case of identically null virtual mass force
        noalias(virtual_mass_force) = ZeroVector(3);
        return;
    }

    else {
        const double volume                     = TBaseElement::CalculateVolume();
        const double delta_t_inv                = 1 / r_current_process_info[DELTA_TIME];
        const array_1d<double, 3>& fluid_acc    = GetGeometry()[0].FastGetSolutionStepValue(FLUID_ACCEL_PROJECTED);
        const array_1d<double, 3>& particle_acc = delta_t_inv * (GetGeometry()[0].FastGetSolutionStepValue(VELOCITY) - GetGeometry()[0].FastGetSolutionStepValue(VELOCITY, 1));
        array_1d<double, 3> slip_acc;

    if (mFluidModelType == 0){ // fluid velocity is modified as a post-process
        noalias(slip_acc) = fluid_acc / mFluidFraction - particle_acc;
    }

    else {
        noalias(slip_acc) = fluid_acc - particle_acc;
    }

    double virtual_mass_coeff = 0.5; // inviscid case

    if (mVirtualMassForceType == 2 || mVirtualMassForceType == 4) { // Zuber (1964) (moderate values of solid fraction)
        virtual_mass_coeff += 1.5 * (1 - mFluidFraction);
    }

    if (mVirtualMassForceType == 3 || mVirtualMassForceType == 4){ // Odar and Hamilton, 1964
        double acc_number;
        ComputeParticleAccelerationNumber(slip_acc, acc_number);
        virtual_mass_coeff *= 2.1 - 0.132 / (SWIMMING_POW_2(acc_number) + 0.12);
    }

    noalias(virtual_mass_force) = virtual_mass_coeff * mFluidDensity * volume * slip_acc;
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeSaffmanLiftForce(array_1d<double, 3>& lift_force, ProcessInfo& r_current_process_info)
{
    if (mSaffmanForceType == 0 || GetGeometry()[0].IsNot(INSIDE) || GetGeometry()[0].Is(BLOCKED)){ // case of identically null lift force
        noalias(lift_force) = ZeroVector(3);

        return;
    }

    else if (mSaffmanForceType >= 1){
        const double& shear_rate                       = GetGeometry()[0].FastGetSolutionStepValue(SHEAR_RATE_PROJECTED);
        const array_1d<double, 3>& vorticity           = GetGeometry()[0].FastGetSolutionStepValue(FLUID_VORTICITY_PROJECTED);
        array_1d<double, 3> vort_cross_slip_vel;
        SWIMMING_SET_TO_CROSS_OF_FIRST_TWO_3(mSlipVel, vorticity, vort_cross_slip_vel)
        const double vorticity_norm                    = SWIMMING_MODULUS_3(vorticity);

        double lift_coeff;

        if (mSaffmanForceType == 1){ // El Samni, E.A. (1949), see paper by R. K. Clark (1994)
            lift_coeff = ComputeElSamniLiftCoefficient(shear_rate, vorticity_norm, r_current_process_info);
        }

        else if (mSaffmanForceType == 2){ // Mei, 1992 (Re ~ 0.1 - 100)
            double reynolds;
            double reynolds_shear;
            ComputeParticleReynoldsNumber(reynolds);
            const double norm_of_vort = SWIMMING_MODULUS_3(vorticity);
            ComputeParticleRotationReynoldsNumber(norm_of_vort, reynolds_shear);
            lift_coeff = ComputeMeiLiftCoefficient(reynolds, reynolds_shear);
        }

        else {
            std::cout << "The integer value designating the Saffman lift coefficient calculation model" << std::endl;
            std::cout << " (mSaffmanForceType = " << mSaffmanForceType << "), is not supported" << std::endl << std::flush;
            return;
        }

        noalias(lift_force) = lift_coeff * vort_cross_slip_vel; // the direction is given by the vorticity x (- slip_vel) (Jackson, 2000), which is normalized here
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeMagnusLiftForce(array_1d<double, 3>& lift_force, ProcessInfo& r_current_process_info)
{
    if (mMagnusForceType == 0 || GetGeometry()[0].IsNot(INSIDE) || GetGeometry()[0].Is(BLOCKED)){
        noalias(lift_force) = ZeroVector(3);

        return;
    }

    const array_1d<double, 3> slip_rot = 0.5 * GetGeometry()[0].FastGetSolutionStepValue(FLUID_VORTICITY_PROJECTED) - GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
    array_1d<double, 3> slip_rot_cross_slip_vel;
    SWIMMING_SET_TO_CROSS_OF_FIRST_TWO_3(slip_rot, mSlipVel, slip_rot_cross_slip_vel)

    if (mMagnusForceType == 1){ // Rubinow and Keller, 1961 (very small Re)

        lift_force = KRATOS_M_PI * SWIMMING_POW_3(mRadius) *  mFluidDensity * slip_rot_cross_slip_vel;
    }

    else if (mMagnusForceType == 2){ // Oesterle and Bui Dihn, 1998 Re < 140
        const double norm_of_slip_rot = SWIMMING_MODULUS_3(slip_rot);
        double reynolds;
        double rot_reynolds;
        ComputeParticleReynoldsNumber(reynolds);
        ComputeParticleRotationReynoldsNumber(norm_of_slip_rot, rot_reynolds);

        if (reynolds == 0.0 || rot_reynolds == 0.0){
            noalias(lift_force) = ZeroVector(3);
        }

        else {
            const double lift_coeff = 0.45  + (rot_reynolds / reynolds - 0.45) * exp(- 0.05684 * pow(rot_reynolds, 0.4) * pow(reynolds, 0.3));
            noalias(lift_force) = 0.5 *  mFluidDensity * KRATOS_M_PI * SWIMMING_POW_2(mRadius) * lift_coeff * mNormOfSlipVel * slip_rot_cross_slip_vel / norm_of_slip_rot;
        }
    }

    else {
        std::cout << "The integer value designating the magnus lift coefficient calculation model" << std::endl;
        std::cout << " (mMagnusForceType = " << mMagnusForceType << "), is not supported" << std::endl << std::flush;
        return;
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeHydrodynamicTorque(array_1d<double, 3>& hydro_torque, ProcessInfo& r_current_process_info)
{
    if (mHydrodynamicTorqueType == 0 || GetGeometry()[0].IsNot(INSIDE) || GetGeometry()[0].Is(BLOCKED)){
        noalias(hydro_torque) = ZeroVector(3);

        return;
    }

    else if (mHydrodynamicTorqueType == 1){
        const array_1d<double, 3> slip_rot = 0.5 * GetGeometry()[0].FastGetSolutionStepValue(FLUID_VORTICITY_PROJECTED) - GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
        const double norm_of_slip_rot = SWIMMING_MODULUS_3(slip_rot);
        double rot_reynolds;
        ComputeParticleRotationReynoldsNumber(norm_of_slip_rot, rot_reynolds);
        double rotational_coeff;

        if (rot_reynolds > 32){ // Rubinow and Keller, 1961 (Re_rot ~ 32 - 1000)
            rotational_coeff = 12.9 / sqrt(rot_reynolds) + 128.4 / rot_reynolds;
        }

        else { // Rubinow and Keller, 1961 (Re_rot < 32)
            rotational_coeff = 64 * KRATOS_M_PI / rot_reynolds;
        }

        noalias(hydro_torque) = 0.5 *  mFluidDensity * SWIMMING_POW_5(mRadius) * rotational_coeff * norm_of_slip_rot * slip_rot;
    }

    else {
        std::cout << "The integer value designating the Hydrodynamic torque calculation model" << std::endl;
        std::cout << " (mHydrodynamicTorqueType = " << mHydrodynamicTorqueType << "), is not supported" << std::endl << std::flush;
        return;
    }

}


//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

template< class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeBrownianMotionForce(array_1d<double, 3>& brownian_motion_force, ProcessInfo& r_current_process_info)
{
    if (mBrownianMotionType == 0 || GetGeometry()[0].IsNot(INSIDE) || GetGeometry()[0].Is(BLOCKED)){
        noalias(brownian_motion_force) = ZeroVector(3);

        return;
    }

    else {
        noalias(brownian_motion_force) = ZeroVector(3);
        /*const double kT = 4.11e-21;
        double current_time = r_current_process_info[TIME] ;
        double delta_t_inv = 1.0 / (current_time - mLastTimeStep);
        mLastTimeStep = current_time;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(-0.5, 0.5);
        double coeff = std::sqrt(24 * kT * 2 / KRATOS_M_PI * ComputeStokesDragCoefficient() * delta_t_inv);
        brownian_motion_force[0] = coeff * dis(gen);
        brownian_motion_force[1] = coeff * dis(gen);
        brownian_motion_force[2] = coeff * dis(gen);*/
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeParticleReynoldsNumber(double& reynolds)
{
    reynolds = 2 * mRadius * mNormOfSlipVel / mKinematicViscosity;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeParticleRotationReynoldsNumber(double norm_of_slip_rot, double& reynolds)
{
    reynolds = 4 * SWIMMING_POW_2(mRadius) * norm_of_slip_rot / mKinematicViscosity;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeParticleAccelerationNumber(const array_1d<double, 3>& slip_acc, double& acc_number)
{
    acc_number = SWIMMING_POW_3(mNormOfSlipVel) / fabs(2 * mRadius * SWIMMING_INNER_PRODUCT_3(mSlipVel, slip_acc));
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::AdditionalCalculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_current_process_info)
{
    if (rVariable == REYNOLDS_NUMBER){

        if (GetGeometry()[0].IsNot(INSIDE)){
            Output = 0.0;
        }

        else {
            mFluidDensity                           = GetGeometry()[0].FastGetSolutionStepValue(FLUID_DENSITY_PROJECTED);
            mKinematicViscosity                     = GetGeometry()[0].FastGetSolutionStepValue(FLUID_VISCOSITY_PROJECTED);
            mFluidFraction                          = GetGeometry()[0].FastGetSolutionStepValue(FLUID_FRACTION_PROJECTED);
            const array_1d<double, 3>& fluid_vel    = GetGeometry()[0].FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
            const array_1d<double, 3>& particle_vel = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);

            if (mFluidModelType == 0){ // fluid velocity is modified as a post-process
                noalias(mSlipVel) = fluid_vel / mFluidFraction - particle_vel;
            }

            else {
                noalias(mSlipVel) = fluid_vel - particle_vel;
            }

            mNormOfSlipVel = SWIMMING_MODULUS_3(mSlipVel);
            ComputeParticleReynoldsNumber(Output);
        }
    }

    else if (rVariable == DRAG_COEFFICIENT){
        Output = ComputeDragCoefficient(r_current_process_info);
    }
}


//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::ComputeDragCoefficient(const ProcessInfo& r_current_process_info)
{
    double drag_coeff;

    if (mDragForceType == 1){
        drag_coeff = ComputeStokesDragCoefficient();
    }

    else if (mDragForceType == 2){ // formulations of Haider (1989) and Chien (1994)
        drag_coeff = ComputeWeatherfordDragCoefficient(r_current_process_info);
    }

    else if (mDragForceType == 3){ // formulation of Ganser (1993)
        drag_coeff = ComputeGanserDragCoefficient(r_current_process_info);
    }

    else if (mDragForceType == 4){ // formulation of Ishii and Zuber (1979)
        drag_coeff = ComputeIshiiDragCoefficient(r_current_process_info);
    }

    else if (mDragForceType == 5){ // Newton regime (Re ~ 1000 - 250000), formulation of Haider and Levenspiel (1989)
        drag_coeff = ComputeNewtonRegimeDragCoefficient();
    }

    else if (mDragForceType == 6){ // Intermediate regime (Re ~ 0.5 - 1000), formulation of Schiller and Naumann (1933)
        drag_coeff = ComputeIntermediateRegimeDragCoefficient();
    }

    else if (mDragForceType == 7){ // All regimes (Re ~ 0 - 250000), formulation of Haider and Levenspiel (1989)
        drag_coeff = ComputeHaiderDragCoefficient();
    }

    else if (mDragForceType == 8){ // Intermediate regimes (Re ~ 0 - 1000), formulation of Beetstra et al. obtained using lattice-Boltzmann (2007)
        drag_coeff = ComputeBeetstraDragCoefficient();
    }

    else if (mDragForceType == 9){ // Coin-shaped Stokesian
        drag_coeff = 2.0 / KRATOS_M_PI * ComputeStokesDragCoefficient();
    }

    else {
        std::string message;
        std::string first_part ("The integer value designating the drag coefficient calculation model.\n");
        std::string second_part ("(mDragForceType) is not supported.");
        message = first_part + second_part;
        KRATOS_THROW_ERROR(std::invalid_argument, message, mDragForceType);
    }

    return drag_coeff;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::ComputeStokesDragCoefficient()
{
    double drag_coeff = 6.0 * KRATOS_M_PI * mKinematicViscosity * mFluidDensity * mRadius;

    return drag_coeff;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::ComputeWeatherfordDragCoefficient(const ProcessInfo& r_current_process_info)
{
    KRATOS_TRY

    const double particle_density              = GetDensity();
    const array_1d<double, 3>& gravity         = r_current_process_info[GRAVITY];
    const int manually_imposed_drag_law_option = r_current_process_info[MANUALLY_IMPOSED_DRAG_LAW_OPTION];
    const int drag_modifier_type               = r_current_process_info[DRAG_MODIFIER_TYPE];
    const double gel_strength                  = GetGeometry()[0].FastGetSolutionStepValue(YIELD_STRESS);
    const double power_law_n                   = GetGeometry()[0].FastGetSolutionStepValue(POWER_LAW_N);
    const double power_law_K                   = GetGeometry()[0].FastGetSolutionStepValue(POWER_LAW_K);
    const double yield_stress                  = GetGeometry()[0].FastGetSolutionStepValue(YIELD_STRESS);

    int non_newtonian_option = 1;

    if (fabs(power_law_n - 1.0) < 0.00001  ||  fabs(yield_stress) < 0.00001) {
        non_newtonian_option = 0;
    }

    const double initial_drag_force            = r_current_process_info[INIT_DRAG_FORCE];
    const double drag_law_slope                = r_current_process_info[DRAG_LAW_SLOPE];
    const double power_law_tol                 = r_current_process_info[POWER_LAW_TOLERANCE];

    const double area                          = KRATOS_M_PI * SWIMMING_POW_2(mRadius);
    const array_1d<double, 3> weight           = mRealMass * gravity;
    const array_1d<double, 3> buoyancy         = mFluidDensity / particle_density * weight; // hydrostatic case!! (only for Weatherford)

    double shahs_term_vel                      = 0.0;
    double beta                                = 0.0;
    double F0                                  = 0.0;
    double regularization_v                    = 0.02 * mRadius;
    double reynolds;
    double drag_coeff;

    if (!non_newtonian_option){ // Newtonian
        ComputeParticleReynoldsNumber(reynolds);

        if (!non_newtonian_option && reynolds < 0.01){
            reynolds = 0.01;
        }

        CalculateNewtonianDragCoefficient(non_newtonian_option, reynolds, mSphericity, drag_coeff, drag_modifier_type);
        drag_coeff = 0.5 *  mFluidDensity * area * drag_coeff * mNormOfSlipVel;
    }

    else {
        shahs_term_vel = CalculateShahsTerm(power_law_n, power_law_K, power_law_tol, particle_density, mSphericity, drag_modifier_type);

        if (!manually_imposed_drag_law_option){
            F0 = 4.0 * gel_strength * area; //initial value
            beta = (SWIMMING_MODULUS_3(weight) - SWIMMING_MODULUS_3(buoyancy) - F0) / shahs_term_vel; //S
        }

        else {
            F0 = initial_drag_force; //initial value
            beta = drag_law_slope; //slope
        }

        if (mNormOfSlipVel >= regularization_v){
            drag_coeff = (F0 + beta * mNormOfSlipVel) / mNormOfSlipVel;
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
template< class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::CalculateNewtonianDragCoefficient(int non_newtonian_option,
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
template< class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::CalculateDragCoeffFromSphericity(const double reynolds,
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
template< class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::CalculateShahsTerm(double power_law_N,
                                                   double power_law_K,
                                                   double power_law_tol,
                                                   const double& particle_density,
                                                   double sphericity,
                                                   int drag_modifier_type)
{
    if (fabs(power_law_N) < power_law_tol || fabs(power_law_K) < power_law_tol){
        std::cout << "WARNING: Shah's method is being used with Power Law data being zero!!" << std::endl << std::flush;
    }

    double shah_A_i = 1 / (6.9148 * power_law_N * power_law_N - 24.838 * power_law_N + 22.642);
    double shah_B_i = 1 / (-0.5067 * power_law_N * power_law_N + 1.3234 * power_law_N - 0.1744);

    double dimensionless_shah = sqrt(pow(13.08, 2 - power_law_N) * pow(2 * mRadius, power_law_N + 2) * pow( mFluidDensity, power_law_N) * pow(particle_density -  mFluidDensity, 2 - power_law_N) / (pow(2, 2 * (power_law_N - 1)) * power_law_K * power_law_K));
    double reynolds = pow(dimensionless_shah * shah_A_i, shah_B_i);
    double fi_i = CalculateDragCoeffFromSphericity(reynolds, 1.0, drag_modifier_type) / CalculateDragCoeffFromSphericity(reynolds, sphericity, drag_modifier_type);
    dimensionless_shah = sqrt(pow(fi_i, 2 - power_law_N)) * dimensionless_shah;
    reynolds = pow(dimensionless_shah * shah_A_i, shah_B_i);

    double terminal_vel =  pow(pow(2, power_law_N - 1) * power_law_K * reynolds / (pow(2 * mRadius, power_law_N) *  mFluidDensity), 1 / (2 - power_law_N)) ;

    return terminal_vel;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::ComputeGanserDragCoefficient(const ProcessInfo& r_current_process_info)
{
    KRATOS_TRY

    const int isometric_shape                = 1; // TEMPORARY!! yes (1) or no (0); shold be given as data
    const double surface_area                = 4 * KRATOS_M_PI * SWIMMING_POW_2(mRadius); // TEMPORARY!! corresponding to a sphere; should be generalized b taking it as a parameter
    const double surface_area_circular_diam  = sqrt(4.0 * surface_area / KRATOS_M_PI);

    double equiv_reynolds;
    double k_1;
    double k_2;
    double drag_coeff;

    ComputeGanserParameters(isometric_shape, surface_area_circular_diam, k_1, k_2);
    ComputeParticleReynoldsNumber(equiv_reynolds);
    equiv_reynolds *= k_1 * k_2;

    // calculating adimensional drag coefficient

    drag_coeff =  k_2 * (24 * (1 + 0.1118 * pow((equiv_reynolds), 0.6567)) / (equiv_reynolds) + 0.4305 / (1 + 3305 / equiv_reynolds));

    // and then the dimensional drag coefficient

    drag_coeff *= 0.5 *  mFluidDensity * surface_area * mNormOfSlipVel;

    return drag_coeff;

    KRATOS_CATCH("")
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::ComputeIshiiDragCoefficient(const ProcessInfo& r_current_process_info)
{
    double coeff = 0.45;
    double reynolds;
    ComputeParticleReynoldsNumber(reynolds);

    if (reynolds <= 1000){
        coeff = (24 + 2.4 * pow(reynolds, 0.75)) / reynolds;
    }

    double drag_coeff = 0.5 * coeff * KRATOS_M_PI * SWIMMING_POW_2(mRadius);

    return drag_coeff;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::ComputeNewtonRegimeDragCoefficient()
{
    double drag_coeff  = 0.5 * KRATOS_M_PI * SWIMMING_POW_2(mRadius) *  mFluidDensity * mNormOfSlipVel;

    drag_coeff *= 0.44;

    return drag_coeff;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::ComputeIntermediateRegimeDragCoefficient()
{
    double reynolds;
    double drag_coeff  = 0.5 * KRATOS_M_PI * SWIMMING_POW_2(mRadius) *  mFluidDensity * mNormOfSlipVel;

    ComputeParticleReynoldsNumber(reynolds);

    drag_coeff *= 24 / reynolds * (1 + 0.15 * pow(reynolds, 0.687));

    return drag_coeff;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::ComputeHaiderDragCoefficient()
{
    const double sphericity = GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_SPHERICITY);
    double drag_coeff       = 0.5 * KRATOS_M_PI * SWIMMING_POW_2(mRadius) *  mFluidDensity * mNormOfSlipVel;

    double A = exp(2.3288 - 6.4581 * sphericity + 2.4486 * sphericity * sphericity);
    double B = 0.0964 + 0.5565 * sphericity;
    double C = exp(4.905  - 13.8944 * sphericity + 18.4222 * sphericity * sphericity - 10.2599 * sphericity * sphericity * sphericity);
    double D = exp(1.4681 + 12.2584 * sphericity - 20.7322 * sphericity * sphericity + 15.8855 * sphericity * sphericity * sphericity);
    double particle_reynolds;
    ComputeParticleReynoldsNumber(particle_reynolds);

    drag_coeff *= (24.0 * (1.0 + A * pow(particle_reynolds, B))) / particle_reynolds + C * particle_reynolds / (particle_reynolds + D);

    return drag_coeff;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::ComputeBeetstraDragCoefficient()
{
    double drag_coeff;
    double particle_reynolds;
    ComputeParticleReynoldsNumber(particle_reynolds);

    if (particle_reynolds < 1.0){
        drag_coeff = ComputeStokesDragCoefficient();
    }

    else {
        double eps = GetGeometry()[0].FastGetSolutionStepValue(FLUID_FRACTION_PROJECTED);

        if (eps > 0.999){
            eps = 0.9;
        }

        const double eps_s = 1.0 - eps;
        particle_reynolds *= eps;

        double A = 180 + 18 * std::pow(eps, 4) / eps_s * (1 + 1.5 * std::sqrt(eps_s));
        double B = 0.31 * (1.0 / eps + 3 * eps_s * eps + 8.4 * std::pow(particle_reynolds, - 0.343)) / (1.0 + std::pow(10.0, 3 * eps_s) * std::pow(particle_reynolds, 2 * eps - 2.5));
        drag_coeff = KRATOS_M_PI_3 * mKinematicViscosity * mFluidDensity * mRadius * (A * eps_s / eps + B * particle_reynolds);
    }

    return drag_coeff;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
 void SphericSwimmingParticle<TBaseElement>::ComputeGanserParameters(const int isometric_shape, const double dn, double& k_1, double& k_2)
 {
     if (isometric_shape){
         k_1 = 3 / (1 + 2 / sqrt(mSphericity));
     }

     else {
         k_1 = 3 / (0.5 * dn / mRadius + 2 / sqrt(mSphericity));
     }

     k_2 = pow(10.0, 1.8148 * pow(- log10(mSphericity), 0.5743));
 }

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ApplyDragPorosityModification(double& drag_coeff)
{
    if (mPorosityCorrectionType == 0){
        return;
    }

    else if (mPorosityCorrectionType == 1){ // Richardson and Zaki, 1954 (fluid fraction ~ 0.01 - 0.2)
        double reynolds;
        ComputeParticleReynoldsNumber(reynolds);
        double K;

        if (reynolds > 500){
            K = 2.39;
        }

        else if (reynolds > 1){
            K = 4.45 * pow(reynolds, - 0.1);
        }

        else if (reynolds > 0.2){
            K = 4.35 * pow(reynolds, - 0.03);
        }

        else {
            K = 4.65;
        }

        drag_coeff *= pow(mFluidFraction, 1 - 2 * K);
    }

    else if (mPorosityCorrectionType == 2){

    }

}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::ComputeElSamniLiftCoefficient(const double norm_of_shear_rate,
                                                              const double vorticity_norm,
                                                              ProcessInfo& r_current_process_info)
{
    if (vorticity_norm > 0.000000000001 && mNormOfSlipVel > 0.000000000001){
         const double yield_stress   = 0.0; // we are considering a Bingham type fluid
         const double power_law_K    = GetGeometry()[0].FastGetSolutionStepValue(POWER_LAW_K);
         const double power_law_n    = GetGeometry()[0].FastGetSolutionStepValue(POWER_LAW_N);
         const double shear_rate_p   = mNormOfSlipVel / mRadius * (4.5 / power_law_n - 3.5); // graphic model by Unhlherr et al. (fit by Wallis, G.B. and Dobson, J.E., 1973)
         double equivalent_viscosity = yield_stress / shear_rate_p + power_law_K * pow(shear_rate_p, power_law_n - 1);
         const double coeff          = std::max(0.09 * mNormOfSlipVel, 5.82 * sqrt(0.5 * mNormOfSlipVel * equivalent_viscosity /  mFluidDensity));
         const double lift_coeff     = 0.5 * KRATOS_M_PI * SWIMMING_POW_2(mRadius) *  mFluidDensity * coeff * mNormOfSlipVel / vorticity_norm;
         return(lift_coeff);
    }

    else {
        return 0.0;
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
 double SphericSwimmingParticle<TBaseElement>::ComputeMeiLiftCoefficient(const double reynolds, const double reynolds_shear)
 {
     if (reynolds != 0.0 && reynolds_shear != 0.0 ){
         double sqrt_beta = sqrt(0.5 * reynolds_shear / reynolds);
         double C;

         if (reynolds < 40){
             C = (1 - 0.3314 * sqrt_beta) * exp(- 0.1 * reynolds) + 0.3314 * sqrt_beta;
         }

         else {
             C = 0.0524 * sqrt_beta * sqrt(reynolds);
         }

         C *= 4.1126 / sqrt(reynolds_shear);

         double lift_coeff = mFluidDensity * KRATOS_M_PI * SWIMMING_POW_3(mRadius) * C;

         return lift_coeff;
     }

     else {
         return 0.0;
     }
 }

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::CustomInitialize()
{
    TBaseElement::CustomInitialize();
    mHasDragForceNodalVar        = GetGeometry()[0].SolutionStepsDataHas(DRAG_FORCE);
    mHasHydroMomentNodalVar      = GetGeometry()[0].SolutionStepsDataHas(HYDRODYNAMIC_MOMENT);
    mHasVirtualMassForceNodalVar = GetGeometry()[0].SolutionStepsDataHas(VIRTUAL_MASS_FORCE);
    mHasLiftForceNodalVar        = GetGeometry()[0].SolutionStepsDataHas(LIFT_FORCE);
    mSphericity                  = GetGeometry()[0].SolutionStepsDataHas(PARTICLE_SPHERICITY);
    mHasDragCoefficientVar       = GetGeometry()[0].SolutionStepsDataHas(DRAG_COEFFICIENT);
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::AdditionalMemberDeclarationFirstStep(const ProcessInfo& r_process_info)
{
    mCouplingType           = r_process_info[COUPLING_TYPE];
    mBuoyancyForceType      = r_process_info[BUOYANCY_FORCE_TYPE];
    mDragForceType          = r_process_info[DRAG_FORCE_TYPE];
    mVirtualMassForceType   = r_process_info[VIRTUAL_MASS_FORCE_TYPE];
    mSaffmanForceType       = r_process_info[LIFT_FORCE_TYPE];
    mMagnusForceType        = r_process_info[MAGNUS_FORCE_TYPE];
    mFluidModelType         = r_process_info[FLUID_MODEL_TYPE];
    mPorosityCorrectionType = r_process_info[DRAG_POROSITY_CORRECTION_TYPE];
    mHydrodynamicTorqueType = r_process_info[HYDRO_TORQUE_TYPE];
    mBrownianMotionType     = 0;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

template class SphericSwimmingParticle<SphericParticle>; //Explicit Instantiation
template class SphericSwimmingParticle<NanoParticle>; //Explicit Instantiation
}  // namespace Kratos.

