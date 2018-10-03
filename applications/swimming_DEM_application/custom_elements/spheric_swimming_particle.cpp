//   Author: G.Casas, gcasas@cimne.upc.edu $

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "swimming_DEM_application.h"
#include "spheric_swimming_particle.h"
#include "../applications/DEM_application/custom_utilities/GeometryFunctions.h"

namespace Kratos
{

template < class TBaseElement >
SphericSwimmingParticle<TBaseElement>& SphericSwimmingParticle<TBaseElement>::operator=(const SphericSwimmingParticle<TBaseElement>& rOther) {

    TBaseElement::operator=(rOther);

    mNeighbourNodes = rOther.mNeighbourNodes;
    mNeighbourNodesDistances = rOther.mNeighbourNodesDistances;
    mHasHydroMomentNodalVar = rOther.mHasHydroMomentNodalVar;
    mHasDragForceNodalVar = rOther.mHasDragForceNodalVar;
    mHasVirtualMassForceNodalVar = rOther.mHasVirtualMassForceNodalVar;
    mHasBassetForceNodalVar = rOther.mHasBassetForceNodalVar;
    mHasLiftForceNodalVar = rOther.mHasLiftForceNodalVar;
    mHasDragCoefficientVar = rOther.mHasDragCoefficientVar;
    mHasOldAdditionalForceVar = rOther.mHasOldAdditionalForceVar;
    mFirstStep = rOther.mFirstStep;
    mCouplingType = rOther.mCouplingType;
    mBuoyancyForceType = rOther.mBuoyancyForceType;
    mDragForceType = rOther.mDragForceType;
    mVirtualMassForceType = rOther.mVirtualMassForceType;
    mBassetForceType = rOther.mBassetForceType;
    mSaffmanForceType = rOther.mSaffmanForceType;
    mMagnusForceType = rOther.mMagnusForceType;
    mFluidModelType = rOther.mFluidModelType;
    mPorosityCorrectionType = rOther.mPorosityCorrectionType;
    mHydrodynamicTorqueType = rOther.mHydrodynamicTorqueType;
    mBrownianMotionType = rOther.mBrownianMotionType;
    mQuadratureOrder = rOther.mQuadratureOrder;
    mFluidDensity = rOther.mFluidDensity;
    mFluidFraction = rOther.mFluidFraction;
    mKinematicViscosity = rOther.mKinematicViscosity;
    mSphericity = rOther.mSphericity;
    mNormOfSlipVel = rOther.mNormOfSlipVel;
    mLastTimeStep = rOther.mLastTimeStep;
    mInitialTime = rOther.mInitialTime;
    mOldDaitchePresentCoefficient = rOther.mOldDaitchePresentCoefficient;
    mLastVirtualMassAddedMass = rOther.mLastVirtualMassAddedMass;
    mLastBassetForceAddedMass = rOther.mLastBassetForceAddedMass;
    noalias(mSlipVel) = rOther.mSlipVel;
    noalias(mOldBassetTerm) = rOther.mOldBassetTerm;

    return *this;
}


template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeAdditionalForces(array_1d<double, 3>& non_contact_force,
                                                                    array_1d<double, 3>& non_contact_moment,
                                                                    const ProcessInfo& r_current_process_info,
                                                                    const array_1d<double,3>& gravity)
{
    KRATOS_TRY

    if (!mCouplingType){
        TBaseElement::ComputeAdditionalForces(non_contact_force, non_contact_moment, r_current_process_info, gravity);
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

        if (mDragForceType == 11) {
            const array_1d<double, 3>& fluid_vel_laplacian = node.FastGetSolutionStepValue(FLUID_VEL_LAPL_PROJECTED);
            noalias(mSlipVel) -= mRadius * mRadius / 6.0 * fluid_vel_laplacian; // adding Faxen term
        }
    }

    mNormOfSlipVel = SWIMMING_MODULUS_3(mSlipVel);
    array_1d<double, 3> weight = ZeroVector(3);
    array_1d<double, 3> buoyancy;
    array_1d<double, 3> drag_force;
    array_1d<double, 3> virtual_mass_plus_undisturbed_flow_force;
    array_1d<double, 3> saffman_lift_force;
    array_1d<double, 3> magnus_lift_force;
    array_1d<double, 3> brownian_motion_force;
    array_1d<double, 3>& basset_force = node.FastGetSolutionStepValue(BASSET_FORCE);

    // The decomposition of forces that is considered here follows Jackson (The Dynamics of Fluidized Particles, 2000);
    // so that the role of f_n1 therein is played by non_contact_force here
    TBaseElement::ComputeAdditionalForces(weight, non_contact_moment, r_current_process_info, gravity); // Could be adding something else apart from weight
    ComputeBuoyancy(node, buoyancy, gravity, r_current_process_info);
    ComputeDragForce(node, drag_force, r_current_process_info);
    ComputeVirtualMassPlusUndisturbedFlowForce(node, virtual_mass_plus_undisturbed_flow_force, r_current_process_info);
    ComputeSaffmanLiftForce(node, saffman_lift_force, r_current_process_info);
    ComputeMagnusLiftForce(node, magnus_lift_force, r_current_process_info);
    ComputeHydrodynamicTorque(node, non_contact_moment, r_current_process_info);
    ComputeBrownianMotionForce(node, brownian_motion_force, r_current_process_info);
    ComputeBassetForce(node, basset_force, r_current_process_info);

    // Adding all forces except Basset's, since they might get averaged in time in a different way
    noalias(non_contact_force) += drag_force
                                + virtual_mass_plus_undisturbed_flow_force
                                + saffman_lift_force
                                + magnus_lift_force
                                + brownian_motion_force
                                + buoyancy
                                + weight;

    const double force_reduction_coeff = mRealMass / (mRealMass + mLastVirtualMassAddedMass + mLastBassetForceAddedMass);

    array_1d<double, 3> non_contact_force_not_altered = non_contact_force;

    if (mHasOldAdditionalForceVar && !mFirstStep){
        ApplyNumericalAveragingWithOldForces(node, non_contact_force, r_current_process_info);
    }

    UpdateNodalValues(node, non_contact_force_not_altered, non_contact_moment, weight, buoyancy, drag_force, virtual_mass_plus_undisturbed_flow_force, basset_force, saffman_lift_force, magnus_lift_force, force_reduction_coeff, r_current_process_info);
    // The Basset force has a different temporal treatment, so first we apply the scheme to the rest of the forces
    // and then we add the Basset force (minus the term proportional to the current acceleration, which is treated implicitly)
    noalias(non_contact_force) += basset_force;
    non_contact_force *= force_reduction_coeff; //TODO: put noalias here?
    mFirstStep = false;

    KRATOS_CATCH( "" )
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
// Here nodal values are modified to record DEM forces that we want to print. In Kratos this is an exception since nodal values are meant to be modified only outside the element. Here it was not possible.
template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::UpdateNodalValues(NodeType& node,
                                                              const array_1d<double, 3>& non_contact_nor_basset_force,
                                                              const array_1d<double, 3>& non_contact_moment,
                                                              const array_1d<double, 3>& weight,
                                                              const array_1d<double, 3>& buoyancy,
                                                              const array_1d<double, 3>& drag_force,
                                                              const array_1d<double, 3>& virtual_mass_plus_undisturbed_flow_force,
                                                              const array_1d<double, 3>& basset_force,
                                                              const array_1d<double, 3>& saffman_lift_force,
                                                              const array_1d<double, 3>& magnus_lift_force,
                                                              const double& force_reduction_coeff,
                                                              const ProcessInfo& r_current_process_info)
{
    noalias(node.FastGetSolutionStepValue(HYDRODYNAMIC_FORCE))       = force_reduction_coeff * (non_contact_nor_basset_force + basset_force - buoyancy - weight);
    noalias(node.FastGetSolutionStepValue(BUOYANCY))                 = buoyancy;
    array_1d<double, 3>& total_force = node.FastGetSolutionStepValue(TOTAL_FORCES);
    total_force *= force_reduction_coeff;

    if (mHasHydroMomentNodalVar){
        noalias(node.FastGetSolutionStepValue(HYDRODYNAMIC_MOMENT))  = non_contact_moment;
    }

    if (mHasDragForceNodalVar){
        noalias(node.FastGetSolutionStepValue(DRAG_FORCE))           = drag_force;
    }

    if (mHasVirtualMassForceNodalVar){ // This only includes the part proportional to the fluid acceleration (undisturbed flow plus added mass terms), since the particle acceleration is treated implicitly. It is added later by the strategy, which calls Calculate here, once the current acceleration is available
        noalias(node.FastGetSolutionStepValue(VIRTUAL_MASS_FORCE))   = virtual_mass_plus_undisturbed_flow_force;
    }

    if (mHasBassetForceNodalVar){
        noalias(node.FastGetSolutionStepValue(BASSET_FORCE))         = basset_force; // This does not include the current-time contribution of the acceleration, which is treated implicitly. It is added later by the strategy, which calls Calculate here, once the current acceleration is available
    }

    if (mHasOldAdditionalForceVar){
        noalias(node.FastGetSolutionStepValue(ADDITIONAL_FORCE_OLD)) = non_contact_nor_basset_force;
    }

    if (mHasLiftForceNodalVar){
        noalias(node.FastGetSolutionStepValue(LIFT_FORCE))           = saffman_lift_force + magnus_lift_force;
    }

    if (mHasDragCoefficientVar){
        double drag_coefficient = ComputeDragCoefficient(r_current_process_info);
        node.FastGetSolutionStepValue(DRAG_COEFFICIENT)              = drag_coefficient;
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ApplyNumericalAveragingWithOldForces(NodeType& node, array_1d<double, 3>& non_contact_nor_basset_force, const ProcessInfo& r_current_process_info)
{
    noalias(non_contact_nor_basset_force) = 0.5 * (3 * non_contact_nor_basset_force - node.FastGetSolutionStepValue(ADDITIONAL_FORCE_OLD));
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeBuoyancy(NodeType& node, array_1d<double, 3>& buoyancy, const array_1d<double, 3>& gravity, const ProcessInfo& r_current_process_info)
{
    if (mBuoyancyForceType == 0 || node.IsNot(INSIDE) || node.Is(BLOCKED)){ // case of identically null buoyancy
        noalias(buoyancy) = ZeroVector(3);
        return;
    }

    else {

        if (mBuoyancyForceType == 2 || mDragForceType == 2){ // Maxey-Riley form of buoyancy (minus the fluid acceleration term); Weatherford
            noalias(buoyancy) = - gravity * GetFluidMass();
        }

        else {
            const array_1d<double, 3>& pressure_grad = node.FastGetSolutionStepValue(PRESSURE_GRAD_PROJECTED);
            noalias(buoyancy) = - CalculateVolume() * pressure_grad;
        }
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeDragForce(NodeType& node, array_1d<double, 3>& drag_force, const ProcessInfo& r_current_process_info)
{
    if (mDragForceType == 0 || node.IsNot(INSIDE) || node.Is(BLOCKED)){ // case of identically null drag force
        noalias(drag_force) = ZeroVector(3);
        return;
    }

    else { // calculating the 'dimensional drag coefficient', i.e., the factor by which the slip velocity must be multiplied to yield the drag force
        double drag_coeff = ComputeDragCoefficient(r_current_process_info);

        ApplyDragPorosityModification(drag_coeff);

        noalias(drag_force) = drag_coeff * mSlipVel;
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeVirtualMassPlusUndisturbedFlowForce(NodeType& node, array_1d<double, 3>& virtual_mass_plus_undisturbed_flow_force, const ProcessInfo& r_current_process_info)
{
    if (mVirtualMassForceType == 0 || node.IsNot(INSIDE) || node.Is(BLOCKED)){ // case of identically null virtual mass force
        noalias(virtual_mass_plus_undisturbed_flow_force) = ZeroVector(3);
        return;
    }

    else {
        // Unperturbed flow force contribution
        const array_1d<double, 3>& fluid_acc = node.FastGetSolutionStepValue(FLUID_ACCEL_PROJECTED);

        // Virtual mass force contribution
        array_1d<double, 3> slip_acc;

        if (mFluidModelType == 0){ // fluid velocity is modified as a post-process
            const array_1d<double, 3>& particle_acc = 1 / GetMass() * GetForce();
            noalias(slip_acc) = fluid_acc / mFluidFraction - particle_acc;
        }

        else {

            if (mVirtualMassForceType == 12){
                noalias(slip_acc) = node.FastGetSolutionStepValue(FLUID_ACCEL_FOLLOWING_PARTICLE_PROJECTED); // the particle acceleration is treated implicitly through the added_mass
            }

            else {
                noalias(slip_acc) = fluid_acc; // the particle acceleration is treated implicitly through the added_mass
            }

            if (mDragForceType == 11) {
                const array_1d<double, 3>& fluid_vel_laplacian_rate = node.FastGetSolutionStepValue(FLUID_VEL_LAPL_RATE_PROJECTED);
                noalias(slip_acc) -= 0.1 * mRadius * mRadius * fluid_vel_laplacian_rate; // add Faxen term
            }
        }

        // Virtual mass force coefficient
        double virtual_mass_coeff = 0.5; // inviscid case

        if (mVirtualMassForceType == 2 || mVirtualMassForceType == 4) { // Zuber (1964) (moderate values of solid fraction)
            virtual_mass_coeff += 1.5 * (1 - mFluidFraction);
        }

        else if (mVirtualMassForceType == 3 || mVirtualMassForceType == 4){ // Odar and Hamilton, 1964
            double acc_number;
            ComputeParticleAccelerationNumber(slip_acc, acc_number);
            virtual_mass_coeff *= 2.1 - 0.132 / (SWIMMING_POW_2(acc_number) + 0.12);
        }

        const double fluid_mass = GetFluidMass();
        mLastVirtualMassAddedMass = virtual_mass_coeff * fluid_mass;
        noalias(virtual_mass_plus_undisturbed_flow_force) = fluid_mass * (virtual_mass_coeff * slip_acc + fluid_acc);
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::GetDaitcheCoefficient(int order, unsigned int n, unsigned int j, const double last_h_over_h, const int n_steps_per_quad_step)
{
    const int l = (int)(last_h_over_h * n_steps_per_quad_step + 0.5) - 1;

    if (order == 1){
        if (j < n){
            return SphericSwimmingParticle<TBaseElement>::mAjs[n_steps_per_quad_step * j + l];
        }
        else {
            return SphericSwimmingParticle<TBaseElement>::mBns[n_steps_per_quad_step * j + l];
        }
    }

    else if (order == 2){
        if (n > 3){
            if (j < n - 1){
                return SphericSwimmingParticle<TBaseElement>::mAjs[n_steps_per_quad_step * j + l];
            }
            else if (j == n - 1){
                return SphericSwimmingParticle<TBaseElement>::mBns[n_steps_per_quad_step * j + l];
            }
            else {
                return SphericSwimmingParticle<TBaseElement>::mCns[n_steps_per_quad_step * j + l];
            }
        }
        else {
            if (n == 1){ // use formula for the phis of first order
                SphericSwimmingParticle<TBaseElement>::GetDaitcheCoefficient(1, n, j, last_h_over_h, n_steps_per_quad_step);
            }

            else if (n == 2){

                const double sqrt_phi_plus_1 = std::sqrt(1 + last_h_over_h);

                if (j == 0){
                    return 4 * sqrt_phi_plus_1 * (4 * last_h_over_h - 1) / (15 * last_h_over_h);
                }
                else if (j == 1){
                    return 4 * SWIMMING_POW_5(sqrt_phi_plus_1) / (15 * last_h_over_h);
                }
                else {
                    return 2 * sqrt_phi_plus_1 * (3 - 2 * last_h_over_h) / 15;
                }
            }
            else {
                const double sqrt_phi_plus_1 = std::sqrt(1 + last_h_over_h);
                const double sqrt_phi_plus_2 = std::sqrt(2 + last_h_over_h);

                if (j == 0){
                    return 4 * sqrt_phi_plus_1 * (4 * last_h_over_h - 1) / (15 * last_h_over_h);
                }
                else if (j == 1){
                    return 4 * SWIMMING_POW_5(sqrt_phi_plus_1) / (15 * last_h_over_h) + 2 * (4 * SWIMMING_POW_3(sqrt_phi_plus_2) * (3 + 4 * last_h_over_h) - SWIMMING_POW_3(sqrt_phi_plus_1) * (9 + 4 * last_h_over_h)) / 15;
                }
                else if (j == 2){
                    return 4 / 15 * (SWIMMING_POW_3(sqrt_phi_plus_2) * (2 - 4 * last_h_over_h) + sqrt_phi_plus_1 * (4 * last_h_over_h * last_h_over_h + 7 * last_h_over_h - 2));
                    //return 2 * sqrt_phi_plus_1 * (3 - 2 * last_h_over_h) / 15 + 2 * (4 * SWIMMING_POW_3(sqrt_phi_plus_2) * (2 * last_h_over_h - 1) + sqrt_phi_plus_1 * (8 * last_h_over_h * (2 + last_h_over_h) - 7)) / 15;
                }
                else {
                    return 2 * (sqrt_phi_plus_1 * (1 - 3 * last_h_over_h - 4 * last_h_over_h * last_h_over_h) + sqrt_phi_plus_2 * (1 + last_h_over_h + 4 * last_h_over_h * last_h_over_h)) / 15;
                }
            }
        }
    }

    else { // not implemented with substeping yet
        if (n > 6){
            if (j < n - 3){
                return SphericSwimmingParticle<TBaseElement>::mAjs[j];
            }
            else if (j == n - 3){
                return SphericSwimmingParticle<TBaseElement>::mBns[n];
            }
            else if (j == n - 2){
                return SphericSwimmingParticle<TBaseElement>::mCns[n];
            }
            else if (j == n - 1){
                return SphericSwimmingParticle<TBaseElement>::mDns[n];
            }
            else {
                return SphericSwimmingParticle<TBaseElement>::mEns[n];
            }
        }

        else if (n == 2){ // use formula for the betas of second order
            SphericSwimmingParticle<TBaseElement>::GetDaitcheCoefficient(2, n, j, last_h_over_h, n_steps_per_quad_step);
        }

        else if (n == 3){
            long double sqrt_3_over_105 = std::sqrt(static_cast<long double>(3)) / 105;
            if (j == 0){
                return 68 * sqrt_3_over_105;
            }
            else if (j == 1){
                return 90 * sqrt_3_over_105;
            }
            else if (j == 2){
                return 36 * sqrt_3_over_105;
            }
            else {
                return 16 * sqrt_3_over_105;
            }
        }
        else if (n == 4){
            long double sqrt_2_over_315 = std::sqrt(static_cast<long double>(2)) / 315;
            long double OneOver315 = 1.0 / 315;

            if (j == 0){
                return 244 * sqrt_2_over_315;
            }
            else if (j == 1){
                return 1888 * OneOver315 - 976 * sqrt_2_over_315;
            }
            else if (j == 2){
                return - 656 * OneOver315 + 1464 * sqrt_2_over_315;
            }
            else if (j == 3){
                return 1632 - 976 * OneOver315;
            }
            else {
                return - 292 * OneOver315 + 244 * sqrt_2_over_315;
            }
        }
        else if (n == 5){
            long double sqrt_2_over_315 = std::sqrt(static_cast<long double>(2)) / 315;
            long double sqrt_3_over_105 = std::sqrt(static_cast<long double>(3)) / 105;
            long double sqrt_5_over_63  = std::sqrt(static_cast<long double>(5)) / 63;

            if (j == 0){
                return 244 * sqrt_2_over_315;
            }
            else if (j == 1){
                return 362 * sqrt_3_over_105 - 976 * sqrt_2_over_315;
            }
            else if (j == 2){
                return 500 * sqrt_5_over_63 - 1448 * sqrt_3_over_105 + 1464 * sqrt_2_over_315;
            }
            else if (j == 3){
                return - 870 * sqrt_5_over_63 + 2172 * sqrt_3_over_105 - 976 * sqrt_2_over_315;
            }
            else if (j == 4){
                return 660 * sqrt_5_over_63 - 1448 * sqrt_3_over_105 + 244 * sqrt_2_over_315;
            }
            else {
                return - 164 * sqrt_5_over_63 + 362 * sqrt_3_over_105;
            }
        }
        else {
            long double sqrt_2_over_315 = std::sqrt(static_cast<long double>(2)) / 315;
            long double sqrt_3_over_105 = std::sqrt(static_cast<long double>(3)) / 105;
            long double sqrt_6_over_105 = std::sqrt(static_cast<long double>(6)) / 105;
            long double OneOver315 = 1.0 / 315;

            if (j == 0){
                return 244 * sqrt_2_over_315;
            }
            else if (j == 1){
                return 362 * sqrt_3_over_105 - 976 * sqrt_2_over_315;
            }
            else if (j == 2){
                return 5584 * OneOver315 - 1448 * sqrt_3_over_105 + 1464 * sqrt_2_over_315;
            }
            else if (j == 3){
                return 1720 * sqrt_6_over_105 - 22336 * OneOver315 + 2172 * sqrt_3_over_105 - 976 * sqrt_2_over_315;
            }
            else if (j == 4){
                return - 3564 * sqrt_6_over_105 + 33504 * OneOver315 - 1448 * sqrt_3_over_105 + 244 * sqrt_2_over_315;
            }
            else if (j == 5){
                return 2808 * sqrt_6_over_105 - 22336 * OneOver315 + 362 * sqrt_3_over_105;
            }
            else {
                return - 754 * sqrt_6_over_105 + 5584 * OneOver315;
            }
        }
    }
    return 0.0;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::CalculateExplicitFractionalDerivative(NodeType& node, array_1d<double, 3>& fractional_derivative,
                                                                                  double& present_coefficient,
                                                                                  DenseVector<double>& historic_integrands,
                                                                                  const double last_h_over_h,
                                                                                  const int n_steps_per_quad_step)
{
    const int N = historic_integrands.size() - 3;
    const int n = (int)N / 3;
    double fast_fractional_derivative[3] = {0.0};

    for (int j = 0; j < n + 1; j++){
        double coefficient = GetDaitcheCoefficient(mQuadratureOrder, n + 1, j + 1, last_h_over_h, n_steps_per_quad_step);
        for (int i_comp = 0; i_comp < 3; i_comp++){
            unsigned int integrand_component_position = N - 3 * j + i_comp;
            fast_fractional_derivative[i_comp] += coefficient * historic_integrands[integrand_component_position];
        }
    }

    present_coefficient = GetDaitcheCoefficient(mQuadratureOrder, n + 1, 0, last_h_over_h, n_steps_per_quad_step);
    noalias(fractional_derivative) = present_coefficient * (node.FastGetSolutionStepValue(SLIP_VELOCITY) - node.FastGetSolutionStepValue(VELOCITY));
    SWIMMING_ADD_SECOND_TO_FIRST(fractional_derivative, fast_fractional_derivative)
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::Phi(const double x)
{
    if (std::abs(x) < 1e-10){
        return (std::exp(x) - 1) / x;
    }
    else {
        return 1 + 0.5 * x + 1.0 / 6 * x * x;
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::Ki(const double alpha, const double beta, const double time)
{
    return alpha * std::exp(beta * time);
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::AddFdi(const int order, array_1d<double, 3>& F, const double t_win, const double alpha, const double beta, const double phi, const double dt, const DenseVector<double>& historic_integrands, const array_1d<double, 3>& oldest_integrand)
{
    if (order == 1){
        const double beta_dt = beta * dt;
        const double coeff = - alpha / beta * std::exp(beta * (t_win - dt + dt * phi));
        const double coeff_N = 1 - Phi(beta_dt);
        const double coeff_N_plus_1 = std::exp(beta_dt) * (Phi(- beta_dt) - 1);

        F[0] +=  coeff * (coeff_N * historic_integrands[0] + coeff_N_plus_1 * oldest_integrand[0]);
        F[1] +=  coeff * (coeff_N * historic_integrands[1] + coeff_N_plus_1 * oldest_integrand[1]);
        F[2] +=  coeff * (coeff_N * historic_integrands[2] + coeff_N_plus_1 * oldest_integrand[2]);
    }

    else if (order == 2){
        if (false){ // unstable
            const double coeff = 0.5 * alpha / (beta * SWIMMING_POW_2(dt * beta));
            const double exp_1 = std::exp(beta * (t_win - dt + dt * phi));
            const double exp_2 = std::exp(beta * dt);
            const double f20 = oldest_integrand[0];
            const double f21 = oldest_integrand[1];
            const double f22 = oldest_integrand[2];
            const double f10 = historic_integrands[0];
            const double f11 = historic_integrands[1];
            const double f12 = historic_integrands[2];
            const double f00 = historic_integrands[3];
            const double f01 = historic_integrands[4];
            const double f02 = historic_integrands[5];
            F[0] += coeff * exp_1 * ((2 * exp_2 * (dt * beta) - dt * beta - 2)                                                * f00
                                   + (exp_2 * (4 - 4 * dt * beta) + 2 * beta * SWIMMING_POW_2(dt * beta) - 4)                 * f10
                                   + (exp_2 * (2 - 3 * dt * beta + 2 * 2 * beta * SWIMMING_POW_2(dt * beta)) + dt * beta - 2) * f20);
            F[1] += coeff * exp_1 * ((2 * exp_2 * (dt * beta) - dt * beta - 2)                                                * f01
                                   + (exp_2 * (4 - 4 * dt * beta) + 2 * beta * SWIMMING_POW_2(dt * beta) - 4)                 * f11
                                   + (exp_2 * (2 - 3 * dt * beta + 2 * 2 * beta * SWIMMING_POW_2(dt * beta)) + dt * beta - 2) * f21);
            F[2] += coeff * exp_1 * ((2 * exp_2 * (dt * beta) - dt * beta - 2)                                                * f02
                                   + (exp_2 * (4 - 4 * dt * beta) + 2 * beta * SWIMMING_POW_2(dt * beta) - 4)                 * f12
                                   + (exp_2 * (2 - 3 * dt * beta + 2 * 2 * beta * SWIMMING_POW_2(dt * beta)) + dt * beta - 2) * f22);
        }
        else {
            const double t_minus_t2 = t_win + dt * phi;
            const double t_minus_t1 = t_win + dt * phi - dt;
            const double t_minus_t0 = t_win + dt * phi - 2 * dt;
            const double Ki2 = Ki(alpha, beta, t_minus_t2);
            const double Ki1 = Ki(alpha, beta, t_minus_t1);
            const double Ki0 = Ki(alpha, beta, t_minus_t0);
            const double f20 = oldest_integrand[0] * Ki2;
            const double f21 = oldest_integrand[1] * Ki2;
            const double f22 = oldest_integrand[2] * Ki2;
            const double f10 = historic_integrands[0] * Ki1;
            const double f11 = historic_integrands[1] * Ki1;
            const double f12 = historic_integrands[2] * Ki1;
            const double f00 = historic_integrands[3] * Ki0;
            const double f01 = historic_integrands[4] * Ki0;
            const double f02 = historic_integrands[5] * Ki0;
            const double aux_coeff = dt / 12;
            F[0] += aux_coeff * (- f00 + 8 * f10 + 5 * f20);
            F[1] += aux_coeff * (- f01 + 8 * f11 + 5 * f21);
            F[2] += aux_coeff * (- f02 + 8 * f12 + 5 * f22);
        }
    }
    else {
        return;
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::AddFre(array_1d<double, 3>& old_Fi, const double beta, const double dt)
{
    const double exp_coeff = std::exp(beta * dt);
    old_Fi[0] = exp_coeff * old_Fi[0];
    old_Fi[1] = exp_coeff * old_Fi[1];
    old_Fi[2] = exp_coeff * old_Fi[2];
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::AddHinsbergTailContribution(NodeType& node, array_1d<double, 3>& fractional_derivative_of_slip_vel, const int order, const int n_steps_per_quad_step, const double time, const double quadrature_delta_time, const double last_h_over_h, DenseVector<double>& historic_integrands)
{
    DenseVector<double>& hinsberg_tail_contributions = node.GetValue(HINSBERG_TAIL_CONTRIBUTIONS);
    int m = hinsberg_tail_contributions.size() / 3 - 1; // number of exponentials: the last three slots hold the components of the oldest historic integrand
    const double t_win = SphericSwimmingParticle<TBaseElement>::mTimeWindow;

    if (n_steps_per_quad_step * last_h_over_h < 1.5 && 2 * n_steps_per_quad_step * (time - t_win) > quadrature_delta_time){ // first calculation after appending and already later than t_win
        array_1d<double, 3> oldest_integrand;
        oldest_integrand[0] = hinsberg_tail_contributions[3 * m];
        oldest_integrand[1] = hinsberg_tail_contributions[3 * m + 1];
        oldest_integrand[2] = hinsberg_tail_contributions[3 * m + 2];
        const std::vector<double>& Ts = SphericSwimmingParticle<SphericParticle>::mTs;
        const double e = std::exp(1);
        array_1d<double, 3> Fi;

        for (int i = 0; i < m; i++){
            const double ti = Ts[i];
            const double beta = - 0.5 / ti;
            const double alpha = std::sqrt(e / ti);
            Fi[0] = hinsberg_tail_contributions[3 * i];
            Fi[1] = hinsberg_tail_contributions[3 * i + 1];
            Fi[2] = hinsberg_tail_contributions[3 * i + 2];
            AddFre(Fi, beta, quadrature_delta_time);
            AddFdi(order, Fi, t_win, alpha, beta, 1.0, quadrature_delta_time, historic_integrands, oldest_integrand);
            hinsberg_tail_contributions[3 * i]     = Fi[0];
            hinsberg_tail_contributions[3 * i + 1] = Fi[1];
            hinsberg_tail_contributions[3 * i + 2] = Fi[2];
        }
    }

    array_1d<double, 3> F_tail = ZeroVector(3);

    for (int i = 0; i < m; i++){
        double ai = SphericSwimmingParticle<TBaseElement>::mAs[i];
        F_tail[0] += ai * hinsberg_tail_contributions[3 * i];
        F_tail[1] += ai * hinsberg_tail_contributions[3 * i + 1];
        F_tail[2] += ai * hinsberg_tail_contributions[3 * i + 2];
    }

    const double sqrt_delta_time_inv = 1.0 / std::sqrt(quadrature_delta_time); // since the multiplication by sqrt(delta_time) corresponding to the F_win part is done to F_win + F_tail outside
    noalias(fractional_derivative_of_slip_vel) += sqrt_delta_time_inv * F_tail;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::AddHinsbergTailContributionStrict(NodeType& node, array_1d<double, 3>& fractional_derivative_of_slip_vel, const int order, const int n_steps_per_quad_step, const double time, const double quadrature_delta_time, const double last_h_over_h, DenseVector<double>& historic_integrands)
{
    DenseVector<double>& hinsberg_tail_contributions = node.GetValue(HINSBERG_TAIL_CONTRIBUTIONS);
    int m = hinsberg_tail_contributions.size() / 3 - 1; // number of exponentials: the last three slots hold the components of the oldest historic integrand

    if (m < 1){ // trivial, 0-exponentials case (there is no tail contribution)
        return;
    }

    const double t_win = SphericSwimmingParticle<TBaseElement>::mTimeWindow;
    const std::vector<double>& Ts = SphericSwimmingParticle<SphericParticle>::mTs;
    const double e = std::exp(1);
    const double delta_time = quadrature_delta_time / n_steps_per_quad_step;

    if (n_steps_per_quad_step * last_h_over_h < 1.5 && 2 * n_steps_per_quad_step * (time - t_win) > quadrature_delta_time){ // calculation right after last append (but at least later than t = t_win, so there is some tail)
        array_1d<double, 3> oldest_integrand;
        oldest_integrand[0] = hinsberg_tail_contributions[3 * m];
        oldest_integrand[1] = hinsberg_tail_contributions[3 * m + 1];
        oldest_integrand[2] = hinsberg_tail_contributions[3 * m + 2];
        array_1d<double, 3> Fi;

        for (int i = 0; i < m; i++){
            const double ti = Ts[i];
            const double beta = - 0.5 / ti;
            const double alpha = std::sqrt(e / ti);
            Fi[0] = hinsberg_tail_contributions[3 * i];
            Fi[1] = hinsberg_tail_contributions[3 * i + 1];
            Fi[2] = hinsberg_tail_contributions[3 * i + 2];
            AddFre(Fi, beta, delta_time);
            AddFdi(order, Fi, t_win, alpha, beta, last_h_over_h, quadrature_delta_time, historic_integrands, oldest_integrand);
            hinsberg_tail_contributions[3 * i]     = Fi[0];
            hinsberg_tail_contributions[3 * i + 1] = Fi[1];
            hinsberg_tail_contributions[3 * i + 2] = Fi[2];
        }
    }

    else { // intermediate step: the time assigned to the tail has not changed (only an exponential factor is needed to correct fot the changing current time, which affects the approximate kernel)

        for (int i = 0; i < m; i++){
            const double ti = Ts[i];
            const double beta = - 0.5 / ti;
            const double exp_beta_dt = std::exp(beta * delta_time);
            hinsberg_tail_contributions[3 * i]     *= exp_beta_dt;
            hinsberg_tail_contributions[3 * i + 1] *= exp_beta_dt;
            hinsberg_tail_contributions[3 * i + 2] *= exp_beta_dt;
        }
    }

    array_1d<double, 3> F_tail = ZeroVector(3);
    const std::vector<double>& As = SphericSwimmingParticle<TBaseElement>::mAs;

    for (int i = 0; i < m; i++){
        const double ai = As[i];
        F_tail[0] += ai * hinsberg_tail_contributions[3 * i];
        F_tail[1] += ai * hinsberg_tail_contributions[3 * i + 1];
        F_tail[2] += ai * hinsberg_tail_contributions[3 * i + 2];
    }

    const double sqrt_delta_time_inv = 1.0 / std::sqrt(quadrature_delta_time); // since the multiplication by sqrt(delta_time) corresponding to the F_win part is done to F_win + F_tail outside
    noalias(fractional_derivative_of_slip_vel) += sqrt_delta_time_inv * F_tail;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeBassetForce(NodeType& node, array_1d<double, 3>& basset_force, const ProcessInfo& r_current_process_info)
{
    if (mBassetForceType == 0 || node.IsNot(INSIDE) || node.Is(BLOCKED)){ // case of identically null history force
        noalias(basset_force) = ZeroVector(3);
        return;
    }

    else {
        const double basset_force_coeff = 6.0 * mRadius * mRadius * mFluidDensity * std::sqrt(Globals::Pi * mKinematicViscosity);
        const double delta_time = r_current_process_info[DELTA_TIME];
        int n_steps_per_quad_step = r_current_process_info[TIME_STEPS_PER_QUADRATURE_STEP];
        const double quadrature_delta_time = n_steps_per_quad_step * delta_time;

        if (r_current_process_info[TIME_STEPS] >= r_current_process_info[NUMBER_OF_INIT_BASSET_STEPS]){
            DenseVector<double>& historic_integrands = node.GetValue(BASSET_HISTORIC_INTEGRANDS);
            const double time = r_current_process_info[TIME];
            const double latest_quadrature_time_step = time + delta_time - r_current_process_info[LAST_TIME_APPENDING];
            array_1d<double, 3> fractional_derivative_of_slip_vel;
            double present_coefficient;
            const double sqrt_of_quad_h_q = std::sqrt(quadrature_delta_time);
            const double last_h_over_h = latest_quadrature_time_step / quadrature_delta_time;


            CalculateExplicitFractionalDerivative(node,
                                                  fractional_derivative_of_slip_vel,
                                                  present_coefficient,
                                                  historic_integrands,
                                                  last_h_over_h,
                                                  n_steps_per_quad_step);

            if (r_current_process_info[FRAME_OF_REFERENCE_TYPE] >= 1){
                const array_1d<double, 3>& displacement = node.FastGetSolutionStepValue(DISPLACEMENT);
                const array_1d<double, 3>& displacement_old = node.FastGetSolutionStepValue(DISPLACEMENT_OLD);
                array_1d<double, 3> aux = displacement - displacement_old;
                const array_1d<double, 3>& omega_frame = r_current_process_info[ANGULAR_VELOCITY_MOVING_FRAME];
                array_1d<double, 3> correction;
                SWIMMING_SET_TO_CROSS_OF_FIRST_TWO_3(omega_frame, aux, correction);
                noalias(fractional_derivative_of_slip_vel) += present_coefficient * correction;
            }

            if (mBassetForceType == 3){
                AddHinsbergTailContribution(node, fractional_derivative_of_slip_vel, mQuadratureOrder, n_steps_per_quad_step, time, quadrature_delta_time, last_h_over_h, historic_integrands);
            }

            if (mBassetForceType == 4){
                AddHinsbergTailContributionStrict(node, fractional_derivative_of_slip_vel, mQuadratureOrder, n_steps_per_quad_step, time, quadrature_delta_time, last_h_over_h, historic_integrands);
            }

            array_1d<double, 3> basset_term    = fractional_derivative_of_slip_vel;
            array_1d<double, 3> old_basset_term;
            const array_1d<double, 3>& vel     = node.FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3>& old_vel = node.FastGetSolutionStepValue(VELOCITY_OLD);
            old_basset_term = mOldBassetTerm + mOldDaitchePresentCoefficient * (old_vel - vel); // the second term corresponds to the part that was treated implicitly in the last step minus a part that was added but did not correspond to the basset term
            noalias(mOldBassetTerm) = basset_term;

            if (r_current_process_info[FRAME_OF_REFERENCE_TYPE] >= 1){
                array_1d<double, 3>& displacement_old = node.FastGetSolutionStepValue(DISPLACEMENT_OLD);
                array_1d<double, 3>& displacement_old_old = node.FastGetSolutionStepValue(VELOCITY_OLD_OLD);
                array_1d<double, 3> aux = mOldDaitchePresentCoefficient * (displacement_old_old - displacement_old);
                const array_1d<double, 3>& omega_frame = r_current_process_info[ANGULAR_VELOCITY_MOVING_FRAME];
                array_1d<double, 3> correction;
                SWIMMING_SET_TO_CROSS_OF_FIRST_TWO_3(omega_frame, aux, correction);
                noalias(displacement_old_old) = displacement_old;
                noalias(displacement_old) = node.FastGetSolutionStepValue(DISPLACEMENT);
                // correcting the old_basset_term for the velocities at different times being subtracted
                noalias(old_basset_term) += correction;
                noalias(aux) = delta_time * (1.5 * basset_term - 0.5 * old_basset_term);
                SWIMMING_SET_TO_CROSS_OF_FIRST_TWO_3(omega_frame, aux, correction);
                // correcting the basset force for the discretization of the time derivative
                noalias(fractional_derivative_of_slip_vel) += correction;
            }

            noalias(fractional_derivative_of_slip_vel) -= old_basset_term;

            mOldDaitchePresentCoefficient = present_coefficient;
            mLastBassetForceAddedMass = basset_force_coeff * sqrt_of_quad_h_q * present_coefficient;

            noalias(basset_force) = basset_force_coeff * sqrt_of_quad_h_q / delta_time * fractional_derivative_of_slip_vel;
        }

        else {
            basset_force = node.FastGetSolutionStepValue(BASSET_FORCE);
            mOldDaitchePresentCoefficient = 0.0;
            mOldBassetTerm = std::sqrt(quadrature_delta_time) / basset_force_coeff * basset_force;
        }
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeSaffmanLiftForce(NodeType& node, array_1d<double, 3>& lift_force, const ProcessInfo& r_current_process_info)
{
    if (mSaffmanForceType == 0 || node.IsNot(INSIDE) || node.Is(BLOCKED)){ // case of identically null lift force
        noalias(lift_force) = ZeroVector(3);

        return;
    }

    else if (mSaffmanForceType >= 1){
        const double& shear_rate                       = node.FastGetSolutionStepValue(SHEAR_RATE_PROJECTED);
        const array_1d<double, 3>& vorticity           = node.FastGetSolutionStepValue(FLUID_VORTICITY_PROJECTED);
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
template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeMagnusLiftForce(NodeType& node, array_1d<double, 3>& lift_force, const ProcessInfo& r_current_process_info)
{
    if (mMagnusForceType == 0 || node.IsNot(INSIDE) || node.Is(BLOCKED)){
        noalias(lift_force) = ZeroVector(3);
        return;
    }

    const array_1d<double, 3> slip_rot = 0.5 * node.FastGetSolutionStepValue(FLUID_VORTICITY_PROJECTED) - node.FastGetSolutionStepValue(ANGULAR_VELOCITY);
    array_1d<double, 3> slip_rot_cross_slip_vel;
    SWIMMING_SET_TO_CROSS_OF_FIRST_TWO_3(slip_rot, mSlipVel, slip_rot_cross_slip_vel)

    if (mMagnusForceType == 1){ // Rubinow and Keller, 1961 (Re_p < 0.1; nondimensional_slip_rot_vel < 0.1)
        lift_force = Globals::Pi * SWIMMING_POW_3(mRadius) * mFluidDensity * slip_rot_cross_slip_vel;
    }

    else if (mMagnusForceType == 2){ // Oesterle and Bui Dihn, 1998 (Re_p < 140)
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
            noalias(lift_force) = 0.5 *  mFluidDensity * Globals::Pi * SWIMMING_POW_2(mRadius) * lift_coeff * mNormOfSlipVel * slip_rot_cross_slip_vel / norm_of_slip_rot;
        }
    }

    else if (mMagnusForceType == 3){ // Loth, 2008 (Re_p < 2000; nondimensional_slip_rot_vel < 20)
        // calculate as in Rubinow and Keller, 1963
        noalias(lift_force) = Globals::Pi * SWIMMING_POW_3(mRadius) * mFluidDensity * slip_rot_cross_slip_vel;
        // correct coefficient
        double reynolds;
        ComputeParticleReynoldsNumber(reynolds);
        const double nondimensional_slip_rot_vel = ComputeNondimensionalRotVelocity(slip_rot);
        const double coeff = 1 - (0.675 + 0.15 * (1 + std::tanh(0.28 * (nondimensional_slip_rot_vel - 2)))) * std::tanh(0.18 * std::sqrt(reynolds));

        noalias(lift_force) = coeff * lift_force;
    }

    else {
        std::cout << "The integer value designating the magnus lift coefficient calculation model" << std::endl;
        std::cout << " (mMagnusForceType = " << mMagnusForceType << "), is not supported" << std::endl << std::flush;
        return;
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeHydrodynamicTorque(NodeType& node, array_1d<double, 3>& hydro_torque, const ProcessInfo& r_current_process_info)
{
    if (mHydrodynamicTorqueType == 0 || node.IsNot(INSIDE) || node.Is(BLOCKED)){
        noalias(hydro_torque) = ZeroVector(3);
        return;
    }

    else if (mHydrodynamicTorqueType == 1 || mHydrodynamicTorqueType == 2){
        const array_1d<double, 3> slip_rot = 0.5 * node.FastGetSolutionStepValue(FLUID_VORTICITY_PROJECTED) - node.FastGetSolutionStepValue(ANGULAR_VELOCITY);
        const double norm_of_slip_rot = SWIMMING_MODULUS_3(slip_rot);
        double rot_reynolds;
        ComputeParticleRotationReynoldsNumberOverNormOfSlipRot(rot_reynolds);
        double rotational_coeff;

        if (rot_reynolds > 32){ // Rubinow and Keller, 1961 (Re_rot ~ 32 - 1000)
            rotational_coeff = 12.9 * std::sqrt(norm_of_slip_rot * rot_reynolds) + 128.4 / rot_reynolds;
        }

        else { // Rubinow and Keller, 1961 (Re_rot < 32)
            rotational_coeff = 64 * Globals::Pi / rot_reynolds;
        }

        if (mHydrodynamicTorqueType == 2){ // Loth, 2008 (Re_rot < 2000, Re_p < 20)
            rotational_coeff *= 1.0 + 5 / (64 * Globals::Pi) * std::pow(rot_reynolds, 0.6);
        }

        noalias(hydro_torque) = 0.5 * mFluidDensity * SWIMMING_POW_5(mRadius) * rotational_coeff * slip_rot;

    }

    else {
        std::cout << "The integer value designating the Hydrodynamic torque calculation model" << std::endl;
        std::cout << " (mHydrodynamicTorqueType = " << mHydrodynamicTorqueType << "), is not supported" << std::endl << std::flush;
        return;
    }

}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeBrownianMotionForce(NodeType& node, array_1d<double, 3>& brownian_motion_force, const ProcessInfo& r_current_process_info)
{
    if (mBrownianMotionType == 0 || node.IsNot(INSIDE) || node.Is(BLOCKED)){
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
        double coeff = std::sqrt(24 * kT * 2 / Globals::Pi * ComputeStokesDragCoefficient() * delta_t_inv);
        brownian_motion_force[0] = coeff * dis(gen);
        brownian_motion_force[1] = coeff * dis(gen);
        brownian_motion_force[2] = coeff * dis(gen);*/
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeParticleReynoldsNumber(double& reynolds)
{
    reynolds = 2 * mRadius * mNormOfSlipVel / mKinematicViscosity;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputePowerLawParticleReynoldsNumber(double& reynolds,
                                                                                  const double consistency_index,
                                                                                  const double flow_behavior_index)
{
    // This function is consistent with Shah 2007 (doi:10.1016/j.ijmultiphaseﬂow.2006.06.006)
    // int coefficient = use_max_shear_rate ? 3 : 2;
    const double& K = consistency_index;
    const double& n = flow_behavior_index;
    reynolds =  2 * std::pow(mRadius, n) * std::pow(mNormOfSlipVel, 2 - n) * mFluidDensity / K;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::ComputeNondimensionalRotVelocity(const array_1d<double, 3>& slip_rot_velocity)
{
    if (mNormOfSlipVel > 0){
        return 2.0 * mRadius * SWIMMING_MODULUS_3(slip_rot_velocity) / mNormOfSlipVel;
    }

    else {
        return 0.0;
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeParticleRotationReynoldsNumber(double norm_of_slip_rot, double& reynolds)
{
    reynolds = 4 * SWIMMING_POW_2(mRadius) * norm_of_slip_rot / mKinematicViscosity;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeParticleRotationReynoldsNumberOverNormOfSlipRot(double& reynolds)
{
    reynolds = 4 * SWIMMING_POW_2(mRadius) / mKinematicViscosity;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeParticleAccelerationNumber(const array_1d<double, 3>& slip_acc, double& acc_number)
{
    acc_number = SWIMMING_POW_3(mNormOfSlipVel) / fabs(2 * mRadius * SWIMMING_INNER_PRODUCT_3(mSlipVel, slip_acc));
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::AdditionalCalculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_current_process_info)
{
    if (rVariable == REYNOLDS_NUMBER){
        NodeType& node = GetGeometry()[0];

        if (node.IsNot(INSIDE)){
            Output = 0.0;
        }

        else {
            mFluidDensity                           = node.FastGetSolutionStepValue(FLUID_DENSITY_PROJECTED);
            mKinematicViscosity                     = node.FastGetSolutionStepValue(FLUID_VISCOSITY_PROJECTED);
            mFluidFraction                          = node.FastGetSolutionStepValue(FLUID_FRACTION_PROJECTED);
            const array_1d<double, 3>& fluid_vel    = node.FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
            const array_1d<double, 3>& particle_vel = node.FastGetSolutionStepValue(VELOCITY);

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

    else if (rVariable == TIME){
        Output = mInitialTime;
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::Calculate(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & Output, const ProcessInfo& r_process_info)
{
    if (rVariable == VIRTUAL_MASS_FORCE) {
        const array_1d<double, 3 > total_forces = GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES);
        Output -= mLastVirtualMassAddedMass / mRealMass * total_forces;
    }

    else if (rVariable == BASSET_FORCE) {
        const array_1d<double, 3 > total_forces = GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES);
        Output -= mLastBassetForceAddedMass / mRealMass * total_forces;
    }

    else {
        TBaseElement::Calculate(rVariable, Output, r_process_info);
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::ComputeDragCoefficient(const ProcessInfo& r_current_process_info)
{
    double drag_coeff;

    if (mDragForceType == 1 || mDragForceType == 10){
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
        drag_coeff = 2.0 / Globals::Pi * ComputeStokesDragCoefficient();
    }

    else if (mDragForceType == 11){ // Maxey-Riley expression with Faxen correction
        drag_coeff = ComputeStokesDragCoefficient(); // temporary
    }

    else if (mDragForceType == 13){ // Re_p < 1000, Shah et al. (2006) (doi:10.1016/j.ijmultiphaseflow.2006.06.006)
        drag_coeff = ComputeShahDragCoefficient(r_current_process_info);
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
template < class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::ComputeStokesDragCoefficient()
{
    double drag_coeff = 6.0 * Globals::Pi * mKinematicViscosity * mFluidDensity * mRadius;

    return drag_coeff;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::ComputeWeatherfordDragCoefficient(const ProcessInfo& r_current_process_info)
{
    KRATOS_TRY

    const double particle_density              = GetDensity();
    const array_1d<double, 3>& gravity         = r_current_process_info[GRAVITY];
    const int manually_imposed_drag_law_option = r_current_process_info[MANUALLY_IMPOSED_DRAG_LAW_OPTION];
    const int drag_modifier_type               = r_current_process_info[DRAG_MODIFIER_TYPE];

    int non_newtonian_option = 1;
    const double yield_stress                  = GetGeometry()[0].FastGetSolutionStepValue(YIELD_STRESS);
    const double power_law_n                   = GetGeometry()[0].FastGetSolutionStepValue(POWER_LAW_N);

    if (fabs(power_law_n - 1.0) < 0.00001  ||  fabs(yield_stress) < 0.00001) {
        non_newtonian_option = 0;
    }

    const double area                          = Globals::Pi * SWIMMING_POW_2(mRadius);
    const array_1d<double, 3> weight           = mRealMass * gravity;
    const array_1d<double, 3> buoyancy         = mFluidDensity / particle_density * weight; // hydrostatic case!! (only for Weatherford)

    double reynolds;
    double drag_coeff;

    if (!non_newtonian_option){ // Newtonian
        ComputeParticleReynoldsNumber(reynolds);

        if (!non_newtonian_option && reynolds < 0.01){
            reynolds = 0.01;
        }

	// CalculateNewtonianDragCoefficient(non_newtonian_option, reynolds, mSphericity, drag_coeff, drag_modifier_type);
	// drag_coeff *= 0.5 *  mFluidDensity * area * drag_coeff * mNormOfSlipVel;

	//to avoid numerical problems it is better to compute the drag coeff as follows (division by zero and multiplication among very small numbers and big numbers are avoided)
	if(reynolds==0){
	  drag_coeff = 6.0 * Globals::Pi * mKinematicViscosity * mFluidDensity * mRadius * (0.5 *  mFluidDensity * area * mNormOfSlipVel) ; //I use the Stoke coefficient
        }else if (reynolds > 1000){
	  drag_coeff = 0.44 *  (0.5 *  mFluidDensity * area * mNormOfSlipVel);
	}else{
	  drag_coeff = 12.0 * mKinematicViscosity / mRadius * (1.0 + 0.15 * pow(reynolds, 0.687)) * (0.5 *  mFluidDensity * area);
      }
    }

    else {
        const double gel_strength                  = GetGeometry()[0].FastGetSolutionStepValue(YIELD_STRESS);
        const double power_law_K                   = GetGeometry()[0].FastGetSolutionStepValue(POWER_LAW_K);

        double beta                                = 0.0;
        double F0                                  = 0.0;
        double regularization_v                    = 0.02 * mRadius;
        const double power_law_tol                 = r_current_process_info[POWER_LAW_TOLERANCE];

        double shahs_term_vel = CalculateShahsTerm(power_law_n, power_law_K, power_law_tol, particle_density, mSphericity, drag_modifier_type);

        if (!manually_imposed_drag_law_option){
            F0 = 4.0 * gel_strength * area; //initial value
            beta = (SWIMMING_MODULUS_3(weight) - SWIMMING_MODULUS_3(buoyancy) - F0) / shahs_term_vel; //S
        }

        else {
            const double initial_drag_force            = r_current_process_info[INIT_DRAG_FORCE];
            const double drag_law_slope                = r_current_process_info[DRAG_LAW_SLOPE];
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
template < class TBaseElement >
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
template < class TBaseElement >
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
template < class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::CalculateShahsTerm(double power_law_N,
                                                   double power_law_K,
                                                   double power_law_tol,
                                                   const double& particle_density,
                                                   double sphericity,
                                                   int drag_modifier_type)
{
    if (std::abs(power_law_N) < power_law_tol || std::abs(power_law_K) < power_law_tol){
        std::cout << "WARNING: Shah's method is being used with Power Law data being zero!!" << std::endl << std::flush;
    }

    double shah_A_i = 1 / (6.9148 * power_law_N * power_law_N - 24.838 * power_law_N + 22.642);
    double shah_B_i = 1 / (-0.5067 * power_law_N * power_law_N + 1.3234 * power_law_N - 0.1744);

    double dimensionless_shah = std::sqrt(std::pow(13.08, 2 - power_law_N) * std::pow(2 * mRadius, power_law_N + 2) * std::pow( mFluidDensity, power_law_N) * std::pow(particle_density -  mFluidDensity, 2 - power_law_N) / (std::pow(2, 2 * (power_law_N - 1)) * power_law_K * power_law_K));
    double reynolds = std::pow(dimensionless_shah * shah_A_i, shah_B_i);
    double fi_i = CalculateDragCoeffFromSphericity(reynolds, 1.0, drag_modifier_type) / CalculateDragCoeffFromSphericity(reynolds, sphericity, drag_modifier_type);
    dimensionless_shah = std::sqrt(std::pow(fi_i, 2 - power_law_N)) * dimensionless_shah;
    reynolds = std::pow(dimensionless_shah * shah_A_i, shah_B_i);

    double terminal_vel =  std::pow(std::pow(2, power_law_N - 1) * power_law_K * reynolds / (std::pow(2 * mRadius, power_law_N) *  mFluidDensity), 1 / (2 - power_law_N)) ;

    return terminal_vel;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::ComputeGanserDragCoefficient(const ProcessInfo& r_current_process_info)
{
    KRATOS_TRY

    const int isometric_shape                = 1; // TEMPORARY!! yes (1) or no (0); shold be given as data
    const double surface_area                = 4 * Globals::Pi * SWIMMING_POW_2(mRadius); // TEMPORARY!! corresponding to a sphere; should be generalized b taking it as a parameter
    const double surface_area_circular_diam  = std::sqrt(4.0 * surface_area / Globals::Pi);

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
template < class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::ComputeIshiiDragCoefficient(const ProcessInfo& r_current_process_info)
{
    double coeff = 0.45;
    double reynolds;
    ComputeParticleReynoldsNumber(reynolds);

    if (reynolds <= 1000){
        coeff = (24 + 2.4 * pow(reynolds, 0.75)) / reynolds;
    }

    double drag_coeff = 0.5 * coeff * Globals::Pi * SWIMMING_POW_2(mRadius);

    return drag_coeff;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::ComputeNewtonRegimeDragCoefficient()
{
    double drag_coeff  = 0.5 * Globals::Pi * SWIMMING_POW_2(mRadius) * mFluidDensity * mNormOfSlipVel;

    drag_coeff *= 0.44;

    return drag_coeff;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::ComputeIntermediateRegimeDragCoefficient()
{
    double reynolds;
    double drag_coeff  = 0.5 * Globals::Pi * SWIMMING_POW_2(mRadius) * mFluidDensity * mNormOfSlipVel;

    ComputeParticleReynoldsNumber(reynolds);

    drag_coeff *= 24 / reynolds * (1 + 0.15 * pow(reynolds, 0.687));

    return drag_coeff;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::ComputeHaiderDragCoefficient()
{
    const double sphericity = GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_SPHERICITY);
    double drag_coeff       = 0.5 * Globals::Pi * SWIMMING_POW_2(mRadius) * mFluidDensity * mNormOfSlipVel;

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
template < class TBaseElement >
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
        drag_coeff = Globals::Pi / 3.0 * mKinematicViscosity * mFluidDensity * mRadius * (A * eps_s / eps + B * particle_reynolds);
    }

    return drag_coeff;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::ComputeShahDragCoefficient(const ProcessInfo& r_current_process_info, const bool use_shahi_correction)
{
    const double power_law_tol = 0.0001;
    const double K = r_current_process_info[POWER_LAW_K];
    const double n = r_current_process_info[POWER_LAW_N];

    if (std::abs(n) < power_law_tol || std::abs(K) < power_law_tol){
        std::cout << "WARNING: Shah's method is being used with Power Law data being zero (n = 0 or K = 0)!!" << std::endl << std::flush;
    }

    double A =   6.9148 * n * n - 24.838 * n + 22.642;
    double B = - 0.5067 * n * n + 1.3234 * n - 0.1744;

    if (use_shahi_correction){ // from 2016 Shahi (doi: 10.1016/j.minpro.2016.06.002)
        A = 1.5269 * A - 3.9375;
        B =  0.892 * B + 0.0326;
    }

    double reynolds;
    ComputePowerLawParticleReynoldsNumber(reynolds, K, n);
    const double exponents_coeff = 1.0 / (2 - n);
    const double area = Globals::Pi * mRadius * mRadius;
    const double dimensional_coefficient = 0.5 * area * mFluidDensity * mNormOfSlipVel;

    return dimensional_coefficient * std::pow(A, exponents_coeff) * std::pow(reynolds, exponents_coeff * (2 * B - 2));
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::ComputeGanserParameters(const int isometric_shape, const double dn, double& k_1, double& k_2)
{
     if (isometric_shape){
         k_1 = 3 / (1 + 2 / std::sqrt(mSphericity));
     }

     else {
         k_1 = 3 / (0.5 * dn / mRadius + 2 / std::sqrt(mSphericity));
     }

     k_2 = pow(10.0, 1.8148 * pow(- log10(mSphericity), 0.5743));
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
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
template < class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::ComputeElSamniLiftCoefficient(const double norm_of_shear_rate,
                                                              const double vorticity_norm,
                                                              const ProcessInfo& r_current_process_info)
{
    if (vorticity_norm > 0.000000000001 && mNormOfSlipVel > 0.000000000001){
         const double yield_stress   = 0.0; // we are considering a Bingham type fluid
         const double power_law_K = r_current_process_info[POWER_LAW_K];
         const double power_law_n = r_current_process_info[POWER_LAW_N];
         const double shear_rate_p   = mNormOfSlipVel / mRadius * (4.5 / power_law_n - 3.5); // graphic model by Unhlherr et al. (fit by Wallis, G.B. and Dobson, J.E., 1973)
         double equivalent_viscosity = yield_stress / shear_rate_p + power_law_K * pow(shear_rate_p, power_law_n - 1);
         const double coeff          = std::max(0.09 * mNormOfSlipVel, 5.82 * std::sqrt(0.5 * mNormOfSlipVel * equivalent_viscosity /  mFluidDensity));
         const double lift_coeff     = 0.5 * Globals::Pi * SWIMMING_POW_2(mRadius) *  mFluidDensity * coeff * mNormOfSlipVel / vorticity_norm;
         return(lift_coeff);
    }

    else {
        return 0.0;
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
 double SphericSwimmingParticle<TBaseElement>::ComputeMeiLiftCoefficient(const double reynolds, const double reynolds_shear)
 {
     if (reynolds != 0.0 && reynolds_shear != 0.0 ){
         double sqrt_beta = std::sqrt(0.5 * reynolds_shear / reynolds);
         double C;

         if (reynolds < 40){
             C = (1 - 0.3314 * sqrt_beta) * exp(- 0.1 * reynolds) + 0.3314 * sqrt_beta;
         }

         else {
             C = 0.0524 * sqrt_beta * std::sqrt(reynolds);
         }

         C *= 4.1126 / std::sqrt(reynolds_shear);

         double lift_coeff = mFluidDensity * Globals::Pi * SWIMMING_POW_3(mRadius) * C;

         return lift_coeff;
     }

     else {
         return 0.0;
     }
 }

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
double SphericSwimmingParticle<TBaseElement>::GetFluidMass()
{
    return mFluidDensity * CalculateVolume();
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
array_1d<double,3> SphericSwimmingParticle<TBaseElement>::ComputeWeight(const array_1d<double,3>& gravity, const ProcessInfo& r_process_info)
{
    array_1d<double,3> weight = TBaseElement::ComputeWeight(gravity, r_process_info);

    if (r_process_info[FRAME_OF_REFERENCE_TYPE] >= 1){
        AddCentrifugalForces(weight, r_process_info);
        AddCoriolisForces(weight, r_process_info);

        if (r_process_info[FRAME_OF_REFERENCE_TYPE] >= 2){
            AddRelativeAccelerationForces(weight, r_process_info);
            AddEulerForces(weight, r_process_info);
        }
    }

   return weight;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::AddCentrifugalForces(array_1d<double,3>& weight, const ProcessInfo& r_process_info)
{
    const array_1d<double,3>& omega = r_process_info[ANGULAR_VELOCITY_MOVING_FRAME];
    const array_1d<double,3>& coordinates = GetGeometry()[0].Coordinates();
    array_1d<double,3> partial_result;
    array_1d<double,3> centrifugal;
    SWIMMING_SET_TO_CROSS_OF_FIRST_TWO_3(omega, coordinates, partial_result)
    SWIMMING_SET_TO_CROSS_OF_FIRST_TWO_3(omega, partial_result, centrifugal)
    noalias(weight) += (GetFluidMass() - GetMass()) * centrifugal;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::AddCoriolisForces(array_1d<double,3>& weight, const ProcessInfo& r_process_info)
{
    const array_1d<double,3>& omega = r_process_info[ANGULAR_VELOCITY_MOVING_FRAME];
    const array_1d<double,3>& particle_velocity = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& fluid_velocity = GetGeometry()[0].FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
    array_1d<double,3> coriolis_particle;
    array_1d<double,3> coriolis_fluid;
    SWIMMING_SET_TO_CROSS_OF_FIRST_TWO_3(omega, particle_velocity, coriolis_particle)
    SWIMMING_SET_TO_CROSS_OF_FIRST_TWO_3(omega, fluid_velocity, coriolis_fluid)
    const double fluid_mass = GetFluidMass();
    const double particle_mass = GetMass();
    noalias(weight) += 2 * (1.5 * fluid_mass * coriolis_fluid - (particle_mass + 0.5 * fluid_mass) * coriolis_particle);
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::AddRelativeAccelerationForces(array_1d<double,3>& weight, const ProcessInfo& r_process_info)
{
    const array_1d<double,3>& origin_acceleration = r_process_info[ACCELERATION_MOVING_FRAME_ORIGIN];
    noalias(weight) += (GetFluidMass() - GetMass()) * origin_acceleration;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::AddEulerForces(array_1d<double,3>& weight, const ProcessInfo& r_process_info)
{
    const array_1d<double,3>& omega_prime = r_process_info[ANGULAR_ACCELERATION_MOVING_FRAME];
    const array_1d<double,3>& coordinates = GetGeometry()[0].Coordinates();
    array_1d<double,3> euler;
    SWIMMING_SET_TO_CROSS_OF_FIRST_TWO_3(omega_prime, coordinates, euler)
    noalias(weight) += (GetFluidMass() - GetMass()) * euler;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::Initialize(const ProcessInfo& r_process_info)
{
    TBaseElement::Initialize(r_process_info);
    NodeType& node = GetGeometry()[0];
    mFirstStep                   = true;
    mHasDragForceNodalVar        = node.SolutionStepsDataHas(DRAG_FORCE);
    mHasHydroMomentNodalVar      = node.SolutionStepsDataHas(HYDRODYNAMIC_MOMENT);
    mHasVirtualMassForceNodalVar = node.SolutionStepsDataHas(VIRTUAL_MASS_FORCE);
    mHasBassetForceNodalVar      = node.SolutionStepsDataHas(BASSET_FORCE);
    mHasLiftForceNodalVar        = node.SolutionStepsDataHas(LIFT_FORCE);


    if (node.SolutionStepsDataHas(PARTICLE_SPHERICITY)){
        node.FastGetSolutionStepValue(PARTICLE_SPHERICITY) = this->GetProperties()[PARTICLE_SPHERICITY];
        mSphericity = node.FastGetSolutionStepValue(PARTICLE_SPHERICITY); //TODO: remove member var mSphericity from everywhere. Care with the occasions when PARTICLE_SPHERICITY is not added to the nodes!
    }
    else {
        mSphericity = 1.0;
    }

    mHasDragCoefficientVar       = node.SolutionStepsDataHas(DRAG_COEFFICIENT);
    mHasOldAdditionalForceVar    = node.SolutionStepsDataHas(ADDITIONAL_FORCE_OLD);
    mLastVirtualMassAddedMass    = 0.0;
    mLastBassetForceAddedMass    = 0.0;

}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SphericSwimmingParticle<TBaseElement>::MemberDeclarationFirstStep(const ProcessInfo& r_process_info)
{
    TBaseElement::MemberDeclarationFirstStep(r_process_info);
    mCouplingType                 = r_process_info[COUPLING_TYPE];
    mBuoyancyForceType            = r_process_info[BUOYANCY_FORCE_TYPE];
    mDragForceType                = r_process_info[DRAG_FORCE_TYPE];
    mVirtualMassForceType         = r_process_info[VIRTUAL_MASS_FORCE_TYPE];
    mBassetForceType              = r_process_info[BASSET_FORCE_TYPE];
    mSaffmanForceType             = r_process_info[LIFT_FORCE_TYPE];
    mMagnusForceType              = r_process_info[MAGNUS_FORCE_TYPE];
    mFluidModelType               = r_process_info[FLUID_MODEL_TYPE];
    mPorosityCorrectionType       = r_process_info[DRAG_POROSITY_CORRECTION_TYPE];
    mHydrodynamicTorqueType       = r_process_info[HYDRO_TORQUE_TYPE];
    mBrownianMotionType           = 0;
    mInitialTime                  = r_process_info[TIME];
    mQuadratureOrder              = r_process_info[QUADRATURE_ORDER];
    mOldBassetTerm                = ZeroVector(3);
    mOldDaitchePresentCoefficient = 0.0;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
//Explicit Instantiations
template class SphericSwimmingParticle<SphericParticle>;
template class SphericSwimmingParticle<NanoParticle>;
template class SphericSwimmingParticle<AnalyticSphericParticle>;

// Definition ( this probably neds to me moved to the .h file )
template <typename T> std::vector<double> SphericSwimmingParticle<T>::mAjs;
template <typename T> std::vector<double> SphericSwimmingParticle<T>::mBns;
template <typename T> std::vector<double> SphericSwimmingParticle<T>::mCns;
template <typename T> std::vector<double> SphericSwimmingParticle<T>::mDns;
template <typename T> std::vector<double> SphericSwimmingParticle<T>::mEns;
template <typename T> std::vector<double> SphericSwimmingParticle<T>::mAs;
template <typename T> std::vector<double> SphericSwimmingParticle<T>::mTs;
template <typename T> std::vector<double> SphericSwimmingParticle<T>::mAlphas;
template <typename T> std::vector<double> SphericSwimmingParticle<T>::mBetas;
template <typename T> double SphericSwimmingParticle<T>::mTimeWindow;
template <typename T> bool SphericSwimmingParticle<T>::mDaitcheVectorsAreFull;

// Instantiation
#define INSTANTIATE_SPHERIC_SWIMMING(_T)                            \
template std::vector<double> SphericSwimmingParticle<_T>::mAjs;     \
template std::vector<double> SphericSwimmingParticle<_T>::mBns;     \
template std::vector<double> SphericSwimmingParticle<_T>::mCns;     \
template std::vector<double> SphericSwimmingParticle<_T>::mDns;     \
template std::vector<double> SphericSwimmingParticle<_T>::mEns;     \
template std::vector<double> SphericSwimmingParticle<_T>::mAs;      \
template std::vector<double> SphericSwimmingParticle<_T>::mTs;      \
template std::vector<double> SphericSwimmingParticle<_T>::mAlphas;  \
template std::vector<double> SphericSwimmingParticle<_T>::mBetas;   \
template double SphericSwimmingParticle<_T>::mTimeWindow;           \
template bool SphericSwimmingParticle<_T>::mDaitcheVectorsAreFull;

INSTANTIATE_SPHERIC_SWIMMING(SphericParticle)
INSTANTIATE_SPHERIC_SWIMMING(NanoParticle)
INSTANTIATE_SPHERIC_SWIMMING(AnalyticSphericParticle)

#undef INSTANTIATE_SPHERIC_SWIMMING


}  // namespace Kratos.
