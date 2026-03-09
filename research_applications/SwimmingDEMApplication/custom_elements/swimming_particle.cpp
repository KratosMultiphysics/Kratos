//   Author: G.Casas, gcasas@cimne.upc.edu $

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "swimming_DEM_application.h"
#include "swimming_particle.h"
#include "custom_utilities/GeometryFunctions.h"

namespace Kratos
{

template < class TBaseElement >
SwimmingParticle<TBaseElement>& SwimmingParticle<TBaseElement>::operator=(const SwimmingParticle<TBaseElement>& rOther) {

    TBaseElement::operator=(rOther);

    mNeighbourNodes = rOther.mNeighbourNodes;
    mNeighbourNodesDistances = rOther.mNeighbourNodesDistances;
    mFirstStep = rOther.mFirstStep;
    mPorosityCorrectionType = rOther.mPorosityCorrectionType;
    mFluidDensity = rOther.mFluidDensity;
    mKinematicViscosity = rOther.mKinematicViscosity;
    mSphericity = rOther.mSphericity;
    mNormOfSlipVel = rOther.mNormOfSlipVel;
    noalias(mSlipVel) = rOther.mSlipVel;
    mHydrodynamicInteractionLaw = rOther.mHydrodynamicInteractionLaw->Clone();

    return *this;
}

template < class TBaseElement >
void SwimmingParticle<TBaseElement>::CreateHydrodynamicInteractionLaws(const ProcessInfo& r_process_info)
{
    mHydrodynamicInteractionLaw = this->GetProperties()[SDEM_HYDRODYNAMIC_INTERACTION_LAW_POINTER]->Clone();
}

template < class TBaseElement >
void SwimmingParticle<TBaseElement>::ComputeAdditionalForces(array_1d<double, 3>& non_contact_force,
                                                             array_1d<double, 3>& non_contact_moment,
                                                             const ProcessInfo& r_current_process_info,
                                                             const array_1d<double,3>& gravity)
{
    KRATOS_TRY
    NodeType& node = GetGeometry()[0];

    if (!r_current_process_info[COUPLING_TYPE] || node.IsNot(INSIDE) || node.Is(BLOCKED)){
        TBaseElement::ComputeAdditionalForces(non_contact_force, non_contact_moment, r_current_process_info, gravity);
        return;
    }

    mFluidDensity                           = node.FastGetSolutionStepValue(FLUID_DENSITY_PROJECTED);
    mKinematicViscosity                     = node.FastGetSolutionStepValue(FLUID_VISCOSITY_PROJECTED);
    const array_1d<double, 3>& fluid_vel    = node.FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
    const array_1d<double, 3>& particle_vel = node.FastGetSolutionStepValue(VELOCITY);


    noalias(mSlipVel) = fluid_vel - particle_vel;

    mNormOfSlipVel = SWIMMING_MODULUS_3(mSlipVel);
    array_1d<double, 3> weight = ZeroVector(3);
    array_1d<double, 3> buoyancy = ZeroVector(3);
    array_1d<double, 3> drag_force = ZeroVector(3);
    array_1d<double, 3> inviscid_force = ZeroVector(3);
    array_1d<double, 3> history_force = ZeroVector(3);
    array_1d<double, 3> vorticity_induced_lift = ZeroVector(3);
    array_1d<double, 3> rotation_induced_lift = ZeroVector(3);
    array_1d<double, 3> steady_viscous_torque = ZeroVector(3);
    Geometry<Node >& r_geometry = GetGeometry();

    // The decomposition of forces that is considered here follows Jackson (The Dynamics of Fluidized Particles, 2000);
    // so that the role of f_n1 therein is played by non_contact_force here

    TBaseElement::ComputeAdditionalForces(weight, non_contact_moment, r_current_process_info, gravity); // Could be adding something else apart from weight

    mHydrodynamicInteractionLaw->ComputeBuoyancyForce(r_geometry,
                                                      mFluidDensity,
                                                      CalculateVolume(),
                                                      gravity,
                                                      buoyancy,
                                                      r_current_process_info);

    mHydrodynamicInteractionLaw->ComputeDragForce(this,
                                                  mRadius,
                                                  mFluidDensity,
                                                  mKinematicViscosity,
                                                  mSlipVel,
                                                  drag_force,
                                                  r_current_process_info);

    mHydrodynamicInteractionLaw->ComputeInviscidForce(r_geometry,
                                                      mFluidDensity,
                                                      CalculateVolume(),
                                                      inviscid_force,
                                                      r_current_process_info);

    mHydrodynamicInteractionLaw->ComputeHistoryForce(r_geometry,
                                                     mRadius,
                                                     mFluidDensity,
                                                     mKinematicViscosity,
                                                     mSlipVel,
                                                     history_force,
                                                     r_current_process_info);

    mHydrodynamicInteractionLaw->ComputeVorticityInducedLift(r_geometry,
                                                             mRadius,
                                                             mFluidDensity,
                                                             mKinematicViscosity,
                                                             mSlipVel,
                                                             vorticity_induced_lift,
                                                             r_current_process_info);

    mHydrodynamicInteractionLaw->ComputeRotationInducedLift(r_geometry,
                                                            mRadius,
                                                            mFluidDensity,
                                                            mKinematicViscosity,
                                                            mSlipVel,
                                                            rotation_induced_lift,
                                                            r_current_process_info);

    mHydrodynamicInteractionLaw->ComputeSteadyViscousTorque(r_geometry,
                                                            mRadius,
                                                            mFluidDensity,
                                                            mKinematicViscosity,
                                                            mSlipVel,
                                                            steady_viscous_torque,
                                                            r_current_process_info);

    // Adding all fluid-related forces except Basset's, since they might get averaged in time in a different way
    noalias(non_contact_force) += weight
                                + buoyancy
                                + drag_force
                                + inviscid_force
                                + vorticity_induced_lift
                                + rotation_induced_lift;

    // Adding all fluid-related moments
    noalias(non_contact_moment) += steady_viscous_torque;

    const double inviscid_added_mass =  mHydrodynamicInteractionLaw->GetInviscidAddedMass(GetGeometry(),
                                                                                          mFluidDensity,
                                                                                          r_current_process_info);

    const double history_force_added_mass =  mHydrodynamicInteractionLaw->GetHistoryForceAddedMass(GetGeometry(),
                                                                                                   r_current_process_info);

    const double force_reduction_coeff = mRealMass / (mRealMass + inviscid_added_mass + history_force_added_mass);

    array_1d<double, 3> non_contact_force_not_altered = non_contact_force;

    if (node.SolutionStepsDataHas(ADDITIONAL_FORCE_OLD) && !mFirstStep){
        ApplyNumericalAveragingWithOldForces(node, non_contact_force, r_current_process_info);
    }

    UpdateNodalValues(node,
                      non_contact_force_not_altered,
                      non_contact_moment,
                      weight,
                      buoyancy,
                      drag_force,
                      inviscid_force,
                      history_force,
                      vorticity_induced_lift,
                      rotation_induced_lift,
                      force_reduction_coeff,
                      r_current_process_info);
    // The Basset force has a different temporal treatment, so first we apply the scheme to the rest of the forces
    // and then we add the Basset force (minus the term proportional to the current acceleration, which is treated implicitly)
    noalias(non_contact_force) += history_force;
    non_contact_force *= force_reduction_coeff; //TODO: put noalias here?
    mFirstStep = false;

    KRATOS_CATCH( "" )
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
// Here nodal values are modified to record DEM forces that we want to print. In Kratos this is an exception since nodal values are meant to be modified only outside the element. Here it was not possible.
template < class TBaseElement >
void SwimmingParticle<TBaseElement>::UpdateNodalValues(NodeType& node,
                                                       const array_1d<double, 3>& non_contact_nor_history_force,
                                                       const array_1d<double, 3>& non_contact_moment,
                                                       const array_1d<double, 3>& weight,
                                                       const array_1d<double, 3>& buoyancy,
                                                       const array_1d<double, 3>& drag_force,
                                                       const array_1d<double, 3>& inviscid_force,
                                                       const array_1d<double, 3>& history_force,
                                                       const array_1d<double, 3>& vorticity_induced_lift,
                                                       const array_1d<double, 3>& rotation_induced_lift,
                                                       const double& force_reduction_coeff,
                                                       const ProcessInfo& r_current_process_info)
{
    noalias(node.FastGetSolutionStepValue(HYDRODYNAMIC_FORCE))       = force_reduction_coeff * (non_contact_nor_history_force + history_force - buoyancy - weight);
    noalias(node.FastGetSolutionStepValue(BUOYANCY))                 = buoyancy;
    array_1d<double, 3>& total_force = node.FastGetSolutionStepValue(TOTAL_FORCES);
    total_force *= force_reduction_coeff;

    if (node.SolutionStepsDataHas(HYDRODYNAMIC_MOMENT)){
        noalias(node.FastGetSolutionStepValue(HYDRODYNAMIC_MOMENT))  = non_contact_moment;
    }

    if (node.SolutionStepsDataHas(DRAG_FORCE)){
        noalias(node.FastGetSolutionStepValue(DRAG_FORCE))           = drag_force;
    }

    if (node.SolutionStepsDataHas(VIRTUAL_MASS_FORCE)){ // This only includes the part proportional to the fluid acceleration (undisturbed flow plus added mass terms), since the particle acceleration is treated implicitly. It is added later by the strategy, which calls Calculate here, once the current acceleration is available
        noalias(node.FastGetSolutionStepValue(VIRTUAL_MASS_FORCE))   = inviscid_force;
    }

    if (node.SolutionStepsDataHas(BASSET_FORCE)){
        noalias(node.FastGetSolutionStepValue(BASSET_FORCE))         = history_force; // This does not include the current-time contribution of the acceleration, which is treated implicitly. It is added later by the strategy, which calls Calculate here, once the current acceleration is available
    }

    if (node.SolutionStepsDataHas(ADDITIONAL_FORCE_OLD)){
        noalias(node.FastGetSolutionStepValue(ADDITIONAL_FORCE_OLD)) = non_contact_nor_history_force;
    }

    if (node.SolutionStepsDataHas(LIFT_FORCE)){
        noalias(node.FastGetSolutionStepValue(LIFT_FORCE))           = vorticity_induced_lift + rotation_induced_lift;
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SwimmingParticle<TBaseElement>::ApplyNumericalAveragingWithOldForces(NodeType& node, array_1d<double, 3>& non_contact_nor_history_force, const ProcessInfo& r_current_process_info)
{
    noalias(non_contact_nor_history_force) = 0.5 * (3 * non_contact_nor_history_force - node.FastGetSolutionStepValue(ADDITIONAL_FORCE_OLD));
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SwimmingParticle<TBaseElement>::AdditionalCalculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_current_process_info)
{
    if (rVariable == REYNOLDS_NUMBER){
        NodeType& node = GetGeometry()[0];

        if (node.IsNot(INSIDE)){
            Output = 0.0;
        }

        else {
            mFluidDensity                           = node.FastGetSolutionStepValue(FLUID_DENSITY_PROJECTED);
            mKinematicViscosity                     = node.FastGetSolutionStepValue(FLUID_VISCOSITY_PROJECTED);
            const array_1d<double, 3>& fluid_vel    = node.FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
            const array_1d<double, 3>& particle_vel = node.FastGetSolutionStepValue(VELOCITY);
            noalias(mSlipVel) = fluid_vel - particle_vel;
            mNormOfSlipVel = SWIMMING_MODULUS_3(mSlipVel);
            mHydrodynamicInteractionLaw->ComputeParticleReynoldsNumber(mRadius,
                                                                       mKinematicViscosity,
                                                                       mNormOfSlipVel);
        }
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SwimmingParticle<TBaseElement>::Calculate(const Variable<array_1d<double, 3 > >& rVariable,
                                               array_1d<double, 3 > & Output,
                                               const ProcessInfo& r_current_process_info)
{
    if (rVariable == VIRTUAL_MASS_FORCE) {
        const array_1d<double, 3 > total_forces = GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES);
        const double virtual_mass_coeff = mHydrodynamicInteractionLaw->GetInviscidAddedMass(GetGeometry(),
                                                                                            mFluidDensity,
                                                                                            r_current_process_info);
        Output -= virtual_mass_coeff / mRealMass * total_forces;
    }

    else if (rVariable == BASSET_FORCE) {
        const array_1d<double, 3 > total_forces = GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES);
        const double history_force_added_mass =  mHydrodynamicInteractionLaw->GetHistoryForceAddedMass(GetGeometry(),
                                                                                                       r_current_process_info);
        Output -= history_force_added_mass / mRealMass * total_forces;
    }

    else {
        TBaseElement::Calculate(rVariable, Output, r_current_process_info);
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
double SwimmingParticle<TBaseElement>::CalculateDragCoeffFromSphericity(const double reynolds,
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
template < class TBaseElement >
void SwimmingParticle<TBaseElement>::ApplyDragPorosityModification(double& drag_coeff)
{
    if (mPorosityCorrectionType == 0){
        return;
    }

    else if (mPorosityCorrectionType == 1){ // Richardson and Zaki, 1954 (fluid fraction ~ 0.01 - 0.2)
        const double reynolds = mHydrodynamicInteractionLaw->ComputeParticleReynoldsNumber(mRadius,
                                                                                           mKinematicViscosity,
                                                                                           mNormOfSlipVel);
        double K;

        if (reynolds > 500){
            K = 2.39;
        }

        else if (reynolds > 1){
            K = 4.45 * std::pow(reynolds, - 0.1);
        }

        else if (reynolds > 0.2){
            K = 4.35 * std::pow(reynolds, - 0.03);
        }

        else {
            K = 4.65;
        }
        const double fluid_volume_fraction = GetGeometry()[0].FastGetSolutionStepValue(FLUID_FRACTION_PROJECTED);
        drag_coeff *= std::pow(fluid_volume_fraction, 1 - 2 * K);
    }

}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
double SwimmingParticle<TBaseElement>::GetFluidMass()
{
    return mFluidDensity * CalculateVolume();
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
array_1d<double,3> SwimmingParticle<TBaseElement>::ComputeWeight(const array_1d<double,3>& gravity, const ProcessInfo& r_process_info)
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
void SwimmingParticle<TBaseElement>::AddCentrifugalForces(array_1d<double,3>& weight, const ProcessInfo& r_process_info)
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
void SwimmingParticle<TBaseElement>::AddCoriolisForces(array_1d<double,3>& weight, const ProcessInfo& r_process_info)
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
void SwimmingParticle<TBaseElement>::AddRelativeAccelerationForces(array_1d<double,3>& weight, const ProcessInfo& r_process_info)
{
    const array_1d<double,3>& origin_acceleration = r_process_info[ACCELERATION_MOVING_FRAME_ORIGIN];
    noalias(weight) += (GetFluidMass() - GetMass()) * origin_acceleration;
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template< class TBaseElement >
void SwimmingParticle<TBaseElement>::AddEulerForces(array_1d<double,3>& weight, const ProcessInfo& r_process_info)
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
void SwimmingParticle<TBaseElement>::Initialize(const ProcessInfo& r_process_info)
{
    TBaseElement::Initialize(r_process_info);
    CreateHydrodynamicInteractionLaws(r_process_info);

    NodeType& node = GetGeometry()[0];
    mFirstStep = true;

    if (node.SolutionStepsDataHas(PARTICLE_SPHERICITY)){
        node.FastGetSolutionStepValue(PARTICLE_SPHERICITY) = this->GetProperties()[PARTICLE_SPHERICITY];
        mSphericity = node.FastGetSolutionStepValue(PARTICLE_SPHERICITY); //TODO: remove member var mSphericity from everywhere. Care with the occasions when PARTICLE_SPHERICITY is not added to the nodes!
    }
    else {
        mSphericity = 1.0;
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
template < class TBaseElement >
void SwimmingParticle<TBaseElement>::MemberDeclarationFirstStep(const ProcessInfo& r_process_info)
{
    TBaseElement::MemberDeclarationFirstStep(r_process_info);
    mPorosityCorrectionType       = r_process_info[DRAG_POROSITY_CORRECTION_TYPE];
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
//Explicit Instantiations
template class SwimmingParticle<SphericParticle>;
template class SwimmingParticle<NanoParticle>;
template class SwimmingParticle<AnalyticSphericParticle>;

}  // namespace Kratos.
