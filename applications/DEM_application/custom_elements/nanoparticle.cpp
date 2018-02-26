// Authors: Miguel Angel Celigueta, maceli@cimne.upc.edu
//          Guillermo Casas, gcasas@cimne.upc.edu

// System includes
#include <string>
#include <iostream>
#include <iomanip> // to improve std::cout precision

// Project includes
#include "nanoparticle.h"

namespace Kratos {

NanoParticle::~NanoParticle(){}

NanoParticle& NanoParticle::operator=(NanoParticle const& rOther) {
    
    SphericParticle::operator=(rOther);
        
    mThicknessOverRadius = rOther.mThicknessOverRadius;
    mInteractionRadius = rOther.mInteractionRadius;

    return *this;
}

void NanoParticle::Initialize(const ProcessInfo& r_process_info) {   
    SphericParticle::Initialize(r_process_info);
    double added_mass_coefficient = 1.0;
    SetMass(added_mass_coefficient * GetDensity() * CalculateVolume());
    this->SetInteractionRadius(2.5 * GetRadius());
    this->SetSearchRadius(3 * GetRadius());
}

void NanoParticle::ComputeAdditionalForces(array_1d<double, 3>& additionally_applied_force,
                             array_1d<double, 3>& additionally_applied_moment,
                             const ProcessInfo& r_current_process_info,
                             const array_1d<double,3>& gravity)
{
    KRATOS_TRY

    array_1d<double, 3> brownian_motion_force; brownian_motion_force.clear();
    array_1d<double, 3> van_der_waals_force; van_der_waals_force.clear();
    array_1d<double, 3> double_layer_force; double_layer_force.clear();
    noalias(additionally_applied_force) += brownian_motion_force + van_der_waals_force + double_layer_force;

    //Now add the contribution of base class function (gravity or other forces added in upper levels):
    SphericParticle::ComputeAdditionalForces(additionally_applied_force, additionally_applied_moment, r_current_process_info, gravity);
    KRATOS_CATCH( "" )
}

void NanoParticle::MemberDeclarationFirstStep(const ProcessInfo& r_process_info) {
    SphericParticle::MemberDeclarationFirstStep(r_process_info);
}

double NanoParticle::GetInteractionRadius(const int radius_index)
{
    return mInteractionRadius;
}

void NanoParticle::SetInteractionRadius(const double radius, const int radius_index)
{
    assert(radius >= GetRadius());
    mInteractionRadius = radius;
}

void NanoParticle::SetDefaultRadiiHierarchy(const double radius)
{
    SetRadius(radius);
    SetInteractionRadius(2.5 * radius);
    SetSearchRadius(3 * radius); // overwriting that established by the strategy
}

double NanoParticle::CalculateVolume()
{
    const double radius = this->GetRadius();
    return Globals::Pi * radius * radius * radius * mThicknessOverRadius;
}

    
} // namespace Kratos
