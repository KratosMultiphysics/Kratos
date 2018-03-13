// Author: Guillermo Casas (gcasas@cimne.upc.edu)

#include "../../../DEM_application/custom_constitutive/DEM_D_Bentonite_Colloid_CL.h"
#include "../../../DEM_application/custom_elements/spheric_particle.h"
#include "bentonite_force_based_inlet.h"

namespace Kratos {

Bentonite_Force_Based_Inlet::Bentonite_Force_Based_Inlet(ModelPart& inlet_modelpart, array_1d<double, 3> injection_force):
    BaseClass(inlet_modelpart, injection_force)
{
    // we want to keep only the direction of the input force (mInjectionForce), since its modulus
    // has to adapt to the prevailing conditions.
    mCationConcentration = inlet_modelpart[CATION_CONCENTRATION];
    const double mod_force = DEM_MODULUS_3(injection_force);
    assert(mod_force > 0.0);
    injection_force = injection_force / mod_force;
    noalias(mInjectionForce) = injection_force;
}

void Bentonite_Force_Based_Inlet::InitializeStep(ModelPart& r_receiver_model_part)
{
    mCationConcentration = GetInletModelPart()[CATION_CONCENTRATION];
    ModelPart& mp = GetInletModelPart();

    for (ElementIteratorType elem_it = mp.ElementsBegin(); elem_it != mp.ElementsEnd(); ++elem_it){
        UpdateInjectionForce((*(elem_it.base())).get());
    }

    for (ElementIteratorType elem_it = r_receiver_model_part.ElementsBegin(); elem_it != r_receiver_model_part.ElementsEnd(); ++elem_it){
        if (elem_it->Is(NEW_ENTITY)){
            UpdateInjectionForce((*(elem_it.base())).get());
        }
    }
}

void Bentonite_Force_Based_Inlet::FixInjectionConditions(Element* p_element, Element* p_injector_element)
{
    UpdateInjectionForce(p_element);
}

void Bentonite_Force_Based_Inlet::UpdateInjectionForce(Element* p_element)
{
    Node<3>& node = p_element->GetGeometry()[0];
    node.FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE) = GetInjectionForce(p_element);
}

void Bentonite_Force_Based_Inlet::FixInjectorConditions(Element* p_element)
{    //AddRandomPerpendicularComponentToGivenVector(mInjectionForce, 60); // the max angle should be an INPUT

    Node<3>& node = p_element->GetGeometry()[0];
    node.FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE) = GetInjectionForce(p_element);
}

array_1d<double, 3> Bentonite_Force_Based_Inlet::GetInjectionForce(Element* p_element)
{
    DEM_D_Bentonite_Colloid* p_law = dynamic_cast<DEM_D_Bentonite_Colloid*>(dynamic_cast<SphericParticle*>(p_element)->GetConstitutiveLawPointer().get());
    const double normal_force_modulus = fabs(p_law->CalculateNormalForce(1e-7, mCationConcentration));
    array_1d<double, 3 > unitary_vector = mInjectionForce;
    //AddRandomPerpendicularComponentToGivenVector(unitary_vector, 60); // the max angle should be an INPUT
    return normal_force_modulus * unitary_vector;
}

} // namespace Kratos
