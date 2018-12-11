// Author: Guillermo Casas (gcasas@cimne.upc.edu)

#include "force_based_inlet.h"
#include "custom_elements/spheric_particle.h"


namespace Kratos {

DEM_Force_Based_Inlet::DEM_Force_Based_Inlet(ModelPart& inlet_modelpart, array_1d<double, 3> injection_force):
               DEM_Inlet(inlet_modelpart), mInjectionForce(injection_force)
{}

void DEM_Force_Based_Inlet::RemoveInjectionConditions(Element &element)
{
    Node<3>& node = element.GetGeometry()[0];
    element.Set(NEW_ENTITY, 0);
    node.Set(NEW_ENTITY, 0);
    node.pGetDof(VELOCITY_X)->FreeDof();
    node.pGetDof(VELOCITY_Y)->FreeDof();
    node.pGetDof(VELOCITY_Z)->FreeDof();
    node.pGetDof(ANGULAR_VELOCITY_X)->FreeDof();
    node.pGetDof(ANGULAR_VELOCITY_Y)->FreeDof();
    node.pGetDof(ANGULAR_VELOCITY_Z)->FreeDof();
    noalias(node.FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE)) = ZeroVector();
}

void DEM_Force_Based_Inlet::FixInjectionConditions(Element* p_element, Element* p_injector_element)
{
    // the injector velocity is not relevant here
    static_cast<void>(p_injector_element);

    Node<3>& node = p_element->GetGeometry()[0];
    node.FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE) = GetInjectionForce(p_element);

//    SphericParticle* p_spheric_particle = dynamic_cast<SphericParticle*>(p_element);
//    p_spheric_particle->SetInteractionRadius(p_spheric_particle->GetRadius());
}

void DEM_Force_Based_Inlet::FixInjectorConditions(Element* p_element)
{
    Node<3>& node = p_element->GetGeometry()[0];
    node.FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE) = GetInjectionForce(p_element);
}

array_1d<double, 3> DEM_Force_Based_Inlet::GetInjectionForce(Element* p_element)
{
    return mInjectionForce;
}

} // namespace Kratos
