// Author: Guillermo Casas (gcasas@cimne.upc.edu)

#include "force_based_inlet.h"


namespace Kratos {

DEM_Force_Based_Inlet::DEM_Force_Based_Inlet(ModelPart& inlet_modelpart, array_1d<double, 3> injection_force) :
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

void DEM_Force_Based_Inlet::FixInjectionConditions(Element* p_element)
{
    Node<3>& node = p_element->GetGeometry()[0];
    node.FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE) = GetInjectionForce();
    KRATOS_WATCH(GetInjectionForce())

            //try
//    p_element->Set(NEW_ENTITY, 0);
//    node.Set(NEW_ENTITY, 0);
}

array_1d<double, 3> DEM_Force_Based_Inlet::GetInjectionForce()
{
    return mInjectionForce;
}

} // namespace Kratos
