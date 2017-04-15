// Author: Guillermo Casas (gcasas@cimne.upc.edu)

#include "force_based_inlet.h"


namespace Kratos {

DEM_Force_Based_Inlet::DEM_Force_Based_Inlet(ModelPart& inlet_modelpart) :
                       DEM_Inlet(inlet_modelpart){}

void DEM_Force_Based_Inlet::FixInjectionConditions(Element* p_element)
{

}


} // namespace Kratos
