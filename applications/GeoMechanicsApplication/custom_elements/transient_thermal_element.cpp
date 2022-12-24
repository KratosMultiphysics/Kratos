// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//

// Application includes
#include "custom_elements/transient_thermal_element.hpp"

namespace Kratos
{
    // ============================================================================================
	// ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    Element::Pointer TransientThermalElement<TDim, TNumNodes>::Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new TransientThermalElement(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    Element::Pointer TransientThermalElement<TDim, TNumNodes>::Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new TransientThermalElement(NewId, pGeom, pProperties));
    }



    // ============================================================================================
    // ============================================================================================
    template class TransientThermalElement<2,3>;
    template class TransientThermalElement<2,4>;
    template class TransientThermalElement<3,4>;
    template class TransientThermalElement<3,8>;
    //
    template class TransientThermalElement<2,6>;
    template class TransientThermalElement<2,8>;
    template class TransientThermalElement<2,9>;
    template class TransientThermalElement<3,10>;
    template class TransientThermalElement<3,20>;
    template class TransientThermalElement<3,27>;

} // Namespace Kratos
