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
//                   John van Esch
//                   Gennady Markelov
//

#include "custom_elements/transient_thermal_element.h"

namespace Kratos {

template <unsigned int TDim, unsigned int TNumNodes>
int TransientThermalElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    CheckDomainSize();
    CheckHasSolutionStepsDataFor(TEMPERATURE);
    CheckHasSolutionStepsDataFor(DT_TEMPERATURE);
    CheckHasDofsFor(TEMPERATURE);

    VerifyProperty(DENSITY_WATER);
    VerifyProperty(POROSITY);
    VerifyProperty(SATURATION);
    VerifyProperty(DENSITY_SOLID);
    VerifyProperty(SPECIFIC_HEAT_CAPACITY_WATER);
    VerifyProperty(SPECIFIC_HEAT_CAPACITY_SOLID);
    VerifyProperty(THERMAL_CONDUCTIVITY_WATER);
    VerifyProperty(THERMAL_CONDUCTIVITY_SOLID_XX);
    VerifyProperty(THERMAL_CONDUCTIVITY_SOLID_YY);
    VerifyProperty(THERMAL_CONDUCTIVITY_SOLID_XY);
    VerifyProperty(LONGITUDINAL_DISPERSIVITY);
    VerifyProperty(TRANSVERSE_DISPERSIVITY);
    VerifyProperty(SOLID_COMPRESSIBILITY);

    const GeometryType& rGeom = GetGeometry();

    if (TDim == 2) {
        auto pos = std::find_if(rGeom.begin(), rGeom.end(),
                                [](const auto& node) { return node.Z() != 0.0; });
        if (pos != rGeom.end()) {
            KRATOS_ERROR << " Node with non-zero Z coordinate found. Id: " << pos->Id()
                         << std::endl;
        }
    }

    if (TDim > 2) {
        VerifyProperty(THERMAL_CONDUCTIVITY_SOLID_ZZ);
        VerifyProperty(THERMAL_CONDUCTIVITY_SOLID_YZ);
        VerifyProperty(THERMAL_CONDUCTIVITY_SOLID_XZ);
    }

    KRATOS_CATCH("")

    return 0;
}

template class TransientThermalElement<2, 3>;
template class TransientThermalElement<2, 4>;
template class TransientThermalElement<2, 6>;
template class TransientThermalElement<2, 8>;
template class TransientThermalElement<2, 9>;
template class TransientThermalElement<2, 10>;
template class TransientThermalElement<2, 15>;
template class TransientThermalElement<3, 4>;
template class TransientThermalElement<3, 8>;
template class TransientThermalElement<3, 10>;
template class TransientThermalElement<3, 20>;
template class TransientThermalElement<3, 27>;

} // Namespace Kratos
