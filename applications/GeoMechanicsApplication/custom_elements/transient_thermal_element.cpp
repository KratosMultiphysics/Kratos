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


    //----------------------------------------------------------------------------------------
    template< unsigned int TDim, unsigned int TNumNodes >
    int TransientThermalElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY
        //KRATOS_INFO("0-TransientThermalElement::Check()") << this->Id() << std::endl;
        //
        const PropertiesType& Prop = this->GetProperties();
        const GeometryType& Geom = this->GetGeometry();
        //
        if (Geom.DomainSize() < 1.0e-15)
            KRATOS_ERROR << "DomainSize < 1.0e-15 for the element " << this->Id() << std::endl;
        //
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            if (Geom[i].SolutionStepsDataHas(TEMPERATURE) == false)
                KRATOS_ERROR << "missing variable TEMPERATURE on node " << Geom[i].Id() << std::endl;
            //
            if (Geom[i].HasDofFor(TEMPERATURE) == false)
                KRATOS_ERROR << "missing variable TEMPERATURE on node " << Geom[i].Id() << std::endl;
        }
        //
        // Verify properties
        if (Prop.Has(DENSITY_WATER) == false || Prop[DENSITY_WATER] < 0.0)
            KRATOS_ERROR << "DENSITY_WATER does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(POROSITY) == false || Prop[POROSITY] < 0.0 || Prop[POROSITY] > 1.0)
            KRATOS_ERROR << "POROSITY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(SATURATION) == false || Prop[SATURATION] < 0.0 || Prop[SATURATION] > 1.0)
            KRATOS_ERROR << "SATURATION does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(DENSITY_SOLID) == false || Prop[DENSITY_SOLID] < 0.0)
            KRATOS_ERROR << "DENSITY_SOLID does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(HEAT_CAPACITY_WATER) == false || Prop[HEAT_CAPACITY_WATER] < 0.0)
            KRATOS_ERROR << "HEAT_CAPACITY_WATER does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(HEAT_CAPACITY_SOLID) == false || Prop[HEAT_CAPACITY_SOLID] < 0.0)
            KRATOS_ERROR << "HEAT_CAPACITY_SOLID does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(THERMAL_CONDUCTIVITY_WATER) == false || Prop[THERMAL_CONDUCTIVITY_WATER] < 0.0)
            KRATOS_ERROR << "THERMAL_CONDUCTIVITY_WATER does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(THERMAL_CONDUCTIVITY_SOLID_XX) == false || Prop[THERMAL_CONDUCTIVITY_SOLID_XX] < 0.0)
            KRATOS_ERROR << "THERMAL_CONDUCTIVITY_SOLID_XX does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(THERMAL_CONDUCTIVITY_SOLID_YY) == false || Prop[THERMAL_CONDUCTIVITY_SOLID_YY] < 0.0)
            KRATOS_ERROR << "THERMAL_CONDUCTIVITY_SOLID_YY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(THERMAL_CONDUCTIVITY_SOLID_XY) == false || Prop[THERMAL_CONDUCTIVITY_SOLID_XY] < 0.0)
            KRATOS_ERROR << "THERMAL_CONDUCTIVITY_SOLID_XY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(THERMAL_CONDUCTIVITY_SOLID_YX) == false || Prop[THERMAL_CONDUCTIVITY_SOLID_YX] < 0.0)
            KRATOS_ERROR << "THERMAL_CONDUCTIVITY_SOLID_YX does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(LONGITUDINAL_DISPERSIVITY) == false || Prop[LONGITUDINAL_DISPERSIVITY] < 0.0)
            KRATOS_ERROR << "LONGITUDINAL_DISPERSIVITY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(TRANSVERSE_DISPERSIVITY) == false || Prop[TRANSVERSE_DISPERSIVITY] < 0.0)
            KRATOS_ERROR << "TRANSVERSE_DISPERSIVITY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(SOLID_COMPRESSIBILITY) == false || Prop[SOLID_COMPRESSIBILITY] < 0.0)
            KRATOS_ERROR << "SOLID_COMPRESSIBILITY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (TDim == 2) {
            // If this is a 2D problem, nodes must be in XY plane
            for (unsigned int i = 0; i < TNumNodes; ++i) {
                if (Geom[i].Z() != 0.0)
                    KRATOS_ERROR << " Node with non-zero Z coordinate found. Id: " << Geom[i].Id() << std::endl;
            }
        }
        //
        return 0;
        //
        KRATOS_CATCH("");
    }


    //----------------------------------------------------------------------------------------------------
    template class TransientThermalElement<2, 3>;
    template class TransientThermalElement<2, 4>;
    template class TransientThermalElement<3, 4>;
    template class TransientThermalElement<3, 8>;
    //
    template class TransientThermalElement<2, 6>;
    template class TransientThermalElement<2, 8>;
    template class TransientThermalElement<2, 9>;
    template class TransientThermalElement<3, 10>;
    template class TransientThermalElement<3, 20>;
    template class TransientThermalElement<3, 27>;

} // Namespace Kratos
