// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "custom_processes/total_structural_mass_process.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
void TotalStructuralMassProcess::Execute()
{
    KRATOS_TRY
    
    // We initialize the total mass
    double total_mass  = 0.0;
    
    const std::size_t dimension = mrThisModelPart.GetProcessInfo()[DOMAIN_SIZE];
    
    // Now we iterate over the elements to calculate the total mass
    ElementsArrayType& elements_array = mrThisModelPart.Elements();
    
    #pragma omp parallel for reduction(+:total_mass)
    for(int i = 0; i < static_cast<int>(elements_array.size()); ++i){
        const auto it_elem = elements_array.begin() + i;
        
        // We get the condition geometry
        const GeometryType& r_this_geometry = it_elem->GetGeometry();
        const std::size_t local_space_dimension = r_this_geometry.LocalSpaceDimension();
        
        // We get the values from the condition
        const Properties& this_properties = it_elem->GetProperties();
        const double density = this_properties[DENSITY];
        
        if (local_space_dimension == 1) { // BEAM CASE
            const double area = this_properties[CROSS_AREA];
            total_mass += density * area * r_this_geometry.Length();
        } else if (local_space_dimension == 2 && dimension == 3) { // SHELL-MEMBRANE
            const double thickness = this_properties[THICKNESS];
            total_mass += density * thickness * r_this_geometry.Area();
        } else { // SOLID
            const double thickness = (dimension == 2) ? (this_properties.Has(THICKNESS)) ? this_properties[THICKNESS] : 1.0 : 1.0;
            const double volume = (dimension == 2) ? r_this_geometry.Area() : r_this_geometry.Volume();
            total_mass += density * thickness * volume;
        }
    }
    
    std::cout << "The total mass of the system is :" << total_mass/1000.0 << " Tn\n" << "Check variable NODAL_MASS in the process info in order to acces to it in any moment" << std::endl;
    
    mrThisModelPart.GetProcessInfo()[NODAL_MASS] = total_mass;
    
    KRATOS_CATCH("")
} // class TotalStructuralMassProcess
} // namespace Kratos
