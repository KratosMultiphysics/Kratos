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

double GetFromProperty(const Properties& rThisProperties, Variable<double>& rVariable)
{
    // The purpose of this function is to avoid silent allocation of memory in case
    // the requested variable does not exist in the Properties!
    if (rThisProperties.Has(rVariable))
        return rThisProperties[rVariable];
    else
        return 0.0;
}

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
        GeometryType& r_this_geometry = it_elem->GetGeometry();
        const std::size_t local_space_dimension = r_this_geometry.LocalSpaceDimension();
        const std::size_t number_of_nodes = r_this_geometry.size();

        // We copy the current coordinates and move the coordinates to the initial configuration
        std::vector<array_1d<double, 3>> current_coordinates(number_of_nodes);
        for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node) {
            noalias(current_coordinates[i_node]) = r_this_geometry[i_node].Coordinates();
            noalias(r_this_geometry[i_node].Coordinates()) = r_this_geometry[i_node].GetInitialPosition().Coordinates();
        }

        // We get the values from the condition
        const Properties& this_properties = it_elem->GetProperties();

        if (local_space_dimension == 0) { // POINT MASSES
            if (it_elem->Has(NODAL_MASS))
                total_mass += it_elem->GetValue(NODAL_MASS);
        } else if (local_space_dimension == 1) { // BEAM CASE
            const double density = GetFromProperty(this_properties,DENSITY);
            const double area = GetFromProperty(this_properties,CROSS_AREA);
            total_mass += density * area * r_this_geometry.Length();
        } else if (local_space_dimension == 2 && dimension == 3) { // SHELL-MEMBRANE
            const double area = r_this_geometry.Area();
            if (this_properties.Has(SHELL_ORTHOTROPIC_LAYERS)) { // composite material
                const auto orthotropic_layers = this_properties[SHELL_ORTHOTROPIC_LAYERS];
                for (std::size_t i=0; i<orthotropic_layers.size1(); ++i)
                    total_mass += orthotropic_layers(i,0) * orthotropic_layers(i,2) * area; // thickness*density*area
            } else {
                const double thickness = GetFromProperty(this_properties,THICKNESS);
                const double density = GetFromProperty(this_properties,DENSITY);
                total_mass += density * thickness * area;
            }
        } else { // SOLID
            const double thickness = (dimension == 2) ? (this_properties.Has(THICKNESS)) ? this_properties[THICKNESS] : 1.0 : 1.0;
            const double volume = (dimension == 2) ? r_this_geometry.Area() : r_this_geometry.Volume();
            const double density = GetFromProperty(this_properties,DENSITY);
            total_mass += density * thickness * volume;
        }

        // We restore the current configuration
        for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node)
            noalias(r_this_geometry[i_node].Coordinates()) = current_coordinates[i_node];
    }

    std::stringstream info_stream;
    info_stream << "Total Mass of ModelPart \"" << mrThisModelPart.Name() << "\"";

    KRATOS_INFO(info_stream.str()) << total_mass << std::endl;
    KRATOS_INFO("Hint")  << "Check variable NODAL_MASS in the process info in "
                         << "order to access to it in any moment" << std::endl;

    mrThisModelPart.GetProcessInfo()[NODAL_MASS] = total_mass;

    KRATOS_CATCH("")
} // class TotalStructuralMassProcess
} // namespace Kratos
