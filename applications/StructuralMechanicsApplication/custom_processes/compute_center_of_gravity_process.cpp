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
#include "custom_processes/compute_center_of_gravity_process.h"
#include "structural_mechanics_application_variables.h"
//#include "geometries/geometry.h"
namespace Kratos
{

/*
Hi Salman,
you can use functions from the class Geometry (kratos/includes/geometry.h)
such as Center() or Area() for your computations
I sugest to compute the mass of an element and use the coordinates of the center
for calculating the COG
*/

namespace {
double GetFromProperty(const Properties& rThisProperties, Variable<double>& rVariable)
{
    // The purpose of this function is to avoid silent allocation of memory in case
    // the requested variable does not exist in the Properties!
    if (rThisProperties.Has(rVariable))
        return rThisProperties[rVariable];
    else
        return 0.0;
}
}


void ComputeCenterOfGravityProcess::Execute()
{
    KRATOS_TRY

    // We initialize the total mass and working variables
    /*beam_holder, pointmass_holder etc are variables to hold intermediate calculations
    where we multiply mass of each element to cog of the same element*/
    double total_mass  = 0.0;
    array_1d<double, 3> cog;
    array_1d<double, 3> beam_elem_holder;
    array_1d<double, 3> pointmass_holder;
    array_1d<double, 3> solid_elem_holder;
    array_1d<double, 3> mass_coor_product;
    array_1d<double, 3> shell_elem_holder;
    array_1d<double, 3> orthotropic_shell_holder;

    //initialization the working variables
    for (int i_i=0;i_i<3;i_i++)
    {
        cog[i_i]=0.0;
        beam_elem_holder[i_i]=0.0;
        pointmass_holder[i_i]=0.0;
        solid_elem_holder[i_i]=0.0;
        shell_elem_holder[i_i]=0.0;
        mass_coor_product[i_i]=0.0;
        orthotropic_shell_holder[i_i]=0.0;
    }

    const std::size_t dimension = mrThisModelPart.GetProcessInfo()[DOMAIN_SIZE];

    // Now we iterate over the elements to calculate the total mass
    ElementsArrayType& elements_array = mrThisModelPart.Elements();

    //#pragma omp parallel for reduction(+:total_mass)
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


        if (local_space_dimension == 0) // POINT MASSES
        {
            if (it_elem->Has(NODAL_MASS))
            {
                total_mass += it_elem->GetValue(NODAL_MASS);
                Point center_pointmass=r_this_geometry.Center();
                pointmass_holder[0]+= (it_elem->GetValue(NODAL_MASS))*center_pointmass.X();
                pointmass_holder[1]+= (it_elem->GetValue(NODAL_MASS))*center_pointmass.Y();
                pointmass_holder[2]+= (it_elem->GetValue(NODAL_MASS))*center_pointmass.Z();
            }
        } else if (local_space_dimension == 1) { // BEAM CASE
            const double density = GetFromProperty(this_properties,DENSITY);
            const double area = GetFromProperty(this_properties,CROSS_AREA);
            const double individual_mass=density * area * r_this_geometry.Length();
            total_mass += individual_mass;
            Point center_beam_elem=r_this_geometry.Center();
            beam_elem_holder[0] += individual_mass * center_beam_elem.X();
            beam_elem_holder[1] += individual_mass * center_beam_elem.Y();
            beam_elem_holder[2] += individual_mass * center_beam_elem.Z();

        } else if (local_space_dimension == 2 && dimension == 3) { // SHELL-MEMBRANE

            const double area = r_this_geometry.Area();
            Point center_shell_elem = r_this_geometry.Center();
            //center will be same for both orthotropic or isotropic shell element
            if (this_properties.Has(SHELL_ORTHOTROPIC_LAYERS)) { // composite material
                double orthotropic_mass=0.0;  // will collect mass of all layers of one element
                const auto orthotropic_layers = this_properties[SHELL_ORTHOTROPIC_LAYERS];
                for (std::size_t i=0; i<orthotropic_layers.size1(); ++i)
                    orthotropic_mass += orthotropic_layers(i,0) * orthotropic_layers(i,2) * area; // thickness*density*area
                total_mass += orthotropic_mass;
                orthotropic_shell_holder[0] += orthotropic_mass * center_shell_elem.X();
                orthotropic_shell_holder[1] += orthotropic_mass * center_shell_elem.Y();
                orthotropic_shell_holder[2] += orthotropic_mass * center_shell_elem.Z();
            } else {
                const double thickness = GetFromProperty(this_properties,THICKNESS);
                const double density = GetFromProperty(this_properties,DENSITY);
                total_mass += density * thickness * area;
                shell_elem_holder[0] += density * thickness * area * center_shell_elem.X();
                shell_elem_holder[1] += density * thickness * area * center_shell_elem.Y();
                shell_elem_holder[2] += density * thickness * area * center_shell_elem.Z();
            }
        } else { // SOLID

            const double thickness = (dimension == 2) ? (this_properties.Has(THICKNESS)) ? this_properties[THICKNESS] : 1.0 : 1.0;
            const double volume = (dimension == 2) ? r_this_geometry.Area() : r_this_geometry.Volume();
            const double density = GetFromProperty(this_properties,DENSITY);
            total_mass += density * thickness * volume;
            Point center_solid_elem = r_this_geometry.Center();
            solid_elem_holder[0] += density * thickness * volume * center_solid_elem.X();
            solid_elem_holder[1] += density * thickness * volume * center_solid_elem.Y();
            solid_elem_holder[2] += density * thickness * volume * center_solid_elem.Z();

        }

        // We restore the current configuration
        for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node)
            noalias(r_this_geometry[i_node].Coordinates()) = current_coordinates[i_node];

    }

    mass_coor_product[0] = beam_elem_holder[0] + pointmass_holder[0] + solid_elem_holder[0] + shell_elem_holder[0] + orthotropic_shell_holder[0];
    mass_coor_product[1] = beam_elem_holder[1] + pointmass_holder[1] + solid_elem_holder[1] + shell_elem_holder[1] + orthotropic_shell_holder[1];
    mass_coor_product[2] = beam_elem_holder[2] + pointmass_holder[2] + solid_elem_holder[2] + shell_elem_holder[2] + orthotropic_shell_holder[2];
    cog[0] = (mass_coor_product[0]) / total_mass;
    cog[1] = (mass_coor_product[1]) / total_mass;
    cog[2] = (mass_coor_product[2]) / total_mass;

    std::stringstream info_stream;
    info_stream << "Center of Gravity of ModelPart \"" << mrThisModelPart.Name() << "\"";

    KRATOS_INFO(info_stream.str()) << cog << std::endl;
    KRATOS_INFO("Hint")  << "Check variable CENTER_OF_GRAVITY in the process info in "
                         << "order to access to it in any moment" << std::endl;



    mrThisModelPart.GetProcessInfo()[CENTER_OF_GRAVITY] = cog;

    KRATOS_CATCH("")
} // class ComputeCenterOfGravityProcess
} // namespace Kratos
