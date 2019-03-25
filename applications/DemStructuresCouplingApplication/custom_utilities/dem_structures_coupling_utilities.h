/*
 * Author: Miguel Angel Celigueta
 *
 *  maceli@cimne.upc.edu
 */

#ifndef KRATOS_STRUCTURES_DEM_COUPLING_UTILITIES_H
#define KRATOS_STRUCTURES_DEM_COUPLING_UTILITIES_H
// /* External includes */

// System includes

// Project includes
#include "includes/variables.h"

/* System includes */
#include <limits>
#include <iostream>
#include <iomanip>
#include <fstream>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "../../DEMApplication/custom_conditions/RigidFace.h"
#include "../../DEMApplication/DEM_application_variables.h"


namespace Kratos {

    class DemStructuresCouplingUtilities {

    public:

    typedef ModelPart::NodesContainerType::ContainerType::iterator NodesIteratorType;

    KRATOS_CLASS_POINTER_DEFINITION(DemStructuresCouplingUtilities);

    DemStructuresCouplingUtilities(){}

    virtual ~DemStructuresCouplingUtilities(){}

    //***************************************************************************************************************
    //***************************************************************************************************************

    void TransferStructuresSkinToDem(ModelPart& r_source_model_part, ModelPart& r_destination_model_part, Properties::Pointer props) {

        std::string error = CheckProvidedProperties(props);

        if (error != "all_ok") KRATOS_ERROR << "The Dem Walls ModelPart has no valid Properties. Missing " << error << " . Exiting." << std::endl;

        r_destination_model_part.Conditions().Sort();
        int id = 1;

        if (r_destination_model_part.Conditions().size()) id = (r_destination_model_part.ConditionsEnd()-1)->Id() + 1;

        ModelPart::ConditionsContainerType& source_conditions = r_source_model_part.Conditions();

        // Adding conditions
        for (unsigned int i = 0; i < source_conditions.size(); i++) {
            ModelPart::ConditionsContainerType::iterator it = r_source_model_part.ConditionsBegin() + i;
            Geometry< Node<3> >::Pointer p_geometry =  it->pGetGeometry();
            Condition::Pointer cond = Condition::Pointer(new RigidFace3D(id, p_geometry, props));
            cond->Set(DEMFlags::STICKY, true);
            r_destination_model_part.AddCondition(cond); //TODO: add all of them in a single sentence! AddConditions. Use a temporary PointerVector as a list (not std::vector!).
            id++;
        }

        // Adding nodes
        r_destination_model_part.AddNodes(r_source_model_part.NodesBegin(), r_source_model_part.NodesEnd());
    }

    std::string CheckProvidedProperties(Properties::Pointer props) {
        std::vector<Variable<double> > list_of_variables_double_to_check = {FRICTION, WALL_COHESION, SEVERITY_OF_WEAR, IMPACT_WEAR_SEVERITY, BRINELL_HARDNESS, YOUNG_MODULUS, POISSON_RATIO};
        std::vector<Variable<bool> > list_of_variables_bool_to_check = {COMPUTE_WEAR};
        for (int i=0; i<(int)list_of_variables_double_to_check.size(); i++) {
            if(!props->Has(list_of_variables_double_to_check[i])) return list_of_variables_double_to_check[i].Name();
        }
        for (int i=0; i<(int)list_of_variables_bool_to_check.size(); i++) {
            if(!props->Has(list_of_variables_bool_to_check[i])) return list_of_variables_bool_to_check[i].Name();
        }
        return "all_ok";
    }

    void ComputeSandProduction(ModelPart& dem_model_part, ModelPart& skin_model_part) {

        const std::string filename = "sand_production_graph.txt";
        std::ifstream ifile(filename.c_str());
        static bool first_time_entered = true;
        if ((bool) ifile && first_time_entered) {
            std::remove("sand_production_graph.txt");
            first_time_entered = false;
        }

        ModelPart::ElementsContainerType& pElements = dem_model_part.GetCommunicator().LocalMesh().Elements();
        double current_total_mass_in_grams = 0.0;

        for (unsigned int k = 0; k < pElements.size(); k++) {

            ModelPart::ElementsContainerType::iterator it = pElements.ptr_begin() + k;
            Element* raw_p_element = &(*it);
            SphericParticle* p_sphere = dynamic_cast<SphericParticle*>(raw_p_element);
            const double particle_radius = p_sphere->GetRadius();
            const double particle_density = p_sphere->GetDensity();
            current_total_mass_in_grams += (4.0/3.0) * Globals::Pi * particle_density * particle_radius * particle_radius * particle_radius * 1000.0;
        }
        static const double initial_total_mass_in_grams = current_total_mass_in_grams;
        const double cumulative_sand_mass_in_grams = initial_total_mass_in_grams - current_total_mass_in_grams;

        ModelPart::NodesContainerType::iterator node_begin = skin_model_part.NodesBegin();
        const double face_pressure_in_psi = node_begin->FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE) * 0.000145;

        static std::ofstream sand_prod_file("sand_production_graph.txt", std::ios_base::out | std::ios_base::app);
        sand_prod_file << face_pressure_in_psi << " " << cumulative_sand_mass_in_grams << '\n';
    }

    virtual std::string Info() const
    {
        return "";
    }

    /// Print information about this object.

    virtual void PrintInfo(std::ostream& rOStream) const
    {
    }

    /// Print object's data.

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    protected:

    private:

    DemStructuresCouplingUtilities & operator=(DemStructuresCouplingUtilities const& rOther);

}; // Class DemStructuresCouplingUtilities

}  // namespace Python.

#endif // KRATOS_STRUCTURES_DEM_COUPLING_UTILITIES_H
