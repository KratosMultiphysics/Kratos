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

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_conditions/RigidFace.h"
#include "custom_conditions/RigidEdge.h"
#include "DEM_application_variables.h"
#include "dem_structures_coupling_application_variables.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos
{
class DemStructuresCouplingUtilities
{
public:
typedef ModelPart::NodesContainerType::ContainerType::iterator NodesIteratorType;

KRATOS_CLASS_POINTER_DEFINITION(DemStructuresCouplingUtilities);

/// Default constructor
DemStructuresCouplingUtilities(){}

/// Destructor
virtual ~DemStructuresCouplingUtilities(){}

//***************************************************************************************************************
//***************************************************************************************************************

void TransferStructuresSkinToDem(ModelPart& r_source_model_part, ModelPart& r_destination_model_part, Properties::Pointer props) {

    std::string error = CheckProvidedProperties(props);
    const int dimension = r_source_model_part.GetProcessInfo()[DOMAIN_SIZE];

    if (error != "all_ok") KRATOS_ERROR << "The Dem Walls ModelPart has no valid Properties. Missing " << error << " . Exiting." << std::endl;

    r_destination_model_part.Conditions().Sort();
    int id = 1;

    if (r_destination_model_part.Conditions().size()) id = (r_destination_model_part.ConditionsEnd()-1)->Id() + 1;

    ModelPart::ConditionsContainerType& source_conditions = r_source_model_part.Conditions();

    // Adding conditions
    for (unsigned int i = 0; i < source_conditions.size(); i++) {
        ModelPart::ConditionsContainerType::iterator it = r_source_model_part.ConditionsBegin() + i;
        Geometry< Node<3> >::Pointer p_geometry =  it->pGetGeometry();
        Condition::Pointer cond;
        if (dimension == 2) {
            cond = Condition::Pointer(new RigidEdge3D(id, p_geometry, props));            
        } else {
            cond = Condition::Pointer(new RigidFace3D(id, p_geometry, props));  
        }

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

void SmoothLoadTrasferredToFem(ModelPart& r_model_part, const double portion_of_the_force_which_is_new) {
    #pragma omp parallel for
    for (int i=0; i<(int)r_model_part.Nodes().size(); i++) {
        auto node_it = r_model_part.NodesBegin() + i;
        array_1d<double, 3> averaged_force;
        array_1d<double, 3>& node_dem_load = node_it->FastGetSolutionStepValue(DEM_SURFACE_LOAD);
        noalias(averaged_force) = portion_of_the_force_which_is_new * node_dem_load + (1.0 - portion_of_the_force_which_is_new) * node_it->FastGetSolutionStepValue(DEM_SURFACE_LOAD, 1);
        noalias(node_dem_load) = averaged_force;
    }
}

void ComputeSandProduction(ModelPart& dem_model_part, ModelPart& outer_walls_model_part, const double time) {

    const std::string sand_prod_filename = "sand_production_graph.txt";
    static std::ofstream ofs_sand_prod_file;
    static bool first_time_entered = true;
    if (first_time_entered) {
        ofs_sand_prod_file.open(sand_prod_filename, std::ofstream::out | std::ofstream::trunc);
        first_time_entered = false;
    }

    ModelPart::ElementsContainerType& pElements = dem_model_part.GetCommunicator().LocalMesh().Elements();
    double current_total_mass_in_grams = 0.0;

    for (unsigned int k = 0; k < pElements.size(); k++) {

        ModelPart::ElementsContainerType::iterator it = pElements.ptr_begin() + k;
        Element* raw_p_element = &(*it);
        SphericParticle* p_sphere = dynamic_cast<SphericParticle*>(raw_p_element);
        if (p_sphere->Is(ISOLATED)) continue;
        const double particle_density = p_sphere->GetDensity();
        const double particle_volume = p_sphere->CalculateVolume();
        current_total_mass_in_grams += particle_volume * particle_density * 1.0e3;
    }
    static const double initial_total_mass_in_grams = current_total_mass_in_grams;
    const double cumulative_sand_mass_in_grams = initial_total_mass_in_grams - current_total_mass_in_grams;

    //ModelPart::ConditionsContainerType::iterator condition_begin = outer_walls_model_part.ConditionsBegin();
    //const double face_pressure_in_psi = condition_begin->GetValue(POSITIVE_FACE_PRESSURE) * 0.000145;
    ProcessInfo& r_process_info = dem_model_part.GetProcessInfo();
    const double Pascals_to_psi_factor = 0.000145;
    const double face_pressure_in_psi = fabs(r_process_info[TARGET_STRESS_Z]) * Pascals_to_psi_factor;

    static std::ofstream sand_prod_file("sand_production_graph.txt", std::ios_base::out | std::ios_base::app);
    sand_prod_file << time << " " << face_pressure_in_psi << " " << cumulative_sand_mass_in_grams << '\n';
    sand_prod_file.flush();
}

void MarkBrokenSpheres(ModelPart& dem_model_part) {

    ModelPart::ElementsContainerType& pElements = dem_model_part.GetCommunicator().LocalMesh().Elements();

    for (unsigned int k = 0; k < pElements.size(); k++) {

        ModelPart::ElementsContainerType::iterator it = pElements.ptr_begin() + k;
        Element* raw_p_element = &(*it);
        SphericContinuumParticle* p_sphere = dynamic_cast<SphericContinuumParticle*>(raw_p_element);
        if (p_sphere->Is(ISOLATED)) continue;
        bool go_to_next_particle = false;

        for (unsigned int i = 0; i < p_sphere->mContinuumInitialNeighborsSize; i++) {

            if (!p_sphere->mIniNeighbourFailureId[i]) {
                go_to_next_particle = true;
                break;
            }
        }
        if (go_to_next_particle) continue;
        else p_sphere->Set(ISOLATED, true);
    }
}

void ComputeSandProductionWithDepthFirstSearchNonRecursiveImplementation(ModelPart& dem_model_part, ModelPart& outer_walls_model_part, const double time) {

    const std::string sand_prod_filename = "sand_production_graph_with_chunks_non_recursive.txt";
    static std::ofstream ofs_sand_prod_file;
    const std::string granulometry_distr_filename = "granulometry_distribution.txt";
    static std::ofstream ofs_granulometry_distr_file;
    static bool first_time_entered = true;
    if (first_time_entered) {
        ofs_sand_prod_file.open(sand_prod_filename, std::ofstream::out | std::ofstream::trunc);
        ofs_granulometry_distr_file.open(granulometry_distr_filename, std::ofstream::out | std::ofstream::trunc);
        first_time_entered = false;
    }

    ModelPart::ElementsContainerType& pElements = dem_model_part.GetCommunicator().LocalMesh().Elements();

    std::vector<double> chunks_masses;

    for (unsigned int k = 0; k < pElements.size(); k++) {
        ModelPart::ElementsContainerType::iterator it = pElements.ptr_begin() + k;
        it->Set(VISITED, false);
    }

    std::vector<SphericContinuumParticle*> stack_of_particles_to_check;

    for (unsigned int k = 0; k < pElements.size(); k++) {
        ModelPart::ElementsContainerType::iterator it = pElements.ptr_begin() + k;
        Element* raw_p_element = &(*it);
        SphericContinuumParticle* p_sphere = dynamic_cast<SphericContinuumParticle*>(raw_p_element);
        double this_chunk_mass = 0.0;
        stack_of_particles_to_check.push_back(p_sphere);
        while (stack_of_particles_to_check.size()) {
            SphericContinuumParticle* current_particle = stack_of_particles_to_check.back();
            stack_of_particles_to_check.pop_back();
            if (current_particle->Is(VISITED)) continue;
            const double particle_density = current_particle->GetDensity();
            const double particle_volume = current_particle->CalculateVolume();
            this_chunk_mass += particle_volume * particle_density * 1.0e3;

            current_particle->Set(VISITED, true);

            for (size_t i = 0; i < current_particle->mContinuumInitialNeighborsSize; i++) {
                SphericParticle* p_neighbour_sphere = current_particle->mNeighbourElements[i];
		        if (p_neighbour_sphere == NULL) continue;

                if (p_neighbour_sphere->Is(VISITED)) continue; //not necessary, but saves increasing and decreasing stack_of_particles_to_check's size
                if (current_particle->mIniNeighbourFailureId[i]) continue;

                auto existing_element_it = dem_model_part.GetMesh(0).Elements().find(p_neighbour_sphere->Id());
                if (existing_element_it == dem_model_part.GetMesh(0).ElementsEnd()) continue;

                SphericContinuumParticle* p_neigh_cont_sphere = dynamic_cast<SphericContinuumParticle*>(p_neighbour_sphere);
                stack_of_particles_to_check.push_back(p_neigh_cont_sphere);
            }
        }
        if (this_chunk_mass) chunks_masses.push_back(this_chunk_mass);
    }

    const double max_mass_of_a_single_chunck = *std::max_element(chunks_masses.begin(), chunks_masses.end());
    const double current_total_mass_in_grams = max_mass_of_a_single_chunck;
    static const double initial_total_mass_in_grams = current_total_mass_in_grams;
    const double cumulative_sand_mass_in_grams = initial_total_mass_in_grams - current_total_mass_in_grams;

    ProcessInfo& r_process_info = dem_model_part.GetProcessInfo();
    const double Pascals_to_psi_factor = 0.000145;
    const double face_pressure_in_psi = fabs(r_process_info[TARGET_STRESS_Z]) * Pascals_to_psi_factor;

    ofs_sand_prod_file << time << " " << face_pressure_in_psi << " " << cumulative_sand_mass_in_grams << '\n';
    ofs_sand_prod_file.flush();

    unsigned int number_of_time_steps_between_granulometry_prints = 1000;
    static unsigned int printing_counter = 0;
    if (printing_counter == number_of_time_steps_between_granulometry_prints) {
        ofs_granulometry_distr_file << time;
        for (unsigned int k = 0; k < chunks_masses.size(); k++) ofs_granulometry_distr_file << " " << chunks_masses[k];
        ofs_granulometry_distr_file << '\n';
        printing_counter = 0;
    }
    printing_counter++;
    ofs_granulometry_distr_file.flush();
}

void ComputeSandProductionWithDepthFirstSearch(ModelPart& dem_model_part, ModelPart& outer_walls_model_part, const double time) {

    const std::string filename = "sand_production_graph_with_chunks.txt";
    std::ifstream ifile(filename.c_str());
    static bool first_time_entered = true;
    if ((bool) ifile && first_time_entered) {
        std::remove(filename.c_str());
        first_time_entered = false;
    }

    ModelPart::ElementsContainerType& pElements = dem_model_part.GetCommunicator().LocalMesh().Elements();

    std::vector<double> chunks_masses;

    for (unsigned int k = 0; k < pElements.size(); k++) {
        ModelPart::ElementsContainerType::iterator it = pElements.ptr_begin() + k;
        it->Set(VISITED, false);
    }

    for (unsigned int k = 0; k < pElements.size(); k++) {
        ModelPart::ElementsContainerType::iterator it = pElements.ptr_begin() + k;
        Element* raw_p_element = &(*it);
        SphericContinuumParticle* p_sphere = dynamic_cast<SphericContinuumParticle*>(raw_p_element);
        double this_chunk_mass = 0.0;
        if( it->IsNot(VISITED) ) {
            DepthFirstSearchVisit(p_sphere, this_chunk_mass);
            chunks_masses.push_back(this_chunk_mass);
        }
    }

    const double max_mass_of_a_single_chunck = *std::max_element(chunks_masses.begin(), chunks_masses.end());
    const double current_total_mass_in_grams = max_mass_of_a_single_chunck;
    static const double initial_total_mass_in_grams = current_total_mass_in_grams;
    const double cumulative_sand_mass_in_grams = initial_total_mass_in_grams - current_total_mass_in_grams;

    ModelPart::ConditionsContainerType::iterator condition_begin = outer_walls_model_part.ConditionsBegin();
    const double Pascals_to_psi_factor = 0.000145;
    const double face_pressure_in_psi = condition_begin->GetValue(POSITIVE_FACE_PRESSURE) * Pascals_to_psi_factor;

    static std::ofstream sand_prod_file(filename, std::ios_base::out | std::ios_base::app);
    sand_prod_file << time << " " << face_pressure_in_psi << " " << cumulative_sand_mass_in_grams << '\n';
    sand_prod_file.flush();
}

void DepthFirstSearchVisit(SphericContinuumParticle* p_sphere, double& this_chunk_mass) {
    p_sphere->Set(VISITED, true);
    const double particle_radius = p_sphere->GetRadius();
    const double particle_density = p_sphere->GetDensity();
    this_chunk_mass += (4.0/3.0) * Globals::Pi * particle_density * particle_radius * particle_radius * particle_radius * 1000.0;
    for (size_t i=0; i<p_sphere->mContinuumInitialNeighborsSize; i++) {
        SphericParticle* p_neighbour_sphere = p_sphere->mNeighbourElements[i];
        if (p_neighbour_sphere==NULL) continue;
        if (p_sphere->mIniNeighbourFailureId[i]) continue;
        if (p_neighbour_sphere->IsNot(VISITED)) {
            SphericContinuumParticle* p_neigh_cont_sphere = dynamic_cast<SphericContinuumParticle*>(p_neighbour_sphere);
            DepthFirstSearchVisit(p_neigh_cont_sphere, this_chunk_mass);
        }
    }
}

void ComputeTriaxialSandProduction(ModelPart& dem_model_part, ModelPart& outer_walls_model_part_1, ModelPart& outer_walls_model_part_2, const double time) {

    const std::string filename = "sand_production_graph.txt";
    std::ifstream ifile(filename.c_str());
    static bool first_time_entered = true;
    if ((bool) ifile && first_time_entered) {
        std::remove(filename.c_str());
        first_time_entered = false;
    }

    ModelPart::ElementsContainerType& pElements = dem_model_part.GetCommunicator().LocalMesh().Elements();
    double current_total_mass_in_grams = 0.0;

    for (unsigned int k = 0; k < pElements.size(); k++) {

        ModelPart::ElementsContainerType::iterator it = pElements.ptr_begin() + k;
        Element* raw_p_element = &(*it);
        SphericParticle* p_sphere = dynamic_cast<SphericParticle*>(raw_p_element);
        if (p_sphere->Is(ISOLATED)) continue;
        const double particle_radius = p_sphere->GetRadius();
        const double particle_density = p_sphere->GetDensity();
        current_total_mass_in_grams += (4.0/3.0) * Globals::Pi * particle_density * particle_radius * particle_radius * particle_radius * 1000.0;
    }
    static const double initial_total_mass_in_grams = current_total_mass_in_grams;
    const double cumulative_sand_mass_in_grams = initial_total_mass_in_grams - current_total_mass_in_grams;

    ModelPart::ConditionsContainerType::iterator condition_begin_1 = outer_walls_model_part_1.ConditionsBegin();
    ModelPart::ConditionsContainerType::iterator condition_begin_2 = outer_walls_model_part_2.ConditionsBegin();

    const double Pascals_to_psi_factor = 0.000145;
    const double face_pressure_in_psi = (condition_begin_1->GetValue(POSITIVE_FACE_PRESSURE) +
                                         condition_begin_2->GetValue(POSITIVE_FACE_PRESSURE) +
                                         3.45e6) * Pascals_to_psi_factor * 0.33333333333333; // 3.45e6 is the sigma_z constant pressure

    static std::ofstream sand_prod_file(filename, std::ios_base::out | std::ios_base::app);
    sand_prod_file << time << " " << face_pressure_in_psi << " " << cumulative_sand_mass_in_grams << '\n';
    sand_prod_file.flush();
}

//***************************************************************************************************************
//***************************************************************************************************************

/// Turn back information as a stemplate<class T, std::size_t dim> tring.

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

/// Assignment operator
DemStructuresCouplingUtilities & operator=(DemStructuresCouplingUtilities const& rOther);


///@}

}; // Class DemStructuresCouplingUtilities

}  // namespace Python.

#endif // KRATOS_STRUCTURES_DEM_COUPLING_UTILITIES_H
