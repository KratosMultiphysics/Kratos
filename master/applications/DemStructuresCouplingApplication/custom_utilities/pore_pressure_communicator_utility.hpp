/*
 * Author: Miguel Angel Celigueta
 *
 *  maceli@cimne.upc.edu
 */

#ifndef PORE_PRESSURE_COMMUNICATOR_UTILITY_H
#define PORE_PRESSURE_COMMUNICATOR_UTILITY_H

#include "includes/variables.h"
#include <limits>
#include <iostream>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "includes/define.h"
#include "includes/condition.h"
#include "includes/model_part.h"
#include "dem_structures_coupling_application_variables.h"
#include "../../PoromechanicsApplication/poromechanics_application_variables.h"
#include "../../DEMApplication/DEM_application_variables.h"

#include "utilities/binbased_fast_point_locator.h"

namespace Kratos {

    class PorePressureCommunicatorUtility {

        public:

        typedef ModelPart::NodesContainerType::ContainerType::iterator NodesIteratorType;

        KRATOS_CLASS_POINTER_DEFINITION(PorePressureCommunicatorUtility);

        PorePressureCommunicatorUtility(ModelPart& r_source_model_part, ModelPart& destination_model_part):mrSourceModelPart(r_source_model_part), mrDestinationModelPart(destination_model_part) {
            Check();
        }

        virtual ~PorePressureCommunicatorUtility() {
            delete mpSearchStructure;
        }

        void Initialize() {
            mpSearchStructure = new BinBasedFastPointLocator<2>(mrSourceModelPart);
            mpSearchStructure->UpdateSearchDatabase();
        }

        void ComputeForceOnParticlesDueToPorePressureGradient() {

            KRATOS_TRY
            
            const int max_results = 10000;

            #pragma omp parallel
            {
                Vector  N;
                typename BinBasedFastPointLocator<2>::ResultContainerType results(max_results);
                typename BinBasedFastPointLocator<2>::ResultIteratorType results_begin = results.begin();

                int property_id = 0;
                for (int i = 0; i < (int)mrDestinationModelPart.Elements().size(); i++) {
                    const auto elem_it = mrDestinationModelPart.ElementsBegin() + i;
                    property_id = elem_it->GetProperties().Id();
                    break;
                }
                const double porosity = mrDestinationModelPart.GetProperties(property_id)[POROSITY];
                const double particle_volume_to_voronoi_volume_factor = 1.0 / (1.0 - porosity);

                #pragma omp for
                for (int i = 0; i < (int)mrDestinationModelPart.Elements().size(); i++) {
                    const auto elem_it = mrDestinationModelPart.ElementsBegin() + i;
                    auto& central_node = elem_it->GetGeometry()[0];
                    const auto& particle_coordinates = central_node.Coordinates();
                    SphericParticle* particle = dynamic_cast<SphericParticle*>(&*elem_it);
                    const double particle_volume = particle->CalculateVolume();
                    const double particle_mass = particle->GetMass();

                    bool is_found = false;
                    Element::Pointer shared_p_element;
                    is_found = mpSearchStructure->FindPointOnMesh(particle_coordinates, N, shared_p_element, results_begin, max_results, 0.0);
                    double gravity = 9.81;
                    const double well_radius = 0.075;
                    if ((particle_coordinates[0] * particle_coordinates[0] + particle_coordinates[1] * particle_coordinates[1]) < well_radius * well_radius) {
                        noalias(central_node.FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE)) = -1.0 * particle_coordinates * particle_mass * gravity / well_radius;
                    } else if (is_found) {
                        const auto& geom = shared_p_element->GetGeometry();
                        array_1d<double, 3> interpolated_gradient_of_pore_pressure = ZeroVector(3);
                        for (size_t j = 0; j < geom.size(); j++) {
                            noalias(interpolated_gradient_of_pore_pressure) += N[j] * geom[j].FastGetSolutionStepValue(LIQUID_PRESSURE_GRADIENT);
                        }
                        noalias(central_node.FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE)) = -1.0 * interpolated_gradient_of_pore_pressure * particle_volume * particle_volume_to_voronoi_volume_factor;
                    }
                }
            }

            KRATOS_CATCH("")
        }

        virtual std::string Info() const { return "";}

        virtual void PrintInfo(std::ostream& rOStream) const {}

        virtual void PrintData(std::ostream& rOStream) const {}

        private:

        ModelPart& mrSourceModelPart;
        ModelPart& mrDestinationModelPart;
        BinBasedFastPointLocator<2>* mpSearchStructure;

        PorePressureCommunicatorUtility& operator= (PorePressureCommunicatorUtility const& rOther);

        void Check() {

        }

    }; // class PorePressureCommunicatorUtility

} // namespace Kratos

#endif // PORE_PRESSURE_COMMUNICATOR_UTILITY_H
