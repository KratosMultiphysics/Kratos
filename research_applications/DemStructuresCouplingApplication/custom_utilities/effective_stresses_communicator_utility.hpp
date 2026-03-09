/*
 * Author: Miguel Angel Celigueta
 *
 *  maceli@cimne.upc.edu
 */

#ifndef EFFECTIVE_STRESSES_COMMUNICATOR_UTILITY_H
#define EFFECTIVE_STRESSES_COMMUNICATOR_UTILITY_H

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

    class EffectiveStressesCommunicatorUtility {

        public:

        typedef ModelPart::NodesContainerType::ContainerType::iterator NodesIteratorType;

        KRATOS_CLASS_POINTER_DEFINITION(EffectiveStressesCommunicatorUtility);

        EffectiveStressesCommunicatorUtility(ModelPart& r_source_model_part, ModelPart& destination_model_part):mrSourceModelPart(r_source_model_part), mrDestinationModelPart(destination_model_part) {
            Check();
        }

        virtual ~EffectiveStressesCommunicatorUtility() {
            delete mpSearchStructure;
        }

        void Initialize() {
            mpSearchStructure = new BinBasedFastPointLocator<2>(mrSourceModelPart);
            mpSearchStructure->UpdateSearchDatabase();
            InitializeVariablesToZero();
        }

        void CopyWallCurrentEffectiveStressesToOldEffectiveStresses() {

            KRATOS_TRY

            #pragma omp parallel for
            for (int i=0; i<(int)mrDestinationModelPart.Nodes().size(); i++) {
                auto node_it = mrDestinationModelPart.NodesBegin() + i;
                node_it->SetValue(OLD_RADIAL_NORMAL_STRESS_COMPONENT, node_it->GetValue(RADIAL_NORMAL_STRESS_COMPONENT));
            }

            KRATOS_CATCH("")
        }

        void CommunicateCurrentRadialEffectiveStressesToDemWalls() {

            KRATOS_TRY
            const int max_results = 10000;

            #pragma omp parallel
            {
                Vector N;
                typename BinBasedFastPointLocator<2>::ResultContainerType results(max_results);
                typename BinBasedFastPointLocator<2>::ResultIteratorType results_begin = results.begin();
                Vector unitary_radial_vector = ZeroVector(3);

                #pragma omp for
                for (int i = 0; i < (int)mrDestinationModelPart.Nodes().size(); i++) {
                    auto node_it = mrDestinationModelPart.NodesBegin() + i;
                    auto& particle_coordinates = node_it->Coordinates();
                    const double norm = MathUtils<double>::Norm(particle_coordinates);
                    if (norm > std::numeric_limits<double>::epsilon()) {
                        const double inv_norm = 1.0 / norm;
                        unitary_radial_vector[0] = particle_coordinates[0] * inv_norm;
                        unitary_radial_vector[1] = particle_coordinates[1] * inv_norm;
                        unitary_radial_vector[2] = particle_coordinates[2] * inv_norm;
                    }
                    else {
                        continue; //In the DEM walls modelpart there is a node which represents the center of gravity of the solid...
                        //KRATOS_ERROR<<"Radial stresses cannot be transferred to a point which is in the origin of coordinates."<<std::endl;
                    }

                    bool is_found = false;
                    Element::Pointer shared_p_element;
                    is_found = mpSearchStructure->FindPointOnMesh(particle_coordinates, N, shared_p_element, results_begin, max_results, 0.0);
                    if (is_found) {
                        Matrix interpolated_effective_stress_tensor = ZeroMatrix(3, 3);
                        const auto& geom = shared_p_element->GetGeometry();
                        for (size_t j = 0; j < geom.size(); j++) {
                            const Matrix& tempM = geom[j].FastGetSolutionStepValue(NODAL_EFFECTIVE_STRESS_TENSOR);
                            noalias(interpolated_effective_stress_tensor) += N[j] * tempM;
                        }
                        Vector tempV = prod(interpolated_effective_stress_tensor, unitary_radial_vector);
                        double radial_stress = MathUtils<double>::Dot(tempV, unitary_radial_vector);
                        node_it->SetValue(RADIAL_NORMAL_STRESS_COMPONENT, radial_stress);
                    }
                }
            }

            KRATOS_CATCH("")
        }

        void CommunicateGivenRadialEffectiveStressesToDemWalls(const Matrix stress_tensor) {

            KRATOS_TRY

            #pragma omp parallel
            {
                Vector unitary_radial_vector = ZeroVector(2);

                #pragma omp for
                for (int i=0; i<(int)mrDestinationModelPart.Nodes().size(); i++) {
                    auto node_it = mrDestinationModelPart.NodesBegin() + i;
                    auto& particle_coordinates = node_it->Coordinates();
                    const double norm = MathUtils<double>::Norm(particle_coordinates);
                    if (norm > std::numeric_limits<double>::epsilon()) {
                        const double inv_norm = 1.0 / norm;
                        unitary_radial_vector[0] = particle_coordinates[0] * inv_norm;
                        unitary_radial_vector[1] = particle_coordinates[1] * inv_norm;
                    }
                    else {
                        continue; //In the DEM walls modelpart there is a node which represents the center of gravity of the solid...
                        //KRATOS_ERROR<<"Radial stresses cannot be transferred to a point which is in the origin of coordinates."<<std::endl;
                    }

                    Vector tempV = prod(stress_tensor, unitary_radial_vector);
                    double radial_stress = MathUtils<double>::Dot(tempV, unitary_radial_vector);
                    node_it->SetValue(RADIAL_NORMAL_STRESS_COMPONENT, radial_stress);
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

        EffectiveStressesCommunicatorUtility& operator= (EffectiveStressesCommunicatorUtility const& rOther);

        void Check() {

        }

        void InitializeVariablesToZero() {
            KRATOS_TRY

            #pragma omp parallel for
            for (int i=0; i<(int)mrDestinationModelPart.Nodes().size(); i++) {
                auto node_it = mrDestinationModelPart.NodesBegin() + i;
                node_it->SetValue(OLD_RADIAL_NORMAL_STRESS_COMPONENT, 0.0);
                node_it->SetValue(RADIAL_NORMAL_STRESS_COMPONENT, 0.0);
            }

            KRATOS_CATCH("")
        }

    }; // class EffectiveStressesCommunicatorUtility

} // namespace Kratos

#endif // EFFECTIVE_STRESSES_COMMUNICATOR_UTILITY_H
