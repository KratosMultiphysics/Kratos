/*
 * Author: Salva Latorre and Ignasi Pouplana
 *
 *  latorre@cimne.upc.edu
 *  ipouplana@cimne.upc.edu
 */

#ifndef INTERPOLATE_STRUCTURAL_SOLUTION_FOR_DEM_UTILITY_H
#define INTERPOLATE_STRUCTURAL_SOLUTION_FOR_DEM_UTILITY_H

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

namespace Kratos {

    class InterpolateStructuralSolutionForDEM {

        public:

        typedef ModelPart::NodesContainerType::ContainerType::iterator NodesIteratorType;

        KRATOS_CLASS_POINTER_DEFINITION(InterpolateStructuralSolutionForDEM);

        InterpolateStructuralSolutionForDEM() {}

        virtual ~InterpolateStructuralSolutionForDEM() {}

        void SaveStructuralSolution(ModelPart& r_structural_model_part) {

            KRATOS_TRY

            const int NNodes = static_cast<int>(r_structural_model_part.Nodes().size());
            ModelPart::NodesContainerType::iterator node_begin = r_structural_model_part.NodesBegin();

            #pragma omp parallel for
            for (int i = 0; i < NNodes; i++) {

                ModelPart::NodesContainerType::iterator itNode = node_begin + i;

                array_1d<double,3>& r_current_velocity = itNode->FastGetSolutionStepValue(BACKUP_LAST_STRUCTURAL_VELOCITY);
                noalias(r_current_velocity) = itNode->FastGetSolutionStepValue(VELOCITY);

                array_1d<double,3>& r_current_displacement = itNode->FastGetSolutionStepValue(BACKUP_LAST_STRUCTURAL_DISPLACEMENT);
                noalias(r_current_displacement) = itNode->FastGetSolutionStepValue(DISPLACEMENT);

                array_1d<double,3>& r_smoothed_velocity = itNode->FastGetSolutionStepValue(SMOOTHED_STRUCTURAL_VELOCITY);
                noalias(r_smoothed_velocity) = 1.0/3.0 * (itNode->FastGetSolutionStepValue(VELOCITY) + 2.0 * itNode->FastGetSolutionStepValue(SMOOTHED_STRUCTURAL_VELOCITY, 1));

            }

            KRATOS_CATCH("")
        }

        void InterpolateStructuralSolution(ModelPart& r_structural_model_part, const double fem_delta_time, const double fem_time, const double dem_delta_time, const double dem_time) {

            KRATOS_TRY

            const double previous_fem_time = fem_time - fem_delta_time;
            const double time_factor = (dem_time - previous_fem_time) / fem_delta_time;
            const double previous_time_factor = (dem_time - dem_delta_time - previous_fem_time) / fem_delta_time;

            const int NNodes = static_cast<int>(r_structural_model_part.Nodes().size());
            ModelPart::NodesContainerType::iterator node_begin = r_structural_model_part.NodesBegin();


            #pragma omp parallel for
            for (int i = 0; i < NNodes; i++) {

                ModelPart::NodesContainerType::iterator it_node = node_begin + i;

                noalias(it_node->Coordinates()) = it_node->GetInitialPosition().Coordinates()
                                                + it_node->FastGetSolutionStepValue(DISPLACEMENT,1)
                                                + (it_node->FastGetSolutionStepValue(BACKUP_LAST_STRUCTURAL_DISPLACEMENT) - it_node->FastGetSolutionStepValue(DISPLACEMENT,1)) * time_factor;

                array_1d<double,3>& r_velocity = it_node->FastGetSolutionStepValue(VELOCITY);
                const array_1d<double,3>& previous_velocity = it_node->FastGetSolutionStepValue(SMOOTHED_STRUCTURAL_VELOCITY,1);
                noalias(r_velocity) = previous_velocity + (it_node->FastGetSolutionStepValue(SMOOTHED_STRUCTURAL_VELOCITY) - previous_velocity) * time_factor;

                array_1d<double,3>& r_displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);
                noalias(r_displacement) = it_node->Coordinates() - it_node->GetInitialPosition().Coordinates();

                array_1d<double, 3> previous_coordinates;
                noalias(previous_coordinates) = it_node->GetInitialPosition().Coordinates()
                                                + it_node->FastGetSolutionStepValue(DISPLACEMENT,1)
                                                + (it_node->FastGetSolutionStepValue(BACKUP_LAST_STRUCTURAL_DISPLACEMENT) - it_node->FastGetSolutionStepValue(DISPLACEMENT,1)) * previous_time_factor;

                array_1d<double,3>& delta_displacement = it_node->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
                noalias(delta_displacement) = it_node->Coordinates() - previous_coordinates;
            }

            KRATOS_CATCH("")
        }

        void RestoreStructuralSolution(ModelPart& r_structural_model_part) {

            KRATOS_TRY

            const int NNodes = static_cast<int>(r_structural_model_part.Nodes().size());
            ModelPart::NodesContainerType::iterator node_begin = r_structural_model_part.NodesBegin();

            #pragma omp parallel for
            for (int i = 0; i < NNodes; i++) {

                ModelPart::NodesContainerType::iterator it_node = node_begin + i;

                array_1d<double,3>& r_velocity = it_node->FastGetSolutionStepValue(VELOCITY);
                noalias(r_velocity) = it_node->FastGetSolutionStepValue(BACKUP_LAST_STRUCTURAL_VELOCITY);

                array_1d<double,3>& r_displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);
                noalias(r_displacement) = it_node->FastGetSolutionStepValue(BACKUP_LAST_STRUCTURAL_DISPLACEMENT);

                noalias(it_node->Coordinates()) = it_node->GetInitialPosition().Coordinates() + it_node->FastGetSolutionStepValue(DISPLACEMENT);
            }

            KRATOS_CATCH("")
        }

        // void SaveStructuralSolution(ModelPart& r_structural_model_part) {

        //     KRATOS_TRY

        //     const int NNodes = static_cast<int>(r_structural_model_part.Nodes().size());
        //     ModelPart::NodesContainerType::iterator node_begin = r_structural_model_part.NodesBegin();

        //     #pragma omp parallel for
        //     for (int i = 0; i < NNodes; i++) {

        //         ModelPart::NodesContainerType::iterator itNode = node_begin + i;

        //         array_1d<double,3>& r_current_velocity = itNode->FastGetSolutionStepValue(BACKUP_LAST_STRUCTURAL_VELOCITY);
        //         noalias(r_current_velocity) = itNode->FastGetSolutionStepValue(VELOCITY);

        //         array_1d<double,3>& r_current_displacement = itNode->FastGetSolutionStepValue(BACKUP_LAST_STRUCTURAL_DISPLACEMENT);
        //         noalias(r_current_displacement) = itNode->FastGetSolutionStepValue(DISPLACEMENT);

        //         array_1d<double,3>& r_smoothed_velocity = itNode->FastGetSolutionStepValue(SMOOTHED_STRUCTURAL_VELOCITY);
        //         noalias(r_smoothed_velocity) = 1.0/3.0 * (itNode->FastGetSolutionStepValue(VELOCITY) + 2.0 * itNode->FastGetSolutionStepValue(SMOOTHED_STRUCTURAL_VELOCITY, 1));

        //     }

        //     KRATOS_CATCH("")
        // }

        // void InterpolateStructuralSolution(ModelPart& r_structural_model_part, const double fem_delta_time, const double fem_time, const double dem_delta_time, const double dem_time) {

        //     KRATOS_TRY

        //     const double previous_fem_time = fem_time - fem_delta_time;
        //     const double time_factor = (dem_time - previous_fem_time) / fem_delta_time;
        //     const double previous_time_factor = (dem_time - dem_delta_time - previous_fem_time) / fem_delta_time;

        //     const int NNodes = static_cast<int>(r_structural_model_part.Nodes().size());
        //     ModelPart::NodesContainerType::iterator node_begin = r_structural_model_part.NodesBegin();


        //     #pragma omp parallel for
        //     for (int i = 0; i < NNodes; i++) {

        //         ModelPart::NodesContainerType::iterator it_node = node_begin + i;

        //         noalias(it_node->Coordinates()) = it_node->GetInitialPosition().Coordinates()
        //                                         + it_node->FastGetSolutionStepValue(DISPLACEMENT,1)
        //                                         + (it_node->FastGetSolutionStepValue(BACKUP_LAST_STRUCTURAL_DISPLACEMENT) - it_node->FastGetSolutionStepValue(DISPLACEMENT,1)) * time_factor;

        //         array_1d<double,3>& r_velocity = it_node->FastGetSolutionStepValue(VELOCITY);
        //         const array_1d<double,3>& previous_velocity = it_node->FastGetSolutionStepValue(SMOOTHED_STRUCTURAL_VELOCITY,1);
        //         noalias(r_velocity) = previous_velocity + (it_node->FastGetSolutionStepValue(SMOOTHED_STRUCTURAL_VELOCITY) - previous_velocity) * time_factor;

        //         array_1d<double,3>& r_displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);
        //         noalias(r_displacement) = it_node->Coordinates() - it_node->GetInitialPosition().Coordinates();

        //         array_1d<double, 3> previous_coordinates;
        //         noalias(previous_coordinates) = it_node->GetInitialPosition().Coordinates()
        //                                         + it_node->FastGetSolutionStepValue(DISPLACEMENT,1)
        //                                         + (it_node->FastGetSolutionStepValue(BACKUP_LAST_STRUCTURAL_DISPLACEMENT) - it_node->FastGetSolutionStepValue(DISPLACEMENT,1)) * previous_time_factor;

        //         array_1d<double,3>& delta_displacement = it_node->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
        //         noalias(delta_displacement) = it_node->Coordinates() - previous_coordinates;
        //     }

        //     KRATOS_CATCH("")
        // }

        // void RestoreStructuralSolution(ModelPart& r_structural_model_part) {

        //     KRATOS_TRY

        //     const int NNodes = static_cast<int>(r_structural_model_part.Nodes().size());
        //     ModelPart::NodesContainerType::iterator node_begin = r_structural_model_part.NodesBegin();

        //     #pragma omp parallel for
        //     for (int i = 0; i < NNodes; i++) {

        //         ModelPart::NodesContainerType::iterator it_node = node_begin + i;

        //         array_1d<double,3>& r_velocity = it_node->FastGetSolutionStepValue(VELOCITY);
        //         noalias(r_velocity) = it_node->FastGetSolutionStepValue(BACKUP_LAST_STRUCTURAL_VELOCITY);

        //         array_1d<double,3>& r_displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);
        //         noalias(r_displacement) = it_node->FastGetSolutionStepValue(BACKUP_LAST_STRUCTURAL_DISPLACEMENT);

        //         noalias(it_node->Coordinates()) = it_node->GetInitialPosition().Coordinates() + it_node->FastGetSolutionStepValue(DISPLACEMENT);
        //     }

        //     KRATOS_CATCH("")
        // }

        virtual std::string Info() const { return "";}

        virtual void PrintInfo(std::ostream& rOStream) const {}

        virtual void PrintData(std::ostream& rOStream) const {}

        private:

        InterpolateStructuralSolutionForDEM& operator= (InterpolateStructuralSolutionForDEM const& rOther);

    }; // class InterpolateStructuralSolutionForDEM

} // namespace Kratos

#endif // INTERPOLATE_STRUCTURAL_SOLUTION_FOR_DEM_UTILITY_H
