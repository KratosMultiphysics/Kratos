/*
 * Author: Salva Latorre and Ignasi Pouplana
 *
 *  latorre@cimne.upc.edu
 *  ipouplana@cimne.upc.edu
 */

#ifndef COMPUTE_DEM_FACE_LOAD_UTILITY_H
#define COMPUTE_DEM_FACE_LOAD_UTILITY_H

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
    
    class ComputeDEMFaceLoadUtility {

        public:
        
        typedef ModelPart::NodesContainerType::ContainerType::iterator NodesIteratorType;

        KRATOS_CLASS_POINTER_DEFINITION(ComputeDEMFaceLoadUtility);

        ComputeDEMFaceLoadUtility() {}

        virtual ~ComputeDEMFaceLoadUtility() {}
        
        void ClearDEMFaceLoads(ModelPart& r_structural_skin_model_part) {
            
            KRATOS_TRY

            ModelPart::NodesContainerType::iterator node_it_begin = r_structural_skin_model_part.NodesBegin();
            ModelPart::NodesContainerType::iterator node_it_end = r_structural_skin_model_part.NodesEnd();
            // TODO: OMP
            for (auto node_it = node_it_begin; node_it != node_it_end; ++node_it) {

                array_1d<double, 3>& node_rhs = node_it->FastGetSolutionStepValue(DEM_SURFACE_LOAD);
                node_rhs = ZeroVector(3);                
            }
            
            KRATOS_CATCH("")
        }
        
        void CalculateDEMFaceLoads(ModelPart& r_structural_skin_model_part, const double DEM_delta_time, const double FEM_delta_time) {
            
            KRATOS_TRY

            ModelPart::NodesContainerType::iterator node_it_begin = r_structural_skin_model_part.NodesBegin();
            ModelPart::NodesContainerType::iterator node_it_end   = r_structural_skin_model_part.NodesEnd();
            
            static bool nodal_area_already_computed = false;
            
            if (!nodal_area_already_computed) {
                // TODO: OMP
                for (ModelPart::NodesContainerType::iterator node_it = node_it_begin; node_it != node_it_end; ++node_it) {
                    double& node_area = node_it->GetSolutionStepValue(DEM_NODAL_AREA);
                    node_area = 0.0;
                }

                ModelPart::ConditionsContainerType& source_conditions = r_structural_skin_model_part.Conditions();

                for (unsigned int i = 0; i < source_conditions.size(); i++) {
                    ModelPart::ConditionsContainerType::iterator it = r_structural_skin_model_part.ConditionsBegin() + i;
                    Condition::GeometryType& geometry =  it->GetGeometry();
                    double Element_Area = geometry.Area();

                    for (unsigned int i = 0; i < geometry.size(); i++) { //talking about each of the three nodes of the condition
                        double& node_area = geometry[i].FastGetSolutionStepValue(DEM_NODAL_AREA);
                        node_area += 0.333333333333333 * Element_Area; //TODO: ONLY FOR TRIANGLE... Generalize for 3 or 4 nodes
                    }
                }
                
                nodal_area_already_computed = true;
            }
            // TODO: OMP
            for (auto node_it = node_it_begin; node_it != node_it_end; ++node_it) {

                double& nodal_area = node_it->FastGetSolutionStepValue(DEM_NODAL_AREA);

                if (nodal_area && FEM_delta_time) {
                    node_it->FastGetSolutionStepValue(DEM_SURFACE_LOAD) += node_it->FastGetSolutionStepValue(CONTACT_FORCES) * DEM_delta_time / (nodal_area * FEM_delta_time);
                }                
            }
            
            KRATOS_CATCH("")
        }

        virtual std::string Info() const { return "";}

        virtual void PrintInfo(std::ostream& rOStream) const {}

        virtual void PrintData(std::ostream& rOStream) const {}

        private:

        ComputeDEMFaceLoadUtility& operator= (ComputeDEMFaceLoadUtility const& rOther);

    }; // class ComputeDEMFaceLoadUtility

} // namespace Kratos

#endif // COMPUTE_DEM_FACE_LOAD_UTILITY_H
