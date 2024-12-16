/*
 * Author: Miguel Angel Celigueta
 *
 *  maceli@cimne.upc.edu
 */

#ifndef PERMEABILITY_TENSOR_COMMUNICATOR_UTILITY_H
#define PERMEABILITY_TENSOR_COMMUNICATOR_UTILITY_H

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

// Project includes
#include "spatial_containers/dem_search.h"
#include "utilities/openmp_utils.h"

// Configures
#include "../../DEMApplication/custom_utilities/discrete_particle_configure.h"
#include "../../DEMApplication/custom_utilities/geometrical_object_configure.h"
#include "../../DEMApplication/custom_utilities/node_configure.h"
#include "../../DEMApplication/custom_utilities/omp_dem_search.h"

// Search
#include "spatial_containers/bins_dynamic_objects.h"
#include "spatial_containers/bins_dynamic.h"
#include "custom_search/bins_dynamic_objects_periodic.h"

namespace Kratos {

    class KRATOS_API(DEM_STRUCTURES_COUPLING_APPLICATION) PermeabilityTensorCommunicatorUtility {

        public:

        typedef ModelPart::NodesContainerType::ContainerType::iterator  NodesIteratorType;
        typedef SpatialSearch                                           SearchType;
        typedef SearchType::ElementsContainerType                       ElementsContainerType;
        typedef SearchType::ElementsContainerType::ContainerType        ContainerType;
        typedef SearchType::NodesContainerType                          NodesContainerType;
        typedef ContainerType::value_type                               PointerType;
        typedef ContainerType::iterator                                 IteratorType;
        typedef ElementsContainerType::iterator                         ElementIteratorType;
        typedef SearchType::ElementsContainerType::ContainerType        ResultContainerType;
        typedef DiscreteParticleConfigure<3>                  ElementConfigureType;
        typedef BinsObjectDynamic<ElementConfigureType>               BinsType;
        typedef std::unique_ptr<BinsType>                             BinsUniquePointerType;
        typedef ElementsContainerType::ContainerType              ResultElementsContainerType;
        typedef SpatialSearch::DistanceType DistanceType;

        KRATOS_CLASS_POINTER_DEFINITION(PermeabilityTensorCommunicatorUtility);

        PermeabilityTensorCommunicatorUtility(ModelPart& r_source_model_part, ModelPart& destination_model_part):mrDEMModelPart(r_source_model_part), mrFEMModelPart(destination_model_part) {
            Check();
        }

        virtual ~PermeabilityTensorCommunicatorUtility() {
            delete mpSearchStructure;
        }

        void Initialize() {
            mpSearchStructure = new BinBasedFastPointLocator<2>(mrDEMModelPart);
            mpSearchStructure->UpdateSearchDatabase();
        }

        void Diagonalize(const BoundedMatrix<double, 3, 3> A, BoundedMatrix<double, 3, 3>& Q, BoundedMatrix<double, 3, 3>& D);
        void EigenVectors(const BoundedMatrix<double, 3, 3>& A,
                                    BoundedMatrix<double, 3, 3>& vectors,
                                    double zero_tolerance =1e-14,
                                    int max_iterations = 50);

        BoundedMatrix<double, 3, 3> Transpose(BoundedMatrix<double, 3, 3> A);

        void TrasferUpdatedPermeabilityTensor();

        virtual std::string Info() const { return "";}
        virtual void PrintInfo(std::ostream& rOStream) const {}
        virtual void PrintData(std::ostream& rOStream) const {}

        private:

        ModelPart& mrDEMModelPart;
        ModelPart& mrFEMModelPart;
        BinBasedFastPointLocator<2>* mpSearchStructure;

        PermeabilityTensorCommunicatorUtility& operator= (PermeabilityTensorCommunicatorUtility const& rOther);

        void Check() {}

    }; // class PermeabilityTensorCommunicatorUtility

} // namespace Kratos

#endif // PERMEABILITY_TENSOR_COMMUNICATOR_UTILITY_H
