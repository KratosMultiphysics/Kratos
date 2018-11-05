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
#include "includes/model_part.h"
#include "../../DEM_application/custom_conditions/RigidFace.h"
#include "../../DEM_application/DEM_application_variables.h"

namespace Kratos {
    
    class ComputeDEMFaceLoadUtility {

        public:
        
        typedef ModelPart::NodesContainerType::ContainerType::iterator NodesIteratorType;

        KRATOS_CLASS_POINTER_DEFINITION(ComputeDEMFaceLoadUtility);

        ComputeDEMFaceLoadUtility() {}

        virtual ~ComputeDEMFaceLoadUtility() {}

        virtual std::string Info() const {
            return "";
        }

        virtual void PrintInfo(std::ostream& rOStream) const {}

        virtual void PrintData(std::ostream& rOStream) const {}

        private:

        ComputeDEMFaceLoadUtility& operator= (ComputeDEMFaceLoadUtility const& rOther);

    }; // class ComputeDEMFaceLoadUtility

} // namespace Kratos

#endif // COMPUTE_DEM_FACE_LOAD_UTILITY_H
