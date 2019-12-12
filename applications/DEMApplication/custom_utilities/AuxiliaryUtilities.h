#pragma once
// Author: Miguel Angel Celigueta maceli@cimne.upc.edu

#include <pybind11/pybind11.h>
#include "includes/define.h"
#include "includes/model_part.h"
#include "pre_utilities.h"


namespace Kratos {

    class AuxiliaryUtilities {

    public:

        typedef ModelPart::ElementsContainerType                         ElementsArrayType;
        typedef ModelPart::NodesContainerType::ContainerType             NodesContainerType;
        typedef GlobalPointersVector<Element>                            ParticleWeakVectorType;
        typedef GlobalPointersVector<Element>::iterator                  ParticleWeakIteratorType;

        KRATOS_CLASS_POINTER_DEFINITION(AuxiliaryUtilities);

        AuxiliaryUtilities() {};
        virtual ~AuxiliaryUtilities() {};

        double ComputeAverageZStressFor2D(ModelPart& rSpheresModelPart) {

            ElementsArrayType& pElements = rSpheresModelPart.GetCommunicator().LocalMesh().Elements();
            double sub_total = 0.0;
            double average_value = 0.0;

            #pragma omp parallel for
            for (int k = 0; k < (int)pElements.size(); k++) {

                ElementsArrayType::iterator it = pElements.ptr_begin() + k;
                Element* p_element = &(*it);
                SphericContinuumParticle* p_sphere = dynamic_cast<SphericContinuumParticle*>(p_element);

                double z_tensor_value = (*p_sphere->mSymmStressTensor)(2,2);
                sub_total += z_tensor_value;
            }
            average_value = sub_total/(int)pElements.size();
            return average_value;
        }
    };
}
