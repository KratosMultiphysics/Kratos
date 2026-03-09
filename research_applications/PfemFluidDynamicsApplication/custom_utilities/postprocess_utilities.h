#ifndef POSTPROCESS_UTILITIES_H
#define POSTPROCESS_UTILITIES_H

#include "utilities/timer.h"
#include "includes/define.h"
#include "includes/variables.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Kratos
{

class PostProcessUtilities
{

public:
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::NodesContainerType NodesContainerType;

    KRATOS_CLASS_POINTER_DEFINITION(PostProcessUtilities);

    /// Default constructor.

    PostProcessUtilities(){};

    /// Destructor.

    virtual ~PostProcessUtilities(){};

    void RebuildPostProcessModelPart(ModelPart &r_post_model_part, ModelPart &r_main_model_part) {

        // Cleaning the Output Model Part
        r_post_model_part.Elements().clear();
        r_post_model_part.Nodes().clear();

        // Adding nodes
        for (size_t i = 0; i < r_main_model_part.NumberOfNodes(); i++) {
            auto node = r_main_model_part.NodesBegin() + i;
            r_post_model_part.AddNode(*(node.base()));
        }

        // Adding elements
        PointerVector<Element> elements_to_be_added;
        for (ModelPart::SubModelPartsContainerType::iterator i_smp = r_main_model_part.SubModelPartsBegin();
            i_smp != r_main_model_part.SubModelPartsEnd(); ++i_smp) {

            // Skipping the Computing Model Part
            if (!((i_smp->Is(ACTIVE) && i_smp->Is(SOLID)) || (i_smp->Is(ACTIVE) && i_smp->Is(FLUID)))) {
                if (i_smp->NumberOfElements()) {
                    ModelPart &sub_model_part = *i_smp;
                    for (size_t i = 0; i < sub_model_part.NumberOfElements(); i++) {

                        auto elem = sub_model_part.ElementsBegin() + i;
                        elements_to_be_added.push_back(*(elem.base()));
                    }
                }
            }
        }

        for (auto &elem : elements_to_be_added) {
            elem.Set(ACTIVE, true);
        }
        r_post_model_part.AddElements(elements_to_be_added.begin(), elements_to_be_added.end());
    }

protected:
}; // Class PostProcessUtilities

} // namespace Kratos.

#endif // POSTPROCESS_UTILITIES_H
