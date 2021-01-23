//-------------------------------------------------------------
//         ___  __           ___ _      _    _
//  KRATOS| _ \/ _|___ _ __ | __| |_  _(_)__| |
//        |  _/  _/ -_) '  \| _|| | || | / _` |
//        |_| |_| \___|_|_|_|_| |_|\_,_|_\__,_|DYNAMICS
//
//  BSD License:    PfemFluidDynamicsApplication/license.txt
//
//  Main authors:   Massimiliano Zecchetto
//  Collaborators:
//
//-------------------------------------------------------------
//

#if !defined(UPDATE_THERMAL_MODEL_PART_PROCESS)
#define UPDATE_THERMAL_MODEL_PART_PROCESS

#include <algorithm>
#include <iostream>
#include <string>
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "utilities/variable_utils.h"

namespace Kratos {

/**
 * This method fills rDestinationModelPart using data obtained from rOriginModelPart.
 * The elements of rDestinationModelPart part use the same connectivity (and id)
 * as in rOriginModelPart but their type is determined by rReferenceElement.
 * Both ModelParts will share the same nodes.
 * @class UpdateThermalModelPartProcess
 * @ingroup PfemFluidDynamicsApplication
 * @author Massimiliano Zecchetto
 */
class UpdateThermalModelPartProcess : public Process {
   public:
    typedef Element BaseType;

    KRATOS_CLASS_POINTER_DEFINITION(UpdateThermalModelPartProcess);

    /// Constructor
    UpdateThermalModelPartProcess(ModelPart& origin_model_part, ModelPart& destination_model_part,
                                  ModelPart& computing_model_part, unsigned int DomainSize)
        : rOriginModelPart(origin_model_part),
          rDestinationModelPart(destination_model_part),
          rComputingModelPart(computing_model_part) {
        rDomainSize = DomainSize;
        if (rDomainSize == 2) {
            ReferenceElement = "EulerianConvDiff2D";
        } else {
            ReferenceElement = "EulerianConvDiff3D";
        }
    }

    /// Destructor.
    ~UpdateThermalModelPartProcess() override {}

    void operator()() { Execute(); }

    void Execute() override {
        KRATOS_TRY;

        this->ResetDestinationModelPart();
        this->CopyNodes();
        this->DuplicateElements();
        this->BuildThermalComputingDomain();

        KRATOS_CATCH("");
    }

    void ExecuteInitialize() override {}

    void ExecuteInitializeSolutionStep() override {}

   protected:
    ModelPart& rOriginModelPart;
    ModelPart& rDestinationModelPart;
    ModelPart& rComputingModelPart;
    std::string ReferenceElement;
    unsigned int rDomainSize;

   private:
    void ResetDestinationModelPart() const {
        rOriginModelPart.RemoveNodesFromAllLevels(TO_ERASE);
        VariableUtils().SetFlag(TO_ERASE, true, rDestinationModelPart.Nodes());
        VariableUtils().SetFlag(TO_ERASE, true, rDestinationModelPart.Elements());
        rDestinationModelPart.RemoveNodesFromAllLevels(TO_ERASE);
        rDestinationModelPart.RemoveElementsFromAllLevels(TO_ERASE);
        VariableUtils().SetFlag(TO_ERASE, false, rOriginModelPart.Nodes());
    }

    void CopyNodes() const {
        // Copy nodes to the rDestinationModelPart
        rDestinationModelPart.AddNodes(rOriginModelPart.NodesBegin(), rOriginModelPart.NodesEnd());

        // Copy nodes to the rDestinationModelPart's sub model parts
        for (auto i_part = rOriginModelPart.SubModelPartsBegin(); i_part != rOriginModelPart.SubModelPartsEnd();
             ++i_part) {
            if (!rDestinationModelPart.HasSubModelPart(i_part->Name())) {
                rDestinationModelPart.CreateSubModelPart(i_part->Name());
            }
            ModelPart& destination_part = rDestinationModelPart.GetSubModelPart(i_part->Name());
            destination_part.AddNodes(i_part->NodesBegin(), i_part->NodesEnd());
        }
    }

    void DuplicateElements() const {
        const Element& mReferenceElement = KratosComponents<Element>::Get(ReferenceElement);
        for (ModelPart::SubModelPartIterator i_mp = rOriginModelPart.SubModelPartsBegin();
             i_mp != rOriginModelPart.SubModelPartsEnd(); i_mp++) {
            if (i_mp->NumberOfElements()) {
                ModelPart& destination_part = rDestinationModelPart.GetSubModelPart(i_mp->Name());
                ModelPart::ElementsContainerType temp_elements;
                temp_elements.reserve(i_mp->NumberOfElements());

                // Skip the ComputingModelPart of the PfemFluidModelPart
                if ((i_mp->Is(SOLID) && i_mp->IsNot(ACTIVE)) || (i_mp->Is(FLUID) && i_mp->IsNot(ACTIVE)) ||
                    (i_mp->Is(BOUNDARY) && i_mp->Is(RIGID))) {
                    for (auto it_elem = i_mp->ElementsBegin(); it_elem != i_mp->ElementsEnd(); ++it_elem) {
                        Properties::Pointer properties = it_elem->pGetProperties();

                        Element::Pointer p_element =
                            mReferenceElement.Create(it_elem->Id(), it_elem->GetGeometry(), properties);
                        temp_elements.push_back(p_element);
                    }
                    rDestinationModelPart.AddElements(temp_elements.begin(), temp_elements.end());
                    destination_part.AddElements(temp_elements.begin(), temp_elements.end());
                }
            }
        }
    }

    void BuildThermalComputingDomain() const {
        // Copy nodes
        rComputingModelPart.AddNodes(rDestinationModelPart.NodesBegin(), rDestinationModelPart.NodesEnd());

        // Copy elements
        std::vector<ModelPart::IndexType> ids;
        ids.reserve(rDestinationModelPart.Elements().size());
        const auto& it_elem_begin = rDestinationModelPart.ElementsBegin();

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rDestinationModelPart.Elements().size()); i++) {
            auto it_elem = it_elem_begin + i;
            #pragma omp critical
            { ids.push_back(it_elem->Id()); }
        }

        rComputingModelPart.AddElements(ids, 0);
    }

};  // Class UpdateThermalModelPartProcess

/// input stream function
inline std::istream& operator>>(std::istream& rIStream, UpdateThermalModelPartProcess& rThis);

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream, const UpdateThermalModelPartProcess& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}  // namespace Kratos.

#endif /* UPDATE_THERMAL_MODEL_PART_PROCESS defined */