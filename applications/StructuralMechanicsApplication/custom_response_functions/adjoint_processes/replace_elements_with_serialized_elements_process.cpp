// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Mahmoud Sesa
//

// System includes

// External includes

// Project includes
#include "utilities/compare_elements_and_conditions_utility.h"
#include "replace_elements_with_serialized_elements_process.h"
#include "structural_mechanics_application_variables.h"
#include "custom_response_functions/adjoint_elements/adjoint_finite_difference_base_element.h"
#include "custom_response_functions/adjoint_elements/adjoint_finite_difference_shell_element.h"
#include "custom_response_functions/adjoint_elements/adjoint_finite_difference_cr_beam_element_3D2N.h"
#include "custom_response_functions/adjoint_elements/adjoint_finite_difference_truss_element_3D2N.h"
#include "custom_response_functions/adjoint_elements/adjoint_finite_difference_truss_element_linear_3D2N.h"
#include "utilities/compare_elements_and_conditions_utility.h"
#include "custom_elements/shell_thin_element_3D3N.hpp"

namespace Kratos
{
    void ReplaceElementsWithSerializedElementsProcess::Execute()
    {
        KRATOS_TRY

        auto it = mrMainModelPart.ElementsBegin();
        this->TypeCast(*it);
        //TODO Mahmoud: use OpenMP
        // for(int i=0; i<static_cast<int>(mrMainModelPart.NumberOfElements()); ++i)
        // {
        //     auto it = mrMainModelPart.ElementsBegin() + i;

        //     // pGetElement index removed because it causes segmentation fault
        //     auto p_it_loaded_element = mrLoadedModelPart.pGetElement(it->Id());

        //     this->TypeCast(*it);

        //     // KratosComponents.GetRegisteredName(p_it_loaded_element)
        //     // get adjoint type name from table
        //     // GEt KratosComponents.registeredElement
        //     //use that to do the dynamic cast
        //     // std::string current_name;
        //     // CompareElementsAndConditionsUtility::GetRegisteredName(*it, current_name);
        //     // Element rReferenceElement = KratosComponents<Element>::Get(current_name);

        //     // // done like that temporarily to compile the code
        //     // auto p_adjoint_element = dynamic_pointer_cast<AdjointFiniteDifferencingBaseElement<ShellThinElement3D3N>>(*it.base());
        //     // if (p_adjoint_element != nullptr)
        //     // {
        //     //     p_adjoint_element->SetPrimalElement(p_it_loaded_element);
        //     // }
        // }

        KRATOS_CATCH("")
    }

    template <class TPrimalElement>
    void ReplaceElementsWithSerializedElementsProcess::TypeCast(TPrimalElement Element)
    {
        for(int i=0; i<static_cast<int>(mrMainModelPart.NumberOfElements()); ++i)
        {
            auto it = mrMainModelPart.ElementsBegin() + i;

            // pGetElement index removed because it causes segmentation fault
            auto p_it_loaded_element = mrLoadedModelPart.pGetElement(it->Id());

            // done like that temporarily to compile the code
            auto p_adjoint_element = dynamic_pointer_cast<AdjointFiniteDifferencingBaseElement<TPrimalElement>>(*it.base());
            if (p_adjoint_element != nullptr)
            {
                p_adjoint_element->SetPrimalElement(p_it_loaded_element);
            }
        }
    }
}