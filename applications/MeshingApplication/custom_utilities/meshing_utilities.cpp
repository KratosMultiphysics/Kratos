// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "includes/key_hash.h"
#include "custom_utilities/meshing_utilities.h"

namespace Kratos
{
namespace MeshingUtilities
{
void KRATOS_API(MESHING_APPLICATION) RecursiveEnsureModelPartOwnsProperties(
    ModelPart& rModelPart,
    const bool RemovePreviousProperties
    )
{
    // First we do in this model part
    EnsureModelPartOwnsProperties(rModelPart, RemovePreviousProperties);

    // Now we do in submodelparts
    for (auto& r_sub_model_part : rModelPart.SubModelParts()) {
        RecursiveEnsureModelPartOwnsProperties(r_sub_model_part, RemovePreviousProperties);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void KRATOS_API(MESHING_APPLICATION) EnsureModelPartOwnsProperties(
    ModelPart& rModelPart,
    const bool RemovePreviousProperties
    )
{
    // First we clear the properties if we want so
    if (RemovePreviousProperties) {
        rModelPart.GetMesh(0).pProperties()->clear();
    }

    // The list of properties
    std::unordered_set<Properties::Pointer, IndexedObjecPointertHasher<Properties::Pointer>, IndexedObjectPointerComparator<Properties::Pointer>> list_of_properties;

    // Iterating over the elements
    auto& r_elements_array = rModelPart.Elements();
    const auto it_elem_begin= r_elements_array.begin();
    const int number_of_elements = static_cast<int>(r_elements_array.size());

    // Iterating over the conditions
    auto& r_conditions_array = rModelPart.Conditions();
    const auto it_cond_begin= r_conditions_array.begin();
    const int number_of_conditions = static_cast<int>(r_conditions_array.size());

    #pragma omp parallel
    {
        // The list of properties
        std::unordered_set<Properties::Pointer, IndexedObjecPointertHasher<Properties::Pointer>, IndexedObjectPointerComparator<Properties::Pointer>> buffer_list_of_properties;

        #pragma omp for schedule(guided, 512) nowait
        for (int i = 0; i < number_of_elements; ++i) {
            auto it_elem = it_elem_begin + i;

            Properties::Pointer p_prop = it_elem->pGetProperties();

            if (buffer_list_of_properties.find(p_prop) == buffer_list_of_properties.end()) {
                buffer_list_of_properties.insert(p_prop);
            }
        }

        #pragma omp for schedule(guided, 512) nowait
        for (int i = 0; i < number_of_conditions; ++i) {
            auto it_cond = it_cond_begin + i;

            Properties::Pointer p_prop = it_cond->pGetProperties();
            if (buffer_list_of_properties.find(p_prop) == buffer_list_of_properties.end()) {
                buffer_list_of_properties.insert(p_prop);
            }
        }

        // Combine buffers together
        #pragma omp critical
        {
            list_of_properties.insert(buffer_list_of_properties.begin(),buffer_list_of_properties.end());
        }
    }

    // Add properties to respective model parts
    for (auto p_prop : list_of_properties) {
        if (!rModelPart.HasProperties(p_prop->Id())) {
            rModelPart.AddProperties(p_prop);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BlockThresholdSizeElements(
    ModelPart& rModelPart,
    Parameters ThisParameters
    )
{
    Parameters default_parameters = Parameters(R"(
    {
        "minimal_size" : 0.1,
        "maximal_size" : 10.0
    })" );

    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // Getting values
    const double minimal_size = ThisParameters["minimal_size"].GetDouble();
    const double maximal_size = ThisParameters["maximal_size"].GetDouble();

    // Compute elements size
    ComputeElementsSize(rModelPart);

    // Iterating over the elements
    ElementsArrayType& r_elements_array = rModelPart.Elements();
    const auto it_elem_begin= r_elements_array.begin();
    const int number_of_elements = static_cast<int>(r_elements_array.size());

    #pragma omp parallel for
    for (int i = 0; i < number_of_elements; ++i) {
        auto it_elem = it_elem_begin + i;

        if (it_elem->IsNot(BLOCKED)) {
            // Getting ELEMENT_H
            const double element_h = it_elem->GetValue(ELEMENT_H);

            // Blocking elements in the threshold sizes
            if (element_h <= minimal_size || element_h >= maximal_size) {
                it_elem->Set(BLOCKED, true);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ComputeElementsSize(ModelPart& rModelPart)
{
    // Iterating over the elements
    ElementsArrayType& r_elements_array = rModelPart.Elements();
    const auto it_elem_begin= r_elements_array.begin();
    const int number_of_elements = static_cast<int>(r_elements_array.size());

    #pragma omp parallel for
    for (int i = 0; i < number_of_elements; ++i) {
        auto it_elem = it_elem_begin + i;
        ComputeElementSize(it_elem);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ComputeElementSize(ElementItType itElement)
{
    auto& this_geometry = itElement->GetGeometry();

    // Here we compute the element size. This process is designed for triangles and tetrahedra, so we only specify for this geometries. Otherwise we take the length (and we throw a warning)
    if (this_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3){ // Triangular elements
        itElement->SetValue(ELEMENT_H, 2.0 * this_geometry.Circumradius());
    } else if(this_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4){ // Tetrahedral elements
        itElement->SetValue(ELEMENT_H,std::pow(12.0 * this_geometry.Volume()/std::sqrt(2.0), 1.0/3.0));
    } else { // In any othe case just considers the length of the element
        KRATOS_WARNING("MetricErrorProcess") << "This process is designed for tetrahedra (3D) and triangles (2D). Error expected" << std::endl;
        itElement->SetValue(ELEMENT_H, this_geometry.Length());
    }
}

} // namespace MeshingUtilities
} // namespace Kratos
