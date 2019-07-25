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

// External includes

// Project includes
#include "custom_utilities/meshing_utilities.h"

namespace Kratos
{
namespace MeshingUtilities
{
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

/***********************************************************************************/
/***********************************************************************************/

void HexahedraMeshToTetrahedraMesh(
    ModelPart& rModelPart,
    Parameters ThisParameters
    )
{
    Parameters default_parameters = Parameters(R"(
    {
        "element_name"   : "Element3D4N",
        "condition_name" : "SurfaceCondition3D3N"
    })" );

    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // Getting values
    const std::string& r_element_name = ThisParameters["element_name"].GetString();
    const std::string& r_condition_name = ThisParameters["condition_name"].GetString();

    // Reference element
    Element const& r_clone_element = KratosComponents<Element>::Get(r_element_name);
    Condition const& r_clone_condition = KratosComponents<Condition>::Get(r_condition_name);

    // The root model part
    ModelPart& r_root_model_part = rModelPart.GetRootModelPart();

    // Iterating over the elements
    ElementsArrayType& r_elements_array = rModelPart.Elements();
    const auto it_elem_begin= r_elements_array.begin();
    const int number_of_elements = static_cast<int>(r_elements_array.size());
    const SizeType total_number_of_elements = r_root_model_part.Elements().size();

    #pragma omp parallel
    {
        // Array of nodes
        GeometryType::PointsArrayType elements_pnodes;

        // Buffer vector
        ElementsArrayType created_elements_vector;

        #pragma omp for
        for (int i = 0; i < number_of_elements; ++i) {
            auto it_elem = it_elem_begin + i;

            const auto& r_geometry = it_elem->GetGeometry();
            if (r_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8) {
                elements_pnodes.push_back(rModelPart.pGetNode(r_geometry[0].Id()));
                elements_pnodes.push_back(rModelPart.pGetNode(r_geometry[4].Id()));
                elements_pnodes.push_back(rModelPart.pGetNode(r_geometry[5].Id()));
                elements_pnodes.push_back(rModelPart.pGetNode(r_geometry[7].Id()));
                auto p_element_1 = r_clone_element.Create(total_number_of_elements + 1 + i * 5, elements_pnodes, it_elem->pGetProperties());
                created_elements_vector.push_back(p_element_1);
                elements_pnodes.clear();

                elements_pnodes.push_back(rModelPart.pGetNode(r_geometry[0].Id()));
                elements_pnodes.push_back(rModelPart.pGetNode(r_geometry[5].Id()));
                elements_pnodes.push_back(rModelPart.pGetNode(r_geometry[1].Id()));
                elements_pnodes.push_back(rModelPart.pGetNode(r_geometry[2].Id()));
                auto p_element_2 = r_clone_element.Create(total_number_of_elements + 2 + i * 5, elements_pnodes, it_elem->pGetProperties());
                created_elements_vector.push_back(p_element_2);
                elements_pnodes.clear();

                elements_pnodes.push_back(rModelPart.pGetNode(r_geometry[0].Id()));
                elements_pnodes.push_back(rModelPart.pGetNode(r_geometry[7].Id()));
                elements_pnodes.push_back(rModelPart.pGetNode(r_geometry[5].Id()));
                elements_pnodes.push_back(rModelPart.pGetNode(r_geometry[2].Id()));
                auto p_element_3 = r_clone_element.Create(total_number_of_elements + 3 + i * 5, elements_pnodes, it_elem->pGetProperties());
                created_elements_vector.push_back(p_element_3);
                elements_pnodes.clear();

                elements_pnodes.push_back(rModelPart.pGetNode(r_geometry[0].Id()));
                elements_pnodes.push_back(rModelPart.pGetNode(r_geometry[7].Id()));
                elements_pnodes.push_back(rModelPart.pGetNode(r_geometry[2].Id()));
                elements_pnodes.push_back(rModelPart.pGetNode(r_geometry[3].Id()));
                auto p_element_4 = r_clone_element.Create(total_number_of_elements + 4 + i * 5, elements_pnodes, it_elem->pGetProperties());
                created_elements_vector.push_back(p_element_4);
                elements_pnodes.clear();

                elements_pnodes.push_back(rModelPart.pGetNode(r_geometry[5].Id()));
                elements_pnodes.push_back(rModelPart.pGetNode(r_geometry[7].Id()));
                elements_pnodes.push_back(rModelPart.pGetNode(r_geometry[6].Id()));
                elements_pnodes.push_back(rModelPart.pGetNode(r_geometry[2].Id()));
                auto p_element_5 = r_clone_element.Create(total_number_of_elements + 5 + i * 5, elements_pnodes, it_elem->pGetProperties());
                created_elements_vector.push_back(p_element_5);
                elements_pnodes.clear();

                it_elem->Set(TO_ERASE);
            }
        }

        // Adding elements to model part
        #pragma omp critical
        {
            rModelPart.AddElements(created_elements_vector.begin(), created_elements_vector.end());
        }
    }

    // Iterating over the conditions
    ConditionsArrayType& r_conditions_array = rModelPart.Conditions();
    const auto it_cond_begin= r_conditions_array.begin();
    const int number_of_conditions = static_cast<int>(r_conditions_array.size());
    const SizeType total_number_of_conditions = r_root_model_part.Conditions().size();

    #pragma omp parallel
    {
        // Array of nodes
        GeometryType::PointsArrayType conditions_pnodes;

        // Buffer vector
        ConditionsArrayType created_conditions_vector;

        #pragma omp for
        for (int i = 0; i < number_of_conditions; ++i) {
            auto it_cond = it_cond_begin + i;

            const auto& r_geometry = it_cond->GetGeometry();
            if (r_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4) {
                conditions_pnodes.push_back(rModelPart.pGetNode(r_geometry[0].Id()));
                conditions_pnodes.push_back(rModelPart.pGetNode(r_geometry[1].Id()));
                conditions_pnodes.push_back(rModelPart.pGetNode(r_geometry[3].Id()));
                auto p_condition_1 = r_clone_condition.Create(total_number_of_conditions + 1 + i * 2, conditions_pnodes, it_cond->pGetProperties());
                created_conditions_vector.push_back(p_condition_1);
                conditions_pnodes.clear();

                conditions_pnodes.push_back(rModelPart.pGetNode(r_geometry[2].Id()));
                conditions_pnodes.push_back(rModelPart.pGetNode(r_geometry[3].Id()));
                conditions_pnodes.push_back(rModelPart.pGetNode(r_geometry[1].Id()));
                auto p_condition_2 = r_clone_condition.Create(total_number_of_conditions + 2 + i * 2, conditions_pnodes, it_cond->pGetProperties());
                created_conditions_vector.push_back(p_condition_2);
                conditions_pnodes.clear();

                it_cond->Set(TO_ERASE);
            }
        }

        // Adding conditions to model part
        #pragma omp critical
        {
            rModelPart.AddConditions(created_conditions_vector.begin(), created_conditions_vector.end());
        }
    }

    // Remove all elements and conditions
    rModelPart.RemoveConditionsFromAllLevels(TO_ERASE);
    rModelPart.RemoveElementsFromAllLevels(TO_ERASE);

    /* Reorder Ids */
    // Iterate over conditions
    auto& r_root_conditions_array = r_root_model_part.Conditions();
    const auto it_root_cond_begin = r_root_conditions_array.begin();
    for(IndexType i = 0; i < r_root_conditions_array.size(); ++i)
        (it_root_cond_begin + i)->SetId(i + 1);

    // Iterate over elements
    auto& r_root_elements_array = r_root_model_part.Elements();
    const auto it_root_elem_begin = r_root_elements_array.begin();
    for(IndexType i = 0; i < r_root_elements_array.size(); ++i)
        (it_root_elem_begin + i)->SetId(i + 1);
}

} // namespace MeshingUtilities
} // namespace Kratos
