// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "custom_processes/assign_nodal_elements_to_nodes.h"
#include "custom_elements/nodal_concentrated_element.h"
#include "structural_mechanics_application_variables.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "utilities/variable_utils.h"

namespace Kratos
{
AssignNodalElementsToNodes::AssignNodalElementsToNodes(
        ModelPart& rThisModelPart,
        Parameters ThisParameters
        ):mrThisModelPart(rThisModelPart),
          mThisParameters(ThisParameters)
{
    KRATOS_TRY

    // We create a list of parameters to partially validate, and we validate them
    Parameters to_validate_parameters = Parameters(R"({})" );
    if (mThisParameters.Has("model_part_name"))
        to_validate_parameters.AddValue("model_part_name", mThisParameters["model_part_name"]);
    if (mThisParameters.Has("rayleigh_damping"))
        to_validate_parameters.AddValue("rayleigh_damping", mThisParameters["rayleigh_damping"]);
    if (mThisParameters.Has("assign_active_flag_node"))
        to_validate_parameters.AddValue("assign_active_flag_node", mThisParameters["assign_active_flag_node"]);
    if (mThisParameters.Has("interval"))
        to_validate_parameters.AddValue("interval", mThisParameters["interval"]);

    Parameters default_parameters = Parameters(R"(
    {
        "model_part_name"                : "",
        "rayleigh_damping"               : false,
        "assign_active_flag_node"        : true,
        "interval"                       : [0.0, 1e30]
    })" );

    to_validate_parameters.ValidateAndAssignDefaults(default_parameters);
    if (mThisParameters.Has("model_part_name"))
        mThisParameters.SetValue("model_part_name", to_validate_parameters["model_part_name"]);
    else
        mThisParameters.AddValue("model_part_name", to_validate_parameters["model_part_name"]);
    if (mThisParameters.Has("rayleigh_damping"))
        mThisParameters.SetValue("rayleigh_damping", to_validate_parameters["rayleigh_damping"]);
    else
        mThisParameters.AddValue("rayleigh_damping", to_validate_parameters["rayleigh_damping"]);
    if (mThisParameters.Has("assign_active_flag_node"))
        mThisParameters.SetValue("assign_active_flag_node", to_validate_parameters["assign_active_flag_node"]);
    else
        mThisParameters.AddValue("assign_active_flag_node", to_validate_parameters["assign_active_flag_node"]);
    if (mThisParameters.Has("interval"))
        mThisParameters.SetValue("interval", to_validate_parameters["interval"]);
    else
        mThisParameters.AddValue("interval", to_validate_parameters["interval"]);

    // List of auxiliar parameters to assign in case not defined
    Parameters auxiliar_parameters = Parameters(R"(
    {
        "nodal_mass"                     : null,
        "nodal_inertia"                  : [null, null, null],
        "nodal_stiffness"                : [null, null, null],
        "nodal_rotational_stiffness"     : [null, null, null],
        "nodal_damping_ratio"            : [null, null, null],
        "nodal_rotational_damping_ratio" : [null, null, null]
    })" );

    if (!mThisParameters.Has("nodal_mass")) mThisParameters.AddValue("nodal_mass", auxiliar_parameters["nodal_mass"]);
    if (!mThisParameters.Has("nodal_inertia")) mThisParameters.AddValue("nodal_inertia", auxiliar_parameters["nodal_inertia"]);
    if (!mThisParameters.Has("nodal_stiffness")) mThisParameters.AddValue("nodal_stiffness", auxiliar_parameters["nodal_stiffness"]);
    if (!mThisParameters.Has("nodal_rotational_stiffness")) mThisParameters.AddValue("nodal_rotational_stiffness", auxiliar_parameters["nodal_rotational_stiffness"]);
    if (!mThisParameters.Has("nodal_damping_ratio")) mThisParameters.AddValue("nodal_damping_ratio", auxiliar_parameters["nodal_damping_ratio"]);
    if (!mThisParameters.Has("nodal_rotational_damping_ratio")) mThisParameters.AddValue("nodal_rotational_damping_ratio", auxiliar_parameters["nodal_rotational_damping_ratio"]);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void AssignNodalElementsToNodes::Execute()
{
    // We execute the different steps of the process
    ExecuteInitialize();
    ExecuteInitializeSolutionStep();
}

/***********************************************************************************/
/***********************************************************************************/

void AssignNodalElementsToNodes::ExecuteInitialize()
{
    KRATOS_TRY

    // Assuming the number of properties is ordered
    ModelPart& r_root_model_part = mrThisModelPart.GetRootModelPart();
    const SizeType number_properties = r_root_model_part.NumberOfProperties();
    Properties::Pointer p_properties = r_root_model_part.pGetProperties(number_properties + 1);

    // Domain size
    const SizeType domain_size = r_root_model_part.GetProcessInfo()[DOMAIN_SIZE];

    // We get the proper model part
    const std::string& model_part_name = mThisParameters["model_part_name"].GetString();
    ModelPart& r_model_part = (model_part_name == "") ? mrThisModelPart : mrThisModelPart.GetSubModelPart(model_part_name);
    r_model_part.AddProperties(p_properties);

    // We assign values for the not null properties
    if (!mThisParameters["nodal_mass"].IsNull())
            p_properties->SetValue(NODAL_MASS, mThisParameters["nodal_mass"].GetDouble());
    if (!mThisParameters["nodal_inertia"][0].IsNull() || !mThisParameters["nodal_inertia"][1].IsNull() || !mThisParameters["nodal_inertia"][2].IsNull()) {
        array_1d<double, 3> nodal_inertia(3, 0.0);
        if (!mThisParameters["nodal_inertia"][0].IsNull())
            nodal_inertia[0] = mThisParameters["nodal_inertia"][0].GetDouble();
        if (!mThisParameters["nodal_inertia"][1].IsNull())
            nodal_inertia[1] = mThisParameters["nodal_inertia"][1].GetDouble();
        if (!mThisParameters["nodal_inertia"][2].IsNull())
            nodal_inertia[2] = mThisParameters["nodal_inertia"][2].GetDouble();

        p_properties->SetValue(NODAL_INERTIA, nodal_inertia);
    }
    if (!mThisParameters["nodal_stiffness"][0].IsNull() || !mThisParameters["nodal_stiffness"][1].IsNull() || !mThisParameters["nodal_stiffness"][2].IsNull()) {
        array_1d<double, 3> nodal_stiffness(3, 0.0);
        if (!mThisParameters["nodal_stiffness"][0].IsNull())
            nodal_stiffness[0] = mThisParameters["nodal_stiffness"][0].GetDouble();
        if (!mThisParameters["nodal_stiffness"][1].IsNull())
            nodal_stiffness[1] = mThisParameters["nodal_stiffness"][1].GetDouble();
        if (!mThisParameters["nodal_stiffness"][2].IsNull())
            nodal_stiffness[2] = mThisParameters["nodal_stiffness"][2].GetDouble();

        p_properties->SetValue(NODAL_STIFFNESS, nodal_stiffness);
    }
    if (!mThisParameters["nodal_rotational_stiffness"][0].IsNull() || !mThisParameters["nodal_rotational_stiffness"][1].IsNull() || !mThisParameters["nodal_rotational_stiffness"][2].IsNull()) {
        array_1d<double, 3> nodal_rotational_stiffness(3, 0.0);
        if (!mThisParameters["nodal_rotational_stiffness"][0].IsNull())
            nodal_rotational_stiffness[0] = mThisParameters["nodal_rotational_stiffness"][0].GetDouble();
        if (!mThisParameters["nodal_rotational_stiffness"][1].IsNull())
            nodal_rotational_stiffness[1] = mThisParameters["nodal_rotational_stiffness"][1].GetDouble();
        if (!mThisParameters["nodal_rotational_stiffness"][2].IsNull())
            nodal_rotational_stiffness[2] = mThisParameters["nodal_rotational_stiffness"][2].GetDouble();

        p_properties->SetValue(NODAL_ROTATIONAL_STIFFNESS, nodal_rotational_stiffness);
    }
    if (!mThisParameters["nodal_damping_ratio"][0].IsNull() || !mThisParameters["nodal_damping_ratio"][1].IsNull() || !mThisParameters["nodal_damping_ratio"][2].IsNull()) {
        array_1d<double, 3> nodal_damping_ratio(3, 0.0);
        if (!mThisParameters["nodal_damping_ratio"][0].IsNull())
            nodal_damping_ratio[0] = mThisParameters["nodal_damping_ratio"][0].GetDouble();
        if (!mThisParameters["nodal_damping_ratio"][1].IsNull())
            nodal_damping_ratio[1] = mThisParameters["nodal_damping_ratio"][1].GetDouble();
        if (!mThisParameters["nodal_damping_ratio"][2].IsNull())
            nodal_damping_ratio[2] = mThisParameters["nodal_damping_ratio"][2].GetDouble();

        p_properties->SetValue(NODAL_DAMPING_RATIO, nodal_damping_ratio);
    }
    if (!mThisParameters["nodal_rotational_damping_ratio"][0].IsNull() || !mThisParameters["nodal_rotational_damping_ratio"][1].IsNull() || !mThisParameters["nodal_rotational_damping_ratio"][2].IsNull()) {
        array_1d<double, 3> nodal_rotational_damping_ratio(3, 0.0);
        if (!mThisParameters["nodal_rotational_damping_ratio"][0].IsNull())
            nodal_rotational_damping_ratio[0] = mThisParameters["nodal_rotational_damping_ratio"][0].GetDouble();
        if (!mThisParameters["nodal_rotational_damping_ratio"][1].IsNull())
            nodal_rotational_damping_ratio[1] = mThisParameters["nodal_rotational_damping_ratio"][1].GetDouble();
        if (!mThisParameters["nodal_rotational_damping_ratio"][2].IsNull())
            nodal_rotational_damping_ratio[2] = mThisParameters["nodal_rotational_damping_ratio"][2].GetDouble();

        p_properties->SetValue(NODAL_ROTATIONAL_DAMPING_RATIO, nodal_rotational_damping_ratio);
    }

    // The number of elements
    const SizeType number_elements = r_root_model_part.NumberOfElements();
    const auto it_elem_begin = r_root_model_part.ElementsBegin();

    // Reorder ids
    #pragma omp parallel for
    for(int i=0; i< static_cast<int>(number_elements); ++i) {
        auto it_elem = it_elem_begin + i;
        it_elem->SetId(i + 1);
    }

    // We get the reference element
    const bool rayleigh_damping = mThisParameters["rayleigh_damping"].GetBool();
    const bool assign_active_flag_node = mThisParameters["assign_active_flag_node"].GetBool();
    PointerVector<NodeType> aux_node_array(1);
    const auto it_node_begin = r_model_part.NodesBegin();
    aux_node_array(0) = *(it_node_begin).base();

    const SizeType number_of_nodes = r_model_part.Nodes().size();
    std::vector<Element::Pointer> auxiliar_elements_vector;

    GeometryType::Pointer p_dummy_geom = GetPointGeometryFromNode(aux_node_array, domain_size);
    const Element& rReferenceElement = NodalConcentratedElement(0, p_dummy_geom, rayleigh_damping, assign_active_flag_node);

    #pragma omp parallel
    {
        // Buffer for new elements if created
        std::vector<Element::Pointer> auxiliar_elements_vector_buffer;

        #pragma omp for
        for(int i=0; i< static_cast<int>(number_of_nodes); ++i) {
            auto it_node = it_node_begin + i;

            PointerVector<NodeType> this_node_array(1);
            this_node_array(0) = *(it_node).base();

            auto p_element = rReferenceElement.Create(number_elements + 1 + i, GetPointGeometryFromNode(this_node_array, domain_size), p_properties);
            auxiliar_elements_vector_buffer.push_back(p_element);

            // Deep copy elemental data and flags
//             p_element->Data() = it_node->Data();
            p_element->Set(Flags(*it_node));
        }

        // Combine buffers together
        #pragma omp critical
        {
            std::move(auxiliar_elements_vector_buffer.begin(),auxiliar_elements_vector_buffer.end(),back_inserter(auxiliar_elements_vector));
        }
    }

    // Adding to the model part
    ElementsArrayType aux_elems;
    aux_elems.GetContainer() = auxiliar_elements_vector;
    r_model_part.AddElements(aux_elems.begin(), aux_elems.end());

    // We Initialize the elements
    InitializeElements(r_model_part);

    // Set the flag ACTIVE
    this->Set(ACTIVE, false);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void AssignNodalElementsToNodes::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    // We get the proper model part
    const std::string& model_part_name = mThisParameters["model_part_name"].GetString();
    ModelPart& r_model_part = (model_part_name == "") ? mrThisModelPart : mrThisModelPart.GetSubModelPart(model_part_name);
    const auto it_elem_begin = r_model_part.Elements().begin();

    // Check the interval
    if (mThisParameters["interval"][0].GetDouble() > 0.0 || mThisParameters["interval"][1].GetDouble() < 1e30) {
        if (this->IsNot(ACTIVE)) {
            // Initialize initial displacement and rotation
            #pragma omp parallel for
            for(int i=0; i< static_cast<int>(r_model_part.Elements().size()); ++i) {
                auto it_elem = it_elem_begin + i;
                if (it_elem->Has(INITIAL_DISPLACEMENT)) {
                    it_elem->SetValue(INITIAL_DISPLACEMENT, it_elem->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT));
                }
                if (it_elem->Has(INITIAL_ROTATION)) {
                    it_elem->SetValue(INITIAL_ROTATION, it_elem->GetGeometry()[0].FastGetSolutionStepValue(ROTATION));
                }
            }
            // Set the flag ACTIVE
            VariableUtils().SetFlag(ACTIVE, true, r_model_part.Elements());
            this->Set(ACTIVE, true);
        }
    } else {
        if (this->Is(ACTIVE)) {
            // Set the flag ACTIVE
            VariableUtils().SetFlag(ACTIVE, false, r_model_part.Elements());
            this->Set(ACTIVE, false);
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void AssignNodalElementsToNodes::InitializeElements(ModelPart& rModelPart)
{
    ElementsArrayType& r_elements_array = rModelPart.Elements();
    const auto it_elem_begin = r_elements_array.begin();
    #pragma omp parallel for
    for(int i=0; i< static_cast<int>(r_elements_array.size()); i++)
        (it_elem_begin + i)->Initialize();

    // Inactive by default
    VariableUtils().SetFlag(ACTIVE, false, rModelPart.Elements());
}

/***********************************************************************************/
/***********************************************************************************/

AssignNodalElementsToNodes::GeometryType::Pointer AssignNodalElementsToNodes::GetPointGeometryFromNode(
    PointerVector<NodeType>& rArrayNodes,
    const SizeType Dimension
    )
{
    if (Dimension == 2) {
        return Kratos::make_shared<Point2D<NodeType>>(rArrayNodes);
    } else {
        return Kratos::make_shared<Point3D<NodeType>>(rArrayNodes);
    }

    return nullptr;
}

// class AssignNodalElementsToNodes
} // namespace Kratos
