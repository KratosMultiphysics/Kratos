// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "custom_processes/assign_nodal_elements_to_nodes_process.h"
#include "custom_elements/nodal_elements/nodal_concentrated_element.h"
#include "structural_mechanics_application_variables.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "utilities/variable_utils.h"
#include "utilities/entities_utilities.h"

namespace Kratos
{

AssignNodalElementsToNodesProcess::AssignNodalElementsToNodesProcess(
    Model& rModel,
    Parameters ThisParameters
    ) : AssignNodalElementsToNodesProcess(rModel.GetModelPart(ThisParameters["model_part_name"].GetString()), ThisParameters)
{

}

/***********************************************************************************/
/***********************************************************************************/

AssignNodalElementsToNodesProcess::AssignNodalElementsToNodesProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ) : mrThisModelPart(rThisModelPart),
        mThisParameters(ThisParameters)
{
    KRATOS_TRY

    // We create a list of parameters to partially validate, and we validate them
    Parameters to_validate_parameters = Parameters(R"({})" );
    if (mThisParameters.Has("rayleigh_damping")) {
        to_validate_parameters.AddValue("rayleigh_damping", mThisParameters["rayleigh_damping"]);
    }
    if (mThisParameters.Has("interval")) {
        to_validate_parameters.AddValue("interval", mThisParameters["interval"]);
    }

    to_validate_parameters.ValidateAndAssignDefaults(GetDefaultParameters());
    if (mThisParameters.Has("rayleigh_damping")) {
        mThisParameters.SetValue("rayleigh_damping", to_validate_parameters["rayleigh_damping"]);
    } else {
        mThisParameters.AddValue("rayleigh_damping", to_validate_parameters["rayleigh_damping"]);
    }
    if (mThisParameters.Has("interval")) {
        mThisParameters.SetValue("interval", to_validate_parameters["interval"]);
    } else {
        mThisParameters.AddValue("interval", to_validate_parameters["interval"]);
    }

    // List of auxiliary parameters to assign in case not defined
    Parameters auxiliary_parameters = Parameters(R"(
    {
        "nodal_mass"                     : null,
        "nodal_inertia"                  : [null, null, null],
        "nodal_stiffness"                : [null, null, null],
        "nodal_rotational_stiffness"     : [null, null, null],
        "nodal_damping_ratio"            : [null, null, null],
        "nodal_rotational_damping_ratio" : [null, null, null]
    })" );

    if (!mThisParameters.Has("nodal_mass")) mThisParameters.AddValue("nodal_mass", auxiliary_parameters["nodal_mass"]);
    if (!mThisParameters.Has("nodal_inertia")) mThisParameters.AddValue("nodal_inertia", auxiliary_parameters["nodal_inertia"]);
    if (!mThisParameters.Has("nodal_stiffness")) mThisParameters.AddValue("nodal_stiffness", auxiliary_parameters["nodal_stiffness"]);
    if (!mThisParameters.Has("nodal_rotational_stiffness")) mThisParameters.AddValue("nodal_rotational_stiffness", auxiliary_parameters["nodal_rotational_stiffness"]);
    if (!mThisParameters.Has("nodal_damping_ratio")) mThisParameters.AddValue("nodal_damping_ratio", auxiliary_parameters["nodal_damping_ratio"]);
    if (!mThisParameters.Has("nodal_rotational_damping_ratio")) mThisParameters.AddValue("nodal_rotational_damping_ratio", auxiliary_parameters["nodal_rotational_damping_ratio"]);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void AssignNodalElementsToNodesProcess::Execute()
{
    // We execute the different steps of the process
    ExecuteInitialize();
    ExecuteInitializeSolutionStep();
}

/***********************************************************************************/
/***********************************************************************************/

void AssignNodalElementsToNodesProcess::ExecuteInitialize()
{
    KRATOS_TRY

    // Assuming the number of properties is ordered
    ModelPart& r_root_model_part = mrThisModelPart.GetRootModelPart();
    const SizeType number_properties = r_root_model_part.NumberOfProperties();
    Properties::Pointer p_properties = r_root_model_part.HasProperties(number_properties + 1) ? r_root_model_part.pGetProperties(number_properties + 1) : r_root_model_part.CreateNewProperties(number_properties + 1);

    // Domain size
    const SizeType domain_size = r_root_model_part.GetProcessInfo()[DOMAIN_SIZE];

    // We assign the properties to the model part
    if (mrThisModelPart.Name() != r_root_model_part.Name()) {
        mrThisModelPart.AddProperties(p_properties);
    }

    // We assign values for the not null properties
    double nodal_mass = 0.0;
    if (!mThisParameters["nodal_mass"].IsNull()) {
        nodal_mass = mThisParameters["nodal_mass"].GetDouble();
        p_properties->SetValue(NODAL_MASS, nodal_mass);
    }
    array_1d<double, 3> nodal_inertia = ZeroVector(3);
    if (!mThisParameters["nodal_inertia"][0].IsNull() || !mThisParameters["nodal_inertia"][1].IsNull() || !mThisParameters["nodal_inertia"][2].IsNull()) {
        if (!mThisParameters["nodal_inertia"][0].IsNull()) {
            nodal_inertia[0] = mThisParameters["nodal_inertia"][0].GetDouble();
        }
        if (!mThisParameters["nodal_inertia"][1].IsNull()) {
            nodal_inertia[1] = mThisParameters["nodal_inertia"][1].GetDouble();
        }
        if (!mThisParameters["nodal_inertia"][2].IsNull()) {
            nodal_inertia[2] = mThisParameters["nodal_inertia"][2].GetDouble();
        }

        p_properties->SetValue(NODAL_INERTIA, nodal_inertia);
    }
    array_1d<double, 3> nodal_stiffness = ZeroVector(3);
    if (!mThisParameters["nodal_stiffness"][0].IsNull() || !mThisParameters["nodal_stiffness"][1].IsNull() || !mThisParameters["nodal_stiffness"][2].IsNull()) {
        if (!mThisParameters["nodal_stiffness"][0].IsNull()) {
            nodal_stiffness[0] = mThisParameters["nodal_stiffness"][0].GetDouble();
        }
        if (!mThisParameters["nodal_stiffness"][1].IsNull()) {
            nodal_stiffness[1] = mThisParameters["nodal_stiffness"][1].GetDouble();
        }
        if (!mThisParameters["nodal_stiffness"][2].IsNull()) {
            nodal_stiffness[2] = mThisParameters["nodal_stiffness"][2].GetDouble();
        }

        p_properties->SetValue(NODAL_DISPLACEMENT_STIFFNESS, nodal_stiffness);
    }
    array_1d<double, 3> nodal_rotational_stiffness = ZeroVector(3);
    if (!mThisParameters["nodal_rotational_stiffness"][0].IsNull() || !mThisParameters["nodal_rotational_stiffness"][1].IsNull() || !mThisParameters["nodal_rotational_stiffness"][2].IsNull()) {
        if (!mThisParameters["nodal_rotational_stiffness"][0].IsNull()) {
            nodal_rotational_stiffness[0] = mThisParameters["nodal_rotational_stiffness"][0].GetDouble();
        }
        if (!mThisParameters["nodal_rotational_stiffness"][1].IsNull()) {
            nodal_rotational_stiffness[1] = mThisParameters["nodal_rotational_stiffness"][1].GetDouble();
        }
        if (!mThisParameters["nodal_rotational_stiffness"][2].IsNull()) {
            nodal_rotational_stiffness[2] = mThisParameters["nodal_rotational_stiffness"][2].GetDouble();
        }

        p_properties->SetValue(NODAL_ROTATIONAL_STIFFNESS, nodal_rotational_stiffness);
    }
    array_1d<double, 3> nodal_damping_ratio = ZeroVector(3);
    if (!mThisParameters["nodal_damping_ratio"][0].IsNull() || !mThisParameters["nodal_damping_ratio"][1].IsNull() || !mThisParameters["nodal_damping_ratio"][2].IsNull()) {
        if (!mThisParameters["nodal_damping_ratio"][0].IsNull()) {
            nodal_damping_ratio[0] = mThisParameters["nodal_damping_ratio"][0].GetDouble();
        }
        if (!mThisParameters["nodal_damping_ratio"][1].IsNull()) {
            nodal_damping_ratio[1] = mThisParameters["nodal_damping_ratio"][1].GetDouble();
        }
        if (!mThisParameters["nodal_damping_ratio"][2].IsNull()) {
            nodal_damping_ratio[2] = mThisParameters["nodal_damping_ratio"][2].GetDouble();
        }

        p_properties->SetValue(NODAL_DAMPING_RATIO, nodal_damping_ratio);
    }
    array_1d<double, 3> nodal_rotational_damping_ratio = ZeroVector(3);
    if (!mThisParameters["nodal_rotational_damping_ratio"][0].IsNull() || !mThisParameters["nodal_rotational_damping_ratio"][1].IsNull() || !mThisParameters["nodal_rotational_damping_ratio"][2].IsNull()) {
        if (!mThisParameters["nodal_rotational_damping_ratio"][0].IsNull()) {
            nodal_rotational_damping_ratio[0] = mThisParameters["nodal_rotational_damping_ratio"][0].GetDouble();
        }
        if (!mThisParameters["nodal_rotational_damping_ratio"][1].IsNull()) {
            nodal_rotational_damping_ratio[1] = mThisParameters["nodal_rotational_damping_ratio"][1].GetDouble();
        }
        if (!mThisParameters["nodal_rotational_damping_ratio"][2].IsNull()) {
            nodal_rotational_damping_ratio[2] = mThisParameters["nodal_rotational_damping_ratio"][2].GetDouble();
        }

        p_properties->SetValue(NODAL_ROTATIONAL_DAMPING_RATIO, nodal_rotational_damping_ratio);
    }

    // The number of elements
    const SizeType number_elements = r_root_model_part.NumberOfElements();
    const auto it_elem_begin = r_root_model_part.ElementsBegin();

    // Reorder ids
    IndexPartition<std::size_t>(number_elements).for_each([&it_elem_begin](std::size_t i) {
        auto it_elem = it_elem_begin + i;
        it_elem->SetId(i + 1);
    });

    // We get the reference element
    const bool rayleigh_damping = mThisParameters["rayleigh_damping"].GetBool();
    if (rayleigh_damping) {
        p_properties->SetValue(CONSIDER_RAYLEIGH_DAMPING, true);
    }
    PointerVector<Node> aux_node_array(1);
    const auto it_node_begin = mrThisModelPart.NodesBegin();
    aux_node_array(0) = *(it_node_begin).base();

    const SizeType number_of_nodes = mrThisModelPart.Nodes().size();

    GeometryType::Pointer p_dummy_geom = GetPointGeometryFromNode(aux_node_array, domain_size);
    const Element& r_reference_element = NodalConcentratedElement(0, p_dummy_geom);

    std::vector<Element::Pointer> auxiliary_elements_vector  = IndexPartition<std::size_t>(number_of_nodes).for_each<AccumReduction<Element::Pointer>>([&](std::size_t i) {
        auto it_node = it_node_begin + i;

        PointerVector<Node> this_node_array(1);
        this_node_array(0) = *(it_node).base();

        auto p_element = r_reference_element.Create(number_elements + 1 + i, GetPointGeometryFromNode(this_node_array, domain_size), p_properties);

        // Deep copy elemental flags
        p_element->Set(Flags(*it_node));

        return p_element;
    });

    // Adding to the model part
    ElementsArrayType aux_elems;
    aux_elems.GetContainer() = auxiliary_elements_vector;
    mrThisModelPart.AddElements(aux_elems.begin(), aux_elems.end());

    // We Initialize the elements
    EntitiesUtilities::InitializeEntities<Element>(mrThisModelPart);

    // Inactive by default
    VariableUtils().SetFlag(ACTIVE, false, mrThisModelPart.Elements());

    // Set the flag ACTIVE
    this->Set(ACTIVE, false);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void AssignNodalElementsToNodesProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    // Check the interval
    if (mThisParameters["interval"][0].GetDouble() > 0.0 || mThisParameters["interval"][1].GetDouble() < 1e30) {
        if (this->IsNot(ACTIVE)) {
            // Initialize initial displacement and rotation
            block_for_each(mrThisModelPart.Elements(), [&](Element& rElement) {
                if (rElement.Has(NODAL_INITIAL_DISPLACEMENT)) {
                    rElement.SetValue(NODAL_INITIAL_DISPLACEMENT, rElement.GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT));
                }
                if (rElement.Has(NODAL_INITIAL_ROTATION)) {
                    rElement.SetValue(NODAL_INITIAL_ROTATION, rElement.GetGeometry()[0].FastGetSolutionStepValue(ROTATION));
                }
            });
            // Set the flag ACTIVE
            VariableUtils().SetFlag(ACTIVE, true, mrThisModelPart.Elements());
            this->Set(ACTIVE, true);
        }
    } else {
        if (this->Is(ACTIVE)) {
            // Set the flag ACTIVE
            VariableUtils().SetFlag(ACTIVE, false, mrThisModelPart.Elements());
            this->Set(ACTIVE, false);
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

AssignNodalElementsToNodesProcess::GeometryType::Pointer AssignNodalElementsToNodesProcess::GetPointGeometryFromNode(
    PointerVector<Node>& rArrayNodes,
    const SizeType Dimension
    )
{
    if (Dimension == 2) {
        return Kratos::make_shared<Point2D<Node>>(rArrayNodes);
    } else {
        return Kratos::make_shared<Point3D<Node>>(rArrayNodes);
    }

    return nullptr;
}

// class AssignNodalElementsToNodesProcess
} // namespace Kratos
