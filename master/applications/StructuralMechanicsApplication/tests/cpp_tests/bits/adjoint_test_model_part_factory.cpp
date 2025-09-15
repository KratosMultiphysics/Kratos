#include "tests/cpp_tests/bits/adjoint_test_model_part_factory.h"

#include "containers/model.h"
#include "structural_mechanics_application_variables.h"
#include "utilities/compare_elements_and_conditions_utility.h"
#include "utilities/variable_utils.h"
#include <string>
#include <unordered_map>

namespace Kratos
{
namespace
{ // cpp internals
namespace adjoint_test_model_part_factory_cpp
{ // unity build unity guard
void AddVariables(ModelPart* pAdjointModelPart, const VariablesList& rCustomVariables = {});
void CopyNodes(ModelPart* pModelPart, const PointerVectorSet<Node, IndexedObject>& rNodes);
void CopyProperties(ModelPart* pModelPart,
                    const PointerVectorSet<Properties, IndexedObject>& rProperties);
void CreateAdjointElements(ModelPart* pAdjointModelPart,
                           const PointerVectorSet<Element, IndexedObject>& rPrimalElements);
void AddDofs(ModelPart* pAdjointModelPart);
void AssignBCs(ModelPart* pAdjointModelPart, const ModelPart& rPrimalModelPart);
void CopySolutionStepData(ModelPart* pDestinationModelPart,
                          const PointerVectorSet<Node, IndexedObject>& rNodes);
} // namespace adjoint_test_model_part_factory_cpp
} // namespace

ModelPart& CreateStructuralMechanicsAdjointTestModelPart(ModelPart* pPrimalModelPart)
{
    using namespace adjoint_test_model_part_factory_cpp;
    const std::string name = pPrimalModelPart->Name() + "(Adjoint)";
    if (pPrimalModelPart->GetModel().HasModelPart(name))
    {
        pPrimalModelPart->GetModel().DeleteModelPart(name);
    }
    ModelPart& adjoint_model_part = pPrimalModelPart->GetModel().CreateModelPart(name);
    AddVariables(&adjoint_model_part, pPrimalModelPart->GetNodalSolutionStepVariablesList());
    CopyNodes(&adjoint_model_part, pPrimalModelPart->Nodes());
    ProcessInfo& r_process_info = adjoint_model_part.GetProcessInfo();
    r_process_info = pPrimalModelPart->GetProcessInfo();
    CopyProperties(&adjoint_model_part, pPrimalModelPart->rProperties());
    CreateAdjointElements(&adjoint_model_part, pPrimalModelPart->Elements());
    adjoint_model_part.SetBufferSize(pPrimalModelPart->GetBufferSize());
    // initialize adjoint time
    adjoint_model_part.CloneTimeStep(r_process_info[TIME] + r_process_info[DELTA_TIME]);
    AddDofs(&adjoint_model_part);
    AssignBCs(&adjoint_model_part, *pPrimalModelPart);
    VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, true,
                                             adjoint_model_part.Nodes());
    CopySolutionStepData(&adjoint_model_part, pPrimalModelPart->Nodes());
    return adjoint_model_part;
}

namespace
{ // cpp internals
namespace adjoint_test_model_part_factory_cpp
{ // unity build unity guard
void AddVariables(ModelPart* pAdjointModelPart, const VariablesList& rCustomVariables)
{
    pAdjointModelPart->GetNodalSolutionStepVariablesList() = rCustomVariables;
    pAdjointModelPart->AddNodalSolutionStepVariable(ADJOINT_DISPLACEMENT);
    pAdjointModelPart->AddNodalSolutionStepVariable(ADJOINT_VECTOR_2);
    pAdjointModelPart->AddNodalSolutionStepVariable(ADJOINT_VECTOR_3);
    pAdjointModelPart->AddNodalSolutionStepVariable(AUX_ADJOINT_VECTOR_1);
    pAdjointModelPart->AddNodalSolutionStepVariable(SHAPE_SENSITIVITY);
}

void CopyNodes(ModelPart* pModelPart, const PointerVectorSet<Node, IndexedObject>& rNodes)
{
    for (const auto& r_node : rNodes)
    {
        auto p_node = pModelPart->CreateNewNode(r_node.Id(), r_node.X0(),
                                                r_node.Y0(), r_node.Z0());
        p_node->Coordinates() = r_node.Coordinates();
    }
}

void CopyProperties(ModelPart* pModelPart,
                    const PointerVectorSet<Properties, IndexedObject>& rProperties)
{
    for (const auto& r_prop : rProperties)
    {
        *pModelPart->CreateNewProperties(r_prop.Id()) = r_prop;
    }
}

const std::string* AdjointName(const std::string& rPrimalName)
{
    // extend this map when adding new adjoint elements for testing
    const static std::unordered_map<std::string, std::string> m = {
        {"TotalLagrangianElement2D3N", "TotalLagrangianAdjointElement2D3N"},
        {"TotalLagrangianElement2D4N", "TotalLagrangianAdjointElement2D4N"},
        {"TotalLagrangianElement2D6N", "TotalLagrangianAdjointElement2D6N"},
        {"TotalLagrangianElement3D8N", "TotalLagrangianAdjointElement3D8N"},
        {"TotalLagrangianElement3D4N", "TotalLagrangianAdjointElement3D4N"}};

    auto search = m.find(rPrimalName);
    return (search != m.end()) ? &search->second : nullptr;
}

const std::string* AdjointName(const Element& rPrimalElement)
{
    std::string primal_name;
    // This is expensive, don't call it in performance critical loops
    CompareElementsAndConditionsUtility::GetRegisteredName(rPrimalElement, primal_name);
    return AdjointName(primal_name);
}

void CreateAdjointElements(ModelPart* pAdjointModelPart,
                           const PointerVectorSet<Element, IndexedObject>& rPrimalElements)
{
    std::vector<std::size_t> node_ids;
    for (const auto& r_elem : rPrimalElements)
    {
        auto p_prop = pAdjointModelPart->pGetProperties(r_elem.GetProperties().Id());
        const auto& r_points = r_elem.GetGeometry().Points();
        node_ids.resize(r_points.size());
        for (std::size_t i = 0; i < r_points.size(); ++i)
        {
            node_ids[i] = r_points[i].Id();
        }
        const std::string* p_adjoint_name = AdjointName(r_elem);
        KRATOS_ERROR_IF_NOT(p_adjoint_name)
            << "Adjoint name not found for element #" << r_elem.Id() << "!";
        pAdjointModelPart->CreateNewElement(*p_adjoint_name, r_elem.Id(), node_ids, p_prop);
    }
}

void AddDofs(ModelPart* pAdjointModelPart)
{
    for (auto& r_node : pAdjointModelPart->Nodes())
    {
        r_node.AddDof(ADJOINT_DISPLACEMENT_X);
        r_node.AddDof(ADJOINT_DISPLACEMENT_Y);
        r_node.AddDof(ADJOINT_DISPLACEMENT_Z);
    }
}

void AssignBCs(ModelPart* pAdjointModelPart, const ModelPart& rPrimalModelPart)
{
    for (auto& r_adjoint_node : pAdjointModelPart->Nodes())
    {
        const auto& r_primal_node = rPrimalModelPart.GetNode(r_adjoint_node.Id());
        if (r_primal_node.IsFixed(DISPLACEMENT_X))
            r_adjoint_node.Fix(ADJOINT_DISPLACEMENT_X);
        if (r_primal_node.IsFixed(DISPLACEMENT_Y))
            r_adjoint_node.Fix(ADJOINT_DISPLACEMENT_Y);
        if (r_primal_node.IsFixed(DISPLACEMENT_Z))
            r_adjoint_node.Fix(ADJOINT_DISPLACEMENT_Z);
    }
}

void CopySolutionStepData(ModelPart* pDestinationModelPart,
                          const PointerVectorSet<Node, IndexedObject>& rNodes)
{
    KRATOS_TRY;
    for (const auto& r_primal_node : rNodes)
    {
        const auto p_variables_list = r_primal_node.pGetVariablesList();
        auto& r_adjoint_node = pDestinationModelPart->GetNode(r_primal_node.Id());
        const std::size_t queue_size = r_primal_node.SolutionStepData().QueueSize();

        for (const auto& r_variable_data : *p_variables_list)
        {
            if (KratosComponents<Variable<array_1d<double, 3>>>::Has(
                    r_variable_data.Name()))
            {
                const auto& r_variable =
                    KratosComponents<Variable<array_1d<double, 3>>>::Get(
                        r_variable_data.Name());
                for (std::size_t i = 0; i < queue_size; ++i)
                    r_adjoint_node.FastGetSolutionStepValue(r_variable, i) =
                        r_primal_node.FastGetSolutionStepValue(r_variable, i);
            }
            else if (KratosComponents<Variable<double>>::Has(r_variable_data.Name()))
            {
                const auto& r_variable =
                    KratosComponents<Variable<double>>::Get(r_variable_data.Name());
                for (std::size_t i = 0; i < queue_size; ++i)
                    r_adjoint_node.FastGetSolutionStepValue(r_variable, i) =
                        r_primal_node.FastGetSolutionStepValue(r_variable, i);
            }
            else
                KRATOS_ERROR << "Unsupported variable: " << r_variable_data.Name()
                             << std::endl;
        }
    }
    KRATOS_CATCH("");
}
} // namespace adjoint_test_model_part_factory_cpp
} // namespace
} // namespace Kratos
