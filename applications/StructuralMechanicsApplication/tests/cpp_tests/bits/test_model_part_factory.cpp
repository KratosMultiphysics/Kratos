#include "tests/cpp_tests/bits/test_model_part_factory.h"

#include "structural_mechanics_application_variables.h"

namespace Kratos
{
namespace
{ // cpp internals
namespace test_model_part_factory_cpp
{ // unity build unity guard
void AddVariables(ModelPart* pModelPart);
void CreateNodes(ModelPart* pModelPart, const Element& rElementPrototype);
void CreatePropertiesAndElements(ModelPart* pModelPart,
                                 const Element& rElementPrototype,
                                 const ConstitutiveLaw& rCLPrototype);
void AddDofs(ModelPart* pModelPart);
} // namespace test_model_part_factory_cpp
} // namespace

ModelPart& CreateStructuralMechanicsTestModelPart(Model* pModel,
                                                  const Element& rElementPrototype,
                                                  const ConstitutiveLaw& rCLPrototype,
                                                  std::function<void(ModelPart*)> CustomFunction)
{
    using namespace test_model_part_factory_cpp;
    const std::string name = "StructuralMechanicsTestModelPart";
    if (pModel->HasModelPart(name))
    {
        pModel->DeleteModelPart(name);
    }
    ModelPart& model_part = pModel->CreateModelPart(name);
    AddVariables(&model_part);
    CreateNodes(&model_part, rElementPrototype);
    model_part.GetProcessInfo()[DOMAIN_SIZE] =
        rElementPrototype.GetGeometry().WorkingSpaceDimension();
    CreatePropertiesAndElements(&model_part, rElementPrototype, rCLPrototype);
    model_part.SetBufferSize(2);
    AddDofs(&model_part);
    CustomFunction(&model_part);
    return model_part;
}

namespace
{ // cpp internals
namespace test_model_part_factory_cpp
{ // unity build unity guard
void AddVariables(ModelPart* pModelPart)
{
    pModelPart->AddNodalSolutionStepVariable(DISPLACEMENT);
    pModelPart->AddNodalSolutionStepVariable(REACTION);
    pModelPart->AddNodalSolutionStepVariable(VELOCITY);
    pModelPart->AddNodalSolutionStepVariable(ACCELERATION);
    pModelPart->AddNodalSolutionStepVariable(DENSITY);
    pModelPart->AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    pModelPart->AddNodalSolutionStepVariable(THICKNESS);
}

void CreateNodes(ModelPart* pModelPart, const Element& rElementPrototype)
{
    Matrix coords;
    rElementPrototype.GetGeometry().PointsLocalCoordinates(coords);
    if (coords.size2() == 2)
    {
        for (std::size_t i = 0; i < coords.size1(); ++i)
        {
            pModelPart->CreateNewNode(i + 1, coords(i, 0), coords(i, 1), 0.0);
        }
    }
    else
    {
        for (std::size_t i = 0; i < coords.size1(); ++i)
        {
            pModelPart->CreateNewNode(i + 1, coords(i, 0), coords(i, 1), coords(i, 2));
        }
    }
}

void CreatePropertiesAndElements(ModelPart* pModelPart,
                                 const Element& rElementPrototype,
                                 const ConstitutiveLaw& rCLPrototype)
{
    // Create properties
    const std::size_t properties_id = 1;
    auto& r_prop = *pModelPart->CreateNewProperties(properties_id);
    r_prop[CONSTITUTIVE_LAW] = rCLPrototype.Clone();
    r_prop[DENSITY] = 1000.0;
    r_prop[YOUNG_MODULUS] = 1400000.0;
    r_prop[POISSON_RATIO] = 0.2;
    r_prop[RAYLEIGH_ALPHA] = 0.0; // 0.02;
    r_prop[RAYLEIGH_BETA] = 0.0;  // 0.03;
    r_prop[VOLUME_ACCELERATION] = ZeroVector(3);
    r_prop[VOLUME_ACCELERATION](1) = 100.0;
    // Create elements
    PointerVector<Node<3>> nodes;
    for (std::size_t i = 0; i < rElementPrototype.GetGeometry().size(); ++i)
    {
        nodes.push_back(pModelPart->pGetNode(i+1));
    }
    auto p_new_elem =
        rElementPrototype.Create(1, nodes, pModelPart->pGetProperties(properties_id));
    pModelPart->AddElement(p_new_elem);
}

void AddDofs(ModelPart* pModelPart)
{
    for (auto& r_node : pModelPart->Nodes())
    {
        r_node.AddDof(DISPLACEMENT_X, REACTION_X);
        r_node.AddDof(DISPLACEMENT_Y, REACTION_Y);
        r_node.AddDof(DISPLACEMENT_Z, REACTION_Z);
    }
}
} // namespace test_model_part_factory_cpp
} // namespace
} // namespace Kratos
