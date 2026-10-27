// License: BSD License (kratos/license.txt)
#include <array>
#include <vector>
#include "tests/cpp_tests/fluid_dynamics_fast_suite.h"
#include "containers/model.h"
#include "custom_constitutive/newtonian_2d_law.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos::Testing {
namespace {
Condition::Pointer GapCondition(ModelPart& part, const bool Connector, double Penalty)
{
    part.AddNodalSolutionStepVariable(VELOCITY);
    part.AddNodalSolutionStepVariable(ACCELERATION);
    part.AddNodalSolutionStepVariable(PRESSURE);
    part.GetProcessInfo()[BDF_COEFFICIENTS] = ZeroVector(3);
    auto material = part.CreateNewProperties(1);
    material->SetValue(DYNAMIC_VISCOSITY, 1.3);
    material->SetValue(CONSTITUTIVE_LAW, Kratos::make_shared<Newtonian2DLaw>());
    auto properties = part.CreateNewProperties(2);
    properties->SetValue(PENALTY_COEFFICIENT, Penalty);
    const double right = Connector?1.0:0.9;
    part.CreateNewNode(1, 0, 0, 0); part.CreateNewNode(2, right, 0, 0);
    part.CreateNewNode(3, right, 1, 0); part.CreateNewNode(4, 0, 1, 0);
    part.CreateNewNode(5, 1, 0, 0); part.CreateNewNode(6, 2, 0, 0);
    part.CreateNewNode(7, 2, 1, 0); part.CreateNewNode(8, 1, 1, 0);
    auto first = part.CreateNewElement("Element2D4N", 1, {1, 2, 3, 4}, material);
    auto second = part.CreateNewElement("Element2D4N", 2, Connector?
        std::vector<ModelPart::IndexType>{2, 6, 7, 3}:std::vector<ModelPart::IndexType>{5, 6, 7, 8}, material);
    for (auto& n:part.Nodes()) {
        n.AddDof(VELOCITY_X); n.AddDof(VELOCITY_Y); n.AddDof(PRESSURE);
        n.FastGetSolutionStepValue(VELOCITY)[0] = 0.1*n.Id();
        n.FastGetSolutionStepValue(VELOCITY)[1] = -0.03*n.Id();
        n.FastGetSolutionStepValue(PRESSURE) = 0.2*n.Id();
    }
    const std::vector<ModelPart::IndexType> ids = Connector?
        std::vector<ModelPart::IndexType>{1, 2, 3, 4, 6, 7}:std::vector<ModelPart::IndexType>{1, 2, 3, 4, 5, 6, 7, 8};
    auto condition = part.CreateNewCondition(Connector?"FluidPatchGapInterfaceCondition2D":"FluidPatchGapCouplingCondition2D8N", 1, ids, properties);
    GlobalPointersVector<Element> parents;
    parents.push_back(GlobalPointer<Element>(first.get())); parents.push_back(GlobalPointer<Element>(second.get()));
    condition->SetValue(NEIGHBOUR_ELEMENTS, parents);
    Matrix segment = ZeroMatrix(2, 3);
    segment(0, 0) = segment(1, 0) = 1; segment(1, 1) = 1;
    if (Connector) { segment(0, 1) = 1; segment(1, 0) = 0.8; segment(1, 1) = 1.3; }
    condition->SetValue(Connector?SURROGATE_BOUNDARY_FACE_COORDINATES:FEM_COUPLING_INTERFACE_COORDINATES, segment);
    condition->Initialize(part.GetProcessInfo());
    return condition;
}
}
KRATOS_TEST_CASE_IN_SUITE(FluidPatchGapCouplingLinearizationAndSharedNodes, FluidDynamicsApplicationFastSuite)
{
    for (const bool connector:{false, true}) {
        for (const double penalty:{0.0, 40.0}) {
            Model model; auto& part = model.CreateModelPart("Fluid");
            auto condition = GapCondition(part, connector, penalty);
            KRATOS_EXPECT_EQ(condition->Check(part.GetProcessInfo()), 0);
            Matrix lhs; Vector rhs;
            condition->CalculateLocalSystem(lhs, rhs, part.GetProcessInfo());
            KRATOS_EXPECT_EQ(lhs.size1(), connector?18:24);
            Matrix weighted = lhs;
            // A pressure-row sign change is needed to compare symmetric
            // Stokes blocks under the q div(u) bulk convention.
            if (penalty>0) {
                for (std::size_t i = 2; i<weighted.size1(); i+=3)
                    for (std::size_t j = 0; j<weighted.size2(); ++j) weighted(i, j)*=-1.0;
            }
            const Matrix error = penalty>0?Matrix(weighted-trans(weighted)):Matrix(weighted+trans(weighted));
            KRATOS_EXPECT_NEAR(norm_frobenius(error), 0.0, 1e-12);
            const double eps = 1e-6;
            for (std::size_t j = 0; j<rhs.size(); ++j) {
                auto& node = condition->GetGeometry()[j/3];
                double& value = j%3==2?node.FastGetSolutionStepValue(PRESSURE):node.FastGetSolutionStepValue(VELOCITY)[j%3];
                value+=eps;
                Matrix perturbed_lhs; Vector perturbed_rhs;
                condition->CalculateLocalSystem(perturbed_lhs, perturbed_rhs, part.GetProcessInfo());
                value-=eps;
                for (std::size_t i = 0; i<rhs.size(); ++i)
                    KRATOS_EXPECT_NEAR((perturbed_rhs[i]-rhs[i])/eps, -lhs(i, j), 1e-7);
            }
            if (!connector) {
                auto fitted = part.CreateNewCondition("FluidPatchCouplingCondition2D8N", 2,
                    {1, 2, 3, 4, 5, 6, 7, 8}, condition->pGetProperties());
                fitted->SetValue(NEIGHBOUR_ELEMENTS, condition->GetValue(NEIGHBOUR_ELEMENTS));
                fitted->SetValue(FEM_COUPLING_INTERFACE_COORDINATES, condition->GetValue(FEM_COUPLING_INTERFACE_COORDINATES));
                KRATOS_EXPECT_EXCEPTION_IS_THROWN(fitted->Check(part.GetProcessInfo()),"Coupling endpoints");
            }
        }
    }
}
} // namespace Kratos::Testing
