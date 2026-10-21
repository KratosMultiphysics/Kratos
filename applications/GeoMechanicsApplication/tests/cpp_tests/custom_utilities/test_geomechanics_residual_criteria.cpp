// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#include "containers/model.h"
#include "custom_utilities/geomechanics_residual_criteria.hpp"
#include "geo_mechanics_application_variables.h"
#include "spaces/default_spaces.h"
#include "test_setup_utilities/model_setup_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <vector>

namespace Kratos::Testing
{

using SparseSpaceType = TDefaultSparseSpace<double>;
using LocalSpaceType  = TDefaultDenseSpace<double>;
using GeoMechanicsResidualCriteriaType = GeoMechanicsResidualCriteria<SparseSpaceType, LocalSpaceType>;

PointerVector<Node> CreateNodesWithWaterPressureDofs(ModelPart& rModelPart, const NodeDefinitionVector& rNodeDefinitions)
{
    auto nodes = ModelSetupUtilities::CreateNodes(rModelPart, rNodeDefinitions);
    for (auto& r_node : nodes) {
        r_node.AddDof(WATER_PRESSURE);
    }
    return nodes;
}

class TestableGeoMechanicsResidualCriteria : public GeoMechanicsResidualCriteriaType
{
public:
    using GeoMechanicsResidualCriteriaType::ComputeActiveDofs;

    const std::vector<int>& GetActiveDofs() const { return this->mActiveDofs; }
};

KRATOS_TEST_CASE_IN_SUITE(GeoMechanicsResidualCriteria_MasksActiveConstraintDofs, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_model_part.GetNodalSolutionStepVariablesList().AddDof(&WATER_PRESSURE);

    CreateNodesWithWaterPressureDofs(
        r_model_part,
        {{1, {0.0, 0.0, 0.0}}, {2, {1.0, 0.0, 0.0}}, {3, {2.0, 0.0, 0.0}}, {4, {3.0, 0.0, 0.0}}});

    auto& r_master_node         = r_model_part.GetNode(1);
    auto& r_slave_node          = r_model_part.GetNode(2);
    auto& r_fixed_node          = r_model_part.GetNode(3);
    auto& r_inactive_slave_node = r_model_part.GetNode(4);

    r_fixed_node.Fix(WATER_PRESSURE);

    r_master_node.pGetDof(WATER_PRESSURE)->SetEquationId(0);
    r_slave_node.pGetDof(WATER_PRESSURE)->SetEquationId(1);
    r_fixed_node.pGetDof(WATER_PRESSURE)->SetEquationId(2);
    r_inactive_slave_node.pGetDof(WATER_PRESSURE)->SetEquationId(3);

    r_model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, r_master_node,
                                                WATER_PRESSURE, r_slave_node, WATER_PRESSURE, 1.0, 0.0);
    auto p_inactive_constraint = r_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 2, r_fixed_node, WATER_PRESSURE, r_inactive_slave_node,
        WATER_PRESSURE, 1.0, 0.0);
    p_inactive_constraint->Set(ACTIVE, false);

    ModelPart::DofsArrayType dof_set;
    dof_set.push_back(r_master_node.pGetDof(WATER_PRESSURE));
    dof_set.push_back(r_slave_node.pGetDof(WATER_PRESSURE));
    dof_set.push_back(r_fixed_node.pGetDof(WATER_PRESSURE));
    dof_set.push_back(r_inactive_slave_node.pGetDof(WATER_PRESSURE));

    TestableGeoMechanicsResidualCriteria criterion;
    criterion.ComputeActiveDofs(r_model_part, dof_set);

    KRATOS_EXPECT_EQ(criterion.GetActiveDofs(), (std::vector<int>{0, 0, 0, 1}));
}

KRATOS_TEST_CASE_IN_SUITE(GeoMechanicsResidualCriteria_ConvergesForRawTiedDofResidual, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_model_part.GetNodalSolutionStepVariablesList().AddDof(&WATER_PRESSURE);

    CreateNodesWithWaterPressureDofs(r_model_part, {{1, {0.0, 0.0, 0.0}}, {2, {1.0, 0.0, 0.0}}});
    auto& r_master_node = r_model_part.GetNode(1);
    auto& r_slave_node  = r_model_part.GetNode(2);
    r_master_node.pGetDof(WATER_PRESSURE)->SetEquationId(0);
    r_slave_node.pGetDof(WATER_PRESSURE)->SetEquationId(1);

    r_model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, r_master_node,
                                                WATER_PRESSURE, r_slave_node, WATER_PRESSURE, 1.0, 0.0);

    ModelPart::DofsArrayType dof_set;
    dof_set.push_back(r_master_node.pGetDof(WATER_PRESSURE));
    dof_set.push_back(r_slave_node.pGetDof(WATER_PRESSURE));

    GeoMechanicsResidualCriteriaType criterion;
    CompressedMatrix                 A;
    Vector                           dx(2);
    Vector                           rb(2);
    rb[0] = 1.0;
    rb[1] = -1.0;

    criterion.InitializeSolutionStep(r_model_part, dof_set, A, dx, rb);
    KRATOS_EXPECT_TRUE(criterion.PostCriteria(r_model_part, dof_set, A, dx, rb));
}

} // namespace Kratos::Testing
