//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Raul Bravo
//

#ifdef KRATOS_USE_FUTURE

// System includes
#include <vector>

/* External includes */

/* Project includes */
#include "containers/model.h"
#include "geometries/line_2d_2.h"
#include "includes/kratos_parameters.h"
#include "testing/testing.h"
#include "future/containers/define_linear_algebra_serial.h"
#include "future/containers/linear_system_tags.h"
#include "future/solving_strategies/schemes/static_scheme.h"
#include "future/solving_strategies/strategies/implicit_strategy_data.h"

/* Application includes */
#include "future/hrom_scheme.h"

namespace Kratos::Testing {
namespace FutureHRomSchemeTestingInternal {

using SchemeType = Future::StaticScheme<Future::SerialLinearAlgebraTraits>;
using HRomSchemeType = Future::HRomScheme<SchemeType>;
using StrategyDataType = Future::ImplicitStrategyData<Future::SerialLinearAlgebraTraits>;

/// Two node element with a unit Laplacian and a uniform source
class DummyLaplacianElement final : public Element
{
public:
    DummyLaplacianElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
    }

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override
    {
        rResult.clear();
        for (const auto& r_node : GetGeometry()) {
            rResult.push_back(r_node.GetDof(TEMPERATURE).EquationId());
        }
    }

    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const override
    {
        rElementalDofList.clear();
        for (const auto& r_node : GetGeometry()) {
            rElementalDofList.push_back(r_node.pGetDof(TEMPERATURE));
        }
    }

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override
    {
        rLeftHandSideMatrix.resize(2, 2, false);
        rRightHandSideVector.resize(2, false);
        rLeftHandSideMatrix(0,0) = 1.0; rLeftHandSideMatrix(0,1) = -1.0;
        rLeftHandSideMatrix(1,0) = -1.0; rLeftHandSideMatrix(1,1) = 1.0;
        rRightHandSideVector[0] = 0.5; rRightHandSideVector[1] = 0.5;
    }

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override
    {
        VectorType rhs;
        CalculateLocalSystem(rLeftHandSideMatrix, rhs, rCurrentProcessInfo);
    }

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override
    {
        MatrixType lhs;
        CalculateLocalSystem(lhs, rRightHandSideVector, rCurrentProcessInfo);
    }
};

/// Three nodes and two elements (1-2 and 2-3)
ModelPart& FillModel(Model& rModel)
{
    auto& r_model_part = rModel.CreateModelPart("main");
    r_model_part.CreateNewProperties(0);
    r_model_part.SetBufferSize(1);
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    for (std::size_t i = 1; i <= 3; ++i) {
        r_model_part.CreateNewNode(i, static_cast<double>(i - 1), 0.0, 0.0)->AddDof(TEMPERATURE);
    }
    auto p_properties = r_model_part.pGetProperties(0);
    for (std::size_t i = 1; i <= 2; ++i) {
        auto p_geometry = Kratos::make_shared<Line2D2<Node>>(r_model_part.pGetNode(i), r_model_part.pGetNode(i+1));
        r_model_part.AddElement(Kratos::make_intrusive<DummyLaplacianElement>(i, p_geometry, p_properties));
    }
    return r_model_part;
}

/// Builds the system with the given scheme and returns the LHS values and the RHS
template<class TSchemeType>
std::pair<std::vector<double>, std::vector<double>> BuildSystem(TSchemeType& rScheme, StrategyDataType& rStrategyData)
{
    auto p_linear_system = rStrategyData.pGetLinearSystem();
    auto& r_lhs = *(p_linear_system->pGetMatrix(Future::LinearSystemTags::SparseMatrixTag::LHS));
    auto& r_rhs = *(p_linear_system->pGetVector(Future::LinearSystemTags::DenseVectorTag::RHS));
    r_lhs.SetValue(0.0);
    r_rhs.SetValue(0.0);
    rScheme.Build(r_lhs, r_rhs);
    std::vector<double> lhs(r_lhs.value_data().begin(), r_lhs.value_data().end());
    std::vector<double> rhs(r_rhs.size());
    for (std::size_t i = 0; i < r_rhs.size(); ++i) {
        rhs[i] = r_rhs[i];
    }
    return {lhs, rhs};
}

Parameters SchemeSettings()
{
    return Parameters(R"({"build_settings" : {"name" : "block_builder"}})");
}

}

KRATOS_TEST_CASE_IN_SUITE(FutureHRomScheme, RomApplicationFastSuite)
{
    using namespace FutureHRomSchemeTestingInternal;

    // Reference: the full order system, and the one of the second element only (the first one deactivated)
    Model model{};
    ModelPart& r_model_part = FillModel(model);
    SchemeType scheme(r_model_part, SchemeSettings());
    StrategyDataType strategy_data;
    scheme.Initialize(strategy_data);
    const auto full_system = BuildSystem(scheme, strategy_data);
    r_model_part.GetElement(1).Set(ACTIVE, false);
    auto second_element_system = BuildSystem(scheme, strategy_data);
    r_model_part.GetElement(1).Set(ACTIVE, true);

    // HROM with all the elements and unit weights: the full order system
    HRomSchemeType hrom_scheme(r_model_part, SchemeSettings());
    StrategyDataType hrom_strategy_data;
    hrom_scheme.Initialize(hrom_strategy_data);
    hrom_scheme.SetHRomEntities({1, 2}, {1.0, 1.0}, {}, {});
    auto hrom_system = BuildSystem(hrom_scheme, hrom_strategy_data);
    KRATOS_EXPECT_VECTOR_NEAR(hrom_system.first, full_system.first, 1.0e-12);
    KRATOS_EXPECT_VECTOR_NEAR(hrom_system.second, full_system.second, 1.0e-12);

    // HROM with the second element only, weight 2 (the first one has a null weight, so it is not assembled):
    // twice the system of the second element. Note that the row of the node 1 stays empty (it is not modified)
    hrom_scheme.SetHRomEntities({1, 2}, {0.0, 2.0}, {}, {});
    KRATOS_EXPECT_EQ(hrom_scheme.NumberOfHRomElements(), 1);
    hrom_system = BuildSystem(hrom_scheme, hrom_strategy_data);
    for (auto& r_value : second_element_system.first) {
        r_value *= 2.0;
    }
    for (auto& r_value : second_element_system.second) {
        r_value *= 2.0;
    }
    KRATOS_EXPECT_VECTOR_NEAR(hrom_system.first, second_element_system.first, 1.0e-12);
    KRATOS_EXPECT_VECTOR_NEAR(hrom_system.second, second_element_system.second, 1.0e-12);
    const std::size_t node_1_eq_id = r_model_part.GetNode(1).GetDof(TEMPERATURE).EquationId();
    KRATOS_EXPECT_DOUBLE_EQ(hrom_system.second[node_1_eq_id], 0.0);
}

} // namespace Kratos::Testing

#endif // KRATOS_USE_FUTURE
