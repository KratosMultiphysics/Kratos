//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Laura Moreno Martinez


// System includes
#include <vector>

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/kratos_components.h"
#include "includes/variables.h"
#include "mpm_application_variables.h"
#include "custom_utilities/mpm_nodal_cauchy_stress_utility.h"
#include "utilities/quadrature_points_utility.h"

namespace Kratos::Testing
{
namespace
{
    void AddMaterialPointElement(
        ModelPart& rMaterialPointModelPart,
        Element& rBackgroundElement,
        const IndexType ElementId,
        const array_1d<double, 3>& rMaterialPointCoordinates,
        const double MaterialPointMass,
        const Vector& rMaterialPointStress)
    {
        auto p_quadrature_point_geometry =
            CreateQuadraturePointsUtility<Node>::CreateFromCoordinates(
                rBackgroundElement.pGetGeometry(), rMaterialPointCoordinates, 1.0);

        const Element& r_mpm_element_prototype =
            KratosComponents<Element>::Get("MPMUpdatedLagrangian2D3N");

        Element::Pointer p_element = r_mpm_element_prototype.Create(
            ElementId,
            p_quadrature_point_geometry,
            rMaterialPointModelPart.pGetProperties(0));

        rMaterialPointModelPart.AddElement(p_element);

        const ProcessInfo& r_current_process_info = rMaterialPointModelPart.GetProcessInfo();
        p_element->SetValuesOnIntegrationPoints(
            MP_MASS, { MaterialPointMass }, r_current_process_info);
        p_element->SetValuesOnIntegrationPoints(
            MP_CAUCHY_STRESS_VECTOR, { rMaterialPointStress }, r_current_process_info);
    }

    void PrepareModelParts(
        ModelPart& rMaterialPointModelPart,
        ModelPart& rGridModelPart)
    {
        rMaterialPointModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
        rGridModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);

        rGridModelPart.AddNodalSolutionStepVariable(CAUCHY_STRESS_VECTOR);
        rGridModelPart.AddNodalSolutionStepVariable(NODAL_MASS);
        rGridModelPart.AddNodalSolutionStepVariable(NODAL_MOMENTUM);
        rGridModelPart.AddNodalSolutionStepVariable(NODAL_INERTIA);

        rGridModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
        rGridModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
        rGridModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);

        // to check that the utility does not crash when encountering nodes that are 
        // not connected to any element
        rGridModelPart.CreateNewNode(4, 2.0, 0.0, 0.0);

        for (auto& r_node : rGridModelPart.Nodes()) {
            Vector& r_nodal_stress =
                r_node.FastGetSolutionStepValue(CAUCHY_STRESS_VECTOR);
            r_nodal_stress = ScalarVector(6, 99.0);
        }

        rGridModelPart.CreateNewProperties(0);
        rMaterialPointModelPart.CreateNewProperties(0);

        auto p_background_element = rGridModelPart.CreateNewElement(
            "Element2D3N", 1, {1, 2, 3}, rGridModelPart.pGetProperties(0));

        Vector mp_cauchy_stress_1 = ZeroVector(3);
        mp_cauchy_stress_1[0] = 10.0;
        mp_cauchy_stress_1[1] = 20.0;
        mp_cauchy_stress_1[2] = 30.0;

        const array_1d<double, 3> mp_coordinates_1{0.2, 0.3, 0.0};
        AddMaterialPointElement( rMaterialPointModelPart, *p_background_element,1,
            mp_coordinates_1,2.0, mp_cauchy_stress_1);

        Vector mp_cauchy_stress_2 = ZeroVector(3);
        mp_cauchy_stress_2[0] = -2.0;
        mp_cauchy_stress_2[1] = 6.0;
        mp_cauchy_stress_2[2] = 14.0;

        const array_1d<double, 3> mp_coordinates_2{0.6, 0.1, 0.0};
        AddMaterialPointElement( rMaterialPointModelPart,*p_background_element,2, 
            mp_coordinates_2, 4.0, mp_cauchy_stress_2);

        for (auto& r_element : rMaterialPointModelPart.Elements()) {
            r_element.AddExplicitContribution(rMaterialPointModelPart.GetProcessInfo());
        }
    }

    void CheckNodalStress(
        const ModelPart& rModelPart,
        const IndexType NodeId,
        const std::vector<double>& rExpectedStress)
    {
        const auto& r_nodal_stress =
            rModelPart.GetNode(NodeId).FastGetSolutionStepValue(CAUCHY_STRESS_VECTOR);

        KRATOS_EXPECT_EQ(r_nodal_stress.size(), rExpectedStress.size());
        KRATOS_EXPECT_VECTOR_NEAR(r_nodal_stress, rExpectedStress, 1.0e-12);
    }
}

/**
* Checks that the nodal Cauchy stress utility resets the nodal values and
* computes the mass and shape-function weighted average of material point stress.
*/
KRATOS_TEST_CASE_IN_SUITE(MPMNodalCauchyStressUtilityWeightedAverage, KratosMPMFastSuite)
{
    Model current_model;
    ModelPart& r_material_point_model_part =
        current_model.CreateModelPart("MPM_Material");
    ModelPart& r_grid_model_part =
        current_model.CreateModelPart("Background_Grid");

    PrepareModelParts(
        r_material_point_model_part,
        r_grid_model_part);

    MPMNodalCauchyStressUtility::CalculateNodalCauchyStress(
        r_material_point_model_part,
        r_grid_model_part);

    CheckNodalStress( r_grid_model_part, 1, {3.4545454545454546, 12.363636363636363, 21.272727272727273});
    CheckNodalStress( r_grid_model_part, 2, {-0.2857142857142857, 8.0, 16.285714285714285});
    CheckNodalStress( r_grid_model_part, 3, {5.2, 14.4, 23.6});
    CheckNodalStress( r_grid_model_part, 4, {0.0, 0.0, 0.0});
}

} // namespace Kratos::Testing
