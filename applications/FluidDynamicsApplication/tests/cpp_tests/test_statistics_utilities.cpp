//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//
//

// Project includes
#include "testing/testing.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"

#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/statistics_record.h"
#include "custom_utilities/statistics_data.h"
#include "custom_utilities/statistics_utilities.h"

namespace Kratos {
namespace Testing  {

namespace Internals {

void TestStatisticsUtilitiesInitializeModelPart(
    ModelPart& rModelPart, double DeltaTime, unsigned int BufferSize) {

    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(PRESSURE);
    rModelPart.SetBufferSize(BufferSize);
    Properties::Pointer p_properties = rModelPart.pGetProperties(0);

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
    std::vector<ModelPart::IndexType> element_nodes{1, 2, 3, 4};
    rModelPart.CreateNewElement("QSVMS3D4N", 1, element_nodes, p_properties);

    // Nodal data
    Element& r_element = *(rModelPart.ElementsBegin());
    Geometry< Node<3> >& r_geometry = r_element.GetGeometry();

    for (unsigned int i = 0; i < 4; i++) {
        Node<3>& r_node = r_geometry[i];
        r_node.FastGetSolutionStepValue(PRESSURE) = 10.0; // * r_node.Id();
        r_node.FastGetSolutionStepValue(VELOCITY_X) = 1.0;//r_node.Id() + 5.0;
        r_node.FastGetSolutionStepValue(VELOCITY_Y) = 2.0;//r_node.Id() + 5.0;
        r_node.FastGetSolutionStepValue(VELOCITY_Z) = 3.0;//r_node.Id() + 5.0;
    }

    // Element data
    r_element.SetValue(C_SMAGORINSKY,0.16);

    // ProcessInfo
    rModelPart.GetProcessInfo().SetValue(OSS_SWITCH,1);

    // Loop starts at 1 because you need one less clone than time steps
    for (unsigned int i = 1; i < BufferSize; i++) {
        rModelPart.CloneTimeStep(i * DeltaTime);

/*        for (unsigned int j = 0; j < 4; j++) {
            Node<3>& r_node = r_geometry[j];
            r_node.FastGetSolutionStepValue(PRESSURE) += float(i);
            r_node.FastGetSolutionStepValue(VELOCITY_Y) = r_node.Id() + i;
        }*/
    }
}

} // namespace internals

KRATOS_TEST_CASE_IN_SUITE(StatisticUtilitiesUsage, FluidDynamicsApplicationFastSuite) {
    ModelPart model_part("TestModelPart");
    Internals::TestStatisticsUtilitiesInitializeModelPart(model_part, 0.1 ,2);

    StatisticsRecord::Pointer p_turbulence_statistics = Kratos::make_shared<StatisticsRecord>();
    auto average_pressure_getter = Kratos::Internals::MakeSamplerAtLocalCoordinate::ValueGetter(PRESSURE);
    StatisticsSampler::Pointer average_pressure = Kratos::make_shared<ScalarAverageSampler>(average_pressure_getter,std::string("p"));
    p_turbulence_statistics->AddResult(average_pressure);
    auto average_velocity_getter = Kratos::Internals::MakeSamplerAtLocalCoordinate::ValueGetter(VELOCITY);
    std::vector<std::string> velocity_tags;
    velocity_tags.push_back(std::string("u"));
    velocity_tags.push_back(std::string("v"));
    velocity_tags.push_back(std::string("w"));
    StatisticsSampler::Pointer average_velocity = StatisticsSampler::Pointer(new VectorAverageSampler<array_1d<double,3>>(average_velocity_getter,3,velocity_tags));
    p_turbulence_statistics->AddResult(average_velocity);

    p_turbulence_statistics->InitializeStorage(model_part.Elements());
    model_part.GetProcessInfo().SetValue(STATISTICS_CONTAINER,p_turbulence_statistics);

    p_turbulence_statistics->SampleIntegrationPointResults(model_part);

    model_part.CloneTimeStep(0.2);
    p_turbulence_statistics->SampleIntegrationPointResults(model_part);

    std::vector<double> expected_output{10.,1.,2.,3.,10.,1.,2.,3.,10.,1.,2.,3.,10.,1.,2.,3.};
    std::vector<double> obtained_output = p_turbulence_statistics->OutputForTest(model_part.Elements());

    //std::cout << "Expected size " << expected_output.size() << " obtained size " << obtained_output.size() << std::endl;
    KRATOS_CHECK_EQUAL(expected_output.size(), obtained_output.size());

    for (unsigned int i = 0; i < expected_output.size(); i++)
    {
        //std::cout << "i: " << i << " expected " << expected_output[i] << " obtained " << obtained_output[i] << std::endl;
        KRATOS_CHECK_NEAR(expected_output[i],obtained_output[i], 1e-12);
    }
}

KRATOS_TEST_CASE_IN_SUITE(StatisticUtilitiesSecondThirdOrder, FluidDynamicsApplicationFastSuite) {
    ModelPart model_part("TestModelPart");
    Internals::TestStatisticsUtilitiesInitializeModelPart(model_part, 0.1 ,2);

    StatisticsRecord::Pointer p_turbulence_statistics = Kratos::make_shared<StatisticsRecord>();
    auto average_pressure_getter = Kratos::Internals::MakeSamplerAtLocalCoordinate::ValueGetter(PRESSURE);
    StatisticsSampler::Pointer average_pressure = Kratos::make_shared<ScalarAverageSampler>(average_pressure_getter,std::string("p"));
    p_turbulence_statistics->AddResult(average_pressure);
    auto average_velocity_getter = Kratos::Internals::MakeSamplerAtLocalCoordinate::ValueGetter(VELOCITY);
    std::vector<std::string> velocity_tags;
    velocity_tags.push_back(std::string("u"));
    velocity_tags.push_back(std::string("v"));
    velocity_tags.push_back(std::string("w"));
    StatisticsSampler::Pointer average_velocity = StatisticsSampler::Pointer(new VectorAverageSampler<array_1d<double,3>>(average_velocity_getter,3,velocity_tags));
    p_turbulence_statistics->AddResult(average_velocity);

    StatisticsSampler::Pointer pressure_correlation = Kratos::make_shared<VarianceSampler>(average_pressure, average_velocity);
    p_turbulence_statistics->AddHigherOrderStatistic(pressure_correlation);
    StatisticsSampler::Pointer reynolds_stresses = Kratos::make_shared<SymmetricVarianceSampler>(average_velocity);
    p_turbulence_statistics->AddHigherOrderStatistic(reynolds_stresses);

    StatisticsSampler::Pointer componentwise_correlation = Kratos::make_shared<ComponentwiseVarianceSampler>(average_pressure,0,average_velocity,0);
    p_turbulence_statistics->AddHigherOrderStatistic(componentwise_correlation);

    StatisticsSampler::Pointer third_order_correlation = Kratos::make_shared<ThirdOrderCorrelationSampler>(
        average_velocity,0,average_velocity,2,average_pressure,0,
        reynolds_stresses, reynolds_stresses->ComponentIndex(0,2), componentwise_correlation,0, pressure_correlation,2);
    p_turbulence_statistics->AddHigherOrderStatistic(third_order_correlation);

    p_turbulence_statistics->InitializeStorage(model_part.Elements());
    model_part.GetProcessInfo().SetValue(STATISTICS_CONTAINER,p_turbulence_statistics);

    for (auto it_node = model_part.NodesBegin(); it_node != model_part.NodesEnd(); ++it_node)
    {
        it_node->FastGetSolutionStepValue(VELOCITY_X) = 0.5;
        it_node->FastGetSolutionStepValue(VELOCITY_Z) = 4.0;
        it_node->FastGetSolutionStepValue(PRESSURE) = 15.0;
    }

    p_turbulence_statistics->SampleIntegrationPointResults(model_part);

    model_part.CloneTimeStep(0.2);
    for (auto it_node = model_part.NodesBegin(); it_node != model_part.NodesEnd(); ++it_node)
    {
        it_node->FastGetSolutionStepValue(VELOCITY_X) = 0.5;
        it_node->FastGetSolutionStepValue(VELOCITY_Z) = 4.0;
        it_node->FastGetSolutionStepValue(PRESSURE) = 10.0;
    }

    p_turbulence_statistics->SampleIntegrationPointResults(model_part);

    model_part.CloneTimeStep(0.2);
    for (auto it_node = model_part.NodesBegin(); it_node != model_part.NodesEnd(); ++it_node)
    {
        it_node->FastGetSolutionStepValue(VELOCITY_X) = 2.0;
        it_node->FastGetSolutionStepValue(VELOCITY_Z) = 1.0;
        it_node->FastGetSolutionStepValue(PRESSURE) = 5.0;
    }

    p_turbulence_statistics->SampleIntegrationPointResults(model_part);

    //p_turbulence_statistics->PrintToFile(model_part, "statistics");

    std::vector<double> expected_output{
    //  p   u  v  w   pu    pv pw   uu  uv  uw    vv vw  ww  pu(cw) uwp
        10.,1.,2.,3., -3.75,0.,7.5, 0.75,0.,-1.5, 0.,0., 3., -3.75, 2.5,
        10.,1.,2.,3., -3.75,0.,7.5, 0.75,0.,-1.5, 0.,0., 3., -3.75, 2.5,
        10.,1.,2.,3., -3.75,0.,7.5, 0.75,0.,-1.5, 0.,0., 3., -3.75, 2.5,
        10.,1.,2.,3., -3.75,0.,7.5, 0.75,0.,-1.5, 0.,0., 3., -3.75, 2.5};
    std::vector<double> obtained_output = p_turbulence_statistics->OutputForTest(model_part.Elements());

    //std::cout << "Expected size " << expected_output.size() << " obtained size " << obtained_output.size() << std::endl;
    KRATOS_CHECK_EQUAL(expected_output.size(), obtained_output.size());

    for (unsigned int i = 0; i < expected_output.size(); i++)
    {
        //std::cout << "i: " << i << " expected " << expected_output[i] << " obtained " << obtained_output[i] << std::endl;
        KRATOS_CHECK_NEAR(expected_output[i],obtained_output[i], 1e-12);
    }
}

}
}