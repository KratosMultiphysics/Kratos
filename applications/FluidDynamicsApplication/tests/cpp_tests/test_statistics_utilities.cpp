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

class MakeSamplerAtLocalCoordinate {
public:
    static std::function<double(const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives)> ValueGetter(Variable<double>& rVariable) {
        return [&rVariable](const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives) -> double {
            KRATOS_DEBUG_ERROR_IF(rGeometry.size() != rShapeFunctions.size()) << "Number of nodes in provided geometry does not match number of shape functions" << std::endl;
            double value = 0.0;
            for (unsigned int i =  0; i < rGeometry.size(); i++) {
                value += rGeometry[i].FastGetSolutionStepValue(rVariable) * rShapeFunctions[i];
            }
            return value;
        };
    }


    static std::function< array_1d<double,3>(const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives) > ValueGetter(Variable<array_1d<double,3>>& rVariable) {
        return [&rVariable](const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives) -> array_1d<double,3> {
            KRATOS_DEBUG_ERROR_IF(rGeometry.size() != rShapeFunctions.size()) << "Number of nodes in provided geometry does not match number of shape functions" << std::endl;
            array_1d<double,3> value = ZeroVector(3);
            for (unsigned int i =  0; i < rGeometry.size(); i++) {
                value += rGeometry[i].FastGetSolutionStepValue(rVariable) * rShapeFunctions[i];
            }
            return value;
        };
    }

    static std::function< Matrix(const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives) > ValueGetter(Variable<Matrix>& rVariable) {
        return [&rVariable](const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives) -> Matrix {
            KRATOS_DEBUG_ERROR_IF(rGeometry.size() != rShapeFunctions.size()) << "Number of nodes in provided geometry does not match number of shape functions" << std::endl;
            Matrix value = ZeroMatrix(3,3);
            for (unsigned int i =  0; i < rGeometry.size(); i++) {
                value += rGeometry[i].FastGetSolutionStepValue(rVariable) * rShapeFunctions[i];
            }
            return value;
        };
    }


    static std::function< array_1d<double,3>(const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives) > GradientGetter(Variable<double>& rVariable) {
        return [&rVariable](const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives) -> array_1d<double,3> {
            KRATOS_DEBUG_ERROR_IF(rGeometry.size() != rShapeDerivatives.size1()) << "Number of nodes in provided geometry does not match number of shape functions" << std::endl;
            array_1d<double,3> gradient = ZeroVector(3);
            for (unsigned int n =  0; n < rGeometry.size(); n++) {
                const auto& value = rGeometry[n].FastGetSolutionStepValue(rVariable);
                for (unsigned int i = 0; i < rShapeDerivatives.size2(); i++)
                    gradient[i] += value * rShapeDerivatives(n,i);
            }
            return gradient;
        };
    }


    static std::function< Matrix(const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives) > GradientGetter(const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, Variable<array_1d<double,3>>& rVariable) {
        return [&rVariable](const Geometry< Node<3> >& rGeometry, const Vector& rShapeFunctions, const Matrix& rShapeDerivatives) -> Matrix {
            KRATOS_DEBUG_ERROR_IF(rGeometry.size() != rShapeDerivatives.size1()) << "Number of nodes in provided geometry does not match number of shape functions" << std::endl;
            Matrix gradient = ZeroMatrix(3,3);
            for (unsigned int n =  0; n < rGeometry.size(); n++) {
                const auto& value = rGeometry[n].FastGetSolutionStepValue(rVariable);
                for (unsigned int i = 0; i < rShapeDerivatives.size2(); i++)
                {
                    for (unsigned int j = 0; j < rShapeDerivatives.size2(); j++)
                        gradient(i,j) += value[i] * rShapeDerivatives(n,j); //dui/dxj
                }
            }
            return gradient;
        };
    }
};

} // namespace internals

KRATOS_TEST_CASE_IN_SUITE(StatisticUtilitiesUsage, FluidDynamicsApplicationFastSuite) {
    ModelPart model_part("TestModelPart");
    Internals::TestStatisticsUtilitiesInitializeModelPart(model_part, 0.1 ,2);

    StatisticsRecord::Pointer p_turbulence_statistics = Kratos::make_shared<StatisticsRecord>();
    auto average_pressure_getter = Internals::MakeSamplerAtLocalCoordinate::ValueGetter(PRESSURE);
    StatisticsSampler::Pointer average_pressure = Kratos::make_shared<ScalarAverageSampler>(average_pressure_getter);
    p_turbulence_statistics->AddResult(average_pressure);
    auto average_velocity_getter = Internals::MakeSamplerAtLocalCoordinate::ValueGetter(VELOCITY);
    auto average_velocity_sampler = VectorAverageSampler<array_1d<double,3>>(average_velocity_getter,3);
    StatisticsSampler::Pointer average_velocity = StatisticsSampler::Pointer(new VectorAverageSampler<array_1d<double,3>>(average_velocity_getter,3));
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


KRATOS_TEST_CASE_IN_SUITE(StatisticUtilitiesVariance, FluidDynamicsApplicationFastSuite) {
    ModelPart model_part("TestModelPart");
    Internals::TestStatisticsUtilitiesInitializeModelPart(model_part, 0.1 ,2);

    StatisticsRecord::Pointer p_turbulence_statistics = Kratos::make_shared<StatisticsRecord>();
    auto average_pressure_getter = Internals::MakeSamplerAtLocalCoordinate::ValueGetter(PRESSURE);
    StatisticsSampler::Pointer average_pressure = Kratos::make_shared<ScalarAverageSampler>(average_pressure_getter);
    p_turbulence_statistics->AddResult(average_pressure);
    auto average_velocity_getter = Internals::MakeSamplerAtLocalCoordinate::ValueGetter(VELOCITY);
    auto average_velocity_sampler = VectorAverageSampler<array_1d<double,3>>(average_velocity_getter,3);
    StatisticsSampler::Pointer average_velocity = StatisticsSampler::Pointer(new VectorAverageSampler<array_1d<double,3>>(average_velocity_getter,3));
    p_turbulence_statistics->AddResult(average_velocity);

    StatisticsSampler::Pointer pressure_correlation = Kratos::make_shared<VarianceSampler>(average_pressure, average_velocity);
    p_turbulence_statistics->AddHigherOrderStatistic(pressure_correlation);
    StatisticsSampler::Pointer reynolds_stresses = Kratos::make_shared<VarianceSampler>(average_velocity, average_velocity);
    p_turbulence_statistics->AddHigherOrderStatistic(reynolds_stresses);

    p_turbulence_statistics->InitializeStorage(model_part.Elements());
    model_part.GetProcessInfo().SetValue(STATISTICS_CONTAINER,p_turbulence_statistics);

    for (auto it_node = model_part.NodesBegin(); it_node != model_part.NodesEnd(); ++it_node)
    {
        it_node->FastGetSolutionStepValue(VELOCITY_X) -= 0.5;
        it_node->FastGetSolutionStepValue(VELOCITY_Z) += 1.0;
        it_node->FastGetSolutionStepValue(PRESSURE) += 5.0;
    }

    p_turbulence_statistics->SampleIntegrationPointResults(model_part);

    model_part.CloneTimeStep(0.2);
    for (auto it_node = model_part.NodesBegin(); it_node != model_part.NodesEnd(); ++it_node)
    {
        it_node->FastGetSolutionStepValue(VELOCITY_X) += 1.0;
        it_node->FastGetSolutionStepValue(VELOCITY_Z) -= 2.0;
        it_node->FastGetSolutionStepValue(PRESSURE) -= 10.0;
    }

    p_turbulence_statistics->SampleIntegrationPointResults(model_part);

    p_turbulence_statistics->PrintToFile(model_part);

    std::vector<double> expected_output{
    //  p   u  v  w   pu  pv pw   uu  uv uw   vu vv vw  wu  wv ww
        10.,1.,2.,3., -5.,0.,10., 0.5,0.,-1., 0.,0.,0., -1.,0.,2.,
        10.,1.,2.,3., -5.,0.,10., 0.5,0.,-1., 0.,0.,0., -1.,0.,2.,
        10.,1.,2.,3., -5.,0.,10., 0.5,0.,-1., 0.,0.,0., -1.,0.,2.,
        10.,1.,2.,3., -5.,0.,10., 0.5,0.,-1., 0.,0.,0., -1.,0.,2.};
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