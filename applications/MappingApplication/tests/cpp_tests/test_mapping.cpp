//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

// Project includes
#include "testing/testing.h"

namespace Kratos {
namespace Testing {


KRATOS_TEST_CASE_IN_SUITE(InterfaceNodeIds, MappingApplicationFastSuite) {

    // TestNodalScalarData nodal_scalar_data;
    // TestNodalVectorData nodal_vector_data;
    // TestElementData element_data;
    // TestPropertiesData properties_data;
    // TestProcessInfoData process_info_data;

    // ModelPart full_model_part("Test Full");

    // constexpr double DeltaTime = 0.1;
    // FluidElementDataTestCompleteModelPart(full_model_part,DeltaTime,2);
    // Element& r_element = *(full_model_part.ElementsBegin());
    // ProcessInfo& r_process_info = full_model_part.GetProcessInfo();

    // nodal_scalar_data.Initialize(r_element,r_process_info);
    // KRATOS_CHECK_EQUAL(nodal_scalar_data.Pressure[0], 11.0);
    // KRATOS_CHECK_EQUAL(nodal_scalar_data.Pressure[1], 21.0);
    // KRATOS_CHECK_EQUAL(nodal_scalar_data.Pressure[2], 31.0);
    // KRATOS_CHECK_EQUAL(nodal_scalar_data.Pressure_OldStep1[0], 10.0);
    // KRATOS_CHECK_EQUAL(nodal_scalar_data.Pressure_OldStep1[1], 20.0);
    // KRATOS_CHECK_EQUAL(nodal_scalar_data.Pressure_OldStep1[2], 30.0);

    // nodal_vector_data.Initialize(r_element,r_process_info);
    // KRATOS_CHECK_EQUAL(nodal_vector_data.Velocity(0,1), 2.0);
    // KRATOS_CHECK_EQUAL(nodal_vector_data.Velocity(1,1), 3.0);
    // KRATOS_CHECK_EQUAL(nodal_vector_data.Velocity(2,1), 4.0);
    // KRATOS_CHECK_EQUAL(nodal_vector_data.Velocity_OldStep1(0,0), 6.0);
    // KRATOS_CHECK_EQUAL(nodal_vector_data.Velocity_OldStep1(1,0), 7.0);
    // KRATOS_CHECK_EQUAL(nodal_vector_data.Velocity_OldStep1(2,0), 8.0);

    // element_data.Initialize(r_element,r_process_info);
    // KRATOS_CHECK_EQUAL(element_data.CSmagorinsky, r_element.GetValue(C_SMAGORINSKY));

    // properties_data.Initialize(r_element,r_process_info);
    // KRATOS_CHECK_EQUAL(properties_data.KinematicViscosity, r_element.GetProperties().GetValue(KINEMATIC_VISCOSITY));

    // process_info_data.Initialize(r_element,r_process_info);
    // KRATOS_CHECK_EQUAL(process_info_data.UseOSS, r_process_info.GetValue(OSS_SWITCH));
    // KRATOS_CHECK_EQUAL(process_info_data.DeltaTime, r_process_info.GetValue(DELTA_TIME));
}

}  // namespace Testing
}  // namespace Kratos