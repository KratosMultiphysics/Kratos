//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Messmer
//
//

// System includes
#include <numeric>

// External includes
#include <queso/kratos_mp_interface/model_part_adapter.hpp>
#include <queso/containers/element_container.hpp>
#include <queso/embedding/brep_operator.h>
#include <queso/utilities/mesh_utilities.h>
#include <queso/quadrature/trimmed_element.hpp>
#include <queso/quadrature/single_element.hpp>

// Project includes
#include "includes/kratos_filesystem.h"
#include "containers/model.h"
#include "integration/integration_point.h"
#include "testing/testing.h"
#include "input_output/stl_io.h"

namespace Kratos {

namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(QuESoTest, KratosExternalLibrariesFastSuite)
{
    // Some typedef's. Note that IntegrationPointType is of 'Kratos::' type
    typedef Kratos::IntegrationPoint<3> IntegrationPointType;
    typedef queso::BoundaryIntegrationPoint BoundaryIntegrationPointType;
    typedef queso::Element<IntegrationPointType, BoundaryIntegrationPointType> ElementType;

    Model model;
    auto& model_part = model.CreateModelPart("Test");

    // Read mode part from STL
    const std::filesystem::path file_name = std::filesystem::path(__FILE__).parent_path().string() + "/data/elephant.stl";

    Parameters settings(R"({
        "open_mode" : "read",
        "new_entity_type" : "element"
    })");

    StlIO stl_read(file_name, settings);
    stl_read.ReadModelPart(model_part);

    // Pass Kratos::ModelPart to QuESo's adapter.
    queso::ModelPartAdapter<ModelPart> triangle_mesh(model_part);

    // Set up QuESo parameters
    queso::Parameters parameters( {queso::Component("b_spline_mesh", true),
                                   queso::Component("lower_bound_xyz", queso::Vector3d{-0.38, -0.56, -0.32}),
                                   queso::Component("upper_bound_xyz", queso::Vector3d{0.38, 0.56, 0.32}),
                                   queso::Component("lower_bound_uvw", queso::Vector3d{-0.38, -0.56, -0.32}),
                                   queso::Component("upper_bound_uvw", queso::Vector3d{0.38, 0.56, 0.32}),
                                   queso::Component("number_of_elements", queso::Vector3i{10, 10, 10})} );


    // Compute volume for reference
    const double initial_volume = queso::MeshUtilities::VolumeOMP(triangle_mesh);
    KRATOS_EXPECT_NEAR(initial_volume,  0.0462012, 1e-4);

    // Set up BrepOperator
    queso::BRepOperator brep_operator(triangle_mesh);

    // This is a very robust and efficient element classifcation scheme, which works reliable for very complex STLs (including flawed models).
    // It groups all elements using a flood-fill algorithm.
    const auto p_el_status = brep_operator.pGetElementClassifications(parameters);
    queso::Mapper mapper(parameters);

    double volume_trimmed_elements = 0.0;
    double volume_full_elements = 0.0;
    double volume_integration_points = 0.0;

    for( IndexType index = 0; index < mapper.NumberOfElements(); ++index){
        // Get element status
        queso::IntersectionStatus status = (*p_el_status)[index];

        // Get bounding box of element
        // The type here is: std::pair<std::array<double,3>, std::array<double,3>>
        const auto bounding_box_xyz = mapper.GetBoundingBoxXYZFromIndex(index);
        const auto bounding_box_uvw = mapper.GetBoundingBoxUVWFromIndex(index);

        // This classifcation function relies on ray tracing. Less robust as the function above, but can be called locally on the bounding box.
        auto status_2 = brep_operator.GetIntersectionState(bounding_box_xyz.first, bounding_box_xyz.second);

        KRATOS_EXPECT_EQ(status, status_2);

        if( status == queso::IntersectionStatus::Inside || status == queso::IntersectionStatus::Trimmed ) {
            // Construct element and check status:
            ElementType new_element(index+1, bounding_box_xyz, bounding_box_uvw);
            bool valid_element = false;

            // Distinguish between trimmed and non-trimmed elements.
            if( status == queso::IntersectionStatus::Trimmed) {
                new_element.SetIsTrimmed(true);
                const double min_vol_element_ratio = 0.0; // This is usually ~0.01-0.001
                const std::size_t num_boundary_triangles = 50;
                auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(bounding_box_xyz.first, bounding_box_xyz.second,
                                                                        min_vol_element_ratio, num_boundary_triangles);
                if( p_trimmed_domain ){ // Don't forget this check.
                    new_element.pSetTrimmedDomain(p_trimmed_domain);

                    const auto& r_intersection_mesh = new_element.pGetTrimmedDomain()->GetTriangleMesh();
                    volume_trimmed_elements += queso::MeshUtilities::Volume(r_intersection_mesh);

                    valid_element = true;
                }

                // If valid solve moment fitting equation
                if( valid_element ){
                    const queso::Vector3i polynomial_order{2, 2, 2};
                    const double moment_fitting_residual = 1e-10;
                    queso::QuadratureTrimmedElement<ElementType>::AssembleIPs(new_element, polynomial_order, moment_fitting_residual, 0);
                    if( new_element.GetIntegrationPoints().size() == 0 ){
                        valid_element = false;
                    }
                }
            }
            else if( status == queso::IntersectionStatus::Inside){
                // Get standard Gauss legendre points
                const queso::Vector3i polynomial_order{2, 2, 2};
                const queso::IntegrationMethodType integration_method = queso::IntegrationMethod::Gauss;
                queso::QuadratureSingleElement<ElementType>::AssembleIPs(new_element, polynomial_order, integration_method);
                valid_element = true;

                volume_full_elements += (bounding_box_xyz.second[0] - bounding_box_xyz.first[0])*
                                        (bounding_box_xyz.second[1] - bounding_box_xyz.first[1])*
                                        (bounding_box_xyz.second[2] - bounding_box_xyz.first[2]);
            }

            if( valid_element ){
                const std::vector<IntegrationPointType>& r_points = new_element.GetIntegrationPoints();
                volume_integration_points += std::accumulate(r_points.begin(), r_points.end(), 0.0, [](double sum, const IntegrationPointType& r_point) {return sum + r_point.Weight(); });
            }
        }
    }

    KRATOS_EXPECT_NEAR(initial_volume, (volume_full_elements+volume_trimmed_elements), 1e-8);
    KRATOS_EXPECT_NEAR(initial_volume, volume_integration_points, 1e-8);
}

} // namespace Testing.
} // namespace Kratos.