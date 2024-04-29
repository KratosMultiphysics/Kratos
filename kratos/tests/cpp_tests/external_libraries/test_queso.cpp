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
#include <queso/kratos_mp_interface/triangle_mesh_kratos.hpp>
#include <queso/containers/triangle_mesh.hpp>
#include <queso/containers/element_container.hpp>
#include <queso/embedding/brep_operator.h>
#include <queso/utilities/mesh_utilities.h>
#include <queso/quadrature/trimmed_element.h>
#include <queso/quadrature/single_element.h>

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "input_output/stl_io.h"

namespace Kratos {

namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(QuESoTest, KratosExternalLibrariesFastSuite)
{
    Model model;
    auto& model_part = model.CreateModelPart("Test");

    // read stl
    std::filesystem::path filename = "/home/manuel/QuESo/examples/stanford_bunny/data/test.stl";

    Parameters settings(R"({
        "open_mode" : "read",
        "new_entity_type" : "element"
    })");

    StlIO stl_read(filename, settings);
    stl_read.ReadModelPart(model_part);

    queso::Parameters parameters( {queso::Component("b_spline_mesh", true),
                                   queso::Component("lower_bound_xyz", queso::Vector3d(-24.0, -43.0, 5.0)),
                                   queso::Component("upper_bound_xyz", queso::Vector3d(85, 46.0, 115)),
                                   queso::Component("lower_bound_uvw", queso::Vector3d(-24.0, -43.0, 5.0)),
                                   queso::Component("upper_bound_uvw", queso::Vector3d(85, 46.0, 115)),
                                   queso::Component("number_of_elements", queso::Vector3i(50, 50, 50))} );

    queso::TriangleMeshKratos<ModelPart> triangle_mesh_kratos(model_part);

    const double initial_volume = queso::MeshUtilities::VolumeOMP(triangle_mesh_kratos);
    KRATOS_EXPECT_NEAR(initial_volume,  279628.299, 0.01);

    queso::BRepOperator brep_operator(triangle_mesh_kratos);
    const auto p_el_status = brep_operator.pGetElementClassifications(parameters);
    queso::Mapper mapper(parameters);

    const double min_vol_element_ratio = 0.0;
    const int num_boundary_triangles = 50;
    const queso::Vector3i polynomial_order(2, 2, 2);
    const double moment_fitting_residual = 1e-10;
    const int echo_level = 0;
    const queso::IntegrationMethodType integration_method = queso::IntegrationMethod::Gauss;
    double volume_points = 0.0;
    queso::ElementContainer element_container(parameters);

    element_container.reserve(mapper.NumberOfElements());

    for( IndexType index = 0; index < mapper.NumberOfElements(); ++index){
        queso::IntersectionStatus status = (*p_el_status)[index];

        if( status == queso::IntersectionStatus::Inside || status == queso::IntersectionStatus::Trimmed ) {
            // Get bounding box of element
            const auto bounding_box_xyz = mapper.GetBoundingBoxXYZFromIndex(index);
            const auto bounding_box_uvw = mapper.GetBoundingBoxUVWFromIndex(index);

            // Mapper can also give you 3d indices
            auto indices = mapper.GetMatrixIndicesFromVectorIndex(index);
            auto index_2 = mapper.GetVectorIndexFromMatrixIndices(indices[0], indices[1], indices[2]);
            // index_2 is now same as index

            // Construct element and check status:
            queso::Unique<queso::Element> p_new_element = queso::MakeUnique<queso::Element>(index+1, bounding_box_xyz, bounding_box_uvw);
            bool valid_element = false;

            // Distinguish between trimmed and non-trimmed elements.
            if( status == queso::IntersectionStatus::Trimmed) {
                p_new_element->SetIsTrimmed(true);
                auto p_trimmed_domain = brep_operator.pGetTrimmedDomain(bounding_box_xyz.first, bounding_box_xyz.second,
                                                                        min_vol_element_ratio, num_boundary_triangles);
                if( p_trimmed_domain ){
                    p_new_element->pSetTrimmedDomain(p_trimmed_domain);
                    valid_element = true;
                }

                // If valid solve moment fitting equation
                if( valid_element ){
                    queso::QuadratureTrimmedElement::AssembleIPs(*p_new_element, polynomial_order, moment_fitting_residual, echo_level);

                    if( p_new_element->GetIntegrationPoints().size() == 0 ){
                        valid_element = false;
                    }
                }
            }
            else if( status == queso::IntersectionStatus::Inside){
                // Get standard gauss legendre points
                queso::QuadratureSingleElement::AssembleIPs(*p_new_element, polynomial_order, integration_method);
                valid_element = true;
            }

            if( valid_element ){
                const std::vector<queso::IntegrationPoint>& r_points = p_new_element->GetIntegrationPoints();
                volume_points += std::accumulate(r_points.begin(), r_points.end(), 0.0, [](double sum, const queso::IntegrationPoint& r_point) {return sum + r_point.GetWeight(); });
                element_container.AddElement(p_new_element);
            }
        }
    }

    // Actually you can move trough the active elements like this
    std::size_t start_id = 5; // Actually start id should be an active element. Not sure if 5 is active.
    std::size_t next_id;
    bool is_found = false;
    bool is_local_end = false;
    auto p_element = element_container.pGetNextElementInX(start_id, next_id, is_found, is_local_end);
    start_id = next_id;
    auto p_element_2 = element_container.pGetNextElementInX(start_id, next_id, is_found, is_local_end);
    start_id = next_id;
    auto p_element_3 = element_container.pGetPreviousElementInX(start_id, next_id, is_found, is_local_end);
    // p_element_3 should now be the same as p_element_2
    // if local_end = true, then you reached the end of the active elements in this particular row.

    KRATOS_EXPECT_NEAR(initial_volume, volume_points, 1e-5);
}

} // namespace Testing.
} // namespace Kratos.