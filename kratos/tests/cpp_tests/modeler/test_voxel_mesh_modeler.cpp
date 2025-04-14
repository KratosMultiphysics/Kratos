//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "includes/kratos_application.h"
#include "includes/kernel.h"

#include "modeler/voxel_mesh_generator_modeler.h"

#include "modeler/surrogate_boundary_modeler.h"

namespace Kratos::Testing {

namespace {
void WriteCubeSkinMeshMdpaFileForVoxelModelerTest()
{
    Kratos::shared_ptr<std::iostream> p_input(new std::stringstream(
        R"input(
	Begin ModelPartData
End ModelPartData
Begin Properties 0
End Properties
Begin Properties 1
End Properties
Begin Nodes
		1     0.03     0.03    0.09
		2     0.03    -0.03    0.09
		3     0.03    -0.03    0.03
		4     0.03     0.03    0.03
		5    -0.03     0.03    0.09
		6    -0.03     0.03    0.03
		7    -0.03    -0.03    0.09
		8    -0.03    -0.03    0.03
End Nodes
Begin Elements Element3D3N
		1     1    1    2    3
		2     1    3    4    1
		3     1    5    1    4
		4     1    4    6    5
		5     1    7    5    6
		6     1    6    8    7
		7     1    7    8    3
		8     1    3    2    7
		9     1    7    2    1
		10    1    1    5    7
		11    1    8    6    4
		12    1    4    3    8
End Elements
Begin SubModelPart workpiece
Begin SubModelPartNodes
		1
		2
		3
		4
		5
		6
		7
		8
End SubModelPartNodes
Begin SubModelPartElements
		1
		2
		3
		4
		5
		6
		7
		8
		9
		10
		11
		12
End SubModelPartElements
End SubModelPart
	)input"));

    Model current_model;

    ModelPartIO model_part_io(p_input);

    ModelPart& model_part_0 =
        current_model.CreateModelPart("tmp_main_model_part");
    model_part_io.ReadModelPart(model_part_0);

    // Create the output .mdpa file
    std::string output_file_name = "cube_skin_mesh";
    std::fstream output_file;
    output_file.open(output_file_name + ".mdpa", std::fstream::out);
    output_file.close();

    // Fill the output .mdpa file
    ModelPartIO* model_part_io_write = new ModelPartIO(output_file_name, IO::WRITE);
    model_part_io_write->WriteModelPart(model_part_0);
}
} // namespace

KRATOS_TEST_CASE_IN_SUITE(VoxelMeshGeneratorModelerKeyPlaneBySize,
                          KratosCoreFastSuite)
{
    using namespace Kratos;

	WriteCubeSkinMeshMdpaFileForVoxelModelerTest();

    Parameters mesher_parameters(R"(
    {
        "output_model_part_name" : "main_model_part",
        "input_model_part_name" : "skin_model_part",
        "mdpa_file_name" : "cube_skin_mesh",
        "key_plane_generator": {
            "Parameters" : {
                "voxel_sizes" : [0.01, 0.025, 0.03],
                "min_point" : [0.01, 0.02, 0.03],
                "max_point" : [0.06, 0.08, 0.08]
            }
        }
    })");

    Model current_model;
    current_model.CreateModelPart("main_model_part");

    // Generate the skin
    current_model.CreateModelPart("skin_model_part");

    // Generating the mesh
    auto voxel_mesher = VoxelMeshGeneratorModeler(current_model, mesher_parameters);
    voxel_mesher.SetupGeometryModel();
    voxel_mesher.SetupModelPart();

    std::vector<double> x_key_planes{0.01, 0.02, 0.03, 0.04, 0.05, 0.06};
    std::vector<double> y_key_planes{0.02, 0.05, 0.08};
    std::vector<double> z_key_planes{0.03, 0.055, 0.08};

    KRATOS_EXPECT_VECTOR_NEAR(x_key_planes, voxel_mesher.GetKeyPlanes(0), 1e-6);
    KRATOS_EXPECT_VECTOR_NEAR(y_key_planes, voxel_mesher.GetKeyPlanes(1), 1e-6);
    KRATOS_EXPECT_VECTOR_NEAR(z_key_planes, voxel_mesher.GetKeyPlanes(2), 1e-6);

}

KRATOS_TEST_CASE_IN_SUITE(XCartesianRayPlaneIntersection, KratosCoreFastSuite)
	{
		Point::Pointer p_point_1=Kratos::make_shared<Point>(.4, 0.00, 0.00);
		Point::Pointer p_point_2=Kratos::make_shared<Point>(.4, 1.00, 0.00);
		Point::Pointer p_point_3=Kratos::make_shared<Point>(.4, 1.00, 1.00);
		Point::Pointer p_point_4=Kratos::make_shared<Point>(.4, 0.00, 1.00);

		Triangle3D3<Point> triangle_1(p_point_1, p_point_2, p_point_3);
		Triangle3D3<Point> triangle_2(p_point_1, p_point_3, p_point_4);

		Kratos::Internals::CartesianRay<Geometry<Point>> ray_1(0, Point(0.00, 0.50, 0.50), Point(1.00, 0.50, 0.50));
		Kratos::Internals::CartesianRay<Geometry<Point>> ray_2(0, Point(0.00, 0.00, 0.00), Point(1.00, 0.00, 0.00));
		Kratos::Internals::CartesianRay<Geometry<Point>> ray_3(0, Point(0.00, 0.20, 0.70), Point(1.00, 0.20, 0.70));

		ray_1.AddIntersection(triangle_1, 1e-9);
		ray_2.AddIntersection(triangle_1, 1e-9);
		ray_3.AddIntersection(triangle_1, 1e-9);

		KRATOS_EXPECT_EQ(ray_1.GetIntersections().size(), 1);
		KRATOS_EXPECT_EQ(ray_2.GetIntersections().size(), 1);
		KRATOS_EXPECT_EQ(ray_3.GetIntersections().size(), 0);

		KRATOS_EXPECT_NEAR(ray_1.GetIntersections()[0].first, .4, 1e-6);
		KRATOS_EXPECT_NEAR(ray_2.GetIntersections()[0].first, .4, 1e-6);

		ray_1.AddIntersection(triangle_2, 1e-9);
		ray_2.AddIntersection(triangle_2, 1e-9);
		ray_3.AddIntersection(triangle_2, 1e-9);

		KRATOS_EXPECT_EQ(ray_1.GetIntersections().size(), 2);
		KRATOS_EXPECT_EQ(ray_2.GetIntersections().size(), 2);
		KRATOS_EXPECT_EQ(ray_3.GetIntersections().size(), 1);

		KRATOS_EXPECT_NEAR(ray_1.GetIntersections()[1].first, .4, 1e-6);
		KRATOS_EXPECT_NEAR(ray_2.GetIntersections()[1].first, .4, 1e-6);
		KRATOS_EXPECT_NEAR(ray_3.GetIntersections()[0].first, .4, 1e-6);

		ray_1.CollapseIntersectionPoints(1e-6);
		ray_2.CollapseIntersectionPoints(1e-6);

		KRATOS_EXPECT_EQ(ray_1.GetIntersections().size(), 1);
		KRATOS_EXPECT_EQ(ray_2.GetIntersections().size(), 0);

		KRATOS_EXPECT_NEAR(ray_1.GetIntersections()[0].first, .4, 1e-6);
		KRATOS_EXPECT_NEAR(ray_2.GetIntersections()[0].first, .4, 1e-6);
	}


	KRATOS_TEST_CASE_IN_SUITE(XCartesianRayLongPlaneIntersection, KratosCoreFastSuite)
	{
		Point::Pointer p_point_1=Kratos::make_shared<Point>(.4, 0.00, 0.00);
		Point::Pointer p_point_2=Kratos::make_shared<Point>(.4, 1.00, 0.00);
		Point::Pointer p_point_3=Kratos::make_shared<Point>(.4, 1.00, 1000.00);
		Point::Pointer p_point_4=Kratos::make_shared<Point>(.4, 0.00, 1000.00);

		Triangle3D3<Point> triangle_1(p_point_1, p_point_2, p_point_3);
		Triangle3D3<Point> triangle_2(p_point_1, p_point_3, p_point_4);

		Kratos::Internals::CartesianRay<Geometry<Point>> ray_1(0, Point(0.00, 0.50, 500.0), Point(1.00, 0.50, 500.0));
		Kratos::Internals::CartesianRay<Geometry<Point>> ray_2(0, Point(0.00, 0.00, 0.00), Point(1.00, 0.00, 0.00));
		Kratos::Internals::CartesianRay<Geometry<Point>> ray_3(0, Point(0.00, 0.20, 700.0), Point(1.00, 0.20, 700.0));

		ray_1.AddIntersection(triangle_1, 1e-9);
		ray_2.AddIntersection(triangle_1, 1e-9);
		ray_3.AddIntersection(triangle_1, 1e-9);

		KRATOS_EXPECT_EQ(ray_1.GetIntersections().size(), 1);
		KRATOS_EXPECT_EQ(ray_2.GetIntersections().size(), 1);
		KRATOS_EXPECT_EQ(ray_3.GetIntersections().size(), 0);

		KRATOS_EXPECT_NEAR(ray_1.GetIntersections()[0].first, .4, 1e-6);
		KRATOS_EXPECT_NEAR(ray_2.GetIntersections()[0].first, .4, 1e-6);

		ray_1.AddIntersection(triangle_2, 1e-9);
		ray_2.AddIntersection(triangle_2, 1e-9);
		ray_3.AddIntersection(triangle_2, 1e-9);

		KRATOS_EXPECT_EQ(ray_1.GetIntersections().size(), 2);
		KRATOS_EXPECT_EQ(ray_2.GetIntersections().size(), 2);
		KRATOS_EXPECT_EQ(ray_3.GetIntersections().size(), 1);

		KRATOS_EXPECT_NEAR(ray_1.GetIntersections()[1].first, .4, 1e-6);
		KRATOS_EXPECT_NEAR(ray_2.GetIntersections()[1].first, .4, 1e-6);
		KRATOS_EXPECT_NEAR(ray_3.GetIntersections()[0].first, .4, 1e-6);

		ray_1.CollapseIntersectionPoints(1e-6);
		ray_2.CollapseIntersectionPoints(1e-6);

		KRATOS_EXPECT_EQ(ray_1.GetIntersections().size(), 1);
		KRATOS_EXPECT_EQ(ray_2.GetIntersections().size(), 0);

		KRATOS_EXPECT_NEAR(ray_1.GetIntersections()[0].first, .4, 1e-6);
		KRATOS_EXPECT_NEAR(ray_2.GetIntersections()[0].first, .4, 1e-6);
	}


	KRATOS_TEST_CASE_IN_SUITE(XCartesianRayPlaneMarkIntersectedIntervals, KratosCoreFastSuite)
	{
		Point::Pointer p_point_1=Kratos::make_shared<Point>(.4, 0.00, 0.00);
		Point::Pointer p_point_2=Kratos::make_shared<Point>(.4, 1.00, 0.00);
		Point::Pointer p_point_3=Kratos::make_shared<Point>(.4, 1.00, 1.00);
		Point::Pointer p_point_4=Kratos::make_shared<Point>(.4, 0.00, 1.00);

		Triangle3D3<Point> triangle_1(p_point_1, p_point_2, p_point_3);
		Triangle3D3<Point> triangle_2(p_point_1, p_point_3, p_point_4);

		Kratos::Internals::CartesianRay<Geometry<Point>> ray_1(0, Point(0.00, 0.50, 0.50), Point(1.00, 0.50, 0.50));

		ray_1.AddIntersection(triangle_1, 1e-9);
		ray_1.AddIntersection(triangle_2, 1e-9);

		ray_1.CollapseIntersectionPoints(1e-6);

		std::vector<double> colors;
		std::vector<double> coordinates_1{0.05, 0.15, 0.25, 0.30, 0.45, 0.50, 0.60, 0.70, 0.80, 0.90};
		ray_1.MarkIntersectedIntervals(coordinates_1, -1, 1, colors, 1e-6);

		for(std::size_t i = 0 ; i < colors.size() ; i++){
			if(i == 3){
				KRATOS_EXPECT_NEAR(colors[i], -1, 1e-6);
			}
			else {
				KRATOS_EXPECT_NEAR(colors[i], 1, 1e-6);
			}
		}

		std::vector<double> coordinates_2{0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90};
		ray_1.MarkIntersectedIntervals(coordinates_2, -1, 1, colors, 1e-6);

		for(std::size_t i = 0 ; i < colors.size() ; i++){
			if(i == 3){
				KRATOS_EXPECT_NEAR(colors[i], -1, 1e-6);
			}
			else {
				KRATOS_EXPECT_NEAR(colors[i], 1, 1e-6);
			}
		}

	}


	KRATOS_TEST_CASE_IN_SUITE(CartesianMeshColorsSetCoordinates, KratosCoreFastSuite)
	{
		Kratos::Internals::CartesianMeshColors mesh_colors;

		array_1d<std::vector<double>, 3> coordinates;
		for(std::size_t i = 0 ; i < 3 ; i++){
			coordinates[i].resize(6);
		}

		Point min_point(1.00, 2.00, 3.00);
		for(std::size_t i = 0 ; i < coordinates.size() ; i++){
			const double min_coordinate_i = min_point[i];
			for(std::size_t j = 0 ; j < coordinates[i].size() ; j++)
				coordinates[i][j] = j * 2.0 + min_coordinate_i;
		}

		mesh_colors.SetCoordinates(coordinates[0], coordinates[1], coordinates[2]);

		for (std::size_t k = 0; k <= 5; k++) {
			for (std::size_t j = 0; j <= 5; j++) {
				for (std::size_t i = 0; i <= 5; i++) {
					double x = 2.00*i + 1.00;
					double y = 2.00*j + 2.00;
					double z = 2.00*k + 3.00;
					Point point = mesh_colors.GetPoint(i,j,k);
                	KRATOS_EXPECT_NEAR(point.X(), x, 1e-6);
                	KRATOS_EXPECT_NEAR(point.Y(), y, 1e-6);
                	KRATOS_EXPECT_NEAR(point.Z(), z, 1e-6);
				}
            }
		}
	}

	KRATOS_TEST_CASE_IN_SUITE(CartesianMeshColorsInitialRays, KratosCoreFastSuite)
	{
		std::size_t size = 4;
		Kratos::Internals::CartesianMeshColors mesh_colors;

		array_1d<std::vector<double>, 3> coordinates;
		for(std::size_t i = 0 ; i < 3 ; i++){
			coordinates[i].resize(size);
		}

		Point min_point(1.00, 2.00, 3.00);
		for(std::size_t i = 0 ; i <  3 ; i++){
			const double min_coordinate_i = min_point[i];
			for(std::size_t j = 0 ; j < size ; j++)
				coordinates[i][j] = j * 2.0 + min_coordinate_i;
		}

		mesh_colors.SetCoordinates(coordinates[0], coordinates[1], coordinates[2]);

		array_1d<std::size_t, 3> min_position;
		array_1d<std::size_t, 3> max_position;

		for(std::size_t i = 0 ; i < 3 ; i++){
			min_position[i] = 0;
			max_position[i] = size;
		}

		mesh_colors.InitializeRays(min_position, max_position, "nodes");

		for (std::size_t i = 0; i < size; i++) {
			for (std::size_t j = 0; j < size; j++) {
				double x = 2.00*i + 1.00;
				double y = 2.00*j + 2.00;
				double z = 2.00*(size - 1.00) + 3.00;
				auto& ray = mesh_colors.GetXYRay(i,j);
				Point& point = ray.GetPoint1();
				KRATOS_EXPECT_NEAR(point.X(), x, 1e-6);
				KRATOS_EXPECT_NEAR(point.Y(), y, 1e-6);
				KRATOS_EXPECT_NEAR(point.Z(), min_point.Z(), 1e-6);

				point = ray.GetPoint2();
				KRATOS_EXPECT_NEAR(point.X(), x, 1e-6);
				KRATOS_EXPECT_NEAR(point.Y(), y, 1e-6);
				KRATOS_EXPECT_NEAR(point.Z(), z, 1e-6);
           }
		}
	}

	KRATOS_TEST_CASE_IN_SUITE(CartesianMeshColorsAddGeometry, KratosCoreFastSuite)
	{

        Model current_model;

		// Generate the skin
		ModelPart &skin_model_part = current_model.CreateModelPart("Skin");

		std::size_t size = 5;
		double cube_size = static_cast<double>(size - 1);

		skin_model_part.AddNodalSolutionStepVariable(VELOCITY);
		skin_model_part.CreateNewNode(901, 0.0, 0.0, 0.0);
		skin_model_part.CreateNewNode(902, cube_size, 0.0, 0.0);
		skin_model_part.CreateNewNode(903, cube_size, cube_size, 0.0);
		skin_model_part.CreateNewNode(904, 0.0, 0.0, cube_size);
		Properties::Pointer p_properties(new Properties(0));
		skin_model_part.CreateNewElement("Element3D3N", 901, { 901,902,903 }, p_properties);
		skin_model_part.CreateNewElement("Element3D3N", 902, { 901,904,903 }, p_properties);
		skin_model_part.CreateNewElement("Element3D3N", 903, { 902,903,904 }, p_properties);
		skin_model_part.CreateNewElement("Element3D3N", 904, { 901,902,904 }, p_properties);
		Kratos::Internals::CartesianMeshColors mesh_colors;

		array_1d<std::vector<double>, 3> coordinates;
		for(std::size_t i = 0 ; i < 3 ; i++){
			coordinates[i].resize(size);
		}

		for(std::size_t i = 0 ; i <  3 ; i++){
			for(std::size_t j = 0 ; j < size ; j++)
				coordinates[i][j] = static_cast<double>(j);
		}

		mesh_colors.SetCoordinates(coordinates[0], coordinates[1], coordinates[2]);

		mesh_colors.ExtendBoundingBox(skin_model_part.Nodes(), 1.0e-6);

		array_1d<std::size_t, 3> min_position;
		array_1d<std::size_t, 3> max_position;

		for(std::size_t i = 0 ; i < 3 ; i++){
			min_position[i] = 0;
			max_position[i] = size;
		}

		mesh_colors.InitializeRays(min_position, max_position, "nodes");

        for(auto& element : skin_model_part.Elements())
        {
            Element::GeometryType& r_geometry = element.GetGeometry();
			mesh_colors.AddGeometry(r_geometry, true);
        }

		for (std::size_t i = 0; i < size; i++) {
			for (std::size_t j = 0; j < size; j++) {
				auto& ray = mesh_colors.GetXYRay(i,j);
				if(i>=j){
					KRATOS_EXPECT_EQ(ray.GetIntersections().size(), 2);
				}
				else {
					KRATOS_EXPECT_EQ(ray.GetIntersections().size(), 0);
				}
           }
		}
	}


	KRATOS_TEST_CASE_IN_SUITE(CartesianMeshColorsXPlaneFaceColoring, KratosCoreFastSuite)
	{

        Model current_model;

		// Generate the skin
		ModelPart &skin_model_part = current_model.CreateModelPart("Skin");

		skin_model_part.CreateNewNode(1, .4, 0.00, 0.00);
		skin_model_part.CreateNewNode(2, .4, 1.00, 0.00);
		skin_model_part.CreateNewNode(3, .4, 1.00, 1.00);
		skin_model_part.CreateNewNode(4, .4, 0.00, 1.00);
		Properties::Pointer p_properties(new Properties(0));
		skin_model_part.CreateNewElement("Element3D3N", 901, { 1,2,3 }, p_properties);
		skin_model_part.CreateNewElement("Element3D3N", 902, { 1,4,3 }, p_properties);

		std::size_t size = 5;
		Kratos::Internals::CartesianMeshColors mesh_colors;

		array_1d<std::vector<double>, 3> coordinates;
		for(std::size_t i = 0 ; i < 3 ; i++){
			coordinates[i].resize(size);
		}

		for(std::size_t i = 0 ; i <  3 ; i++){
			for(std::size_t j = 0 ; j < size ; j++)
				coordinates[i][j] = static_cast<double>(j)/static_cast<double>(size);
		}

		array_1d<std::size_t, 3> min_position;
		array_1d<std::size_t, 3> max_position;

		for(std::size_t i = 0 ; i < 3 ; i++){
			min_position[i] = 0;
			max_position[i] = size-1;
		}

		mesh_colors.SetCoordinates(coordinates[0], coordinates[1], coordinates[2]);
		mesh_colors.InitializeRays(min_position, max_position, "face_of_elements");
		mesh_colors.SetAllColors(1);

        for(auto& element : skin_model_part.Elements())
        {
            Element::GeometryType& r_geometry = element.GetGeometry();
			mesh_colors.AddGeometry(r_geometry, false);
        }

		for (std::size_t k = 0; k < size - 1; k++) {
			for (std::size_t j = 0; j < size - 1; j++) {
				for (std::size_t i = 0; i < size -1; i++) {
					if(coordinates[0][i] >= .4){
						mesh_colors.GetElementalColor(i,j,k) = -2;
					}
				}
			}
		}

		mesh_colors.CalculateElementalFaceColors(min_position, max_position, -1, 1, -2);


		for (std::size_t k = 0; k < size - 1; k++) {
			for (std::size_t j = 0; j < size - 1; j++) {
				for (std::size_t i = 0; i < size -1; i++) {
					auto colors = mesh_colors.GetElementalFaceColor(i,j,k);
					if(i == 2){
						KRATOS_EXPECT_NEAR(colors[0], -1, 1e-6);
					}
					else{
						KRATOS_EXPECT_NEAR(colors[0], 1, 1e-6);
					}
					KRATOS_EXPECT_NEAR(colors[1], 1, 1e-6);
					KRATOS_EXPECT_NEAR(colors[2], 1, 1e-6);
					KRATOS_EXPECT_NEAR(colors[3], 1, 1e-6);
					KRATOS_EXPECT_NEAR(colors[4], 1, 1e-6);
					KRATOS_EXPECT_NEAR(colors[5], 1, 1e-6);
				}
            }
		}

	}

	KRATOS_TEST_CASE_IN_SUITE(CartesianMeshColorsYPlaneFaceColoring, KratosCoreFastSuite)
	{
		Node::Pointer p_point_1=Kratos::make_intrusive<Node>(1, 0.00, .4, 0.00);
		Node::Pointer p_point_2=Kratos::make_intrusive<Node>(2, 1.00, .4, 0.00);
		Node::Pointer p_point_3=Kratos::make_intrusive<Node>(3, 1.00, .4, 1.00);
		Node::Pointer p_point_4=Kratos::make_intrusive<Node>(4, 0.00, .4, 1.00);

		Triangle3D3<Node> triangle_1(p_point_1, p_point_2, p_point_3);
		Triangle3D3<Node> triangle_2(p_point_1, p_point_3, p_point_4);

		std::size_t size = 5;
		Kratos::Internals::CartesianMeshColors mesh_colors;

		array_1d<std::vector<double>, 3> coordinates;
		for(std::size_t i = 0 ; i < 3 ; i++){
			coordinates[i].resize(size);
		}

		for(std::size_t i = 0 ; i <  3 ; i++){
			for(std::size_t j = 0 ; j < size ; j++)
				coordinates[i][j] = static_cast<double>(j)/static_cast<double>(size);
		}

		array_1d<std::size_t, 3> min_position;
		array_1d<std::size_t, 3> max_position;

		for(std::size_t i = 0 ; i < 3 ; i++){
			min_position[i] = 0;
			max_position[i] = size-1;
		}

		mesh_colors.SetCoordinates(coordinates[0], coordinates[1], coordinates[2]);
		mesh_colors.InitializeRays(min_position, max_position, "face_of_elements");
		mesh_colors.SetAllColors(1);
		mesh_colors.AddGeometry(triangle_1,false);
		mesh_colors.AddGeometry(triangle_2, false);

		for (std::size_t k = 0; k < size - 1; k++) {
			for (std::size_t j = 0; j < size - 1; j++) {
				for (std::size_t i = 0; i < size -1; i++) {
					if(coordinates[1][j] >= .4){
						mesh_colors.GetElementalColor(i,j,k) = -2;
					}
				}
			}
		}

		mesh_colors.CalculateElementalFaceColors(min_position, max_position, -1, 1, -2);


		for (std::size_t k = 0; k < size - 1; k++) {
			for (std::size_t j = 0; j < size - 1; j++) {
				for (std::size_t i = 0; i < size -1; i++) {
					auto colors = mesh_colors.GetElementalFaceColor(i,j,k);
					if(j == 2){
						KRATOS_EXPECT_NEAR(colors[1], -1, 1e-6);
					}
					else{
						KRATOS_EXPECT_NEAR(colors[1], 1, 1e-6);
					}
					KRATOS_EXPECT_NEAR(colors[0], 1, 1e-6);
					KRATOS_EXPECT_NEAR(colors[2], 1, 1e-6);
					KRATOS_EXPECT_NEAR(colors[3], 1, 1e-6);
					KRATOS_EXPECT_NEAR(colors[4], 1, 1e-6);
					KRATOS_EXPECT_NEAR(colors[5], 1, 1e-6);
				}
            }
		}

	}

	KRATOS_TEST_CASE_IN_SUITE(CartesianMeshColorsZPlaneFaceColoring, KratosCoreFastSuite)
	{
		Node::Pointer p_point_1=Kratos::make_intrusive<Node>(1, 0.00, 0.00, .4);
		Node::Pointer p_point_2=Kratos::make_intrusive<Node>(2, 1.00, 0.00, .4);
		Node::Pointer p_point_3=Kratos::make_intrusive<Node>(3, 1.00, 1.00, .4);
		Node::Pointer p_point_4=Kratos::make_intrusive<Node>(4, 0.00, 1.00, .4);

		Triangle3D3<Node> triangle_1(p_point_1, p_point_2, p_point_3);
		Triangle3D3<Node> triangle_2(p_point_1, p_point_3, p_point_4);

		std::size_t size = 5;
		Kratos::Internals::CartesianMeshColors mesh_colors;

		array_1d<std::vector<double>, 3> coordinates;
		for(std::size_t i = 0 ; i < 3 ; i++){
			coordinates[i].resize(size);
		}

		for(std::size_t i = 0 ; i <  3 ; i++){
			for(std::size_t j = 0 ; j < size ; j++)
				coordinates[i][j] = static_cast<double>(j)/static_cast<double>(size);
		}

		array_1d<std::size_t, 3> min_position;
		array_1d<std::size_t, 3> max_position;

		for(std::size_t i = 0 ; i < 3 ; i++){
			min_position[i] = 0;
			max_position[i] = size-1;
		}

		mesh_colors.SetCoordinates(coordinates[0], coordinates[1], coordinates[2]);
		mesh_colors.InitializeRays(min_position, max_position, "face_of_elements");
		mesh_colors.SetAllColors(1);
		mesh_colors.AddGeometry(triangle_1,false);
		mesh_colors.AddGeometry(triangle_2, false);

		for (std::size_t k = 0; k < size - 1; k++) {
			for (std::size_t j = 0; j < size - 1; j++) {
				for (std::size_t i = 0; i < size -1; i++) {
					if(coordinates[2][k] >= .4){
						mesh_colors.GetElementalColor(i,j,k) = -2;
					}
				}
			}
		}

		mesh_colors.CalculateElementalFaceColors(min_position, max_position, -1, 1, -2);


		for (std::size_t k = 0; k < size - 1; k++) {
			for (std::size_t j = 0; j < size - 1; j++) {
				for (std::size_t i = 0; i < size -1; i++) {
					auto colors = mesh_colors.GetElementalFaceColor(i,j,k);
					if(k == 2){
						KRATOS_EXPECT_NEAR(colors[2], -1, 1e-6);
					}
					else{
						KRATOS_EXPECT_NEAR(colors[2], 1, 1e-6);
					}
					KRATOS_EXPECT_NEAR(colors[0], 1, 1e-6);
					KRATOS_EXPECT_NEAR(colors[1], 1, 1e-6);
					KRATOS_EXPECT_NEAR(colors[3], 1, 1e-6);
					KRATOS_EXPECT_NEAR(colors[4], 1, 1e-6);
					KRATOS_EXPECT_NEAR(colors[5], 1, 1e-6);
				}
            }
		}

	}

	KRATOS_TEST_CASE_IN_SUITE(VoxelMesherWithVectorDistances, KratosCoreFastSuite)
	{
		using namespace Kratos;

		WriteCubeSkinMeshMdpaFileForVoxelModelerTest();
		std::cout << "Cube writen" << std::endl;

		Parameters mesher_parameters(R"(
			{
				"output_model_part_name" : "main_model_part",
				"input_model_part_name" : "skin_model_part",
				"mdpa_file_name" : "cube_skin_mesh",
				"key_plane_generator": {
					"Parameters" : {
						"voxel_sizes" : [0.025, 0.025, 0.025],
						"min_point" : [-0.05, -0.05, 0],
						"max_point" : [0.05, 0.05, 0.1]
					}
				},
				"coloring_settings_list": [
				{
					"type" : "cells_in_touch",
					"model_part_name": "skin_model_part.workpiece",
					"color": 14
				},
				{
					"type" : "cells_with_inside_center",
					"model_part_name": "skin_model_part.workpiece",
					"color": 14
				}
				],
				"entities_generator_list": [
				{
					"type" : "elements_with_cell_color",
					"model_part_name": "main_model_part.workpiece",
					"color": 14,
					"properties_id": 1
				} 
				]
			})");
			
		Model current_model;
		current_model.CreateModelPart("main_model_part");

		// Generate the skin
		current_model.CreateModelPart("skin_model_part");

		std::cout << "Modelpart created" << std::endl;

		// Generating the mesh
		auto voxel_mesher = SurrogateBoundaryModeler(current_model, mesher_parameters);
		voxel_mesher.SetupGeometryModel();
		voxel_mesher.SetupModelPart();

		std::cout << "SurrogateBoundaryModeler created" << std::endl;

		//voxel_mesher.FindDistanceToSkin();
		//voxel_mesher.ApplyColoringToNodes();

		voxel_mesher.ComputeSurrogateBoundary();

		std::cout << "Distances and colors computed" << std::endl;
		
		auto& nodalSBdata = voxel_mesher.GetSurrogateBoundaryData();

		for (SurrogateBoundaryModeler::SurrogateBoundaryNode& node : nodalSBdata) 
		{
			if (node.IsActive()) 
			{
				std::cout << *node.GetNodePtr() 
						  << " \n  Vector distance to skin: " << node.GetVectorDistance()
						  << " \n  Signed distance to skin: " << node.GetSignedDistance() 
						  << " \n  Is inside: " << node.IsInside();
			}
		}
	}

} // namespace Kratos::Testing
