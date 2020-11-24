//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//					 Ruben Zorrilla
//

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/checks.h"
// #include "includes/gid_io.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/tetrahedra_3d_4.h"
#include "processes/coarse_voxel_mesh_generator_process.h"

namespace Kratos {
  namespace Testing {


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

		KRATOS_CHECK_EQUAL(ray_1.GetIntersections().size(), 1);
		KRATOS_CHECK_EQUAL(ray_2.GetIntersections().size(), 1);
		KRATOS_CHECK_EQUAL(ray_3.GetIntersections().size(), 0);

		KRATOS_CHECK_NEAR(ray_1.GetIntersections()[0].first, .4, 1e-6);
		KRATOS_CHECK_NEAR(ray_2.GetIntersections()[0].first, .4, 1e-6);

		ray_1.AddIntersection(triangle_2, 1e-9);
		ray_2.AddIntersection(triangle_2, 1e-9);
		ray_3.AddIntersection(triangle_2, 1e-9);

		KRATOS_CHECK_EQUAL(ray_1.GetIntersections().size(), 2);
		KRATOS_CHECK_EQUAL(ray_2.GetIntersections().size(), 2);
		KRATOS_CHECK_EQUAL(ray_3.GetIntersections().size(), 1);

		KRATOS_CHECK_NEAR(ray_1.GetIntersections()[1].first, .4, 1e-6);
		KRATOS_CHECK_NEAR(ray_2.GetIntersections()[1].first, .4, 1e-6);
		KRATOS_CHECK_NEAR(ray_3.GetIntersections()[0].first, .4, 1e-6);

		ray_1.CollapseIntersectionPoints(1e-6);
		ray_2.CollapseIntersectionPoints(1e-6);

		KRATOS_CHECK_EQUAL(ray_1.GetIntersections().size(), 1);
		KRATOS_CHECK_EQUAL(ray_2.GetIntersections().size(), 0);
		
		KRATOS_CHECK_NEAR(ray_1.GetIntersections()[0].first, .4, 1e-6);
		KRATOS_CHECK_NEAR(ray_2.GetIntersections()[0].first, .4, 1e-6);
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

		KRATOS_CHECK_EQUAL(ray_1.GetIntersections().size(), 1);
		KRATOS_CHECK_EQUAL(ray_2.GetIntersections().size(), 1);
		KRATOS_CHECK_EQUAL(ray_3.GetIntersections().size(), 0);

		KRATOS_CHECK_NEAR(ray_1.GetIntersections()[0].first, .4, 1e-6);
		KRATOS_CHECK_NEAR(ray_2.GetIntersections()[0].first, .4, 1e-6);

		ray_1.AddIntersection(triangle_2, 1e-9);
		ray_2.AddIntersection(triangle_2, 1e-9);
		ray_3.AddIntersection(triangle_2, 1e-9);

		KRATOS_CHECK_EQUAL(ray_1.GetIntersections().size(), 2);
		KRATOS_CHECK_EQUAL(ray_2.GetIntersections().size(), 2);
		KRATOS_CHECK_EQUAL(ray_3.GetIntersections().size(), 1);

		KRATOS_CHECK_NEAR(ray_1.GetIntersections()[1].first, .4, 1e-6);
		KRATOS_CHECK_NEAR(ray_2.GetIntersections()[1].first, .4, 1e-6);
		KRATOS_CHECK_NEAR(ray_3.GetIntersections()[0].first, .4, 1e-6);

		ray_1.CollapseIntersectionPoints(1e-6);
		ray_2.CollapseIntersectionPoints(1e-6);

		KRATOS_CHECK_EQUAL(ray_1.GetIntersections().size(), 1);
		KRATOS_CHECK_EQUAL(ray_2.GetIntersections().size(), 0);
		
		KRATOS_CHECK_NEAR(ray_1.GetIntersections()[0].first, .4, 1e-6);
		KRATOS_CHECK_NEAR(ray_2.GetIntersections()[0].first, .4, 1e-6);
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
				KRATOS_CHECK_NEAR(colors[i], -1, 1e-6);
			}
			else {
				KRATOS_CHECK_NEAR(colors[i], 1, 1e-6);
			}
		}
		
		std::vector<double> coordinates_2{0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90};
		ray_1.MarkIntersectedIntervals(coordinates_2, -1, 1, colors, 1e-6);

		for(std::size_t i = 0 ; i < colors.size() ; i++){
			if(i == 3){
				KRATOS_CHECK_NEAR(colors[i], -1, 1e-6);
			}
			else {
				KRATOS_CHECK_NEAR(colors[i], 1, 1e-6);
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
                	KRATOS_CHECK_NEAR(point.X(), x, 1e-6);
                	KRATOS_CHECK_NEAR(point.Y(), y, 1e-6);
                	KRATOS_CHECK_NEAR(point.Z(), z, 1e-6);
				}
            }
		}
	}

	KRATOS_TEST_CASE_IN_SUITE(CartesianMeshCalculateCenterOfElementPosition, KratosCoreFastSuite)
	{
		Kratos::Internals::CartesianMeshColors mesh_colors;
		
		array_1d<std::vector<double>, 3> coordinates;
		coordinates[0] = {1.00,2.5,5.0};
		coordinates[1] =  {0.00, 1.00};
		coordinates[2] =  {2.00, 3.00};
	
		mesh_colors.SetCoordinates(coordinates[0], coordinates[1], coordinates[2]);

		KRATOS_CHECK_EQUAL(mesh_colors.CalculateCenterOfElementPosition(0.5, 0), 0);
		KRATOS_CHECK_EQUAL(mesh_colors.CalculateCenterOfElementPosition(0.99999999, 0), 0);
		KRATOS_CHECK_EQUAL(mesh_colors.CalculateCenterOfElementPosition(1.1, 0), 0);
		KRATOS_CHECK_EQUAL(mesh_colors.CalculateCenterOfElementPosition(1.74999999, 0), 0);
		KRATOS_CHECK_EQUAL(mesh_colors.CalculateCenterOfElementPosition(1.75, 0), 0);
		KRATOS_CHECK_EQUAL(mesh_colors.CalculateCenterOfElementPosition(1.751, 0), 1);
		KRATOS_CHECK_EQUAL(mesh_colors.CalculateCenterOfElementPosition(2.499999, 0), 1);
		KRATOS_CHECK_EQUAL(mesh_colors.CalculateCenterOfElementPosition(2.5, 0), 1);
		KRATOS_CHECK_EQUAL(mesh_colors.CalculateCenterOfElementPosition(2.51, 0), 1);
		KRATOS_CHECK_EQUAL(mesh_colors.CalculateCenterOfElementPosition(4.999, 0), 1);
		KRATOS_CHECK_EQUAL(mesh_colors.CalculateCenterOfElementPosition(5.0, 0), 1);
		KRATOS_CHECK_EQUAL(mesh_colors.CalculateCenterOfElementPosition(10., 0), 1);
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
				KRATOS_CHECK_NEAR(point.X(), x, 1e-6);
				KRATOS_CHECK_NEAR(point.Y(), y, 1e-6);
				KRATOS_CHECK_NEAR(point.Z(), min_point.Z(), 1e-6);
 				
				point = ray.GetPoint2();
				KRATOS_CHECK_NEAR(point.X(), x, 1e-6);
				KRATOS_CHECK_NEAR(point.Y(), y, 1e-6);
				KRATOS_CHECK_NEAR(point.Z(), z, 1e-6);
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
					KRATOS_CHECK_EQUAL(ray.GetIntersections().size(), 2);
				}
				else {
					KRATOS_CHECK_EQUAL(ray.GetIntersections().size(), 0);
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
						KRATOS_CHECK_NEAR(colors[0], -1, 1e-6);
					}
					else{
						KRATOS_CHECK_NEAR(colors[0], 1, 1e-6);
					}
					KRATOS_CHECK_NEAR(colors[1], 1, 1e-6);
					KRATOS_CHECK_NEAR(colors[2], 1, 1e-6);
					KRATOS_CHECK_NEAR(colors[3], 1, 1e-6);
					KRATOS_CHECK_NEAR(colors[4], 1, 1e-6);
					KRATOS_CHECK_NEAR(colors[5], 1, 1e-6);
				}
            }
		}

	}

	KRATOS_TEST_CASE_IN_SUITE(CartesianMeshColorsYPlaneFaceColoring, KratosCoreFastSuite)
	{
		Node<3>::Pointer p_point_1=Kratos::make_intrusive<Node<3>>(1, 0.00, .4, 0.00);
		Node<3>::Pointer p_point_2=Kratos::make_intrusive<Node<3>>(2, 1.00, .4, 0.00);
		Node<3>::Pointer p_point_3=Kratos::make_intrusive<Node<3>>(3, 1.00, .4, 1.00);
		Node<3>::Pointer p_point_4=Kratos::make_intrusive<Node<3>>(4, 0.00, .4, 1.00);

		Triangle3D3<Node<3>> triangle_1(p_point_1, p_point_2, p_point_3);
		Triangle3D3<Node<3>> triangle_2(p_point_1, p_point_3, p_point_4);

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
						KRATOS_CHECK_NEAR(colors[1], -1, 1e-6);
					}
					else{
						KRATOS_CHECK_NEAR(colors[1], 1, 1e-6);
					}
					KRATOS_CHECK_NEAR(colors[0], 1, 1e-6);
					KRATOS_CHECK_NEAR(colors[2], 1, 1e-6);
					KRATOS_CHECK_NEAR(colors[3], 1, 1e-6);
					KRATOS_CHECK_NEAR(colors[4], 1, 1e-6);
					KRATOS_CHECK_NEAR(colors[5], 1, 1e-6);
				}
            }
		}

	}

	KRATOS_TEST_CASE_IN_SUITE(CartesianMeshColorsZPlaneFaceColoring, KratosCoreFastSuite)
	{
		Node<3>::Pointer p_point_1=Kratos::make_intrusive<Node<3>>(1, 0.00, 0.00, .4);
		Node<3>::Pointer p_point_2=Kratos::make_intrusive<Node<3>>(2, 1.00, 0.00, .4);
		Node<3>::Pointer p_point_3=Kratos::make_intrusive<Node<3>>(3, 1.00, 1.00, .4);
		Node<3>::Pointer p_point_4=Kratos::make_intrusive<Node<3>>(4, 0.00, 1.00, .4);

		Triangle3D3<Node<3>> triangle_1(p_point_1, p_point_2, p_point_3);
		Triangle3D3<Node<3>> triangle_2(p_point_1, p_point_3, p_point_4);

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
						KRATOS_CHECK_NEAR(colors[2], -1, 1e-6);
					}
					else{
						KRATOS_CHECK_NEAR(colors[2], 1, 1e-6);
					}
					KRATOS_CHECK_NEAR(colors[0], 1, 1e-6);
					KRATOS_CHECK_NEAR(colors[1], 1, 1e-6);
					KRATOS_CHECK_NEAR(colors[3], 1, 1e-6);
					KRATOS_CHECK_NEAR(colors[4], 1, 1e-6);
					KRATOS_CHECK_NEAR(colors[5], 1, 1e-6);
				}
            }
		}

	}

	KRATOS_TEST_CASE_IN_SUITE(VoxelMeshGeneratorProcessNodesPositions, KratosCoreFastSuite)
	{
		Parameters mesher_parameters(R"(
		{
			"number_of_divisions":   [5,5,5],
			"element_name":     "Element3D4N"
		})");

        Model current_model;
		ModelPart &volume_part = current_model.CreateModelPart("Volume");

		// Generate the skin
		ModelPart &skin_model_part = current_model.CreateModelPart("Skin");


		// Generating the mesh
		VoxelMeshGeneratorProcess(Point{1.00, 2.00, 3.00}, Point{11.00, 12.00, 13.00}, volume_part, skin_model_part, mesher_parameters).Execute();

		auto i_node = volume_part.NodesBegin();
		for (std::size_t k = 0; k <= 5; k++) {
			for (std::size_t j = 0; j <= 5; j++) {
				for (std::size_t i = 0; i <= 5; i++) {
					auto& node = *i_node++;
					double x = 2.00*i + 1.00;
					double y = 2.00*j + 2.00;
					double z = 2.00*k + 3.00;
                	KRATOS_CHECK_NEAR(node.X(), x, 1e-6);
                	KRATOS_CHECK_NEAR(node.Y(), y, 1e-6);
                	KRATOS_CHECK_NEAR(node.Z(), z, 1e-6);
				}
            }
		}
	}


	KRATOS_TEST_CASE_IN_SUITE(VoxelMeshGeneratorProcessXYZPositions, KratosCoreFastSuite)
	{
		Parameters mesher_parameters(R"(
		{
			"number_of_divisions":   [5,5,5],
			"element_name":     "Element3D4N"
		})");

        Model current_model;
		ModelPart &volume_part = current_model.CreateModelPart("Volume");

		// Generate the skin
		ModelPart &skin_model_part = current_model.CreateModelPart("Skin");


		// Generating the mesh
		VoxelMeshGeneratorProcess({1.00, 3.00, 5.00, 7.00, 9.00, 11.00}, {2.00, 4.00, 6.00, 8.00, 10.00, 12.00}, {3.00, 5.00, 7.00, 9.00, 11.00, 13.00}, volume_part, skin_model_part, mesher_parameters).Execute();

		auto i_node = volume_part.NodesBegin();
		for (std::size_t k = 0; k <= 5; k++) {
			for (std::size_t j = 0; j <= 5; j++) {
				for (std::size_t i = 0; i <= 5; i++) {
					auto& node = *i_node++;
					double x = 2.00*i + 1.00;
					double y = 2.00*j + 2.00;
					double z = 2.00*k + 3.00;
                	KRATOS_CHECK_NEAR(node.X(), x, 1e-6);
                	KRATOS_CHECK_NEAR(node.Y(), y, 1e-6);
                	KRATOS_CHECK_NEAR(node.Z(), z, 1e-6);
				}
            }
		}
	}

	KRATOS_TEST_CASE_IN_SUITE(VoxelMeshGeneratorProcessRectilinearCoordinatesOutput, KratosCoreFastSuite)
	{
		Parameters mesher_parameters(R"(
		{
			"number_of_divisions":   [5,5,5],
			"element_name":     "Element3D4N",
			"output": "rectilinear_coordinates"
		})");

        Model current_model;
		ModelPart &volume_part = current_model.CreateModelPart("Volume");

		// Generate the skin
		ModelPart &skin_model_part = current_model.CreateModelPart("Skin");


		// Generating the mesh
		VoxelMeshGeneratorProcess({1.00, 3.00, 5.00, 7.00, 9.00, 11.00}, {2.00, 4.00, 6.00, 8.00, 10.00, 12.00}, {3.00, 5.00, 7.00, 9.00, 11.00, 13.00}, volume_part, skin_model_part, mesher_parameters).Execute();

		KRATOS_CHECK_EQUAL(volume_part.NumberOfNodes(), 0);
		KRATOS_CHECK_EQUAL(volume_part.NumberOfElements(), 0);
	
		for (std::size_t i = 0; i <= 5; i++) {
			double x = 2.00*i + 1.00;
			double y = 2.00*i + 2.00;
			double z = 2.00*i + 3.00;
			KRATOS_CHECK_NEAR(volume_part.GetValue(RECTILINEAR_X_COORDINATES)[i], x, 1e-6);
			KRATOS_CHECK_NEAR(volume_part.GetValue(RECTILINEAR_Y_COORDINATES)[i], y, 1e-6);
			KRATOS_CHECK_NEAR(volume_part.GetValue(RECTILINEAR_Z_COORDINATES)[i], z, 1e-6);
		}
	}

	KRATOS_TEST_CASE_IN_SUITE(VoxelMeshGeneratorProcessCentersPositions, KratosCoreFastSuite)
	{
		Parameters mesher_parameters(R"(
		{
			"number_of_divisions":   [10,10,10],
			"element_name":     "Element3D4N",
			"entities_to_generate": "center_of_elements"
		})");

        Model current_model;
		ModelPart &volume_part = current_model.CreateModelPart("Volume");

		// Generate the skin
		ModelPart &skin_model_part = current_model.CreateModelPart("Skin");


		// Generating the mesh
		VoxelMeshGeneratorProcess(Point{1.00, 2.00, 3.00}, Point{11.00, 12.00, 13.00}, volume_part, skin_model_part, mesher_parameters).Execute();

		auto i_node = volume_part.NodesBegin();
		for (std::size_t k = 0; k < 10; k++) {
			for (std::size_t j = 0; j < 10; j++) {
				for (std::size_t i = 0; i < 10; i++) {
					auto& node = *i_node++;
					double x = i+1.50;
					double y = j+2.50;
					double z = k+3.50;
                	KRATOS_CHECK_NEAR(node.X(), x, 1e-6);
                	KRATOS_CHECK_NEAR(node.Y(), y, 1e-6);
                	KRATOS_CHECK_NEAR(node.Z(), z, 1e-6);
				}
            }
		}
	}


	KRATOS_TEST_CASE_IN_SUITE(VoxelMeshGeneratorProcessTetrahedraElementColoring, KratosCoreFastSuite)
	{
		Parameters mesher_parameters(R"(
		{
			"number_of_divisions":   [10,10,10],
			"element_name":     "Element3D4N",
			"coloring_settings_list": [
				{
					"model_part_name": "SkinPart",
					"inside_color": -1,
					"outside_color": 1,
					"apply_outside_color": true,
					"coloring_entities": "center_of_elements"
				}
			]
		})");

        Model current_model;
		ModelPart &volume_part = current_model.CreateModelPart("Volume");
		volume_part.AddNodalSolutionStepVariable(DISTANCE);

		// Generate the skin
		ModelPart &skin_model_part = current_model.CreateModelPart("Skin");
        ModelPart &skin_part = skin_model_part.CreateSubModelPart("SkinPart");

		skin_part.AddNodalSolutionStepVariable(VELOCITY);
		skin_part.CreateNewNode(901, 2.0, 2.0, 2.0);
		skin_part.CreateNewNode(902, 6.0, 2.0, 2.0);
		skin_part.CreateNewNode(903, 4.0, 6.0, 2.0);
		skin_part.CreateNewNode(904, 4.0, 4.0, 7.0);
		Properties::Pointer p_properties(new Properties(0));
		skin_part.CreateNewElement("Element3D3N", 901, { 901,902,903 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 902, { 901,904,903 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 903, { 902,903,904 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 904, { 901,902,904 }, p_properties);

		// Generating the mesh
		VoxelMeshGeneratorProcess(Point{0.00, 0.00, 0.00}, Point{10.00, 10.00, 10.00}, volume_part, skin_model_part, mesher_parameters).Execute();
		// Compute distance
		// VoxelMeshColoringProcess(Point{0.00, 0.00, 0.00}, Point{10.00, 10.00, 10.00}, volume_part, skin_model_part, mesher_parameters).Execute();


		Tetrahedra3D4<Node<3>> tetrahedra(skin_part.pGetNode(901), skin_part.pGetNode(902), skin_part.pGetNode(903), skin_part.pGetNode(904));

        Point dummy(0.00,0.00,0.00);
        for(auto& element : volume_part.Elements()){
            if(tetrahedra.IsInside(element.GetGeometry().Center(),dummy)){
                KRATOS_CHECK_NEAR(element.GetValue(DISTANCE), -1.00, 1e-6);
            }
		}
	}


	KRATOS_TEST_CASE_IN_SUITE(VoxelMeshGeneratorProcessTetrahedraOnlyInsideElements, KratosCoreFastSuite)
	{
		Parameters mesher_parameters(R"(
		{
			"number_of_divisions":   [10,10,10],
			"element_name":     "Element3D4N",
			"output" : "inside_elements",
			"coloring_settings_list": [
				{
					"model_part_name": "SkinPart",
					"inside_color": -1,
					"outside_color": 1,
					"apply_outside_color": true,
					"coloring_entities": "center_of_elements"
				}
			]
		})");

        Model current_model;
		ModelPart &volume_part = current_model.CreateModelPart("Volume");
		volume_part.AddNodalSolutionStepVariable(DISTANCE);

		// Generate the skin
		ModelPart &skin_model_part = current_model.CreateModelPart("Skin");
        ModelPart &skin_part = skin_model_part.CreateSubModelPart("SkinPart");

		skin_part.AddNodalSolutionStepVariable(VELOCITY);
		skin_part.CreateNewNode(901, 2.0, 2.0, 2.0);
		skin_part.CreateNewNode(902, 6.0, 2.0, 2.0);
		skin_part.CreateNewNode(903, 4.0, 6.0, 2.0);
		skin_part.CreateNewNode(904, 4.0, 4.0, 7.0);
		Properties::Pointer p_properties(new Properties(0));
		skin_part.CreateNewElement("Element3D3N", 901, { 901,902,903 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 902, { 901,904,903 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 903, { 902,903,904 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 904, { 901,902,904 }, p_properties);

		// Generating the mesh
		VoxelMeshGeneratorProcess(Point{0.00, 0.00, 0.00}, Point{10.00, 10.00, 10.00}, volume_part, skin_model_part, mesher_parameters).Execute();
		// Compute distance
		// VoxelMeshColoringProcess(Point{0.00, 0.00, 0.00}, Point{10.00, 10.00, 10.00}, volume_part, skin_model_part, mesher_parameters).Execute();


		Tetrahedra3D4<Node<3>> tetrahedra(skin_part.pGetNode(901), skin_part.pGetNode(902), skin_part.pGetNode(903), skin_part.pGetNode(904));

        Point dummy(0.00,0.00,0.00);
        for(auto& element : volume_part.Elements()){
			KRATOS_CHECK(tetrahedra.IsInside(element.GetGeometry().Center(),dummy));
			KRATOS_CHECK_NEAR(element.GetValue(DISTANCE), -1.00, 1e-6);
		}
	}


	KRATOS_TEST_CASE_IN_SUITE(VoxelMeshGeneratorProcessCubeCenterOfElementsColoring, KratosCoreFastSuite)
	{
		Parameters mesher_parameters(R"(
		{
			"number_of_divisions":   [10,10,10],
			"element_name":     "Element3D4N",
			"entities_to_generate": "center_of_elements",
			"coloring_settings_list": [
				{
					"model_part_name": "SkinPart",
					"inside_color": -1,
					"outside_color": 1,
					"apply_outside_color": true,
					"coloring_entities": "center_of_elements"
				}
			]
		})");

        Model current_model;
		ModelPart &volume_part = current_model.CreateModelPart("Volume");
		volume_part.AddNodalSolutionStepVariable(DISTANCE);

		// Generate the skin
		ModelPart &skin_model_part = current_model.CreateModelPart("Skin");
        ModelPart &skin_part = skin_model_part.CreateSubModelPart("SkinPart");

		skin_part.AddNodalSolutionStepVariable(VELOCITY);
		skin_part.CreateNewNode(901, 2.0, 2.0, 2.0);
		skin_part.CreateNewNode(902, 5.9, 2.0, 2.0);
		skin_part.CreateNewNode(903, 5.9, 5.9, 2.0);
		skin_part.CreateNewNode(904, 2.0, 5.9, 2.0);
		skin_part.CreateNewNode(905, 2.0, 2.0, 7.0);
		skin_part.CreateNewNode(906, 5.9, 2.0, 7.0);
		skin_part.CreateNewNode(907, 5.9, 5.9, 7.0);
		skin_part.CreateNewNode(908, 2.0, 5.9, 7.0);
		Properties::Pointer p_properties(new Properties(0));
		skin_part.CreateNewElement("Element3D3N", 901, { 901,902,903 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 902, { 901,904,903 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 903, { 901,902,906 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 904, { 901,905,906 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 905, { 902,903,907 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 906, { 902,907,906 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 907, { 903,904,908 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 908, { 903,907,908 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 909, { 904,901,905 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 910, { 904,905,908 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 911, { 905,906,907 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 912, { 905,908,907 }, p_properties);

		// Generating the mesh
		VoxelMeshGeneratorProcess(Point{0.00, 0.00, 0.00}, Point{10.00, 10.00, 10.00}, volume_part, skin_model_part, mesher_parameters).Execute();

		auto i_node = volume_part.NodesBegin();
		for (std::size_t k = 0; k <= 4; k++) {
			for (std::size_t j = 0; j <= 4; j++) {
				for (std::size_t i = 0; i <= 4; i++) {
					auto& node = *i_node++;
					if((node.X() > 2.00) && (node.X() < 6.00) && (node.Y() > 2.00) && (node.Y() < 6.00) && (node.Z() > 2.00) && (node.Z() < 7.00)){
               			KRATOS_CHECK_NEAR(node.GetSolutionStepValue(DISTANCE), -1.00, 1e-6);
					}
					else{
               			KRATOS_CHECK_NEAR(node.GetSolutionStepValue(DISTANCE), 1.00, 1e-6);
					}
				}
            }
		}

	}

	KRATOS_TEST_CASE_IN_SUITE(VoxelMeshGeneratorProcessOutsideCubeCenterOfElementsColoring, KratosCoreFastSuite)
	{
		Parameters mesher_parameters(R"(
		{
			"number_of_divisions":   [10,10,10],
			"element_name":     "Element3D4N",
			"entities_to_generate": "center_of_elements",
			"coloring_settings_list": [
				{
					"model_part_name": "SkinPart",
					"inside_color": -1,
					"outside_color": 1,
					"apply_outside_color": true,
					"coloring_entities": "center_of_elements"
				}
			]
		})");

        Model current_model;
		ModelPart &volume_part = current_model.CreateModelPart("Volume");
		volume_part.AddNodalSolutionStepVariable(DISTANCE);

		// Generate the skin
		ModelPart &skin_model_part = current_model.CreateModelPart("Skin");
        ModelPart &skin_part = skin_model_part.CreateSubModelPart("SkinPart");

		skin_part.AddNodalSolutionStepVariable(VELOCITY);
		skin_part.CreateNewNode(901, 2.0, 2.0, 2.0);
		skin_part.CreateNewNode(902, 5.9, 2.0, 2.0);
		skin_part.CreateNewNode(903, 5.9, 5.9, 2.0);
		skin_part.CreateNewNode(904, 2.0, 5.9, 2.0);
		skin_part.CreateNewNode(905, 2.0, 2.0, 7.0);
		skin_part.CreateNewNode(906, 5.9, 2.0, 7.0);
		skin_part.CreateNewNode(907, 5.9, 5.9, 7.0);
		skin_part.CreateNewNode(908, 2.0, 5.9, 7.0);
		Properties::Pointer p_properties(new Properties(0));
		skin_part.CreateNewElement("Element3D3N", 901, { 901,902,903 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 902, { 901,904,903 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 903, { 901,902,906 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 904, { 901,905,906 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 905, { 902,903,907 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 906, { 902,907,906 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 907, { 903,904,908 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 908, { 903,907,908 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 909, { 904,901,905 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 910, { 904,905,908 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 911, { 905,906,907 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 912, { 905,908,907 }, p_properties);

		// Generating the mesh
		VoxelMeshGeneratorProcess(Point{2.50, 2.50, 2.50}, Point{5.00, 5.00, 5.00}, volume_part, skin_model_part, mesher_parameters).Execute();

		auto i_node = volume_part.NodesBegin();
		for (std::size_t k = 0; k < 10; k++) {
			for (std::size_t j = 0; j < 10; j++) {
				for (std::size_t i = 0; i < 10; i++) {
              		KRATOS_CHECK_NEAR(i_node->GetSolutionStepValue(DISTANCE), -1.00, 1e-6);
				}
            }
		}

	}


	KRATOS_TEST_CASE_IN_SUITE(VoxelMeshGeneratorProcessCubeCenterOfElementsRectilinearColors, KratosCoreFastSuite)
	{
		Parameters mesher_parameters(R"(
		{
			"number_of_divisions":   [10,10,10],
			"element_name":     "Element3D4N",
			"entities_to_generate": "center_of_elements",
			"output" : "rectilinear_coordinates",
			"output_filename" : "Cube.vtr",
			"coloring_settings_list": [
				{
					"model_part_name": "SkinPart",
					"inside_color": -1,
					"outside_color": 1,
					"apply_outside_color": true,
					"coloring_entities": "center_of_elements"
				}
			]
		})");

        Model current_model;
		ModelPart &volume_part = current_model.CreateModelPart("Volume");
		volume_part.AddNodalSolutionStepVariable(DISTANCE);

		// Generate the skin
		ModelPart &skin_model_part = current_model.CreateModelPart("Skin");
        ModelPart &skin_part = skin_model_part.CreateSubModelPart("SkinPart");

		skin_part.AddNodalSolutionStepVariable(VELOCITY);
		skin_part.CreateNewNode(901, 2.0, 2.0, 2.0);
		skin_part.CreateNewNode(902, 5.9, 2.0, 2.0);
		skin_part.CreateNewNode(903, 5.9, 5.9, 2.0);
		skin_part.CreateNewNode(904, 2.0, 5.9, 2.0);
		skin_part.CreateNewNode(905, 2.0, 2.0, 7.0);
		skin_part.CreateNewNode(906, 5.9, 2.0, 7.0);
		skin_part.CreateNewNode(907, 5.9, 5.9, 7.0);
		skin_part.CreateNewNode(908, 2.0, 5.9, 7.0);
		Properties::Pointer p_properties(new Properties(0));
		skin_part.CreateNewElement("Element3D3N", 901, { 901,902,903 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 902, { 901,904,903 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 903, { 901,902,906 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 904, { 901,905,906 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 905, { 902,903,907 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 906, { 902,907,906 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 907, { 903,904,908 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 908, { 903,907,908 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 909, { 904,901,905 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 910, { 904,905,908 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 911, { 905,906,907 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 912, { 905,908,907 }, p_properties);

		// Generating the mesh
		VoxelMeshGeneratorProcess(Point{0.00, 0.00, 0.00}, Point{10.00, 10.00, 10.00}, volume_part, skin_model_part, mesher_parameters).Execute();

		auto& x_coordinates = volume_part.GetValue(RECTILINEAR_X_COORDINATES);
		auto& y_coordinates = volume_part.GetValue(RECTILINEAR_Y_COORDINATES);
		auto& z_coordinates = volume_part.GetValue(RECTILINEAR_Z_COORDINATES);
		auto& colors = volume_part.GetValue(COLORS);

		KRATOS_CHECK_EQUAL(colors.size(), 1000);

		auto i_color = colors.begin();
		for (std::size_t k = 0; k < 10; k++) {
			for (std::size_t j = 0; j < 10; j++) {
				for (std::size_t i = 0; i < 10; i++) {
					if((x_coordinates[i] > 2.00) && (x_coordinates[i] < 6.00) && (y_coordinates[j]  > 2.00) && (y_coordinates[j] < 6.00) && (z_coordinates[k] > 2.00) && (z_coordinates[k] < 7.00)){
               			KRATOS_CHECK_NEAR(*i_color++, -1.00, 1e-6);
					}
					else{
               			KRATOS_CHECK_NEAR(*i_color++, 1.00, 1e-6);
					}
				}
            }
		}

	}

	KRATOS_TEST_CASE_IN_SUITE(VoxelMeshGeneratorProcessCubeRayCollapseBug, KratosCoreFastSuite)
	{
		Parameters mesher_parameters(R"(
		{
			"number_of_divisions":   [3,3,3],
			"element_name":     "Element3D4N",
			"entities_to_generate": "center_of_elements",
			"output" : "rectilinear_coordinates",
			"output_filename" : "CubeBug.vtr",
			"coloring_settings_list": [
				{
					"model_part_name": "SkinPart",
					"inside_color": -1,
					"outside_color": 1,
					"apply_outside_color": true,
					"coloring_entities": "center_of_elements"
				}
			]
		})");

        Model current_model;
		ModelPart &volume_part = current_model.CreateModelPart("Volume");
		volume_part.AddNodalSolutionStepVariable(DISTANCE);

		// Generate the skin
		ModelPart &skin_model_part = current_model.CreateModelPart("Skin");
        ModelPart &skin_part = skin_model_part.CreateSubModelPart("SkinPart");

		skin_part.AddNodalSolutionStepVariable(VELOCITY);
		skin_part.CreateNewNode(469, 0.018, 0.008, 0.01 );
		skin_part.CreateNewNode(470, 0.018, 0.003, 0.01 );
		skin_part.CreateNewNode(471, 0.018, 0.003, 0.009);
		skin_part.CreateNewNode(472, 0.018, 0.008, 0.009);
		skin_part.CreateNewNode(473, 0.015, 0.008, 0.01 );
		skin_part.CreateNewNode(474, 0.015, 0.008, 0.009);
		skin_part.CreateNewNode(475, 0.015, 0.003, 0.01 );
		skin_part.CreateNewNode(476, 0.015, 0.003, 0.009);
		Properties::Pointer p_properties(new Properties(0));
		skin_part.CreateNewElement("Element3D3N", 929, { 469,470,471 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 930, { 471,472,469 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 931, { 473,469,472 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 932, { 472,474,473 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 933, { 475,473,474 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 934, { 474,476,475 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 935, { 471,470,475 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 936, { 475,476,471 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 937, { 471,476,474 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 938, { 474,472,471 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 939, { 469,473,475 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 940, { 475,470,469 }, p_properties);



		// Generating the mesh
		VoxelMeshGeneratorProcess({0.0165, 0.017, 0.0175, 0.018}, { 0.0035, 0.004, 0.0045, 0.005}, {0.0095, 0.01}, volume_part, skin_model_part, mesher_parameters).Execute();

		auto& colors = volume_part.GetValue(COLORS);

		KRATOS_CHECK_EQUAL(colors.size(), 9);

		for(auto color : colors){
			KRATOS_CHECK_EQUAL(color, -1);
		}
	}

	KRATOS_TEST_CASE_IN_SUITE(VoxelMeshGeneratorProcessCubeFaceOfElementsRectilinearColors, KratosCoreFastSuite)
	{
		Parameters mesher_parameters(R"(
		{
			"number_of_divisions":   [10,10,10],
			"element_name":     "Element3D4N",
			"entities_to_generate": "center_of_elements",
			"output" : "rectilinear_coordinates",
				"output_filename" : "TestFace.vtr",
			"coloring_settings_list": [
				{
					"model_part_name": "SkinPart",
					"inside_color": -1,
					"outside_color": 1,
					"apply_outside_color": true,
					"coloring_entities": "center_of_elements"
				},
				{
					"model_part_name": "SkinPart",
					"inside_color": -1,
					"outside_color": 1,
					"interface_color":-2,
					"apply_outside_color": false,
					"coloring_entities": "face_of_elements"
				}
			]
		})");

        Model current_model;
		ModelPart &volume_part = current_model.CreateModelPart("Volume");
		volume_part.AddNodalSolutionStepVariable(DISTANCE);

		// Generate the skin
		ModelPart &skin_model_part = current_model.CreateModelPart("Skin");
        ModelPart &skin_part = skin_model_part.CreateSubModelPart("SkinPart");

		skin_part.AddNodalSolutionStepVariable(VELOCITY);
		skin_part.CreateNewNode(901, 2.0, 2.1, 1.9);
		skin_part.CreateNewNode(902, 5.9, 2.1, 1.9);
		skin_part.CreateNewNode(903, 5.9, 5.9, 1.9);
		skin_part.CreateNewNode(904, 2.0, 5.9, 1.9);
		skin_part.CreateNewNode(905, 2.0, 2.1, 7.3);
		skin_part.CreateNewNode(906, 5.9, 2.1, 7.3);
		skin_part.CreateNewNode(907, 5.9, 5.9, 7.3);
		skin_part.CreateNewNode(908, 2.0, 5.9, 7.3);
		Properties::Pointer p_properties(new Properties(0));
		skin_part.CreateNewElement("Element3D3N", 901, { 901,902,903 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 902, { 901,904,903 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 903, { 901,902,906 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 904, { 901,905,906 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 905, { 902,903,907 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 906, { 902,907,906 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 907, { 903,904,908 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 908, { 903,907,908 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 909, { 904,901,905 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 910, { 904,905,908 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 911, { 905,906,907 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 912, { 905,908,907 }, p_properties);

		// Generating the mesh
		VoxelMeshGeneratorProcess(Point{0.00, 0.00, 0.00}, Point{10.00, 10.00, 10.00}, volume_part, skin_model_part, mesher_parameters).Execute();

		auto& x_coordinates = volume_part.GetValue(RECTILINEAR_X_COORDINATES);
		auto& y_coordinates = volume_part.GetValue(RECTILINEAR_Y_COORDINATES);
		auto& z_coordinates = volume_part.GetValue(RECTILINEAR_Z_COORDINATES);
		auto& colors = volume_part.GetValue(VOXEL_FACE_COLORS);

		KRATOS_CHECK_EQUAL(colors.size1(), 1000);
		KRATOS_CHECK_EQUAL(colors.size2(), 6);

		std::size_t index = 0;
		for (std::size_t k = 0; k < 10; k++) {
			for (std::size_t j = 0; j < 10; j++) {
				for (std::size_t i = 0; i < 10; i++) {
					if((i >= 1) && (i < 9) && (y_coordinates[j] > 2.1) && (y_coordinates[j] < 5.9) && (z_coordinates[k] > 1.9) && (z_coordinates[k] < 7.3)){
						if((x_coordinates[i-1] < 2.00) && (x_coordinates[i] > 2.00)){
							KRATOS_CHECK_NEAR(colors(index, 0), -2.00, 1e-6);
							KRATOS_CHECK_NEAR(colors(index, 3), 1.00, 1e-6);
						}
						else if((x_coordinates[i] < 5.90) && (x_coordinates[i+1] > 5.90)){
							KRATOS_CHECK_NEAR(colors(index, 0), 1.00, 1e-6);
							KRATOS_CHECK_NEAR(colors(index, 3), -2.00, 1e-6);
						}
						else{
							KRATOS_CHECK_NEAR(colors(index, 0), 1.00, 1e-6);
							KRATOS_CHECK_NEAR(colors(index, 3), 1.00, 1e-6);
						}
					}
					
					index++;
				}
            }
		}
	}

	KRATOS_TEST_CASE_IN_SUITE(VoxelMeshGeneratorProcessCubeFaceOfElementsCoarseRectilinearColors, KratosCoreFastSuite)
	{
		Parameters mesher_parameters(R"(
		{
			"number_of_divisions":   [10,10,10],
			"element_name":     "Element3D4N",
			"entities_to_generate": "center_of_elements",
			"output" : "rectilinear_coordinates",
				"output_filename" : "TestFace.vtr",
			"coloring_settings_list": [
				{
					"model_part_name": "SkinPart",
					"inside_color": -1,
					"outside_color": 1,
					"apply_outside_color": true,
					"coloring_entities": "center_of_elements"
				},
				{
					"model_part_name": "SkinPart",
					"inside_color": -1,
					"outside_color": 1,
					"interface_color": -2,
					"apply_outside_color": false,
					"coloring_entities": "face_of_elements"
				}
			]
		})");

        Model current_model;
		ModelPart &volume_part = current_model.CreateModelPart("Volume");
		volume_part.AddNodalSolutionStepVariable(DISTANCE);

		// Generate the skin
		ModelPart &skin_model_part = current_model.CreateModelPart("Skin");
        ModelPart &skin_part = skin_model_part.CreateSubModelPart("SkinPart");

		skin_part.AddNodalSolutionStepVariable(VELOCITY);
		skin_part.CreateNewNode(901, 2.0, 2.1, 1.9);
		skin_part.CreateNewNode(902, 5.9, 2.1, 1.9);
		skin_part.CreateNewNode(903, 5.9, 5.9, 1.9);
		skin_part.CreateNewNode(904, 2.0, 5.9, 1.9);
		skin_part.CreateNewNode(905, 2.0, 2.1, 7.3);
		skin_part.CreateNewNode(906, 5.9, 2.1, 7.3);
		skin_part.CreateNewNode(907, 5.9, 5.9, 7.3);
		skin_part.CreateNewNode(908, 2.0, 5.9, 7.3);
		Properties::Pointer p_properties(new Properties(0));
		skin_part.CreateNewElement("Element3D3N", 901, { 901,902,903 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 902, { 901,904,903 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 903, { 901,902,906 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 904, { 901,905,906 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 905, { 902,903,907 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 906, { 902,907,906 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 907, { 903,904,908 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 908, { 903,907,908 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 909, { 904,901,905 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 910, { 904,905,908 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 911, { 905,906,907 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 912, { 905,908,907 }, p_properties);

		// Generating the mesh
		CoarseVoxelMeshGeneratorProcess(Point{0.00, 0.00, 0.00}, Point{10.00, 10.00, 10.00}, volume_part, skin_model_part, mesher_parameters).Execute();

		auto& colors = volume_part.GetValue(VOXEL_FACE_COLORS);

		KRATOS_CHECK_EQUAL(colors.size1(), 27);
		KRATOS_CHECK_EQUAL(colors.size2(), 6);

		std::size_t index = 0;
		for (std::size_t k = 0; k < 3; k++) {
			for (std::size_t j = 0; j < 3; j++) {
				for (std::size_t i = 0; i < 3; i++) {
					if((i == 1) && (j == 1) && (k == 1)){
						KRATOS_CHECK_NEAR(colors(index, 0), -2.00, 1e-6);
						KRATOS_CHECK_NEAR(colors(index, 1), -2.00, 1e-6);
						KRATOS_CHECK_NEAR(colors(index, 2), -2.00, 1e-6);
						KRATOS_CHECK_NEAR(colors(index, 3), -2.00, 1e-6);
						KRATOS_CHECK_NEAR(colors(index, 4), -2.00, 1e-6);
						KRATOS_CHECK_NEAR(colors(index, 5), -2.00, 1e-6);
					}
					else{
						KRATOS_CHECK_NEAR(colors(index, 0), 1.00, 1e-6);
						KRATOS_CHECK_NEAR(colors(index, 1), 1.00, 1e-6);
						KRATOS_CHECK_NEAR(colors(index, 2), 1.00, 1e-6);
						KRATOS_CHECK_NEAR(colors(index, 3), 1.00, 1e-6);
						KRATOS_CHECK_NEAR(colors(index, 4), 1.00, 1e-6);
						KRATOS_CHECK_NEAR(colors(index, 5), 1.00, 1e-6);
					}
					
					index++;
				}
            }
		}
	}

	KRATOS_TEST_CASE_IN_SUITE(VoxelMeshGeneratorProcessCubeCoarseMesh, KratosCoreFastSuite)
	{
		Parameters mesher_parameters(R"(
		{
			"number_of_divisions":   [10,10,10],
			"element_name":     "Element3D4N",
			"mesh_type": "coarse",
			"entities_to_generate": "center_of_elements",
			"coloring_settings_list": [
				{
					"model_part_name": "SkinPart",
					"inside_color": -1,
					"outside_color": 1,
					"apply_outside_color": true,
					"coloring_entities": "center_of_elements"
				}
			]
		})");

        Model current_model;
		ModelPart &volume_part = current_model.CreateModelPart("Volume");
		volume_part.AddNodalSolutionStepVariable(DISTANCE);

		// Generate the skin
		ModelPart &skin_model_part = current_model.CreateModelPart("Skin");
        ModelPart &skin_part = skin_model_part.CreateSubModelPart("SkinPart");

		skin_part.AddNodalSolutionStepVariable(VELOCITY);
		skin_part.CreateNewNode(901, 2.0, 2.0, 2.0);
		skin_part.CreateNewNode(902, 5.9, 2.0, 2.0);
		skin_part.CreateNewNode(903, 5.9, 5.9, 2.0);
		skin_part.CreateNewNode(904, 2.0, 5.9, 2.0);
		skin_part.CreateNewNode(905, 2.0, 2.0, 7.0);
		skin_part.CreateNewNode(906, 5.9, 2.0, 7.0);
		skin_part.CreateNewNode(907, 5.9, 5.9, 7.0);
		skin_part.CreateNewNode(908, 2.0, 5.9, 7.0);
		Properties::Pointer p_properties(new Properties(0));
		skin_part.CreateNewElement("Element3D3N", 901, { 901,902,903 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 902, { 901,904,903 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 903, { 901,902,906 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 904, { 901,905,906 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 905, { 902,903,907 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 906, { 902,907,906 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 907, { 903,904,908 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 908, { 903,907,908 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 909, { 904,901,905 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 910, { 904,905,908 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 911, { 905,906,907 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 912, { 905,908,907 }, p_properties);

		// Generating the mesh
		CoarseVoxelMeshGeneratorProcess(Point{0.00, 0.00, 0.00}, Point{10.00, 10.00, 10.00}, volume_part, skin_model_part, mesher_parameters).Execute();

		KRATOS_CHECK_EQUAL(volume_part.NumberOfNodes(), 4*4*4);

		std::array<double,4> xy{0.00, 2.00, 6.00, 10.00};
		std::array<double,4> z{0.00, 2.00, 7.00, 10.00};

		auto i_node = volume_part.NodesBegin();
		for (std::size_t k = 0; k < 4; k++) {
			for (std::size_t j = 0; j < 4; j++) {
				for (std::size_t i = 0; i < 4; i++) {
					auto& node = *i_node++;
                	KRATOS_CHECK_NEAR(node.X(), xy[i], 1e-6);
                	KRATOS_CHECK_NEAR(node.Y(), xy[j], 1e-6);
                	KRATOS_CHECK_NEAR(node.Z(), z[k], 1e-6);
				}
            }
		}

	}


	KRATOS_TEST_CASE_IN_SUITE(VoxelMeshGeneratorProcessCubeCoarseMeshRectilinearOutput, KratosCoreFastSuite)
	{
		Parameters mesher_parameters(R"(
		{
			"number_of_divisions":   [10,10,10],
			"element_name":     "Element3D4N",
			"mesh_type": "coarse",
			"entities_to_generate": "center_of_elements",
			"output" : "rectilinear_coordinates",
			"coloring_settings_list": [
				{
					"model_part_name": "SkinPart",
					"inside_color": -1,
					"outside_color": 1,
					"apply_outside_color": true,
					"coloring_entities": "center_of_elements"
				}
			]
		})");

        Model current_model;
		ModelPart &volume_part = current_model.CreateModelPart("Volume");
		volume_part.AddNodalSolutionStepVariable(DISTANCE);

		// Generate the skin
		ModelPart &skin_model_part = current_model.CreateModelPart("Skin");
        ModelPart &skin_part = skin_model_part.CreateSubModelPart("SkinPart");

		skin_part.AddNodalSolutionStepVariable(VELOCITY);
		skin_part.CreateNewNode(901, 2.0, 2.0, 2.0);
		skin_part.CreateNewNode(902, 5.9, 2.0, 2.0);
		skin_part.CreateNewNode(903, 5.9, 5.9, 2.0);
		skin_part.CreateNewNode(904, 2.0, 5.9, 2.0);
		skin_part.CreateNewNode(905, 2.0, 2.0, 7.0);
		skin_part.CreateNewNode(906, 5.9, 2.0, 7.0);
		skin_part.CreateNewNode(907, 5.9, 5.9, 7.0);
		skin_part.CreateNewNode(908, 2.0, 5.9, 7.0);
		Properties::Pointer p_properties(new Properties(0));
		skin_part.CreateNewElement("Element3D3N", 901, { 901,902,903 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 902, { 901,904,903 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 903, { 901,902,906 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 904, { 901,905,906 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 905, { 902,903,907 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 906, { 902,907,906 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 907, { 903,904,908 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 908, { 903,907,908 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 909, { 904,901,905 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 910, { 904,905,908 }, p_properties);

		skin_part.CreateNewElement("Element3D3N", 911, { 905,906,907 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 912, { 905,908,907 }, p_properties);

		// Generating the mesh
		CoarseVoxelMeshGeneratorProcess(Point{0.00, 0.00, 0.00}, Point{10.00, 10.00, 10.00}, volume_part, skin_model_part, mesher_parameters).Execute();

		std::array<double,4> xy{0.00, 2.00, 6.00, 10.00};
		std::array<double,4> z{0.00, 2.00, 7.00, 10.00};

		auto& x_coordinates = volume_part.GetValue(RECTILINEAR_X_COORDINATES);
		auto& y_coordinates = volume_part.GetValue(RECTILINEAR_Y_COORDINATES);
		auto& z_coordinates = volume_part.GetValue(RECTILINEAR_Z_COORDINATES);
		auto& colors = volume_part.GetValue(COLORS);

		KRATOS_CHECK_EQUAL(x_coordinates.size(), 4);
		KRATOS_CHECK_EQUAL(y_coordinates.size(), 4);
		KRATOS_CHECK_EQUAL(z_coordinates.size(), 4);
		KRATOS_CHECK_EQUAL(colors.size(), 27);

		for(int i = 0 ; i < 4 ; i++){
			KRATOS_CHECK_EQUAL(x_coordinates[i], xy[i]);
			KRATOS_CHECK_EQUAL(y_coordinates[i], xy[i]);
			KRATOS_CHECK_EQUAL(z_coordinates[i], z[i]);
		}

		auto i_color = colors.begin();
		for(int k = 0 ; k < 3 ; k++){
			for(int j = 0 ; j < 3 ; j++){
				for(int i = 0 ; i < 3 ; i++){
					if((i==1) && (j==1) && (k == 1)){
						KRATOS_CHECK_EQUAL(*i_color++, -1);
					}
					else {
						KRATOS_CHECK_EQUAL(*i_color++, 1);
					}
				}
			}
		}
	}


	KRATOS_TEST_CASE_IN_SUITE(VoxelMeshGeneratorProcessTetrahedraElementNodesColoring, KratosCoreFastSuite)
	{
		Parameters mesher_parameters(R"(
		{
			"number_of_divisions":   [10,10,10],
			"element_name":     "Element3D4N",
			"entities_to_generate": "center_of_elements",
			"coloring_settings_list": [
				{
					"model_part_name": "SkinPart",
					"inside_color": -1,
					"outside_color": 1,
					"apply_outside_color": true,
					"coloring_entities": "nodes"
				}
			]
		})");

        Model current_model;
		ModelPart &volume_part = current_model.CreateModelPart("Volume");
		volume_part.AddNodalSolutionStepVariable(VELOCITY);
		volume_part.AddNodalSolutionStepVariable(DISTANCE);
		volume_part.AddNodalSolutionStepVariable(EMBEDDED_VELOCITY);

		// Generate the skin
		ModelPart &skin_model_part = current_model.CreateModelPart("Skin");
        ModelPart &skin_part = skin_model_part.CreateSubModelPart("SkinPart");

		skin_part.AddNodalSolutionStepVariable(VELOCITY);
		skin_part.CreateNewNode(901, 2.0, 2.0, 2.0);
		skin_part.CreateNewNode(902, 6.0, 2.0, 2.0);
		skin_part.CreateNewNode(903, 4.0, 6.0, 2.0);
		skin_part.CreateNewNode(904, 4.0, 4.0, 7.0);
		Properties::Pointer p_properties(new Properties(0));
		skin_part.CreateNewElement("Element3D3N", 901, { 901,902,903 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 902, { 901,904,903 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 903, { 902,903,904 }, p_properties);
		skin_part.CreateNewElement("Element3D3N", 904, { 901,902,904 }, p_properties);

		// Generating the mesh
		VoxelMeshGeneratorProcess(Point{0.00, 0.00, 0.00}, Point{10.00, 10.00, 10.00}, volume_part, skin_model_part, mesher_parameters).Execute();
		// Compute distance
		// VoxelMeshColoringProcess(Point{0.00, 0.00, 0.00}, Point{10.00, 10.00, 10.00}, volume_part, skin_model_part, mesher_parameters).Execute();


		Tetrahedra3D4<Node<3>> tetrahedra(skin_part.pGetNode(901), skin_part.pGetNode(902), skin_part.pGetNode(903), skin_part.pGetNode(904));

        Point dummy(0.00,0.00,0.00);
        for(auto& node : volume_part.Nodes()){
            if(!tetrahedra.IsInside(node,dummy)){
				// KRATOS_WATCH(node.GetSolutionStepValue(DISTANCE));
                KRATOS_CHECK_NEAR(node.GetSolutionStepValue(DISTANCE), -1.00, 1e-6);
            }
		}
	}


}
}  // namespace Kratos.
