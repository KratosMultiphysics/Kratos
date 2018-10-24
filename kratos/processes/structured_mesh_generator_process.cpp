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
//
//

// System includes
#include <vector>
// External includes

// Project includes
#include "processes/structured_mesh_generator_process.h"
#include "geometries/geometry.h"
#include "geometries/point.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "includes/checks.h"


namespace Kratos
{
    StructuredMeshGeneratorProcess::StructuredMeshGeneratorProcess(const GeometryType& rGeometry, ModelPart& rOutputModelPart, Parameters& TheParameters)
		: Process()
		, mrGeometry(rGeometry)
        , mrOutputModelPart(rOutputModelPart) {

		Parameters default_parameters(R"(
            {
	            "create_skin_sub_model_part": true,
	            "start_node_id":1,
                "start_element_id":1,
                "start_condition_id":1,
                "number_of_divisions":1,
                "elements_properties_id":0,
                "conditions_properties_id":0,
                "element_name": "PLEASE SPECIFY IT",
                "condition_name": "PLEASE SPECIFY IT"
            }  )");

		TheParameters["element_name"]; // Should be given by caller! if not thorws an error

		TheParameters.ValidateAndAssignDefaults(default_parameters);

		mStartNodeId = TheParameters["start_node_id"].GetInt();
		mStartElementId = TheParameters["start_element_id"].GetInt();
		mStartConditionId = TheParameters["start_condition_id"].GetInt();

        mNumberOfDivisions = TheParameters["number_of_divisions"].GetInt();
		mElementPropertiesId = TheParameters["elements_properties_id"].GetInt();
		mConditiongPropertiesId = TheParameters["conditions_properties_id"].GetInt();
		mElementName = TheParameters["element_name"].GetString();
		mConditionName = TheParameters["condition_name"].GetString();
        mCreateSkinSubModelPart = TheParameters["create_skin_sub_model_part"].GetBool();

        PerformChecks();

    }

	StructuredMeshGeneratorProcess::~StructuredMeshGeneratorProcess() {

	}

	void StructuredMeshGeneratorProcess::Execute() {
        if (mCreateSkinSubModelPart)
            mrOutputModelPart.CreateSubModelPart("Skin");
		if (mrGeometry.LocalSpaceDimension() == 2)
			Generate2DMesh();
		else if (mrGeometry.LocalSpaceDimension() == 3)
            Generate3DMesh();
		else
			KRATOS_ERROR << "Not supported geometry is given" << std::endl;

	}

	std::string StructuredMeshGeneratorProcess::Info() const {
		return "StructuredMeshGeneratorProcess";
	}

	void StructuredMeshGeneratorProcess::PrintInfo(std::ostream& rOStream) const {
		rOStream << Info();
	}

	void StructuredMeshGeneratorProcess::PrintData(std::ostream& rOStream) const {

	}

	void StructuredMeshGeneratorProcess::Generate2DMesh() {
		Point min_point(1.00, 1.00, 1.00);
		Point max_point(-1.00, -1.00, -1.00);
		GetLocalCoordinatesRange(min_point, max_point);

		GenerateNodes2D(min_point, max_point);

		GenerateTriangularElements();
	}

	void StructuredMeshGeneratorProcess::Generate3DMesh() {
		Point min_point(1.00, 1.00, 1.00);
		Point max_point(-1.00, -1.00, -1.00);
		GetLocalCoordinatesRange(min_point, max_point);

		GenerateNodes3D(min_point, max_point);

		GenerateTetrahedraElements();
	}

	void StructuredMeshGeneratorProcess::GenerateNodes2D(Point const& rMinPoint, Point const& rMaxPoint) {
		GeometryType::CoordinatesArrayType local_element_size = rMaxPoint - rMinPoint;
		local_element_size /= mNumberOfDivisions;
		//const std::size_t local_space_dimension = mrGeometry.LocalSpaceDimension();
		Point local_coordinates = rMinPoint;
		Point global_coordinates = ZeroVector(3);
		std::size_t node_id = mStartNodeId;

		for (std::size_t j = 0; j <= mNumberOfDivisions; j++) {
			for (std::size_t i = 0; i <= mNumberOfDivisions; i++) {
				local_coordinates[0] = rMinPoint[0] + (i * local_element_size[0]);
				local_coordinates[1] = rMinPoint[1] + (j * local_element_size[1]);
				mrGeometry.GlobalCoordinates(global_coordinates, local_coordinates);
				if (mCreateSkinSubModelPart && (
                    (i == 0) || (i == mNumberOfDivisions) || (j == 0) || (j == mNumberOfDivisions)))  // Is on skin
					mrOutputModelPart.GetSubModelPart("Skin").CreateNewNode(node_id++, global_coordinates[0],
                                                                                       global_coordinates[1],
                                                                                       global_coordinates[2]);
				else
					mrOutputModelPart.CreateNewNode(node_id++, global_coordinates[0],
                                                               global_coordinates[1],
                                                               global_coordinates[2]);
			}
		}
	}

	void StructuredMeshGeneratorProcess::GenerateNodes3D(Point const& rMinPoint, Point const& rMaxPoint) {
		GeometryType::CoordinatesArrayType local_element_size = rMaxPoint - rMinPoint;
		local_element_size /= mNumberOfDivisions;
		Point local_coordinates = rMinPoint;
		Point global_coordinates = ZeroVector(3);
		std::size_t node_id = mStartNodeId;

		for (std::size_t k = 0; k <= mNumberOfDivisions; k++) {
			for (std::size_t j = 0; j <= mNumberOfDivisions; j++) {
				for (std::size_t i = 0; i <= mNumberOfDivisions; i++) {
					local_coordinates[0] = rMinPoint[0] + (i * local_element_size[0]);
					local_coordinates[1] = rMinPoint[1] + (j * local_element_size[1]);
					local_coordinates[2] = rMinPoint[2] + (k * local_element_size[2]);
					mrGeometry.GlobalCoordinates(global_coordinates, local_coordinates);
					if (mCreateSkinSubModelPart && (
                        (i == 0) || (i == mNumberOfDivisions) || (j == 0) ||
                        (j == mNumberOfDivisions) || (k == 0) || (k == mNumberOfDivisions)))  // Is on skin
						mrOutputModelPart.GetSubModelPart("Skin").CreateNewNode(node_id++, global_coordinates[0],
                                                                                           global_coordinates[1],
                                                                                           global_coordinates[2]);
					else
						mrOutputModelPart.CreateNewNode(node_id++, global_coordinates[0],
                                                                   global_coordinates[1],
                                                                   global_coordinates[2]);
				}
			}
		}
	}

	void StructuredMeshGeneratorProcess::GenerateTriangularElements() {
		std::size_t element_id = mStartElementId;

		Properties::Pointer p_properties = mrOutputModelPart.pGetProperties(mElementPropertiesId);
		std::vector<ModelPart::IndexType> element_connectivity(3);

		for (std::size_t j = 0; j < mNumberOfDivisions; j++) {
			for (std::size_t i = 0; i < mNumberOfDivisions; i++) {
				element_connectivity = { GetNodeId(i,j,0), GetNodeId(i + 1,j + 1,0), GetNodeId(i + 1,j,0) };
				mrOutputModelPart.CreateNewElement(mElementName, element_id++, element_connectivity, p_properties);

				element_connectivity = { GetNodeId(i,j,0), GetNodeId(i,j + 1,0), GetNodeId(i + 1,j + 1,0) };
				mrOutputModelPart.CreateNewElement(mElementName, element_id++, element_connectivity, p_properties);
			}
		}
	}

	void StructuredMeshGeneratorProcess::GenerateTetrahedraElements() {
		Properties::Pointer p_properties = mrOutputModelPart.pGetProperties(mElementPropertiesId);

		for (std::size_t k = 0; k < mNumberOfDivisions; k++) {
			for (std::size_t j = 0; j < mNumberOfDivisions; j++) {
				for (std::size_t i = 0; i < mNumberOfDivisions; i++) {
					CreateCellTetrahedra(i, j, k, p_properties);
				}
			}
		}
	}

	void  StructuredMeshGeneratorProcess::CreateCellTetrahedra(std::size_t I, std::size_t J, std::size_t K, Properties::Pointer pProperties) {
		using point_in_cell_position_type = std::array<std::size_t, 3>;
		using tetrahedra_connectivity_in_cell_type = std::array<std::size_t, 4>;
		constexpr std::size_t number_of_cases = 6;
		constexpr point_in_cell_position_type cell_points[8] = { {{ 0,0,0 }},{{ 1,0,0 }},{{ 1,1,0 }},{{ 0,1,0 }},
																 {{ 0,0,1 }},{{ 1,0,1 }},{{ 1,1,1 }},{{ 0,1,1 }}  };

		constexpr tetrahedra_connectivity_in_cell_type connectivity_cases[number_of_cases] = { {{ 0,3,6,2 }},{{ 3,6,7,0 }},{{ 4,7,6,0 }},
																							   {{ 0,4,5,6 }},{{ 0,1,2,6 }},{{ 1,5,6,0 }} };
		std::vector<ModelPart::IndexType> element_connectivity(4);

		for (std::size_t i_case = 0; i_case < number_of_cases; i_case++) {
			auto connectivity = connectivity_cases[i_case];
			for (std::size_t i_position = 0; i_position < 4; i_position++)
			{
				auto& cell_point = cell_points[connectivity[i_position]];
				element_connectivity[i_position] = GetNodeId(I + cell_point[0], J + cell_point[1], K + cell_point[2]);
			}
			mrOutputModelPart.CreateNewElement(mElementName, mStartElementId++, element_connectivity, pProperties);
		}
	}

	std::size_t StructuredMeshGeneratorProcess::GetNodeId(std::size_t I, std::size_t J, std::size_t K) {
		return mStartNodeId + (K * (mNumberOfDivisions + 1) * (mNumberOfDivisions + 1)) + (J * (mNumberOfDivisions + 1)) + I;
	}

    void StructuredMeshGeneratorProcess::GetLocalCoordinatesRange(Point& rMinPoint, Point& rMaxPoint) {
		const std::size_t local_space_dimension = mrGeometry.LocalSpaceDimension();
		Matrix geometry_points_local_coordinates;
		mrGeometry.PointsLocalCoordinates(geometry_points_local_coordinates);

		const std::size_t number_of_points = mrGeometry.size();

		for (std::size_t i_point = 0; i_point < number_of_points; i_point++) {
			for (std::size_t i_dimension = 0; i_dimension < local_space_dimension; i_dimension++) {
				rMinPoint[i_dimension] = std::min(rMinPoint[i_dimension], geometry_points_local_coordinates(i_point, i_dimension));
				rMaxPoint[i_dimension] = std::max(rMaxPoint[i_dimension], geometry_points_local_coordinates(i_point, i_dimension));
			}
		}
	}

    int StructuredMeshGeneratorProcess::Check()
    {
        KRATOS_TRY

        KRATOS_CHECK(CheckDomainGeometry());
        KRATOS_CHECK(KratosComponents<Element>::Has(mElementName));

        if ((mrGeometry.GetGeometryType() != GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4) &
            (mrGeometry.GetGeometryType() != GeometryData::KratosGeometryType::Kratos_Hexahedra3D8))
            KRATOS_ERROR << "An unsupported geometry was given. Only Quadrilateral2D4 and Hexahedra3D8 are supported and given geometry is : " << mrGeometry << std::endl;

        KRATOS_CHECK_NOT_EQUAL(mNumberOfDivisions, 0);

        return 0;

        KRATOS_CATCH("")
    }

    bool StructuredMeshGeneratorProcess::CheckDomainGeometry() {
        if (mrGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4) {
            return CheckDomainGeometryConnectivityForQuadrilateral2D4();
        }

        else if (mrGeometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8) {
            return CheckDomainGeometryConnectivityForHexahedra3D8();
        }

        return true;
    }

    bool StructuredMeshGeneratorProcess::CheckDomainGeometryConnectivityForQuadrilateral2D4() {
        using triangle_connectivity_in_cell_type = std::array<std::size_t, 3>;
        constexpr std::size_t number_of_cases = 2;

        std::vector<std::array<double, 3> > cell_points;
        for (std::size_t i = 0; i < 4; ++i){
            std::array<double, 3> coordinates{{mrGeometry[i][0], mrGeometry[i][1], mrGeometry[i][2]}};
            cell_points.push_back(coordinates);
        }

        constexpr triangle_connectivity_in_cell_type connectivity_cases[number_of_cases] = { {{ 0,2,1 }},{{ 0,3,2 }} };

        std::vector<Point::Pointer> my_points(3);
        double min_area = 1.0;

        for (std::size_t i_case = 0; i_case < number_of_cases; i_case++) {
            auto connectivity = connectivity_cases[i_case];
            for (std::size_t i_position = 0; i_position < 3; i_position++)
            {
                auto& cell_point = cell_points[connectivity[i_position]];
                Point::Pointer pPi(new Point(cell_point[0], cell_point[1], cell_point[2]));
                my_points[i_position] = pPi;
            }

            Triangle2D3<Point > trial_triangle(my_points[0], my_points[1], my_points[2]);
            min_area = std::min(min_area, trial_triangle.DomainSize());
        }

        bool all_triangles_have_positive_area = min_area > 0.0;
        return all_triangles_have_positive_area;
    }


    bool StructuredMeshGeneratorProcess::CheckDomainGeometryConnectivityForHexahedra3D8() {
        using tetrahedra_connectivity_in_cell_type = std::array<std::size_t, 4>;
        constexpr std::size_t number_of_cases = 6;

        std::vector<std::array<double, 3> > cell_points;
        for (std::size_t i = 0; i < 8; ++i){
            std::array<double, 3> coordinates{{mrGeometry[i][0], mrGeometry[i][1], mrGeometry[i][2]}};
            cell_points.push_back(coordinates);
        }

        constexpr tetrahedra_connectivity_in_cell_type connectivity_cases[number_of_cases] = { {{ 0,3,6,2 }},{{ 3,6,7,0 }},{{ 4,7,6,0 }},
                                                                                               {{ 0,4,5,6 }},{{ 0,1,2,6 }},{{ 1,5,6,0 }} };
        std::vector<Point::Pointer> my_points(4);
        double min_volume = 1.0;

        for (std::size_t i_case = 0; i_case < number_of_cases; i_case++) {
            auto connectivity = connectivity_cases[i_case];
            for (std::size_t i_position = 0; i_position < 4; i_position++)
            {
                auto& cell_point = cell_points[connectivity[i_position]];
                Point::Pointer pPi(new Point(cell_point[0], cell_point[1], cell_point[2]));
                my_points[i_position] = pPi;
            }

            Tetrahedra3D4<Point > trial_tetra(my_points[0], my_points[1], my_points[2], my_points[3]);
            min_volume = std::min(min_volume, trial_tetra.DomainSize());
        }

        bool all_tetrahedra_have_positive_volume = min_volume > 0.0;
        return all_tetrahedra_have_positive_volume;
    }



}  // namespace Kratos.
