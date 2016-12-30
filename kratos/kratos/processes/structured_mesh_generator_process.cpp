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
#include "includes/checks.h"


namespace Kratos
{
	StructuredMeshGeneratorProcess::StructuredMeshGeneratorProcess(GeometryType& rGeometry, ModelPart& rOutputModelPart, Parameters TheParameters)
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
		mCrateSkinSubModelPart = TheParameters["create_skin_sub_model_part"].GetBool();

		KRATOS_CHECK(KratosComponents<Element>::Has(mElementName));
		Element const& r_element = KratosComponents<Element>::Get(mElementName);

		// I cannot do this test because the GetGeometryType is not a constant method! Pooyan.
		//if (r_element.GetGeometry().GetGeometryType() != GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4)
		//	KRATOS_ERROR << "Un supported geometry is given. Only Quadrilateral2D4 is supported and given geometry is : " << rGeometry << std::endl;

		if (rGeometry.GetGeometryType() != GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4)
			KRATOS_ERROR << "Un supported geometry is given. Only Quadrilateral2D4 is supported and given geometry is : " << rGeometry << std::endl;

		KRATOS_CHECK_NOT_EQUAL(mNumberOfDivisions, 0);
	}

	StructuredMeshGeneratorProcess::~StructuredMeshGeneratorProcess() {

	}

	void StructuredMeshGeneratorProcess::Execute() {
		if (mCrateSkinSubModelPart)
			mrOutputModelPart.CreateSubModelPart("Skin");
		if (mrGeometry.LocalSpaceDimension() == 2)
			Generate2DMesh();
		//else if (mrGeometry.LocalSpaceDimension() == 3)
		//	Generate3DMesh();
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
		Point<3> min_point(1.00, 1.00, 1.00);
		Point<3> max_point(-1.00, -1.00, -1.00);
		GetLocalCoordinatesRange(min_point, max_point);

		GenerateNodes2D(min_point, max_point);

		GenerateTriangularElements();
	}

	void StructuredMeshGeneratorProcess::GenerateNodes2D(Point<3> const& rMinPoint, Point<3> const& rMaxPoint) {
		GeometryType::CoordinatesArrayType local_element_size = rMaxPoint - rMinPoint;
		local_element_size /= mNumberOfDivisions;
		//const std::size_t local_space_dimension = mrGeometry.LocalSpaceDimension();
		Point<3> local_coordinates = rMinPoint;
		Point<3> global_coordinates = ZeroVector(3);
		std::size_t node_id = mStartNodeId;

		for (std::size_t i = 0; i <= mNumberOfDivisions; i++) {
			for (std::size_t j = 0; j <= mNumberOfDivisions; j++) {
				local_coordinates[0] = rMinPoint[0] + (i * local_element_size[0]);
				local_coordinates[1] = rMinPoint[1] + (j * local_element_size[1]);
				mrGeometry.GlobalCoordinates(global_coordinates, local_coordinates);
				if((i == 0) || (i == mNumberOfDivisions) || (j == 0) || (j == mNumberOfDivisions))  // Is on skin
					mrOutputModelPart.GetSubModelPart("Skin").CreateNewNode(node_id++, global_coordinates[0], global_coordinates[1], global_coordinates[2]);
				else
					mrOutputModelPart.CreateNewNode(node_id++, global_coordinates[0], global_coordinates[1], global_coordinates[2]);
			}
		}
	}

	void StructuredMeshGeneratorProcess::GenerateTriangularElements() {
		//std::size_t number_of_nodes = (mNumberOfDivisions + 1) * (mNumberOfDivisions + 1);
		//std::size_t first_node_id = mStartNodeId;
		std::size_t element_id = mStartElementId;

		Properties::Pointer p_properties = mrOutputModelPart.pGetProperties(mElementPropertiesId);
		std::vector<ModelPart::IndexType> element_connectivity(3);

		for (std::size_t i = 0; i < mNumberOfDivisions; i++) {
			for (std::size_t j = 0; j < mNumberOfDivisions; j++) {
				element_connectivity = { GetNodeId(i,j), GetNodeId(i + 1,j + 1), GetNodeId(i + 1,j)};
				mrOutputModelPart.CreateNewElement(mElementName, element_id++, element_connectivity, p_properties);
				element_connectivity = { GetNodeId(i,j), GetNodeId(i,j + 1), GetNodeId(i + 1,j + 1) };
				mrOutputModelPart.CreateNewElement(mElementName, element_id++, element_connectivity, p_properties);
			}
		}
	}

	std::size_t StructuredMeshGeneratorProcess::GetNodeId(std::size_t Row, std::size_t Column) {
		return mStartNodeId + (Row * (mNumberOfDivisions + 1)) + Column;
	}

    void StructuredMeshGeneratorProcess::GetLocalCoordinatesRange(Point<3>& rMinPoint, Point<3>& rMaxPoint) {
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



}  // namespace Kratos.
