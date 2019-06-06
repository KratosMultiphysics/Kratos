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
#include "processes/voxel_mesh_generator_process.h"
#include "geometries/geometry.h"
#include "geometries/point.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/hexahedra_3d_8.h"
#include "includes/checks.h"
#include "processes/voxel_mesh_coloring_process.h"



namespace Kratos
{
    VoxelMeshGeneratorProcess::VoxelMeshGeneratorProcess(Point const& MinPoint, Point const& MaxPoint,
        ModelPart& rVolumePart,
        ModelPart& rSkinPart, Parameters& TheParameters)
		: Process()
        , mGeometry(Point(MinPoint[0], MinPoint[1], MinPoint[2]),
                    Point( MaxPoint[0], MinPoint[1], MinPoint[2]),
                    Point( MaxPoint[0],  MaxPoint[1], MinPoint[2]),
                    Point(MinPoint[0],  MaxPoint[1], MinPoint[2]),
                    Point(MinPoint[0], MinPoint[1],  MaxPoint[2]),
                    Point( MaxPoint[0], MinPoint[1],  MaxPoint[2]),
                    Point( MaxPoint[0],  MaxPoint[1],  MaxPoint[2]),
                    Point(MinPoint[0],  MaxPoint[1],  MaxPoint[2]))
		, mMinPoint(MinPoint)
        , mMaxPoint(MaxPoint)
        , mrVolumePart(rVolumePart), mrSkinPart(rSkinPart), mFindIntersectedObjectsProcess(rVolumePart, rSkinPart) {

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
                "condition_name": "PLEASE SPECIFY IT",
				"coloring_settings_list": [],
				"entities_to_generate": "elements"
            }  )");

		TheParameters["element_name"]; // Should be given by caller! if not thorws an error

		TheParameters.ValidateAndAssignDefaults(default_parameters);

		mStartNodeId = TheParameters["start_node_id"].GetInt();
		mStartElementId = TheParameters["start_element_id"].GetInt();
		mStartConditionId = TheParameters["start_condition_id"].GetInt();

        mNumberOfDivisions[0] = TheParameters["number_of_divisions"].GetInt();
        mNumberOfDivisions[1] = TheParameters["number_of_divisions"].GetInt();
        mNumberOfDivisions[2] = TheParameters["number_of_divisions"].GetInt();
		mElementPropertiesId = TheParameters["elements_properties_id"].GetInt();
		mConditiongPropertiesId = TheParameters["conditions_properties_id"].GetInt();
		mElementName = TheParameters["element_name"].GetString();
		mConditionName = TheParameters["condition_name"].GetString();
        mCreateSkinSubModelPart = TheParameters["create_skin_sub_model_part"].GetBool();
        mCellSizes = mMaxPoint - mMinPoint;
		mColoringParameters = TheParameters["coloring_settings_list"];
        for(int i = 0 ; i < 3 ; i++)
            mCellSizes[i] /= mNumberOfDivisions[i];
		mEntitiesToGenerate=TheParameters["entities_to_generate"].GetString();

        Check();
    }

	VoxelMeshGeneratorProcess::~VoxelMeshGeneratorProcess() {

	}

	void VoxelMeshGeneratorProcess::Execute() {

		if(mEntitiesToGenerate == "center_of_elements")
			GenerateCenterOfElements(mMinPoint, mMaxPoint);

		if(mEntitiesToGenerate == "nodes" || mEntitiesToGenerate == "elements"){
			if (mCreateSkinSubModelPart)
				if(!mrVolumePart.HasSubModelPart("Skin"))
					mrVolumePart.CreateSubModelPart("Skin");

			GenerateNodes3D(mMinPoint, mMaxPoint);
		}

		if(mEntitiesToGenerate == "elements")
        	Generate3DMesh();

		auto number_of_divisions = mNumberOfDivisions;
		auto min_point = mMinPoint;
		auto max_point = mMaxPoint;


		// if(mEntitiesToGenerate == "center_of_elements")
		// 	for(std::size_t i = 0 ; i < 3 ; i++){
		// 		number_of_divisions[i] -= 1;
		// 		min_point[i] += 0.5 * mCellSizes[i];
		// 		max_point[i] -= 0.5 * mCellSizes[i];
		// 	}
			
		for(auto item : mColoringParameters){
			if(mEntitiesToGenerate != "elements")
				KRATOS_ERROR_IF(item["coloring_entities"].GetString() == "elements") << "The coloring entities is set to element but there are no elements generated."
																					 << " Please set the entities_to_generate to 'elements' or set coloring entities to 'nodes'" << std::endl;
			std::string model_part_name = item["model_part_name"].GetString();
			if(model_part_name == mrSkinPart.Name())
				VoxelMeshColoringProcess(min_point, max_point, number_of_divisions, mrVolumePart, mrSkinPart, item).Execute();
			else {
				ModelPart& skin_part = mrSkinPart.GetSubModelPart(model_part_name);
				VoxelMeshColoringProcess(min_point, max_point, number_of_divisions, mrVolumePart, skin_part, item).Execute();
			}
		}

	}

	std::string VoxelMeshGeneratorProcess::Info() const {
		return "VoxelMeshGeneratorProcess";
	}

	void VoxelMeshGeneratorProcess::PrintInfo(std::ostream& rOStream) const {
		rOStream << Info();
	}

	void VoxelMeshGeneratorProcess::PrintData(std::ostream& rOStream) const {

	}

	void VoxelMeshGeneratorProcess::Generate3DMesh() {
        if(!mrVolumePart.HasProperties(mElementPropertiesId))
            mrVolumePart.CreateNewProperties(mElementPropertiesId);

        Properties::Pointer p_properties = mrVolumePart.pGetProperties(mElementPropertiesId);

        std::size_t cell_index = 0;
		Geometry<Node<3>>::PointsArrayType points(8);

		for (std::size_t K = 0; K < mNumberOfDivisions[2]; K++) {
			for (std::size_t J = 0; J < mNumberOfDivisions[1]; J++) {
				for (std::size_t i = 0; i < mNumberOfDivisions[0]; i++) {
					//Hexahedra3D8<Node<3>> hexahedra(pGetNode(i, J, K), pGetNode(i, J+1, K), pGetNode(i+1, J+1, K), pGetNode(i+1, J, K), pGetNode(i, J, K+1), pGetNode(i, J+1, K+1), pGetNode(i+1, J+1, K+1), pGetNode(i+1, J, K+1));
					points(0) = pGetNode(i, J, K);
					points(1) = pGetNode(i, J+1, K);
					points(2) = pGetNode(i+1, J+1, K);
					points(3) = pGetNode(i+1, J, K);
					points(4) = pGetNode(i, J, K+1);
					points(5) = pGetNode(i, J+1, K+1);
					points(6) = pGetNode(i+1, J+1, K+1);
					points(7) = pGetNode(i+1, J, K+1);
					mrVolumePart.CreateNewElement("Element3D8N", mStartElementId + cell_index, points, p_properties);
					cell_index++;
				}
			}
		}
	}

	void VoxelMeshGeneratorProcess::GenerateNodes3D(Point const& rMinPoint, Point const& rMaxPoint) {
		GeometryType::CoordinatesArrayType local_element_size = mCellSizes;
		Point local_coordinates = rMinPoint;
		auto global_coordinates = Point{ZeroVector(3)};
		std::size_t node_id = mStartNodeId;
		ModelPart* p_skin_part = nullptr;
		if (mCreateSkinSubModelPart)
			p_skin_part = &(mrVolumePart.GetSubModelPart("Skin"));

		for (std::size_t k = 0; k <= mNumberOfDivisions[2]; k++) {
			for (std::size_t j = 0; j <= mNumberOfDivisions[1]; j++) {
				for (std::size_t i = 0; i <= mNumberOfDivisions[0]; i++) {
					global_coordinates[0] = rMinPoint[0] + (i * local_element_size[0]);
					global_coordinates[1] = rMinPoint[1] + (j * local_element_size[1]);
					global_coordinates[2] = rMinPoint[2] + (k * local_element_size[2]);
					if (mCreateSkinSubModelPart && (
                        (i == 0) || (i == mNumberOfDivisions[0]) || (j == 0) ||
                        (j == mNumberOfDivisions[1]) || (k == 0) || (k == mNumberOfDivisions[2])))  // Is on skin
						p_skin_part->CreateNewNode(node_id++, global_coordinates[0],
                                                                                           global_coordinates[1],
                                                                                           global_coordinates[2]);
					else
						mrVolumePart.CreateNewNode(node_id++, global_coordinates[0],
                                                                   global_coordinates[1],
                                                                   global_coordinates[2]);
				}
			}
		}
	}

	void VoxelMeshGeneratorProcess::GenerateCenterOfElements(Point const& rMinPoint, Point const& rMaxPoint) {
		GeometryType::CoordinatesArrayType local_element_size = mCellSizes;
		Point local_coordinates = rMinPoint;
		auto global_coordinates = Point{ZeroVector(3)};
		std::size_t node_id = mStartNodeId;

		for (std::size_t k = 0; k < mNumberOfDivisions[2]; k++) {
			for (std::size_t j = 0; j < mNumberOfDivisions[1]; j++) {
				for (std::size_t i = 0; i < mNumberOfDivisions[0]; i++) {
					global_coordinates[0] = rMinPoint[0] + ((i + 0.5) * local_element_size[0]);
					global_coordinates[1] = rMinPoint[1] + ((j + 0.5) * local_element_size[1]);
					global_coordinates[2] = rMinPoint[2] + ((k + 0.5) * local_element_size[2]);
					mrVolumePart.CreateNewNode(node_id++, global_coordinates[0],
                                                          global_coordinates[1],
                                                          global_coordinates[2]);
				}
			}
		}
	}

	Node<3>::Pointer VoxelMeshGeneratorProcess::pGetNode(std::size_t I, std::size_t J, std::size_t K) {
		return *(mrVolumePart.NodesBegin() + (K * (mNumberOfDivisions[1] + 1) * (mNumberOfDivisions[0] + 1)) + (J * (mNumberOfDivisions[0] + 1)) + I).base();
	}

    int VoxelMeshGeneratorProcess::Check()
    {
        KRATOS_TRY

		KRATOS_ERROR_IF(mEntitiesToGenerate != "nodes" && mEntitiesToGenerate != "elements" && mEntitiesToGenerate != "center_of_elements") << mEntitiesToGenerate 
			<< " is not accepted as entities_to_generate. The valid options are 'nodes', 'elements' and 'center_of_elements'." << std::endl;

        return 0;

        KRATOS_CATCH("")
    }

}  // namespace Kratos.
