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
#include "geometries/tetrahedra_3d_4.h"
#include "includes/checks.h"


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
                "condition_name": "PLEASE SPECIFY IT"
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
        for(int i = 0 ; i < 3 ; i++)
            mCellSizes[i] /= mNumberOfDivisions[i];

        Check();
    }

	VoxelMeshGeneratorProcess::~VoxelMeshGeneratorProcess() {

	}

	void VoxelMeshGeneratorProcess::Execute() {

        /// Fill container with objects

        Generate3DMesh();
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
		if (mCreateSkinSubModelPart)
			if(!mrVolumePart.HasSubModelPart("Skin"))
				mrVolumePart.CreateSubModelPart("Skin");

		GenerateNodes3D(mMinPoint, mMaxPoint);

        if(!mrVolumePart.HasProperties(mElementPropertiesId))
            mrVolumePart.CreateNewProperties(mElementPropertiesId);

        Properties::Pointer p_properties = mrVolumePart.pGetProperties(mElementPropertiesId);

        std::size_t cell_index = 0;

		for (std::size_t K = 0; K < mNumberOfDivisions[2]; K++) {
			for (std::size_t J = 0; J < mNumberOfDivisions[1]; J++) {
				for (std::size_t i = 0; i < mNumberOfDivisions[0]; i++) {
					std::vector<ModelPart::IndexType> element_connectivity(8);
					element_connectivity[0] = GetNodeId(i, J, K);
					element_connectivity[1] = GetNodeId(i, J+1, K);
					element_connectivity[2] = GetNodeId(i+1, J+1, K);
					element_connectivity[3] = GetNodeId(i+1, J, K);
					element_connectivity[4] = GetNodeId(i, J, K+1);
					element_connectivity[5] = GetNodeId(i, J+1, K+1);
					element_connectivity[6] = GetNodeId(i+1, J+1, K+1);
					element_connectivity[7] = GetNodeId(i+1, J, K+1);
					mrVolumePart.CreateNewElement("Element3D8N", mStartElementId + cell_index, element_connectivity, p_properties);
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

	std::size_t VoxelMeshGeneratorProcess::GetNodeId(std::size_t I, std::size_t J, std::size_t K) {
		return mStartNodeId + (K * (mNumberOfDivisions[1] + 1) * (mNumberOfDivisions[0] + 1)) + (J * (mNumberOfDivisions[0] + 1)) + I;
	}

    int VoxelMeshGeneratorProcess::Check()
    {
        KRATOS_TRY

        return 0;

        KRATOS_CATCH("")
    }

}  // namespace Kratos.
