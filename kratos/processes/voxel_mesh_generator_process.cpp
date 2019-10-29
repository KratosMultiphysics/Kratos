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
#include "utilities/timer.h"



namespace Kratos
{
    VoxelMeshGeneratorProcess::VoxelMeshGeneratorProcess(Point const& MinPoint, Point const& MaxPoint,
        ModelPart& rVolumePart,
        ModelPart& rSkinPart, Parameters& TheParameters)
		: Process()
        , mMinPoint(MinPoint)
        , mMaxPoint(MaxPoint)
        , mrVolumePart(rVolumePart), mrSkinPart(rSkinPart), mCoarseMeshType(false) {

		Parameters default_parameters(R"(
            {
	            "create_skin_sub_model_part": true,
	            "start_node_id":1,
                "start_element_id":1,
                "start_condition_id":1,
                "number_of_divisions":[1,1,1],
                "elements_properties_id":0,
                "conditions_properties_id":0,
                "element_name": "PLEASE SPECIFY IT",
                "condition_name": "PLEASE SPECIFY IT",
				"coloring_settings_list": [],
				"entities_to_generate": "elements",
				"output" : "mesh",
				"mesh_type": "uniform",
				"output_filename" : ""
            }  )");

		TheParameters["element_name"]; // Should be given by caller! if not thorws an error

		TheParameters.ValidateAndAssignDefaults(default_parameters);

		mStartNodeId = TheParameters["start_node_id"].GetInt();
		mStartElementId = TheParameters["start_element_id"].GetInt();
		mStartConditionId = TheParameters["start_condition_id"].GetInt();

        mNumberOfDivisions[0] = TheParameters["number_of_divisions"].GetArrayItem(0).GetInt();
        mNumberOfDivisions[1] = TheParameters["number_of_divisions"].GetArrayItem(1).GetInt();
        mNumberOfDivisions[2] = TheParameters["number_of_divisions"].GetArrayItem(2).GetInt();
		mElementPropertiesId = TheParameters["elements_properties_id"].GetInt();
		mConditiongPropertiesId = TheParameters["conditions_properties_id"].GetInt();
		mElementName = TheParameters["element_name"].GetString();
		mConditionName = TheParameters["condition_name"].GetString();
        mCreateSkinSubModelPart = TheParameters["create_skin_sub_model_part"].GetBool();
		mOutputFilename = TheParameters["output_filename"].GetString();
        mCellSizes = mMaxPoint - mMinPoint;
		mColoringParameters = TheParameters["coloring_settings_list"];
        for(int i = 0 ; i < 3 ; i++)
            mCellSizes[i] /= mNumberOfDivisions[i];
		mEntitiesToGenerate=TheParameters["entities_to_generate"].GetString();
		if(TheParameters["mesh_type"].GetString() == "coarse")
			mCoarseMeshType=true;

		mOutput = TheParameters["output"].GetString();

		mHasColor = false;

        Check();

		array_1d<std::vector<double>, 3> coordinates;
		for(std::size_t i = 0 ; i < 3 ; i++){
			coordinates[i].resize(mNumberOfDivisions[i]+1);
		}

		for(std::size_t i = 0 ; i < coordinates.size() ; i++){
			const double min_coordinate_i = mMinPoint[i];
			for(std::size_t j = 0 ; j < coordinates[i].size() ; j++)
				coordinates[i][j] = j*mCellSizes[i] + min_coordinate_i;
		}
		
		mColors.SetCoordinates(coordinates[0], coordinates[1], coordinates[2]);

		double Margine = (mMaxPoint[0] - mMinPoint[0]) * 1.0e-2;

		mColors.ExtendBoundingBox(mrSkinPart.Nodes(), Margine);
    }

    VoxelMeshGeneratorProcess::VoxelMeshGeneratorProcess(std::vector<double> const& XCoordinates, std::vector<double> const& YCoordinates, std::vector<double> const& ZCoordinates,
        ModelPart& rVolumePart,
        ModelPart& rSkinPart, Parameters& TheParameters)
		: Process()
        , mMinPoint(XCoordinates.front(), YCoordinates.front(), ZCoordinates.front())
        , mMaxPoint(XCoordinates.back(), YCoordinates.back(), ZCoordinates.back())
        , mrVolumePart(rVolumePart), mrSkinPart(rSkinPart), mCoarseMeshType(false) {
		
		mColors.SetCoordinates(XCoordinates, YCoordinates, ZCoordinates);

		double Margine = (XCoordinates.back() - XCoordinates.front()) * 1.0e-2;

		mColors.ExtendBoundingBox(mrSkinPart.Nodes(), Margine);

		Parameters default_parameters(R"(
            {
	            "create_skin_sub_model_part": true,
	            "start_node_id":1,
                "start_element_id":1,
                "start_condition_id":1,
                "number_of_divisions":[1,1,1],
                "elements_properties_id":0,
                "conditions_properties_id":0,
                "element_name": "PLEASE SPECIFY IT",
                "condition_name": "PLEASE SPECIFY IT",
				"coloring_settings_list": [],
				"entities_to_generate": "elements",
				"mesh_type": "uniform",
				"output_filename" : "",
				"output" : "mesh"
            }  )");

		TheParameters["element_name"]; // Should be given by caller! if not thorws an error

		TheParameters.ValidateAndAssignDefaults(default_parameters);

		mStartNodeId = TheParameters["start_node_id"].GetInt();
		mStartElementId = TheParameters["start_element_id"].GetInt();
		mStartConditionId = TheParameters["start_condition_id"].GetInt();

        mNumberOfDivisions[0] = XCoordinates.size() - 1;
        mNumberOfDivisions[1] = YCoordinates.size() - 1;
        mNumberOfDivisions[2] = ZCoordinates.size() - 1;
		
		mElementPropertiesId = TheParameters["elements_properties_id"].GetInt();
		mConditiongPropertiesId = TheParameters["conditions_properties_id"].GetInt();
		mElementName = TheParameters["element_name"].GetString();
		mConditionName = TheParameters["condition_name"].GetString();
        mCreateSkinSubModelPart = TheParameters["create_skin_sub_model_part"].GetBool();
		mOutputFilename = TheParameters["output_filename"].GetString();
        mCellSizes = mMaxPoint - mMinPoint;
		mColoringParameters = TheParameters["coloring_settings_list"];
        for(int i = 0 ; i < 3 ; i++)
            mCellSizes[i] /= mNumberOfDivisions[i];
		mEntitiesToGenerate=TheParameters["entities_to_generate"].GetString();
		if(TheParameters["mesh_type"].GetString() == "coarse")
			mCoarseMeshType=true;

		mHasColor = false;
		mOutput = TheParameters["output"].GetString();

        Check();

    }

	VoxelMeshGeneratorProcess::~VoxelMeshGeneratorProcess() {

	}

	void VoxelMeshGeneratorProcess::Execute() {

		Timer::Start("Voxel Mesh Coloring");
		for(auto parameters : mColoringParameters){

			mHasColor=true;


			Parameters default_parameters(R"(
				{
					"model_part_name": "PLEASE SPECIFY IT",
					"inside_color": -1,
					"outside_color": 1,
					"apply_outside_color": true,
					"coloring_entities" : "nodes"
				}  )");

			parameters.ValidateAndAssignDefaults(default_parameters);
			
			std::string model_part_name = parameters["model_part_name"].GetString();
			ModelPart& skin_part = (model_part_name == mrSkinPart.Name()) ? mrSkinPart : mrSkinPart.GetSubModelPart(model_part_name);

			ApplyColoring(skin_part, parameters);
		}
		Timer::Stop("Voxel Mesh Coloring");

		if(mOutputFilename != "")			
			mColors.WriteParaViewVTR(mOutputFilename);

		if(mCoarseMeshType){
			Generate3DCoarseMesh();
		}
		else{
			if(mEntitiesToGenerate == "center_of_elements")
				GenerateCenterOfElements();
			else if(mEntitiesToGenerate == "nodes" || mEntitiesToGenerate == "elements"){
				Generate3DMesh();
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

	void VoxelMeshGeneratorProcess::ApplyColoring(ModelPart const& TheSkinModelPart, Parameters parameters) {

		std::string coloring_entities = parameters["coloring_entities"].GetString();

		if(mEntitiesToGenerate != "elements")
			KRATOS_ERROR_IF(coloring_entities == "elements") << "The coloring entities is set to element but there are no elements generated. "
															 << "Please set the entities_to_generate to 'elements' or set coloring entities to 'nodes'" << std::endl;

		double inside_color = parameters["inside_color"].GetDouble();
		double outside_color = parameters["outside_color"].GetDouble();
		bool apply_outside_color = parameters["apply_outside_color"].GetBool();

		if(apply_outside_color)
			mColors.SetAllColors(outside_color);
		
		array_1d< std::size_t, 3 > min_ray_position;
		array_1d< std::size_t, 3 > max_ray_position;
		if(coloring_entities == "center_of_elements"){
			mColors.CalculateMinMaxCenterOfElementPositions(TheSkinModelPart.Nodes(), min_ray_position, max_ray_position);
		}
		else{
			mColors.CalculateMinMaxNodePositions(TheSkinModelPart.Nodes(), min_ray_position, max_ray_position);
		}

		mColors.InitializeRays(min_ray_position, max_ray_position, coloring_entities);

		bool is_nodal = (coloring_entities == "node") ? true : false;
	
		for(auto& element : TheSkinModelPart.Elements())
		{
			Element::GeometryType& r_geometry = element.GetGeometry();
			mColors.AddGeometry(r_geometry, is_nodal);
		}

		if(coloring_entities == "nodes"){
			mColors.CalculateNodalRayColors(min_ray_position, max_ray_position, inside_color, outside_color);
		}
		else if(coloring_entities == "center_of_elements") {
			mColors.CalculateElementalRayColors(min_ray_position, max_ray_position, inside_color, outside_color);
		}
	}

	void VoxelMeshGeneratorProcess::Generate3DMesh() {
        if(!mrVolumePart.HasProperties(mElementPropertiesId))
            mrVolumePart.CreateNewProperties(mElementPropertiesId);

		if (mCreateSkinSubModelPart) {
			if(!mrVolumePart.HasSubModelPart("Skin")) {
				mrVolumePart.CreateSubModelPart("Skin");
			}
		}
		GenerateNodes3D();

		if(mOutput == "rectilinear_coordinates") { // we don't need to explictly generate the elements.
			return;
		}

		if(mEntitiesToGenerate == "elements") {

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
						auto p_element = mrVolumePart.CreateNewElement("Element3D8N", mStartElementId + cell_index, points, p_properties);
						cell_index++;
						if(mHasColor){
							p_element->GetValue(DISTANCE) = mColors.GetElementalColor(i,J,K);
						}
					}
				}
			}
		}
	}

	void VoxelMeshGeneratorProcess::Generate3DCoarseMesh(){
		Timer::Start("Creating Coarse Mesh");

		KRATOS_ERROR_IF(mEntitiesToGenerate != "center_of_elements") << "The coarse mesh can be generated only when entities to generate is center_of_elements" << std::endl;

		std::vector<bool> x_cell_coarse(mNumberOfDivisions[0]+1,false);
		std::vector<bool> y_cell_coarse(mNumberOfDivisions[1]+1,false);
		std::vector<bool> z_cell_coarse(mNumberOfDivisions[2]+1,false);

		x_cell_coarse[0]=true;
		y_cell_coarse[0]=true;
		z_cell_coarse[0]=true;

		x_cell_coarse[mNumberOfDivisions[0]]=true;
		y_cell_coarse[mNumberOfDivisions[1]]=true;
		z_cell_coarse[mNumberOfDivisions[2]]=true;

		for (std::size_t k = 1; k < mNumberOfDivisions[2]; k++) {
			for (std::size_t j = 1; j < mNumberOfDivisions[1]; j++) {
				for (std::size_t i = 1; i < mNumberOfDivisions[0]; i++) {
					double color = mColors.GetElementalColor(i,j,k);
					if(color != mColors.GetElementalColor(i-1,j,k))
						x_cell_coarse[i]=true;
					if(color != mColors.GetElementalColor(i,j-1,k))
						y_cell_coarse[j]=true;
					if(color != mColors.GetElementalColor(i,j,k-1))
						z_cell_coarse[k]=true;
				}
			}
		}

		if(mOutput == "rectilinear_coordinates") { 
			std::vector<double> x_key_planes;
			std::vector<double> y_key_planes;
			std::vector<double> z_key_planes;

			for(std::size_t i = 0 ; i < mNumberOfDivisions[0]; i++){
				if(x_cell_coarse[i])
					x_key_planes.push_back(i*mCellSizes[0]+mMinPoint[0]);
			}

			for(std::size_t i = 0 ; i < mNumberOfDivisions[1]; i++){
				if(y_cell_coarse[i])
					y_key_planes.push_back(i*mCellSizes[1]+mMinPoint[1]);
			}

			for(std::size_t i = 0 ; i < mNumberOfDivisions[2]; i++){
				if(z_cell_coarse[i])
					z_key_planes.push_back(i*mCellSizes[2]+mMinPoint[2]);
			}

			x_key_planes.push_back(mMaxPoint[0]);
			y_key_planes.push_back(mMaxPoint[1]);
			z_key_planes.push_back(mMaxPoint[2]);

			auto& x_coordinates = mrVolumePart.GetValue(RECTILINEAR_X_COORDINATES);
 			auto& y_coordinates = mrVolumePart.GetValue(RECTILINEAR_Y_COORDINATES);
 			auto& z_coordinates = mrVolumePart.GetValue(RECTILINEAR_Z_COORDINATES);

			x_coordinates.resize(x_key_planes.size(), false);
			y_coordinates.resize(y_key_planes.size(), false);
			z_coordinates.resize(z_key_planes.size(), false);

			std::copy(x_key_planes.begin(), x_key_planes.end(), x_coordinates.begin());
			std::copy(y_key_planes.begin(), y_key_planes.end(), y_coordinates.begin());
			std::copy(z_key_planes.begin(), z_key_planes.end(), z_coordinates.begin());

			auto& colors = mrVolumePart.GetValue(COLORS);
			colors.resize(mColors.GetNodalColors().size());
			for(int i = 0 ; i < static_cast<int>(colors.size()) ; i++){
				colors[i] = static_cast<int>(mColors.GetNodalColors()[i]);
			}
			std::size_t index = 0;
			for (std::size_t k = 0; k < mNumberOfDivisions[2]; k++) {
				for (std::size_t j = 0; j < mNumberOfDivisions[1]; j++) {
					for (std::size_t i = 0; i < mNumberOfDivisions[0]; i++) {
						if(x_cell_coarse[i]&y_cell_coarse[j]&z_cell_coarse[k]){
							colors[index++] = mColors.GetElementalColor(i,j,k);
						}
					}
				}
			}

			return;
		}
		
		std::size_t node_id = mStartNodeId;
		for (std::size_t k = 0; k < mNumberOfDivisions[2] + 1; k++) {
			for (std::size_t j = 0; j < mNumberOfDivisions[1] + 1; j++) {
				for (std::size_t i = 0; i < mNumberOfDivisions[0] + 1; i++) {
					if(x_cell_coarse[i]&y_cell_coarse[j]&z_cell_coarse[k]){
						Point global_coordinates = mColors.GetPoint(i,j,k);
						auto p_node = mrVolumePart.CreateNewNode(node_id++, global_coordinates[0],
                                                                   global_coordinates[1],
                                                                   global_coordinates[2]);
						if(mHasColor ){
							if((i == mNumberOfDivisions[0]) || (j == mNumberOfDivisions[1]) || (k == mNumberOfDivisions[2])){
								auto ii = (i== mNumberOfDivisions[0]) ? i - 1 : i;
								auto jj = (j== mNumberOfDivisions[1]) ? j - 1 : j;
								auto kk = (k== mNumberOfDivisions[2]) ? k - 1 : k;
								p_node->GetSolutionStepValue(DISTANCE) = mColors.GetElementalColor(ii,jj,kk);
							}
							else {
								p_node->GetSolutionStepValue(DISTANCE) = mColors.GetElementalColor(i,j,k);
							}
						}
					}
				}
			}
		}

		Timer::Stop("Creating Coarse Mesh");
			
	}

	void VoxelMeshGeneratorProcess::GenerateNodes3D() {
		GeometryType::CoordinatesArrayType local_element_size = mCellSizes;
		Point local_coordinates = mMinPoint;
		auto global_coordinates = Point{ZeroVector(3)};
		std::size_t node_id = mStartNodeId;
		ModelPart* p_skin_part = nullptr;
		if (mCreateSkinSubModelPart)
			p_skin_part = &(mrVolumePart.GetSubModelPart("Skin"));

		if(mOutput == "rectilinear_coordinates") { 
			auto& x_coordinates = mrVolumePart.GetValue(RECTILINEAR_X_COORDINATES);
			auto& y_coordinates = mrVolumePart.GetValue(RECTILINEAR_Y_COORDINATES);
			auto& z_coordinates = mrVolumePart.GetValue(RECTILINEAR_Z_COORDINATES);

			x_coordinates.resize(mColors.GetNodalCoordinates(0).size(), false);
			y_coordinates.resize(mColors.GetNodalCoordinates(1).size(), false);
			z_coordinates.resize(mColors.GetNodalCoordinates(2).size(), false);

			std::copy(mColors.GetNodalCoordinates(0).begin(), mColors.GetNodalCoordinates(0).end(), x_coordinates.begin());
			std::copy(mColors.GetNodalCoordinates(1).begin(), mColors.GetNodalCoordinates(1).end(), y_coordinates.begin());
			std::copy(mColors.GetNodalCoordinates(2).begin(), mColors.GetNodalCoordinates(2).end(), z_coordinates.begin());

			auto& colors = mrVolumePart.GetValue(COLORS);
			colors.resize(mColors.GetElementalColors().size());
			for(int i = 0 ; i < static_cast<int>(colors.size()) ; i++){
				colors[i] = static_cast<int>(mColors.GetElementalColors()[i]);
			}

			return;
		}

		for (std::size_t k = 0; k < mNumberOfDivisions[2] + 1; k++) {
			for (std::size_t j = 0; j < mNumberOfDivisions[1] + 1; j++) {
				for (std::size_t i = 0; i < mNumberOfDivisions[0] + 1; i++) {
					Point global_coordinates = mColors.GetPoint(i,j,k);
					if (mCreateSkinSubModelPart && (
                        (i == 0) || (i == mNumberOfDivisions[0]) || (j == 0) ||
                        (j == mNumberOfDivisions[1]) || (k == 0) || (k == mNumberOfDivisions[2]))) { // Is on skin
						auto p_node = p_skin_part->CreateNewNode(node_id++, global_coordinates[0],
                                                                                           global_coordinates[1],
                                                                                           global_coordinates[2]);
						if(mHasColor){
							p_node->GetSolutionStepValue(DISTANCE) = mColors.GetNodalColor(i,j,k);
						}
					}
					else {
						auto p_node = mrVolumePart.CreateNewNode(node_id++, global_coordinates[0],
                                                                   global_coordinates[1],
                                                                   global_coordinates[2]);
						if(mHasColor){
							p_node->GetSolutionStepValue(DISTANCE) = mColors.GetNodalColor(i,j,k);
						}
					}
				}
			}
		}

	}

	void VoxelMeshGeneratorProcess::GenerateCenterOfElements() {
		GeometryType::CoordinatesArrayType local_element_size = mCellSizes;
		Point local_coordinates = mMinPoint;
		auto global_coordinates = Point{ZeroVector(3)};
		std::size_t node_id = mStartNodeId;

		if(mOutput == "rectilinear_coordinates") { 
			auto& x_coordinates = mrVolumePart.GetValue(RECTILINEAR_X_COORDINATES);
			auto& y_coordinates = mrVolumePart.GetValue(RECTILINEAR_Y_COORDINATES);
			auto& z_coordinates = mrVolumePart.GetValue(RECTILINEAR_Z_COORDINATES);

			x_coordinates.resize(mColors.GetElementCenterCoordinates(0).size(), false);
			y_coordinates.resize(mColors.GetElementCenterCoordinates(1).size(), false);
			z_coordinates.resize(mColors.GetElementCenterCoordinates(2).size(), false);

			std::copy(mColors.GetElementCenterCoordinates(0).begin(), mColors.GetElementCenterCoordinates(0).end(), x_coordinates.begin());
			std::copy(mColors.GetElementCenterCoordinates(1).begin(), mColors.GetElementCenterCoordinates(1).end(), y_coordinates.begin());
			std::copy(mColors.GetElementCenterCoordinates(2).begin(), mColors.GetElementCenterCoordinates(2).end(), z_coordinates.begin());


			return;
		}

		for (std::size_t k = 0; k < mNumberOfDivisions[2]; k++) {
			for (std::size_t j = 0; j < mNumberOfDivisions[1]; j++) {
				for (std::size_t i = 0; i < mNumberOfDivisions[0]; i++) {
					Point global_coordinates = mColors.GetCenterOfElement(i,j,k);
					auto p_node = mrVolumePart.CreateNewNode(node_id++, global_coordinates[0],
                                                          global_coordinates[1],
                                                          global_coordinates[2]);
					if(mHasColor){
						p_node->GetSolutionStepValue(DISTANCE) = mColors.GetElementalColor(i,j,k);
					}
				}
			}
		}

		Timer::Stop("GenerateCenterOfElements");
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
