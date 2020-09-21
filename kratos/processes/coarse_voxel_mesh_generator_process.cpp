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
#include "processes/coarse_voxel_mesh_generator_process.h"
#include "geometries/geometry.h"
#include "geometries/point.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/hexahedra_3d_8.h"
#include "includes/checks.h"
#include "utilities/timer.h"



namespace Kratos
{

    CoarseVoxelMeshGeneratorProcess::CoarseVoxelMeshGeneratorProcess(Point const& MinPoint, Point const& MaxPoint,
        ModelPart& rVolumePart,
        ModelPart& rSkinPart, Parameters& TheParameters)
		: VoxelMeshGeneratorProcess(MinPoint, MaxPoint, rVolumePart, rSkinPart, TheParameters) {
			
		}

    CoarseVoxelMeshGeneratorProcess::CoarseVoxelMeshGeneratorProcess(std::vector<double> const& XCoordinates, std::vector<double> const& YCoordinates, std::vector<double> const& ZCoordinates,
        ModelPart& rVolumePart,
        ModelPart& rSkinPart, Parameters& TheParameters)
		: VoxelMeshGeneratorProcess(XCoordinates, YCoordinates, ZCoordinates, rVolumePart, rSkinPart, TheParameters) {
    }

	CoarseVoxelMeshGeneratorProcess::~CoarseVoxelMeshGeneratorProcess() {

	}

	void CoarseVoxelMeshGeneratorProcess::Execute() {

		Timer::Start("Voxel Mesh Coloring");
		for(auto parameters : mColoringParameters){

			mHasColor=true;


			Parameters default_parameters(R"(
				{
					"model_part_name": "PLEASE SPECIFY IT",
					"inside_color": -1,
					"outside_color": 1,
					"interface_color": 0,
					"apply_outside_color": true,
					"coloring_entities" : "nodes"
				}  )");

			parameters.ValidateAndAssignDefaults(default_parameters);
			
			std::string model_part_name = parameters["model_part_name"].GetString();
			ModelPart& skin_part = (model_part_name == mrSkinPart.Name()) ? mrSkinPart : mrSkinPart.GetSubModelPart(model_part_name);

			ApplyColoring(skin_part, parameters);
		}
		Timer::Stop("Voxel Mesh Coloring");

		Generate3DCoarseMesh();

		if(mOutputFilename != "")			
			mCoarseMeshColors.WriteParaViewVTR(mOutputFilename);
	}

	std::string CoarseVoxelMeshGeneratorProcess::Info() const {
		return "CoarseVoxelMeshGeneratorProcess";
	}

	void CoarseVoxelMeshGeneratorProcess::PrintInfo(std::ostream& rOStream) const {
		rOStream << Info();
	}

	void CoarseVoxelMeshGeneratorProcess::PrintData(std::ostream& rOStream) const {

	}

	void CoarseVoxelMeshGeneratorProcess::Generate3DCoarseMesh(){
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

		for (std::size_t k = 0; k < mNumberOfDivisions[2]; k++) {
			for (std::size_t j = 0; j < mNumberOfDivisions[1]; j++) {
				for (std::size_t i = 0; i < mNumberOfDivisions[0]; i++) {
					double color = mColors.GetElementalColor(i,j,k);
					auto face_color = mColors.GetElementalFaceColor(i,j,k);
					if((i > 0) && (color != mColors.GetElementalColor(i-1,j,k)))
						x_cell_coarse[i]=true;
					if((j > 0) && (color != mColors.GetElementalColor(i,j-1,k)))
						y_cell_coarse[j]=true;
					if((k > 0) && (color != mColors.GetElementalColor(i,j,k-1)))
						z_cell_coarse[k]=true;


					if(i > 0){
						auto& previous_x_face_color = mColors.GetElementalFaceColor(i-1, j,k);
						for(std::size_t i_face = 0 ; i_face < 6 ; i_face++){
							if((i_face != 0) && (i_face != 3)) // assuming that there are no face condition inside a volume
								if(face_color[i_face] != previous_x_face_color[i_face]) 
									x_cell_coarse[i]=true;
						}
					}
						
					if(j > 0){
						auto& previous_y_face_color = mColors.GetElementalFaceColor(i,j-1,k);
						for(std::size_t i_face = 0 ; i_face < 6 ; i_face++){
							if((i_face != 1) && (i_face != 4)) // assuming that there are no face condition inside a volume
								if(face_color[i_face] != previous_y_face_color[i_face]) 
									y_cell_coarse[j]=true;
						}
					}
					
					if(k > 0){
						auto& previous_z_face_color = mColors.GetElementalFaceColor(i,j,k-1);
						for(std::size_t i_face = 0 ; i_face < 6 ; i_face++){
							if((i_face != 2) && (i_face != 5)) // assuming that there are no face condition inside a volume
								if(face_color[i_face] != previous_z_face_color[i_face])
									z_cell_coarse[k]=true;
						}
					}
				}
			}
		}

			std::vector<double> x_key_planes;
			std::vector<double> y_key_planes;
			std::vector<double> z_key_planes;

			std::vector<std::size_t> min_i;
			std::vector<std::size_t> min_j;
			std::vector<std::size_t> min_k;

			std::vector<std::size_t> max_i;
			std::vector<std::size_t> max_j;
			std::vector<std::size_t> max_k;

			for(std::size_t i = 0 ; i < mNumberOfDivisions[0]; i++){
				if(x_cell_coarse[i]){
					x_key_planes.push_back(mColors.GetNodalCoordinates(0)[i]);
					min_i.push_back(i);
				}
				if(x_cell_coarse[i+1]){
					max_i.push_back(i);
				}
			}

			for(std::size_t i = 0 ; i < mNumberOfDivisions[1]; i++){
				if(y_cell_coarse[i]){
					y_key_planes.push_back(mColors.GetNodalCoordinates(1)[i]);
					min_j.push_back(i);
				}
				if(y_cell_coarse[i+1]){
					max_j.push_back(i);
				}
			}

			for(std::size_t i = 0 ; i < mNumberOfDivisions[2]; i++){
				if(z_cell_coarse[i]){
					z_key_planes.push_back(mColors.GetNodalCoordinates(2)[i]);
					min_k.push_back(i);
				}
				if(z_cell_coarse[i+1]){
					max_k.push_back(i);
				}
			}

			x_key_planes.push_back(mColors.GetNodalCoordinates(0).back());
			y_key_planes.push_back(mColors.GetNodalCoordinates(1).back());
			z_key_planes.push_back(mColors.GetNodalCoordinates(2).back());

		mCoarseMeshColors.SetCoordinates(x_key_planes, y_key_planes, z_key_planes);

		auto coarse_x_size=x_key_planes.size() - 1;
		auto coarse_y_size=y_key_planes.size() - 1;
		auto coarse_z_size=z_key_planes.size() - 1;

		for (std::size_t coarse_k = 0; coarse_k < coarse_z_size; coarse_k++) {
			auto k = min_k[coarse_k];
			for (std::size_t coarse_j = 0; coarse_j < coarse_y_size; coarse_j++) {
				auto j = min_j[coarse_j];
				for (std::size_t coarse_i = 0; coarse_i < coarse_x_size; coarse_i++) {
					auto i = min_i[coarse_i];
					mCoarseMeshColors.GetElementalColor(coarse_i, coarse_j, coarse_k) = mColors.GetElementalColor(i,j,k);
					auto& coarse_face_color = mCoarseMeshColors.GetElementalFaceColor(coarse_i, coarse_j, coarse_k);
					auto min_face_color = mColors.GetElementalFaceColor(i,j,k);
					auto max_face_color = mColors.GetElementalFaceColor(max_i[coarse_i],max_j[coarse_j],max_k[coarse_k]);
					coarse_face_color[0] =  min_face_color[0];
					coarse_face_color[1] =  min_face_color[1];
					coarse_face_color[2] =  min_face_color[2];
					coarse_face_color[3] =  max_face_color[3];
					coarse_face_color[4] =  max_face_color[4];
					coarse_face_color[5] =  max_face_color[5];

				}
			}
		}

		if(mOutput == "rectilinear_coordinates") { 

			auto& x_coordinates = mrVolumePart.GetValue(RECTILINEAR_X_COORDINATES);
 			auto& y_coordinates = mrVolumePart.GetValue(RECTILINEAR_Y_COORDINATES);
 			auto& z_coordinates = mrVolumePart.GetValue(RECTILINEAR_Z_COORDINATES);

			x_coordinates.resize(x_key_planes.size(), false);
			y_coordinates.resize(y_key_planes.size(), false);
			z_coordinates.resize(z_key_planes.size(), false);

			std::copy(x_key_planes.begin(), x_key_planes.end(), x_coordinates.begin());
			std::copy(y_key_planes.begin(), y_key_planes.end(), y_coordinates.begin());
			std::copy(z_key_planes.begin(), z_key_planes.end(), z_coordinates.begin());

			const std::size_t size = mCoarseMeshColors.GetElementalColors().size();

			auto& colors = mrVolumePart.GetValue(COLORS);
			colors.resize(size, false);
			std::copy(mCoarseMeshColors.GetElementalColors().begin(), mCoarseMeshColors.GetElementalColors().end(), colors.begin());
			
			auto& face_colors = mrVolumePart.GetValue(VOXEL_FACE_COLORS);
			face_colors.resize(size, 6, false);

			for(std::size_t i_element = 0 ; i_element < size ; i_element++){
				for(std::size_t j_face = 0 ; j_face < 6 ; j_face++){
					face_colors(i_element, j_face) = mCoarseMeshColors.GetElementalFaceColors()[i_element][j_face];
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

    int CoarseVoxelMeshGeneratorProcess::Check()
    {
        KRATOS_TRY

		KRATOS_ERROR_IF(mEntitiesToGenerate != "nodes" && mEntitiesToGenerate != "elements" && mEntitiesToGenerate != "center_of_elements") << mEntitiesToGenerate 
			<< " is not accepted as entities_to_generate. The valid options are 'nodes', 'elements' and 'center_of_elements'." << std::endl;

        return 0;

        KRATOS_CATCH("")
    }

}  // namespace Kratos.
