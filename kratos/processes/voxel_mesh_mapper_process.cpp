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
#include <set>

// External includes

// Project includes
#include "processes/voxel_mesh_mapper_process.h"
#include "geometries/geometry.h"
#include "geometries/point.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/hexahedra_3d_8.h"
#include "includes/checks.h"
#include "utilities/timer.h"
#include "utilities/parallel_utilities.h"



namespace Kratos
{
    VoxelMeshMapperProcess::VoxelMeshMapperProcess(std::string InputFileName, ModelPart& rVolumePart,
        ModelPart& rSkinPart, Parameters& TheParameters) : VoxelMeshGeneratorProcess(rVolumePart, rSkinPart, TheParameters)
	{

		Timer::Start("Reading Input Mesh");
        ReadInputFile(InputFileName);
		Timer::Stop("Reading Input Mesh");

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

		TheParameters.ValidateAndAssignDefaults(default_parameters);

        mMinPoint = mColors.GetMinPoint();
        mMaxPoint = mColors.GetMaxPoint();

        Check();


    }

	VoxelMeshMapperProcess::~VoxelMeshMapperProcess() {

	}

	void VoxelMeshMapperProcess::Execute() {

		Timer::Start("Coloring Input Mesh");
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
		Timer::Stop("Coloring Input Mesh");

        mColors.WriteParaViewVTR("input_mesh_colors.vtr");

		Timer::Start("Mapping The Results");
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

			MapResults(skin_part, parameters);
		}
		Timer::Stop("Mapping The Results");        
	}

	std::string VoxelMeshMapperProcess::Info() const {
		return "VoxelMeshMapperProcess";
	}

	void VoxelMeshMapperProcess::PrintInfo(std::ostream& rOStream) const {
		rOStream << Info();
	}

	void VoxelMeshMapperProcess::PrintData(std::ostream& rOStream) const {

	}

	void VoxelMeshMapperProcess::ReadInputFile(std::string InputFileName) {
        std::ifstream input(InputFileName);
        std::size_t n_x = 0;
        std::size_t n_y = 0;
        std::size_t n_z = 0;

        input >> n_x >> n_y >> n_z;

        const int size = n_x * n_y * n_z;
        KRATOS_WATCH(size);
        Matrix input_data(size, 4); // x y z value

        std::set<double> x_centers;
        std::set<double> y_centers;
        std::set<double> z_centers;
        for(int i = 0 ; i < size ; i++){
            double& x = input_data(i,0);
            double& y = input_data(i,1);
            double& z = input_data(i,2);
            double& result = input_data(i,3);
            input >> x >> y >> z >> result;
            x_centers.insert(x);
            y_centers.insert(y);
            z_centers.insert(z);
        }

        // Passing center coordinates to the nodal coordinates
        std::vector<double> x_planes(x_centers.begin(), x_centers.end());
        std::vector<double> y_planes(y_centers.begin(), y_centers.end());
        std::vector<double> z_planes(z_centers.begin(), z_centers.end());

        CenterToNodalCoordinates(x_planes);
        CenterToNodalCoordinates(y_planes);
        CenterToNodalCoordinates(z_planes);


        mInputMesh.SetCoordinates(x_planes, y_planes, z_planes); 
        mColors.SetCoordinates(x_planes, y_planes, z_planes); 

        // Applying nodal data
        for(int i_node = 0 ; i_node < size ; i_node++){
            std::size_t i = mInputMesh.CalculateCellPosition(input_data(i_node,0), 0);
            std::size_t j = mInputMesh.CalculateCellPosition(input_data(i_node,1), 1);
            std::size_t k = mInputMesh.CalculateCellPosition(input_data(i_node,2), 2);
            KRATOS_ERROR_IF(i>n_x || j > n_y || k > n_z) << "The given point (" << input_data(i_node,0) << "," << input_data(i_node,1) << "," << input_data(i_node,2) << ") is out of mesh" << std::endl;
            mInputMesh.GetElementalColor(i,j,k) = input_data(i_node,3);
            std::cout << i << "," << j << "," << k << ":" << input_data(i_node,3) << std::endl;
        }

        mInputMesh.WriteParaViewVTR(InputFileName + "_input.vtr"); 


	}

    void VoxelMeshMapperProcess::CenterToNodalCoordinates(std::vector<double>& rCoordinates){
        KRATOS_ERROR_IF(rCoordinates.size() < 2) << "At least two center coordinates are needed" << std::endl;
        std::vector<double> nodal_coordinates;

        double half_cell_size = (rCoordinates[1] - rCoordinates[0]) * 0.5;
        nodal_coordinates.push_back(rCoordinates[0] - half_cell_size);

        for(std::size_t i = 0 ; i < rCoordinates.size()-1 ; i++){
            nodal_coordinates.push_back((rCoordinates[i+1]+rCoordinates[i])*0.5);
        }

        half_cell_size = rCoordinates.back() - nodal_coordinates.back();
        nodal_coordinates.push_back(rCoordinates.back() + half_cell_size);
        
        rCoordinates.swap(nodal_coordinates);
    }

    int VoxelMeshMapperProcess::Check()
    {
        KRATOS_TRY

        return 0;

        KRATOS_CATCH("")
    }

    void VoxelMeshMapperProcess::MapResults(ModelPart const& TheModelPart, Parameters parameters){

        double inside_color = parameters["inside_color"].GetDouble();
        std::size_t offset=1;

        const auto& coordinates_x =  mInputMesh.GetElementCenterCoordinates(0);
        const auto& coordinates_y =  mInputMesh.GetElementCenterCoordinates(1);
        const auto& coordinates_z =  mInputMesh.GetElementCenterCoordinates(2);

        const std::size_t n_x =  coordinates_x.size();
        const std::size_t n_y =  coordinates_y.size();
        const std::size_t n_z =  coordinates_z.size();

//        block_for_each(TheModelPart.Nodes(), [&](Node<3>& rNode)
        for(auto& rNode : TheModelPart.Nodes())
        {
            std::size_t i_position = mInputMesh.CalculateCellPosition(rNode.X(), 0);
            std::size_t j_position = mInputMesh.CalculateCellPosition(rNode.Y(), 1);
            std::size_t k_position = mInputMesh.CalculateCellPosition(rNode.Z(), 2);

            double cell_color = mColors.GetElementalColor(i_position,j_position,k_position);

            if(cell_color == inside_color){
                rNode.GetSolutionStepValue(TEMPERATURE) = mInputMesh.GetElementalColor(i_position,j_position,k_position);
            }
            else{
                double min_distance = std::numeric_limits<double>::max();


                std::size_t min_i = (i_position >= offset) ? i_position - offset : 0;
                std::size_t min_j = (j_position >= offset) ? j_position - offset : 0;
                std::size_t min_k = (k_position >= offset) ? k_position - offset : 0;
                std::size_t max_i = (i_position < n_x - offset) ? i_position + offset : n_x;
                std::size_t max_j = (j_position < n_y - offset) ? j_position + offset : n_y;
                std::size_t max_k = (k_position < n_z - offset) ? k_position + offset : n_z;
                std::size_t nearest_i = i_position;
                std::size_t nearest_j = j_position;
                std::size_t nearest_k = k_position;


                Point point_1(coordinates_x[i_position], coordinates_y[j_position], coordinates_z[k_position]);

                for(std::size_t i = min_i ; i < max_i ; i++){
                    for(std::size_t j = min_j ; j < max_j ; j++){
                        for(std::size_t k = min_k ; k < max_k ; k++){
                            if(mColors.GetElementalColor(i,j,k) == inside_color){
                                Point point_2(coordinates_x[i], coordinates_y[j], coordinates_z[k]);
                                double distance = norm_2(point_2 - point_1);
                                if(distance < min_distance){
                                    nearest_i = i;
                                    nearest_j = j;
                                    nearest_k = k;
                                }
                            }
                        }
                    }
                }
                rNode.GetSolutionStepValue(TEMPERATURE) = mInputMesh.GetElementalColor(nearest_i,nearest_j,nearest_k);
            }
        }
        //);

    }

}  // namespace Kratos.
