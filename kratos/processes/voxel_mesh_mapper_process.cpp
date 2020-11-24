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
        int n_x = 0;
        int n_y = 0;
        int n_z = 0;

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

        // Applying nodal data
        for(int i_node = 0 ; i_node < size ; i_node++){
            double& x = input_data(i_node,0);
            double& y = input_data(i_node,1);
            double& z = input_data(i_node,2);
            double& result = input_data(i_node,3);
            std::size_t i = mInputMesh.CalculateCellPosition(input_data(i_node,0), 0);
            std::size_t j = mInputMesh.CalculateCellPosition(input_data(i_node,1), 1);
            std::size_t k = mInputMesh.CalculateCellPosition(input_data(i_node,2), 2);
            KRATOS_ERROR_IF(i>n_x || j > n_y || k > n_z) << "The given point (" << input_data(i_node,0) << "," << input_data(i_node,1) << "," << input_data(i_node,2) << ") is out of mesh" << std::endl;
            mInputMesh.GetElementalColor(i,j,k) = input_data(i_node,3);
            std::cout << i << "," << j << "," << k << ":" << input_data(i_node,3) << std::endl;
        }

        mInputMesh.WriteParaViewVTR(InputFileName + ".vtr"); 


	}

    void VoxelMeshMapperProcess::CenterToNodalCoordinates(std::vector<double>& rCoordinates){
        KRATOS_ERROR_IF(rCoordinates.size() < 2) << "At least two center coordinates are needed" << std::endl;
        std::vector<double> nodal_coordinates;

        double half_cell_size = (rCoordinates[1] - rCoordinates[0]) * 0.5;
        nodal_coordinates.push_back(rCoordinates[0] - half_cell_size);

        for(int i = 0 ; i < rCoordinates.size()-1 ; i++){
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

}  // namespace Kratos.
