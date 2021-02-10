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
#include "input_output/vtr_io.h"



namespace Kratos
{
    VoxelMeshMapperProcess::VoxelMeshMapperProcess(std::string InputFileName, ModelPart& rVolumePart,
        ModelPart& rSkinPart, Parameters& TheParameters) : VoxelMeshGeneratorProcess(rVolumePart, rSkinPart, TheParameters)
	{

        TheParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());

        ConstructFromParameters(TheParameters);

        mMappingParameters = TheParameters["mapping_results_to_variable"];

        KRATOS_WATCH(mMappingParameters);

		Timer::Start("Reading Input Mesh");
        ReadInputVtrFile(InputFileName);
        ReadInputFemFile(InputFileName);
		Timer::Stop("Reading Input Mesh");

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
            if(model_part_name == mrVolumePart.Name()){
                MapResults(mrVolumePart, parameters);
            }
            else if(mrVolumePart.HasSubModelPart(model_part_name)){
                MapResults(mrVolumePart.GetSubModelPart(model_part_name), parameters);
            }
            else{
                KRATOS_WARNING("Mapper") << "The given volume modelpart does not have the sub modelpart \"" << model_part_name << "\"" << std::endl;
            }

			
		}

		Timer::Stop("Mapping The Results");        
	}


	const Parameters VoxelMeshMapperProcess::GetDefaultParameters() const {
		Parameters default_parameters(R"(
            {
                "mapping_results_to_variable": []
            }  )");

		// Getting base class default parameters
        const Parameters base_default_parameters = VoxelMeshGeneratorProcess::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);

		return default_parameters;
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

	void VoxelMeshMapperProcess::ReadInputVtrFile(std::string InputFileName) {
        VtrIO vtr_io(InputFileName + ".vtr");
        std::vector<Kratos::Internals::CartesianMeshColors> multi_block_mesh;

        vtr_io.Read(multi_block_mesh);

        KRATOS_ERROR_IF(multi_block_mesh.size() != 1) << "The voxel mesh mapper only works with one block" << std::endl;

        mInputMesh = multi_block_mesh[0];
        mColors.SetCoordinates(mInputMesh.GetNodalCoordinates(0), mInputMesh.GetNodalCoordinates(1), mInputMesh.GetNodalCoordinates(2)); 

        mInputMesh.WriteParaViewVTR(InputFileName + "_input.vtr"); 
	}

	void VoxelMeshMapperProcess::ReadInputFemFile(std::string InputFileName) {
        std::ifstream input(InputFileName + ".fem");
        std::string line;
        ModelPart* p_current_sub_model_part = nullptr;
        if(!mrVolumePart.HasProperties(0))
            mrVolumePart.CreateNewProperties(0);

        Properties::Pointer p_properties = mrVolumePart.pGetProperties(0);

        while(!input.eof()){
            std::getline(input, line);
            if(line.substr(0,18) == "Begin SubModelPart"){
                std::string sub_model_part_name=line.substr(19);
                sub_model_part_name.erase(sub_model_part_name.begin(), std::find_if(sub_model_part_name.begin(), sub_model_part_name.end(), [](unsigned char ch) {return !std::isspace(ch);}));
                sub_model_part_name.erase(std::find_if(sub_model_part_name.begin(), sub_model_part_name.end(), [](unsigned char ch) {return std::isspace(ch);}), sub_model_part_name.end());
                KRATOS_INFO("Mapper") << "Reading modelpart \"" << sub_model_part_name << "\"" << std::endl;
                if(!mrVolumePart.HasSubModelPart(sub_model_part_name)){
                    p_current_sub_model_part = &mrVolumePart.CreateSubModelPart(sub_model_part_name);
                }
                else{
                    p_current_sub_model_part = &mrVolumePart.GetSubModelPart(sub_model_part_name);
                }
            }
            if(line.substr(0,4) == "GRID"){
                std::size_t index = std::stoi(line.substr(4,12));
                double x = ReadDouble(line.substr(24, 8));
                double y = ReadDouble(line.substr(32, 8));
                double z = ReadDouble(line.substr(40));

                if(x < mInputMesh.GetMinPoint()[0] || x > mInputMesh.GetMaxPoint()[0]){
                    KRATOS_WATCH(line);
                    std::cout << index << " : " << x << " " << y << " " << z << std::endl;
                }

                p_current_sub_model_part->CreateNewNode(index, x, y, z);
            }
            if(line.substr(0,6) == "CTETRA"){
                std::size_t index = std::stoi(line.substr(7,12));
                std::size_t n1 = std::stoi(line.substr(24, 8));
                std::size_t n2 = std::stoi(line.substr(32, 8));
                std::size_t n3 = std::stoi(line.substr(40, 8));
                std::size_t n4 = std::stoi(line.substr(48));

                // KRATOS_INFO("Mapper") << "Reading element " << index << " with nodes " << n1 << "," << n2 << "," << n3 << "," << n4 << std::endl;

                p_current_sub_model_part->CreateNewElement("Element3D4N", index, {n1, n2, n3, n4}, p_properties);
            }
        }
        KRATOS_WATCH(mrVolumePart);
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

        std::vector<std::pair<Variable<double>, std::string>> variables_to_map;

        for(auto mapping_parameters : mMappingParameters){
            std::string map_from = mapping_parameters["from_result"].GetString();
            std::string variable_name = mapping_parameters["to_variable"].GetString();

            KRATOS_ERROR_IF(!mInputMesh.HasElementalData(map_from)) << "The data to map \"" << map_from << "\" is not in the input data" << std::endl;
            if(KratosComponents<Variable<double> >::Has(variable_name))
            {
                auto& the_variable  = KratosComponents<Variable<double> >::Get(variable_name);
                variables_to_map.push_back(std::make_pair(the_variable,map_from));
            }

            KRATOS_WATCH(map_from);
            KRATOS_WATCH(variable_name);
        }
//        block_for_each(TheModelPart.Nodes(), [&](Node<3>& rNode)
        for(auto& rNode : TheModelPart.Nodes())
        {
            std::size_t i_position = mInputMesh.CalculateCellPosition(rNode.X(), 0);
            std::size_t j_position = mInputMesh.CalculateCellPosition(rNode.Y(), 1);
            std::size_t k_position = mInputMesh.CalculateCellPosition(rNode.Z(), 2);

            double cell_color = mColors.GetElementalColor(i_position,j_position,k_position);

            if(cell_color == inside_color){
                for(auto& variable_data : variables_to_map){
                    rNode.GetSolutionStepValue(variable_data.first) = mInputMesh.GetElementalData(variable_data.second , i_position,j_position,k_position);
                }
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
                for(auto& variable_data : variables_to_map){
                    rNode.GetSolutionStepValue(variable_data.first) = mInputMesh.GetElementalData(variable_data.second, nearest_i,nearest_j,nearest_k);
                }
            }
        }
        //);

    }

    double VoxelMeshMapperProcess::ReadDouble(std::string&& Input){
        std::size_t pos = Input.rfind('-');

        if(pos != 0 && pos != std::string::npos){
            Input.replace(pos, 1, "e-");
        }

        pos = Input.rfind('+');
        if(pos != 0 && pos != std::string::npos){
            Input.replace(pos, 1, "e+");
        }

        return std::stod(Input);
    }


}  // namespace Kratos.
