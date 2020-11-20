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
#include <unordered_set>

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
    VoxelMeshMapperProcess::VoxelMeshMapperProcess(std::string InputFileName,
        ModelPart& rSkinPart, Parameters& TheParameters)
		: Process()
        , mrSkinPart(rSkinPart) {

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

        Check();


    }

	VoxelMeshMapperProcess::~VoxelMeshMapperProcess() {

	}

	void VoxelMeshMapperProcess::Execute() {

		Timer::Start("Reading Input Mesh");
		Timer::Stop("Reading Input Mesh");
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

        std::unordered_set<double> x_planes;
        std::unordered_set<double> y_planes;
        std::unordered_set<double> z_planes;
        for(int i = 0 ; i < size ; i++){
            double& x = input_data(i,0);
            double& y = input_data(i,1);
            double& z = input_data(i,2);
            double result = input_data(i,3);
            input >> x >> y >> z >> result;
            x_planes.insert(x);
            y_planes.insert(y);
            z_planes.insert(z);
        }
        KRATOS_WATCH(x_planes.size());
        KRATOS_WATCH(y_planes.size());
        KRATOS_WATCH(z_planes.size());

        mInputMesh.SetCoordinates(std::vector<double>(x_planes.begin(), x_planes.end()), 
                                  std::vector<double>(y_planes.begin(), y_planes.end()),
                                  std::vector<double>(z_planes.begin(), z_planes.end())); 

        // Applying nodal data
        for(int i_node = 0 ; i_node < size ; i_node++){
            std::size_t i = mInputMesh.CalculateNodePosition(input_data(i_node,0), 0);
            std::size_t j = mInputMesh.CalculateNodePosition(input_data(i_node,1), 1);
            std::size_t k = mInputMesh.CalculateNodePosition(input_data(i_node,2), 2);
            mInputMesh.GetNodalColor(i,j,k) = input_data(i_node,3);
        }

        mInputMesh.WriteParaViewVTR(InputFileName + ".vtr"); 


	}


    int VoxelMeshMapperProcess::Check()
    {
        KRATOS_TRY

        return 0;

        KRATOS_CATCH("")
    }

}  // namespace Kratos.
