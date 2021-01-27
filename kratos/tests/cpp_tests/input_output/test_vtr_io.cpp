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

// System includes
#include <sstream>

// External includes


// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "input_output/vtr_io.h"



namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(ReadMeshFromAsciiVTR, KratosCoreFastSuite)
{
std::stringstream* p_input = new std::stringstream(R"input(
<VTKFile type="RectilinearGrid" version="0.1" byte_order="LittleEndian">
<RectilinearGrid WholeExtent="0 10 0 10 0 10">
<Piece Extent="0 10 0 10 0 10">
<Coordinates> 
<DataArray type="Float64" Name="coordinates" NumberOfComponents="1" format="ascii">
0 1 2 3 4 5 6 7 8 9 10 </DataArray> 
<DataArray type="Float64" Name="coordinates" NumberOfComponents="1" format="ascii">
0 1 2 3 4 5 6 7 8 9 10 </DataArray> 
<DataArray type="Float64" Name="coordinates" NumberOfComponents="1" format="ascii">
0 1 2 3 4 5 6 7 8 9 10 </DataArray> 
</Coordinates> 
</Piece>
</RectilinearGrid>
</VTKFile>

)input");  

    // Model current_model;
    // ModelPart& r_model_part = current_model.CreateModelPart("Main");

    VtrIO vtr_io(p_input);
    std::vector<Kratos::Internals::CartesianMeshColors> multi_block_mesh;

    vtr_io.Read(multi_block_mesh);

    KRATOS_CHECK_EQUAL(multi_block_mesh.size(), 1);

    auto& mesh = multi_block_mesh[0];

    std::vector<double> const& coordinates_x = mesh.GetNodalCoordinates(0);
    std::vector<double> const& coordinates_y = mesh.GetNodalCoordinates(1);
    std::vector<double> const& coordinates_z = mesh.GetNodalCoordinates(2);
    std::vector<double> reference_coordinates{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

    KRATOS_CHECK_EQUAL(coordinates_x.size(), 11);
    for(std::size_t i = 0 ; i < coordinates_x.size() ; i++){
        KRATOS_CHECK_EQUAL(coordinates_x[i], reference_coordinates[i]);
    }

    KRATOS_CHECK_EQUAL(coordinates_y.size(), 11);
    for(std::size_t i = 0 ; i < coordinates_y.size() ; i++){
        KRATOS_CHECK_EQUAL(coordinates_y[i], reference_coordinates[i]);
    }

    KRATOS_CHECK_EQUAL(coordinates_z.size(), 11);
    for(std::size_t i = 0 ; i < coordinates_z.size() ; i++){
        KRATOS_CHECK_EQUAL(coordinates_z[i], reference_coordinates[i]);
    }

}

KRATOS_TEST_CASE_IN_SUITE(ReadMeshAndCellDataFromAsciiVTR, KratosCoreFastSuite)
{
std::stringstream* p_input = new std::stringstream(R"input(
<VTKFile type="RectilinearGrid" version="0.1" byte_order="LittleEndian">
	<RectilinearGrid WholeExtent="0 3 0 4 0 5">
		<Piece Extent="0 3 0 4 0 5">
			<Coordinates>
				<DataArray type="Float64" Name="coordinates" NumberOfComponents="1" format="ascii">
0.000000 0.033333 0.066667 0.100000 
				</DataArray>
				<DataArray type="Float64" Name="coordinates" NumberOfComponents="1" format="ascii">
0.000000 0.020000 0.040000 0.060000 0.080000 
				</DataArray>
				<DataArray type="Float64" Name="coordinates" NumberOfComponents="1" format="ascii">
0.000000 0.020000 0.040000 0.060000 0.080000 0.100000 
				</DataArray>
			</Coordinates>
			<CellData Scalars="Temperature">
				<DataArray type="Int64" Name="Matrial" NumberOfComponents="1" format="ascii">
1 1 1 
1 1 1 
1 1 1 
1 1 1 
1 1 1 
1 1 1 
1 1 1 
1 1 1 
1 1 1 
1 1 1 
1 1 1 
1 1 1 
1 1 1 
1 1 1 
1 1 1 
1 1 1 
1 1 1 
1 1 1 
1 1 1 
1 1 1 
				</DataArray>
				<DataArray type="Float64" Name="VOF" NumberOfComponents="1" format="ascii">
0.753880 0.062963 0.753880 0.585248 0.284320 0.585248 0.641427 0.145841 0.641427 0.773637 
0.049395 0.773637 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 
				</DataArray>
				<DataArray type="Float64" Name="Temperature" NumberOfComponents="1" format="ascii">
303.905182 303.906128 303.905182 303.910126 303.911438 303.910126 303.910370 303.911987 303.910370 303.904816 
303.905396 303.904816 298.149994 298.149994 298.149994 298.149994 298.149994 298.149994 298.149994 298.149994 
298.149994 298.149994 298.149994 298.149994 298.149994 298.149994 298.149994 298.149994 298.149994 298.149994 
298.149994 298.149994 298.149994 298.149994 298.149994 298.149994 298.149994 298.149994 298.149994 298.149994 
298.149994 298.149994 298.149994 298.149994 298.149994 298.149994 298.149994 298.149994 298.149994 298.149994 
298.149994 298.149994 298.149994 298.149994 298.149994 298.149994 298.149994 298.149994 298.149994 298.149994 
				</DataArray>
				<DataArray type="Float64" Name="Density" NumberOfComponents="1" format="ascii">
1093.572388 1093.572388 1093.572388 1093.572388 1093.572388 1093.572388 1093.572388 1093.572388 1093.572388 1093.572388 
1093.572388 1093.572388 1.224500 1.224500 1.224500 1.224500 1.224500 1.224500 1.224500 1.224500 
1.224500 1.224500 1.224500 1.224500 1.224500 1.224500 1.224500 1.224500 1.224500 1.224500 
1.224500 1.224500 1.224500 1.224500 1.224500 1.224500 1.224500 1.224500 1.224500 1.224500 
1.224500 1.224500 1.224500 1.224500 1.224500 1.224500 1.224500 1.224500 1.224500 1.224500 
1.224500 1.224500 1.224500 1.224500 1.224500 1.224500 1.224500 1.224500 1.224500 1.224500 
				</DataArray>
			</CellData>
		</Piece>
	</RectilinearGrid>
</VTKFile>
)input");  

    // Model current_model;
    // ModelPart& r_model_part = current_model.CreateModelPart("Main");

    VtrIO vtr_io(p_input);
    std::vector<Kratos::Internals::CartesianMeshColors> multi_block_mesh;

    vtr_io.Read(multi_block_mesh);

    KRATOS_WATCH(multi_block_mesh.size());
    KRATOS_CHECK_EQUAL(multi_block_mesh.size(), 1);

    auto& mesh = multi_block_mesh[0];

    std::vector<double> const& coordinates_x = mesh.GetNodalCoordinates(0);
    std::vector<double> const& coordinates_y = mesh.GetNodalCoordinates(1);
    std::vector<double> const& coordinates_z = mesh.GetNodalCoordinates(2);
    std::vector<double> reference_x_coordinates{0.000000, 0.033333, 0.066667, 0.100000};
    std::vector<double> reference_y_coordinates{0.000000, 0.020000, 0.040000, 0.060000, 0.080000};
    std::vector<double> reference_z_coordinates{0.000000, 0.020000, 0.040000, 0.060000, 0.080000, 0.100000};

    KRATOS_CHECK_EQUAL(coordinates_x.size(), 4);
    for(std::size_t i = 0 ; i < coordinates_x.size() ; i++){
        KRATOS_CHECK_EQUAL(coordinates_x[i], reference_x_coordinates[i]);
    }

    KRATOS_CHECK_EQUAL(coordinates_y.size(), 5);
    for(std::size_t i = 0 ; i < coordinates_y.size() ; i++){
        KRATOS_CHECK_EQUAL(coordinates_y[i], reference_y_coordinates[i]);
    }

    KRATOS_CHECK_EQUAL(coordinates_z.size(), 6);
    for(std::size_t i = 0 ; i < coordinates_z.size() ; i++){
        KRATOS_CHECK_EQUAL(coordinates_z[i], reference_z_coordinates[i]);
    }

    KRATOS_CHECK(mesh.HasElementalData("Matrial"));
    auto& material = mesh.GetElementalData("Matrial");
    KRATOS_CHECK_EQUAL(material.size(), 60);
    for(auto& i_material : material)
        KRATOS_CHECK_EQUAL(i_material, 1);

    KRATOS_CHECK(mesh.HasElementalData("VOF"));
    auto& vof = mesh.GetElementalData("VOF");
    std::vector<double> vof_reference{
        0.753880, 0.062963, 0.753880, 0.585248, 0.284320, 0.585248, 0.641427, 0.145841, 0.641427, 0.773637, 
        0.049395, 0.773637, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
        0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
        0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
        0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
        0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
    KRATOS_CHECK_EQUAL(vof.size(), 60);
    for(std::size_t i = 0 ; i < vof.size() ; i++)
        KRATOS_CHECK_EQUAL(vof[i], vof_reference[i]);

    KRATOS_CHECK(mesh.HasElementalData("Temperature"));
    auto& temperature = mesh.GetElementalData("Temperature");
    std::vector<double> temperature_reference{
        303.905182, 303.906128, 303.905182, 303.910126, 303.911438, 303.910126, 303.910370, 303.911987, 303.910370, 303.904816, 
        303.905396, 303.904816, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 
        298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 
        298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 
        298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 
        298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994, 298.149994};
    
    KRATOS_CHECK_EQUAL(temperature.size(), 60);
    for(std::size_t i = 0 ; i < temperature.size() ; i++)
        KRATOS_CHECK_EQUAL(temperature[i], temperature_reference[i]);

    KRATOS_CHECK(mesh.HasElementalData("Density"));
    auto& density = mesh.GetElementalData("Density");
    std::vector<double> density_reference{
        1093.572388, 1093.572388, 1093.572388, 1093.572388, 1093.572388, 1093.572388, 1093.572388, 1093.572388, 1093.572388, 1093.572388, 
        1093.572388, 1093.572388, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 
        1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 
        1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 
        1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 
        1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500, 1.224500};
    
    KRATOS_CHECK_EQUAL(density.size(), 60);
    for(std::size_t i = 0 ; i < density.size() ; i++)
        KRATOS_CHECK_EQUAL(density[i], density_reference[i]);

}

}
}  // namespace Kratos.
