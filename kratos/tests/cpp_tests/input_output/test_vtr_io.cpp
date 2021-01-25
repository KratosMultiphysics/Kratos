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

}
}  // namespace Kratos.
