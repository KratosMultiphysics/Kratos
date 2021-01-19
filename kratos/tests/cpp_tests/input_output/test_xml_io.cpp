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
#include "input_output/xml_io.h"



namespace Kratos {
namespace Testing {


KRATOS_TEST_CASE_IN_SUITE(ReadXMLDocumentBlock, KratosCoreFastSuite)
{
    std::stringstream* p_input = new std::stringstream(R"input(
    < VTKFile type="RectilinearGrid"      version="0.1 "          byte_order="LittleEndian">
        </VTKFile>

    )input");  

    XmlIO xml_io(p_input);

    xml_io.SetBlockAction("VTKFile", [](XmlIO& TheIO){
        auto& current_block_info = TheIO.GetCurrentBlockInfo();
        KRATOS_CHECK_EQUAL(current_block_info.Name(), "VTKFile");
        KRATOS_CHECK(current_block_info.HasAttribute("type"));
        KRATOS_CHECK(current_block_info.HasAttribute("version"));
        KRATOS_CHECK(current_block_info.HasAttribute("byte_order"));
        KRATOS_CHECK(current_block_info.GetAttribute("type")=="RectilinearGrid");
        KRATOS_CHECK(current_block_info.GetAttribute("version")=="0.1 ");
        KRATOS_CHECK(current_block_info.GetAttribute("byte_order")=="LittleEndian");
    });

    xml_io.Read();
}

KRATOS_TEST_CASE_IN_SUITE(ReadXMLElementBlock, KratosCoreFastSuite)
{
    std::stringstream* p_input = new std::stringstream(R"input(
    < VTKFile type="RectilinearGrid"      version="0.1 "          byte_order="LittleEndian">
        <RectilinearGrid WholeExtent="0 10 0 10 0 10">
        </RectilinearGrid>
    </VTKFile>

    )input");  

    XmlIO xml_io(p_input);

    xml_io.SetBlockAction("VTKFile", [](XmlIO& TheIO){
        auto& current_block_info = TheIO.GetCurrentBlockInfo();
        KRATOS_CHECK_EQUAL(current_block_info.Name(), "VTKFile");
        KRATOS_CHECK(current_block_info.HasAttribute("type"));
        KRATOS_CHECK(current_block_info.HasAttribute("version"));
        KRATOS_CHECK(current_block_info.HasAttribute("byte_order"));
        TheIO.Read();
    });

     xml_io.SetBlockAction("RectilinearGrid", [](XmlIO& TheIO){
        auto& current_block_info = TheIO.GetCurrentBlockInfo();
        KRATOS_CHECK_EQUAL(current_block_info.Name(), "RectilinearGrid");
        KRATOS_CHECK(current_block_info.HasAttribute("WholeExtent"));
        KRATOS_CHECK_IS_FALSE(current_block_info.HasAttribute("type"));
        KRATOS_CHECK_IS_FALSE(current_block_info.HasAttribute("version"));
        KRATOS_CHECK_IS_FALSE(current_block_info.HasAttribute("byte_order"));
    });

   xml_io.Read();
}
// #include "zlib.h"
// size_t UncompressBuffer(char const* compressedData,
//   size_t compressedSize, char* uncompressedData, size_t uncompressedSize)
// {
//   uLongf us = static_cast<uLongf>(uncompressedSize);
//   Bytef* ud = reinterpret_cast<Bytef*>(uncompressedData);
//   const Bytef* cd = reinterpret_cast<const Bytef*>(compressedData);
//   uLong cs = static_cast<uLong>(compressedSize);

//   // Call zlib's uncompress function.
//   if (uncompress(ud, &us, cd, cs) != Z_OK)
//   {
//     KRATOS_ERROR << "Zlib error while uncompressing data.";
//     return 0;
//   }

//   // Make sure the output size matched that expected.
//   if (us != static_cast<uLongf>(uncompressedSize))
//   {
//     KRATOS_ERROR << "Decompression produced incorrect size.\n" <<
//                   "Expected "
//       << uncompressedSize << " and got " << us;
//     return 0;
//   }

//   return static_cast<size_t>(us);
// }

KRATOS_TEST_CASE_IN_SUITE(ReadXMLElementBlockContent, KratosCoreFastSuite)
{
    std::stringstream* p_input = new std::stringstream(R"input(
    <Coordinates> 
    <DataArray type="Float64" Name="coordinates" NumberOfComponents="1" format="ascii">
    0 1 2 3 4 5 6 7 8 9 10 </DataArray> 
    <DataArray type="Float64" Name="coordinates" NumberOfComponents="1" format="ascii">
    0 1 2 3 4 5 6 7 8 9 10 </DataArray> 
    <DataArray type="Float64" Name="coordinates" NumberOfComponents="1" format="ascii">
    0 1 2 3 4 5 6 7 8 9 10 </DataArray> 
    </Coordinates> 

    )input");  

    XmlIO xml_io(p_input);

    xml_io.SetBlockAction("Coordinates", [](XmlIO& TheIO){
        std::vector<double> x_coordinates(11, -1.00);
        TheIO.ReadBlock("DataArray", [&x_coordinates](XmlIO& TheIO){ // Here I am forcing to have DataArray block in the input with the given action.
            TheIO.ReadBlockContent(x_coordinates);
        });
        for(std::size_t i = 0 ; i < x_coordinates.size() ; i++)
            KRATOS_CHECK_EQUAL(x_coordinates[i], i);

        std::vector<double> y_coordinates(11, -1.00);
        TheIO.ReadBlock("DataArray", [&y_coordinates](XmlIO& TheIO){ // Here I am forcing to have DataArray block in the input with the given action.
            TheIO.ReadBlockContent(y_coordinates);
        });
        for(std::size_t i = 0 ; i < y_coordinates.size() ; i++)
            KRATOS_CHECK_EQUAL(y_coordinates[i], i);

        std::vector<double> z_coordinates(11, -1.00);
        TheIO.ReadBlock("DataArray", [&z_coordinates](XmlIO& TheIO){ // Here I am forcing to have DataArray block in the input with the given action.
            TheIO.ReadBlockContent(z_coordinates);
        });
        for(std::size_t i = 0 ; i < z_coordinates.size() ; i++)
            KRATOS_CHECK_EQUAL(z_coordinates[i], i);
    });

   xml_io.Read();
//     int size = 8*8*11*32;
//     std::string input(R"input(AQAAAACAAAAACwAAFQEAAA==eJzFjzEOAVEYhM/DlnqNUlRcQEKhsPeQiAgq0ehEZCUKiVZD5wJ6N1AiL9mV32T+SRbJFl92s997OzP366B2L4jKdCKd5zPHvHXo0QVOx+EH55STeaI7w7lvfchndzLnedtfOfRs/7o9LgyVr/qp/mof27+5zQpD5at+qr/ax/ZXR/PCUPmqn+qv9rH9z8XyzYPwT49gPjrWzzrrmcP7bH9rvXrTJHzj8d0D89F53uuPDj3b328kuYkNeb11IT+Gb9Z53vZXDj3bX7rsamWH6AcXkTMREPKV87ztr5y9XyYukPT2udmm5PVbAPPRZZ7dQ58Qh/9n+zv1wwfdFO87811yRt3NCPnKed7rjw492/8CcpdexA==)input");
//     char output[size]; 
//     auto is_ok = UncompressBuffer(input.c_str(), input.size(), output, size);
//     KRATOS_CHECK(is_ok);
//     std::cout << std::endl << "output : " << output << std::endl; 
}

}
}  // namespace Kratos.
