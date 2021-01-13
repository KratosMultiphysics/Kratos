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

}
}  // namespace Kratos.
