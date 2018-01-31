//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//
//

// Project includes
#include "testing/testing.h"
#include "includes/serializer.h"

namespace Kratos {
    namespace Testing {

        KRATOS_TEST_CASE_IN_SUITE(SerializerMap, KratosCoreFastSuite)
        {
            Serializer serializer;

            const std::string tag_string("TestString");

            std::map<std::string, int> first;

            std::map<std::string, int> loaded;

            first["www"]=10;
            first["eee"]=30;
            first["rrr"]=50;
            first["ttt"]=70;

            serializer.save(tag_string, first);

            serializer.load(tag_string, loaded);

            KRATOS_CHECK_EQUAL(first,loaded);


        }
    } // namespace Testing
}  // namespace Kratos.
