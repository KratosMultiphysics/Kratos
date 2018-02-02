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
#include "includes/node.h"

namespace Kratos {
    namespace Testing {

        /*********************************************************************/
        /* Helper Functions */
        /*********************************************************************/
        template<typename TObjectType>
        void SaveAndLoadObjects(const TObjectType& rObjectToBeSaved, TObjectType& rObjectToBeLoaded)
        {
            Serializer serializer;

            const std::string tag_string("TestString");

            serializer.save(tag_string, rObjectToBeSaved);
            serializer.load(tag_string, rObjectToBeLoaded);
        }

        template<typename TObjectType>
        void TestObjectSerialization(const TObjectType& rObjectToBeSaved, TObjectType& rObjectToBeLoaded)
        {
            SaveAndLoadObjects(rObjectToBeSaved, rObjectToBeLoaded);
            KRATOS_CHECK_EQUAL(rObjectToBeLoaded, rObjectToBeSaved);
        }

        template<typename TObjectType>
        void TestObjectSerializationComponentwise1D(const TObjectType& rObjectToBeSaved, TObjectType& rObjectToBeLoaded)
        {
            SaveAndLoadObjects(rObjectToBeSaved, rObjectToBeLoaded);

            KRATOS_ERROR_IF_NOT(rObjectToBeSaved.size() == rObjectToBeLoaded.size()) 
                << "Object sizes are not consistent" << std::endl;

            for (std::size_t i=0; i< rObjectToBeSaved.size(); ++i)
                KRATOS_CHECK_EQUAL(rObjectToBeLoaded[i], rObjectToBeSaved[i]);
        }

        template<typename TObjectType>
        void TestObjectSerializationComponentwise2D(const TObjectType& rObjectToBeSaved, TObjectType& rObjectToBeLoaded)
        {
            SaveAndLoadObjects(rObjectToBeSaved, rObjectToBeLoaded);   

            KRATOS_ERROR_IF_NOT(rObjectToBeSaved.size1() == rObjectToBeLoaded.size1() &&
                                rObjectToBeSaved.size2() == rObjectToBeLoaded.size2()) 
                << "Object sizes are not consistent" << std::endl;

            for (std::size_t i=0; i<rObjectToBeSaved.size1(); ++i) {
                for (std::size_t j=0; j<rObjectToBeSaved.size2(); ++j) {
                    KRATOS_CHECK_EQUAL(rObjectToBeLoaded(i,j), rObjectToBeSaved(i,j));
                }
            }
        }

        /*********************************************************************/
        /* Testing the Datatypes that for which the 
           "KRATOS_SERIALIZATION_DIRECT_LOAD" macro is used */
        /*********************************************************************/
        KRATOS_TEST_CASE_IN_SUITE(SerializerBool, KratosCoreFastSuite)
        {
            bool object_to_be_saved = true;
            bool object_to_be_loaded;
            
            TestObjectSerialization(object_to_be_saved, object_to_be_loaded);
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerInt, KratosCoreFastSuite)
        {
            int object_to_be_saved = -105;
            int object_to_be_loaded;
            
            TestObjectSerialization(object_to_be_saved, object_to_be_loaded);
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerLong, KratosCoreFastSuite)
        {
            long object_to_be_saved = -1598456605;
            long object_to_be_loaded;
            
            TestObjectSerialization(object_to_be_saved, object_to_be_loaded);
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerDouble, KratosCoreFastSuite)
        {
            double object_to_be_saved = -159845.6605;
            double object_to_be_loaded;
            
            TestObjectSerialization(object_to_be_saved, object_to_be_loaded);
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerUnsignedLong, KratosCoreFastSuite)
        {
            unsigned long object_to_be_saved = 1598456605;
            unsigned long object_to_be_loaded;
            
            TestObjectSerialization(object_to_be_saved, object_to_be_loaded);
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerUnsignedInt, KratosCoreFastSuite)
        {
            unsigned int object_to_be_saved = 15905;
            unsigned int object_to_be_loaded;
            
            TestObjectSerialization(object_to_be_saved, object_to_be_loaded);
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerStdString, KratosCoreFastSuite)
        {
            std::string object_to_be_saved = "MyStringToBeSerialized";
            std::string object_to_be_loaded;
            
            TestObjectSerialization(object_to_be_saved, object_to_be_loaded);
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerMatrix, KratosCoreFastSuite)
        {
            Matrix object_to_be_saved(5,3);
            Matrix object_to_be_loaded;

            for (std::size_t i=0; i<object_to_be_saved.size1(); ++i) {
                for (std::size_t j=0; j<object_to_be_saved.size2(); ++j) {
                    object_to_be_saved(i,j) = i*j - 51.21;
                }
            }
            
            TestObjectSerializationComponentwise2D(object_to_be_saved, object_to_be_loaded);
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerMatrix3, KratosCoreFastSuite)
        {
            using Matrix3 = bounded_matrix<double,3,3>; // see serializer.h
            Matrix3 object_to_be_saved;
            Matrix3 object_to_be_loaded;

            for (std::size_t i=0; i<object_to_be_saved.size1(); ++i) {
                for (std::size_t j=0; j<object_to_be_saved.size2(); ++j) {
                    object_to_be_saved(i,j) = i*j - 51.21;
                }
            }
            
            TestObjectSerializationComponentwise2D(object_to_be_saved, object_to_be_loaded);
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerLongLong, KratosCoreFastSuite)
        {
            long long object_to_be_saved = -1598456546843565605;
            long long object_to_be_loaded;
            
            TestObjectSerialization(object_to_be_saved, object_to_be_loaded);
        }

        /*********************************************************************/
        /* Testing the Datatypes that have a specific save/load implementation */
        /*********************************************************************/

        KRATOS_TEST_CASE_IN_SUITE(SerializerKratosSharedPtr, KratosCoreFastSuite)
        {
            Serializer serializer;

            const std::string tag_string("TestString");
            const std::string tag_string_2("TestString2");

            Node<3>::Pointer p_node = Kratos::make_shared<Node<3>>();
            Point::Pointer p_point = p_node;

            Node<3>::Pointer p_loaded_node;
            Point::Pointer p_loaded_point;

            serializer.save(tag_string, p_point);
            serializer.save(tag_string_2, p_node);

            serializer.load(tag_string, p_loaded_point);
            serializer.load(tag_string_2, p_loaded_node);

            KRATOS_CHECK_EQUAL(p_loaded_node, p_loaded_point);
            KRATOS_CHECK_EQUAL(*p_node, *p_loaded_node);
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerVariable, KratosCoreFastSuite)
        {
            Serializer serializer;

            const std::string tag_string("TestString");

            Variable<double> object_to_be_saved();
            Variable<double> object_to_be_loaded();

            // serializer.save(tag_string, object_to_be_saved);
            // serializer.load(tag_string, object_to_be_loaded);
            
            // KRATOS_CHECK_EQUAL(object_to_be_loaded, object_to_be_saved);
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerStdVector, KratosCoreFastSuite)
        {
            std::vector<double> object_to_be_saved{1.235, 4.66, 1.23, -88.2, -66.1};
            std::vector<double> object_to_be_loaded(3); // initialized smaller to check if resizing works

            TestObjectSerializationComponentwise1D(object_to_be_saved, object_to_be_loaded);
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerUblasVector, KratosCoreFastSuite)
        {
            constexpr std::size_t size_vector = 5;
            Vector object_to_be_saved(size_vector);
            Vector object_to_be_loaded(3); // initialized smaller to check if resizing works

            // initializer list initialization is not available for "Vector", therefore doing it manually
            for (std::size_t i = 0; i < size_vector; ++i)
                 object_to_be_saved[i] = i*i*0.2;

            TestObjectSerializationComponentwise1D(object_to_be_saved, object_to_be_loaded);
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerMap, KratosCoreFastSuite)
        {
            std::map <std::string, double> object_to_be_saved { 
                { "42", -30.556 }, 
                { "3",   10.258 } 
            };
            std::map<std::string, double> object_to_be_loaded;

            TestObjectSerialization(object_to_be_saved, object_to_be_loaded);
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerUnorderedMap, KratosCoreFastSuite)
        {
            std::unordered_map <std::string, double> object_to_be_saved { 
                { "42", -30.556 }, 
                { "3",   10.258 } 
            };
            std::unordered_map<std::string, double> object_to_be_loaded;

            TestObjectSerialization(object_to_be_saved, object_to_be_loaded);
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerArray1D, KratosCoreFastSuite)
        {
            constexpr std::size_t size_object = 5;
            array_1d<double, size_object> object_to_be_saved;
            array_1d<double, size_object> object_to_be_loaded;

            for (std::size_t i = 0; i < size_object; ++i)
                 object_to_be_saved[i] = i*i*0.2;

            TestObjectSerializationComponentwise1D(object_to_be_saved, object_to_be_loaded);
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerPair, KratosCoreFastSuite)
        {
            std::pair<int, double> object_to_be_saved(42, 0.123);
            std::pair<int, double> object_to_be_loaded(42, 0.123);

            TestObjectSerialization(object_to_be_saved, object_to_be_loaded);
        }

    } // namespace Testing
}  // namespace Kratos.
