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
#include "includes/file_serializer.h"
#include "includes/stream_serializer.h"
#include "geometries/point.h"

namespace Kratos {
    namespace Testing {

        /*********************************************************************/
        /* Helper Functions */
        /*********************************************************************/
        template<typename TObjectType>
        void SaveAndLoadObjects(const TObjectType& rObjectToBeSaved, TObjectType& rObjectToBeLoaded)
        {
            StreamSerializer serializer;

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

            KRATOS_CHECK_EQUAL(rObjectToBeLoaded.size(), rObjectToBeSaved.size());

            for (std::size_t i=0; i< rObjectToBeSaved.size(); ++i)
                KRATOS_CHECK_EQUAL(rObjectToBeLoaded[i], rObjectToBeSaved[i]);
        }

        template<typename TObjectType>
        void TestObjectSerializationComponentwise2D(const TObjectType& rObjectToBeSaved, TObjectType& rObjectToBeLoaded)
        {
            SaveAndLoadObjects(rObjectToBeSaved, rObjectToBeLoaded);

            KRATOS_CHECK_EQUAL(rObjectToBeLoaded.size1(), rObjectToBeSaved.size1());
            KRATOS_CHECK_EQUAL(rObjectToBeLoaded.size2(), rObjectToBeSaved.size2());

            for (std::size_t i=0; i<rObjectToBeSaved.size1(); ++i) {
                for (std::size_t j=0; j<rObjectToBeSaved.size2(); ++j) {
                    KRATOS_CHECK_EQUAL(rObjectToBeLoaded(i,j), rObjectToBeSaved(i,j));
                }
            }
        }

        template<typename TObjectType>
        void FillVectorWithValues(TObjectType& rObject)
        {
            for (std::size_t i = 0; i < rObject.size(); ++i)
                 rObject[i] = i*i*0.2 + 5.333;
        }

        template<typename TObjectType>
        void FillMatrixWithValues(TObjectType& rObject)
        {
            for (std::size_t i=0; i<rObject.size1(); ++i) {
                for (std::size_t j=0; j<rObject.size2(); ++j) {
                    rObject(i,j) = i*j - 51.21;
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

            FillMatrixWithValues(object_to_be_saved);

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
            StreamSerializer serializer;

            const std::string tag_string("TestString");
            const std::string tag_string_2("TestString2");

            Point::Pointer p_point = Kratos::make_shared<Point>(-0.25, 3.5, 8.55);
            Kratos::shared_ptr<array_1d<double,3>> p_array = p_point;

            Point::Pointer p_loaded_point;
            Kratos::shared_ptr<array_1d<double,3>> p_loaded_array;

            serializer.save(tag_string, p_array);
            serializer.save(tag_string_2, p_point);

            serializer.load(tag_string, p_loaded_array);
            serializer.load(tag_string_2, p_loaded_point);

            KRATOS_CHECK_EQUAL(*p_point, *p_loaded_point);
            for (std::size_t i=0; i<(*p_array).size(); ++i)
                KRATOS_CHECK_EQUAL((*p_loaded_array)[i], (*p_array)[i]);
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerStdArray, KratosCoreFastSuite)
        {
            using Array5Dbl = std::array<double,5>;
            using Array7Int = std::array<int,7>;

            Array5Dbl object_to_be_saved_1;
            Array5Dbl object_to_be_loaded_1;
            Array7Int object_to_be_saved_2;
            Array7Int object_to_be_loaded_2;

            FillVectorWithValues(object_to_be_saved_1);
            FillVectorWithValues(object_to_be_saved_2);

            TestObjectSerializationComponentwise1D(object_to_be_saved_1, object_to_be_loaded_1);
            TestObjectSerializationComponentwise1D(object_to_be_saved_2, object_to_be_loaded_2);
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerStdVector, KratosCoreFastSuite)
        {
            std::vector<double> object_to_be_saved(5);
            std::vector<double> object_to_be_loaded(3); // initialized smaller to check if resizing works

            FillVectorWithValues(object_to_be_saved);

            TestObjectSerializationComponentwise1D(object_to_be_saved, object_to_be_loaded);
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerUblasVector, KratosCoreFastSuite)
        {
            Vector object_to_be_saved(5);
            Vector object_to_be_loaded(3); // initialized smaller to check if resizing works

            FillVectorWithValues(object_to_be_saved);

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

            FillVectorWithValues(object_to_be_saved);

            TestObjectSerializationComponentwise1D(object_to_be_saved, object_to_be_loaded);
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerPair, KratosCoreFastSuite)
        {
            std::pair<int, double> object_to_be_saved(42, 0.123);
            std::pair<int, double> object_to_be_loaded(42, 0.123);

            TestObjectSerialization(object_to_be_saved, object_to_be_loaded);
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerBoundedVector, KratosCoreFastSuite)
        {
            using Vector5Double = BoundedVector<double,5>;
            using Vector6Int = BoundedVector<int,6>;

            Vector5Double object_to_be_saved_1;
            Vector5Double object_to_be_loaded_1;
            Vector6Int object_to_be_saved_2;
            Vector6Int object_to_be_loaded_2;

            FillVectorWithValues(object_to_be_saved_1);
            FillVectorWithValues(object_to_be_saved_2);

            TestObjectSerializationComponentwise1D(object_to_be_saved_1, object_to_be_loaded_1);
            TestObjectSerializationComponentwise1D(object_to_be_saved_2, object_to_be_loaded_2);
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerBoundedMatrix, KratosCoreFastSuite)
        {
            using Matrix53Double = BoundedMatrix<double,5,3>;
            using Matrix26Int = BoundedMatrix<int,2,6>;

            Matrix53Double object_to_be_saved_1;
            Matrix53Double object_to_be_loaded_1;
            Matrix26Int object_to_be_saved_2;
            Matrix26Int object_to_be_loaded_2;

            FillMatrixWithValues(object_to_be_saved_1);
            FillMatrixWithValues(object_to_be_saved_2);

            TestObjectSerializationComponentwise2D(object_to_be_saved_1, object_to_be_loaded_1);
            TestObjectSerializationComponentwise2D(object_to_be_saved_2, object_to_be_loaded_2);
        }

    } // namespace Testing
}  // namespace Kratos.
