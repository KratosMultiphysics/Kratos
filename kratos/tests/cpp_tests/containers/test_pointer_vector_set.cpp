//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Pooyan Dadvand
//
//

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "containers/model.h"
#include "containers/pointer_vector_set.h"
#include "testing/testing.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(PointerVectorSetCBeginAndCEnd, KratosCoreFastSuite)
{
    PointerVectorSet<const Element, IndexedObject> test_container;
    auto p_element_1 = Kratos::make_intrusive<Element>(1);
    auto p_element_2 = Kratos::make_intrusive<Element>(2);
    auto p_element_3 = Kratos::make_intrusive<Element>(3);
    test_container.insert(p_element_1);
    test_container.insert(p_element_2);
    test_container.insert(p_element_3);

    KRATOS_EXPECT_EQ(test_container.cbegin()->Id(), 1);
    KRATOS_EXPECT_EQ((test_container.cend()-1)->Id(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(PointerVectorSetInsert1, KratosCoreFastSuite)
{
    PointerVectorSet<const Element> test_container;
    auto p_element_1 = Kratos::make_intrusive<Element>(1);
    auto p_element_2 = Kratos::make_intrusive<Element>(2);
    auto p_element_3 = Kratos::make_intrusive<Element>(3);
    auto p_element_4 = Kratos::make_intrusive<Element>(4);
    auto p_element_5 = Kratos::make_intrusive<Element>(5);
    auto p_element_6 = Kratos::make_intrusive<Element>(6);
    auto p_element_3_ptr_copy = Kratos::intrusive_ptr<Element>(p_element_3);
    KRATOS_EXPECT_EQ(&*test_container.insert(p_element_3), &*p_element_3);
    KRATOS_EXPECT_EQ(&*test_container.insert(p_element_2), &*p_element_2);
    KRATOS_EXPECT_EQ(&*test_container.insert(p_element_1), &*p_element_1);
    KRATOS_EXPECT_EQ(&*test_container.insert(p_element_5), &*p_element_5);
    KRATOS_EXPECT_EQ(&*test_container.insert(p_element_6), &*p_element_6);
    KRATOS_EXPECT_EQ(&*test_container.insert(p_element_4), &*p_element_4);
    KRATOS_EXPECT_EQ(&*test_container.insert(p_element_3_ptr_copy), &*p_element_3);

    KRATOS_EXPECT_EQ(test_container.size(), 6);

    auto itr = test_container.begin();
    for (; itr != test_container.end() - 1; ++itr) {
        KRATOS_EXPECT_TRUE(&*(itr) < &*(itr + 1));
    }
}

KRATOS_TEST_CASE_IN_SUITE(PointerVectorSetInsert2, KratosCoreFastSuite)
{
    PointerVectorSet<const Element, IndexedObject> test_container;
    auto p_element_1 = Kratos::make_intrusive<Element>(1);
    auto p_element_2 = Kratos::make_intrusive<Element>(2);
    auto p_element_3 = Kratos::make_intrusive<Element>(3);
    auto p_element_4 = Kratos::make_intrusive<Element>(4);
    auto p_element_5 = Kratos::make_intrusive<Element>(5);
    auto p_element_6 = Kratos::make_intrusive<Element>(6);
    auto p_element_3_copy = Kratos::make_intrusive<Element>(3);
    KRATOS_EXPECT_EQ(&*test_container.insert(p_element_3), &*p_element_3);
    KRATOS_EXPECT_EQ(&*test_container.insert(p_element_2), &*p_element_2);
    KRATOS_EXPECT_EQ(&*test_container.insert(p_element_1), &*p_element_1);
    KRATOS_EXPECT_EQ(&*test_container.insert(p_element_5), &*p_element_5);
    KRATOS_EXPECT_EQ(&*test_container.insert(p_element_6), &*p_element_6);
    KRATOS_EXPECT_EQ(&*test_container.insert(p_element_4), &*p_element_4);
    KRATOS_EXPECT_EQ(&*test_container.insert(p_element_3_copy), &*p_element_3);

    KRATOS_EXPECT_EQ(test_container.size(), 6);

    auto itr = test_container.begin();
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_1);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_2);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_3);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_4);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_5);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_6);
}

KRATOS_TEST_CASE_IN_SUITE(PointerVectorSetInsert3, KratosCoreFastSuite)
{
    PointerVectorSet<const Element, IndexedObject> test_container;
    auto p_element_1 = Kratos::make_intrusive<Element>(1);
    auto p_element_2 = Kratos::make_intrusive<Element>(2);
    auto p_element_3 = Kratos::make_intrusive<Element>(3);
    auto p_element_4 = Kratos::make_intrusive<Element>(4);
    auto p_element_5 = Kratos::make_intrusive<Element>(5);
    auto p_element_6 = Kratos::make_intrusive<Element>(6);
    auto p_element_3_ptr_copy = Kratos::intrusive_ptr<Element>(p_element_3);
    KRATOS_EXPECT_EQ(&*test_container.insert(test_container.end(), p_element_3), &*p_element_3);
    KRATOS_EXPECT_EQ(&*test_container.insert(test_container.end(), p_element_6), &*p_element_6);
    KRATOS_EXPECT_EQ(&*test_container.insert(test_container.end(), p_element_2), &*p_element_2);
    KRATOS_EXPECT_EQ(&*test_container.insert(test_container.begin(), p_element_1), &*p_element_1);
    KRATOS_EXPECT_EQ(&*test_container.insert(test_container.begin(), p_element_5), &*p_element_5);
    KRATOS_EXPECT_EQ(&*test_container.insert(test_container.begin() + 3, p_element_4), &*p_element_4);
    KRATOS_EXPECT_EQ(&*test_container.insert(test_container.begin() + 3, p_element_1), &*p_element_1);
    KRATOS_EXPECT_EQ(&*test_container.insert(test_container.begin() + 3, p_element_3_ptr_copy), &*p_element_3);

    KRATOS_EXPECT_EQ(test_container.size(), 6);

    auto itr = test_container.begin();
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_1);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_2);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_3);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_4);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_5);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_6);
}

KRATOS_TEST_CASE_IN_SUITE(PointerVectorSetInsert4, KratosCoreFastSuite)
{
    PointerVectorSet<const Element, IndexedObject> test_container;
    auto p_element_1 = Kratos::make_intrusive<Element>(1);
    auto p_element_2 = Kratos::make_intrusive<Element>(2);
    auto p_element_3 = Kratos::make_intrusive<Element>(3);
    auto p_element_4 = Kratos::make_intrusive<Element>(4);
    auto p_element_5 = Kratos::make_intrusive<Element>(5);
    auto p_element_6 = Kratos::make_intrusive<Element>(6);
    auto p_element_3_copy = Kratos::make_intrusive<Element>(3);

    std::vector<Element::Pointer> tmp;
    tmp.push_back(p_element_2);
    tmp.push_back(p_element_1);
    tmp.push_back(p_element_4);
    tmp.push_back(p_element_3);
    tmp.push_back(p_element_4);
    tmp.push_back(p_element_6);
    tmp.push_back(p_element_1);
    tmp.push_back(p_element_5);
    tmp.push_back(p_element_3_copy);
    test_container.insert(tmp.begin(), tmp.end());

    KRATOS_EXPECT_EQ(test_container.size(), 6);

    auto itr = test_container.begin();
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_1);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_2);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_3);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_4);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_5);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_6);
}

KRATOS_TEST_CASE_IN_SUITE(PointerVectorSetInsert5, KratosCoreFastSuite)
{
    PointerVectorSet<const Element> test_container;
    auto p_element_1 = Kratos::make_intrusive<Element>(1);
    auto p_element_2 = Kratos::make_intrusive<Element>(2);
    auto p_element_3 = Kratos::make_intrusive<Element>(3);
    auto p_element_4 = Kratos::make_intrusive<Element>(4);
    auto p_element_5 = Kratos::make_intrusive<Element>(5);
    auto p_element_6 = Kratos::make_intrusive<Element>(6);
    auto p_element_7 = Kratos::make_intrusive<Element>(7);
    auto p_element_8 = Kratos::make_intrusive<Element>(8);
    auto p_element_3_ptr_copy = Kratos::intrusive_ptr<Element>(p_element_3);

    std::vector<Element::Pointer> tmp;
    tmp.push_back(p_element_2);
    tmp.push_back(p_element_1);
    tmp.push_back(p_element_5);
    tmp.push_back(p_element_8);
    test_container.insert(tmp.begin(), tmp.end());

    KRATOS_EXPECT_EQ(test_container.size(), 4);
    auto itr = test_container.begin();
    for (; itr != test_container.end() - 1; ++itr) {
        KRATOS_EXPECT_TRUE(&*(itr) < &*(itr + 1));
    }

    tmp.clear();
    tmp.push_back(p_element_3);
    tmp.push_back(p_element_4);
    tmp.push_back(p_element_5);
    tmp.push_back(p_element_3);
    tmp.push_back(p_element_3_ptr_copy);
    tmp.push_back(p_element_1);
    tmp.push_back(p_element_7);
    tmp.push_back(p_element_6);
    test_container.insert(tmp.begin(), tmp.end());

    KRATOS_EXPECT_EQ(test_container.size(), 8);
    itr = test_container.begin();
    for (; itr != test_container.end() - 1; ++itr) {
        KRATOS_EXPECT_TRUE(&*(itr) < &*(itr + 1));
    }
}

KRATOS_TEST_CASE_IN_SUITE(PointerVectorSetInsert6, KratosCoreFastSuite)
{
    PointerVectorSet<const Element, IndexedObject> test_container;
    auto p_element_1 = Kratos::make_intrusive<Element>(1);
    auto p_element_2 = Kratos::make_intrusive<Element>(2);
    auto p_element_3 = Kratos::make_intrusive<Element>(3);
    auto p_element_4 = Kratos::make_intrusive<Element>(4);
    auto p_element_5 = Kratos::make_intrusive<Element>(5);
    auto p_element_6 = Kratos::make_intrusive<Element>(6);
    auto p_element_7 = Kratos::make_intrusive<Element>(7);
    auto p_element_8 = Kratos::make_intrusive<Element>(8);
    auto p_element_3_copy = Kratos::make_intrusive<Element>(3);

    std::vector<Element::Pointer> tmp;
    tmp.push_back(p_element_2);
    tmp.push_back(p_element_1);
    tmp.push_back(p_element_5);
    tmp.push_back(p_element_8);
    test_container.insert(tmp.begin(), tmp.end());

    KRATOS_EXPECT_EQ(test_container.size(), 4);
    auto itr = test_container.begin();
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_1);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_2);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_5);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_8);

    tmp.clear();
    tmp.push_back(p_element_3);
    tmp.push_back(p_element_4);
    tmp.push_back(p_element_5);
    tmp.push_back(p_element_3);
    tmp.push_back(p_element_3_copy);
    tmp.push_back(p_element_1);
    tmp.push_back(p_element_7);
    tmp.push_back(p_element_6);
    test_container.insert(tmp.begin(), tmp.end());

    KRATOS_EXPECT_EQ(test_container.size(), 8);
    itr = test_container.begin();
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_1);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_2);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_3);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_4);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_5);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_6);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_7);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_8);
}

KRATOS_TEST_CASE_IN_SUITE(PointerVectorSetInsert7, KratosCoreFastSuite)
{
    PointerVectorSet<const Element, IndexedObject> test_container_1, test_container_2, test_container_3;
    auto p_element_1 = Kratos::make_intrusive<Element>(1);
    auto p_element_2 = Kratos::make_intrusive<Element>(2);
    auto p_element_3 = Kratos::make_intrusive<Element>(3);
    auto p_element_4 = Kratos::make_intrusive<Element>(4);
    auto p_element_5 = Kratos::make_intrusive<Element>(5);
    auto p_element_6 = Kratos::make_intrusive<Element>(6);
    auto p_element_7 = Kratos::make_intrusive<Element>(7);
    auto p_element_8 = Kratos::make_intrusive<Element>(8);
    auto p_element_3_copy = Kratos::make_intrusive<Element>(3);

    std::vector<Element::Pointer> tmp;
    tmp.push_back(p_element_2);
    tmp.push_back(p_element_1);
    tmp.push_back(p_element_5);
    tmp.push_back(p_element_8);

    test_container_1.insert(tmp.begin(), tmp.end());
    KRATOS_EXPECT_EQ(test_container_1.size(), 4);
    auto itr = test_container_1.begin();
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_1);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_2);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_5);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_8);

    test_container_2.insert(test_container_1);
    KRATOS_EXPECT_EQ(test_container_2.size(), 4);
    itr = test_container_2.begin();
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_1);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_2);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_5);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_8);

    tmp.clear();
    tmp.push_back(p_element_3);
    tmp.push_back(p_element_4);
    tmp.push_back(p_element_5);
    tmp.push_back(p_element_3);
    tmp.push_back(p_element_3_copy);
    tmp.push_back(p_element_1);
    tmp.push_back(p_element_7);
    tmp.push_back(p_element_6);
    test_container_3.insert(tmp.begin(), tmp.end());

    test_container_2.insert(test_container_3);
    KRATOS_EXPECT_EQ(test_container_2.size(), 8);
    itr = test_container_2.begin();
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_1);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_2);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_3);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_4);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_5);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_6);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_7);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_8);
}

KRATOS_TEST_CASE_IN_SUITE(PointerVectorSetInsert8, KratosCoreFastSuite)
{
    PointerVectorSet<const Element, IndexedObject> test_container_1, test_container_2, test_container_3;
    std::vector<Element::Pointer> elements;
    for (IndexType i = 0; i < 50; ++i) {
        elements.push_back(Kratos::make_intrusive<Element>(i + 1));
    }
    auto p_element_4_copy = Kratos::make_intrusive<Element>(4);
    auto p_element_10_copy = Kratos::make_intrusive<Element>(10);

    std::vector<Element::Pointer> tmp;
    tmp.push_back(elements[2]);
    tmp.push_back(elements[1]);
    tmp.push_back(elements[25]);
    tmp.push_back(elements[28]);

    test_container_1.insert(tmp.begin(), tmp.end());
    KRATOS_EXPECT_EQ(test_container_1.size(), 4);
    auto itr = test_container_1.begin();
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[1]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[2]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[25]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[28]);

    test_container_2.insert(test_container_1);
    KRATOS_EXPECT_EQ(test_container_2.size(), 4);
    itr = test_container_2.begin();
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[1]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[2]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[25]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[28]);

    tmp.clear();
    tmp.push_back(elements[3]);
    tmp.push_back(elements[4]);
    tmp.push_back(elements[5]);
    tmp.push_back(elements[3]);
    tmp.push_back(p_element_4_copy);
    tmp.push_back(p_element_10_copy);
    tmp.push_back(elements[7]);
    tmp.push_back(elements[6]);
    test_container_3.insert(tmp.begin(), tmp.end());

    test_container_2.insert(test_container_3);
    KRATOS_EXPECT_EQ(test_container_2.size(), 10);
    itr = test_container_2.begin();
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[1]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[2]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[3]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[4]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[5]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[6]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[7]);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_10_copy);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[25]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[28]);

    tmp.clear();
    tmp.push_back(elements[39]);
    tmp.push_back(elements[29]);
    tmp.push_back(elements[48]);
    test_container_2.insert(tmp.begin(), tmp.end());

    KRATOS_EXPECT_EQ(test_container_2.size(), 13);
    itr = test_container_2.begin();
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[1]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[2]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[3]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[4]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[5]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[6]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[7]);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_10_copy);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[25]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[28]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[29]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[39]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[48]);

    tmp.clear();
    tmp.push_back(elements[28]);
    tmp.push_back(elements[31]);
    tmp.push_back(elements[45]);
    test_container_2.insert(tmp.begin(), tmp.end());

    KRATOS_EXPECT_EQ(test_container_2.size(), 15);
    itr = test_container_2.begin();
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[1]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[2]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[3]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[4]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[5]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[6]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[7]);
    KRATOS_EXPECT_EQ(&*(itr++), &*p_element_10_copy);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[25]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[28]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[29]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[31]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[39]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[45]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[48]);
}

KRATOS_TEST_CASE_IN_SUITE(TestPointerVectorSetInsert9, KratosCoreFastSuite)
{
    PointerVectorSet<Element, IndexedObject, std::less<IndexType>, std::equal_to<IndexType>, std::shared_ptr<Element>> test_container;

    std::vector<std::shared_ptr<Element>> elements;
    for (IndexType i = 0; i < 50; ++i) {
        elements.push_back(Kratos::make_shared<Element>(i + 1));
    }

    std::vector<std::shared_ptr<Element>> tmp;
    tmp.push_back(elements[2]);
    tmp.push_back(elements[1]);
    tmp.push_back(elements[28]);
    tmp.push_back(elements[25]);
    tmp.push_back(elements[28]);
    test_container.insert(tmp.begin(), tmp.end());

    KRATOS_EXPECT_EQ(test_container.size(), 4);
    auto itr = test_container.begin();
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[1]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[2]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[25]);
    KRATOS_EXPECT_EQ(&*(itr++), &*elements[28]);
}

KRATOS_TEST_CASE_IN_SUITE(TestPointerVectorSet, KratosCoreFastSuite)
{
    // create model and model part
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("Main");
    r_model_part.SetBufferSize(3);

    // create 2 elements
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    r_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
    r_model_part.CreateNewNode(5, 0.0, 0.0, 1.0);
    r_model_part.CreateNewNode(6, 1.0, 0.0, 1.0);
    r_model_part.CreateNewNode(7, 0.0, 1.0, 1.0);
    r_model_part.CreateNewNode(8, 0.0, 0.0, 2.0);
    std::vector<ModelPart::IndexType> elNodes1 {1, 2, 3, 4};
    std::vector<ModelPart::IndexType> elNodes2 {5, 6, 7, 8};
    auto p_elem_prop = r_model_part.CreateNewProperties(1);
    r_model_part.CreateNewElement("Element3D4N", 1, elNodes1, p_elem_prop);
    r_model_part.CreateNewElement("Element3D4N", 2, elNodes2, p_elem_prop);
    Element::Pointer p_element1 = r_model_part.pGetElement(1);
    Element::Pointer p_element2 = r_model_part.pGetElement(2);

    // set TO_ERASE = true first element
    p_element1->Set(TO_ERASE, true);
    // remove element 1 from model part
    r_model_part.RemoveElements();

    // list with only element 2
    std::vector<int> ids = {2};
    const auto& container = r_model_part.Elements();

    for(const int id : ids ) {
        const auto it = container.find(id);
        KRATOS_EXPECT_TRUE(it->Id() == 2);
    }
}

}
} // namespace Kratos.
