//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/python_object_cpp_wrapper_utility.h"

namespace Kratos
{
PythonObjectCppWrapperUtility::PythonObjectCppWrapperUtility(ListType& rObjectList)
{
    const std::size_t size_objects = len(rObjectList);
    mListPythonObjects.resize(size_objects);
    for (std::size_t i_object = 0; i_object < size_objects; ++i_object)
        mListPythonObjects[i_object] = pybind11::cast<ObjectType>(rObjectList[i_object]);
}

/***********************************************************************************/
/***********************************************************************************/

PythonObjectCppWrapperUtility::PythonObjectCppWrapperUtility(ObjectType& rObject)
{
    mListPythonObjects.resize(1);
    mListPythonObjects[0] = rObject;
}

/***********************************************************************************/
/***********************************************************************************/

void PythonObjectCppWrapperUtility::AddObject(ObjectType& rObject)
{
    mListPythonObjects.push_back(rObject);
}

/***********************************************************************************/
/***********************************************************************************/

void PythonObjectCppWrapperUtility::AddObjects(ListType& rObjectList)
{
    const std::size_t size_objects = len(rObjectList);
    for (std::size_t i_object = 0; i_object < size_objects; ++i_object)
        mListPythonObjects.push_back(pybind11::cast<ObjectType>(rObjectList[i_object]));
}

/***********************************************************************************/
/***********************************************************************************/

void PythonObjectCppWrapperUtility::Execute(const std::string& rNameMethod)
{
    for (auto& process : mListPythonObjects)
        process.attr(rNameMethod.c_str())();
}
}  // namespace Kratos.
