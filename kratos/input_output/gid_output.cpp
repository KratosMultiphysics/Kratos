//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Janosch Stascheit
//                   Pooyan Dadvand
//

// System includes

// External includes

// Project includes
#include "input_output/gid_output.h"

namespace Kratos
{
GidIOBase& GidIOBase::GetInstance()
{
    if (mpInstance == nullptr) {
        Create();
    }

    return *mpInstance;
}

/***********************************************************************************/
/***********************************************************************************/

void GidIOBase::Create() {
    static GidIOBase sGidIOBase;
    mpInstance = &sGidIOBase;
}

/***********************************************************************************/
/***********************************************************************************/

int GidIOBase::GetData() {
    return this -> data;
}

/***********************************************************************************/
/***********************************************************************************/

void GidIOBase::SetData(int data) {
    this -> data = data;
}

/***********************************************************************************/
/***********************************************************************************/

GidIOBase* GidIOBase::mpInstance = nullptr;

} // namespace Kratos.
