//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Angel Celigueta
//

// System includes

// External includes

// Project includes
#include "includes/gid_io.h"

namespace Kratos
{
    GidIOBase& GidIOBase::GetInstance() {
        if (mpInstance == nullptr) {
            Create();
        }

        return *mpInstance;
    }

    void GidIOBase::Create() {
        static GidIOBase gid_io_base;
        mpInstance = &gid_io_base;
    }

    int GidIOBase::GetData() {
        return this -> data;
    }

    void GidIOBase::SetData(int data) {
        this -> data = data;
    }

    GidIOBase* GidIOBase::mpInstance = nullptr;

    // GidIO default instantiation
    template class GidIO<GidGaussPointsContainer,GidMeshContainer>;
}  // namespace Kratos.


