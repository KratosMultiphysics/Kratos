//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

#ifndef KRATOS_MPI_MANAGER_H_INCLUDED
#define KRATOS_MPI_MANAGER_H_INCLUDED

#include "includes/parallel_environment.h"

namespace Kratos
{

class KRATOS_API(KRATOS_MPI_CORE) MPIManager: public EnvironmentManager
{
public:
    typedef std::unique_ptr<MPIManager> Pointer;

    MPIManager(MPIManager& rOther) = delete;

    ~MPIManager() override;

    static EnvironmentManager::Pointer Create();

    bool IsInitialized() const override;

    bool IsFinalized() const override;
private:
    MPIManager();
};

}

#endif