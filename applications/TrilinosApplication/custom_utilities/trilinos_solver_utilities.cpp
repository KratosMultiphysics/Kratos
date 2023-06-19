//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// External includes
#include "Epetra_MpiComm.h"

// Project includes
#include "trilinos_solver_utilities.h"

namespace Kratos {
namespace TrilinosSolverUtilities {

void SetTeuchosParameters(const Parameters rSettings, Teuchos::ParameterList& rParameterlist)
{
    for (auto it = rSettings.begin(); it != rSettings.end(); ++it) {
        if      (it->IsString()) rParameterlist.set(it.name(), it->GetString());
        else if (it->IsInt())    rParameterlist.set(it.name(), it->GetInt());
        else if (it->IsBool())   rParameterlist.set(it.name(), it->GetBool());
        else if (it->IsDouble()) rParameterlist.set(it.name(), it->GetDouble());
    }
}

MPI_Comm GetMPICommFromEpetraComm(const Epetra_Comm& rEpetraComm)
{
    // see https://github.com/trilinos/Trilinos/issues/10122#issuecomment-1021614956
    const Epetra_MpiComm& r_epetra_mpi_comm = dynamic_cast<const Epetra_MpiComm&>(rEpetraComm); // cannot use static_cast due to virtual inheritance
    return r_epetra_mpi_comm.Comm();
}

}  // namespace TrilinosSolverUtilities.
}  // namespace Kratos.
