//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

#if !defined (KRATOS_TRILINOS_SOLVER_UTILITIES_H_INCLUDED)
#define KRATOS_TRILINOS_SOLVER_UTILITIES_H_INCLUDED

// External includes
#include "Teuchos_ParameterList.hpp"

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"

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

}  // namespace TrilinosSolverUtilities.
}  // namespace Kratos.

#endif // KRATOS_TRILINOS_SOLVER_UTILITIES_H_INCLUDED defined
