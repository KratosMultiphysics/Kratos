//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// Project includes
#include "Epetra_FEVector.h"
#include "trilinos_space.h"

#include "custom_strategies/relaxed_dof_updater.h"

namespace Kratos
{
template <>
void RelaxedDofUpdater<TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>>::Initialize(
    const RelaxedDofUpdater<TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>>::DofsArrayType& rDofSet,
    const RelaxedDofUpdater<TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>>::SystemVectorType& rDx)
{
    int system_size = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>::Size(rDx);
    int number_of_dofs = rDofSet.size();
    std::vector<int> index_array(number_of_dofs);

    // filling the array with the global ids
    unsigned int counter = 0;
    for (typename DofsArrayType::const_iterator i_dof = rDofSet.begin();
         i_dof != rDofSet.end(); ++i_dof)
    {
        int id = i_dof->EquationId();
        if (id < system_size)
        {
            index_array[counter++] = id;
        }
    }

    std::sort(index_array.begin(), index_array.end());
    std::vector<int>::iterator new_end =
        std::unique(index_array.begin(), index_array.end());
    index_array.resize(new_end - index_array.begin());

    int check_size = -1;
    int tot_update_dofs = index_array.size();
    rDx.Comm().SumAll(&tot_update_dofs, &check_size, 1);
    if ((check_size < system_size) && (rDx.Comm().MyPID() == 0))
    {
        std::stringstream msg;
        msg << "Dof count is not correct. There are less dofs then expected." << std::endl;
        msg << "Expected number of active dofs: " << system_size
            << ", dofs found: " << check_size << std::endl;
        KRATOS_ERROR << msg.str();
    }

    // defining a map as needed
    Epetra_Map dof_update_map(-1, index_array.size(), &(*(index_array.begin())),
                              0, rDx.Comm());

    // defining the import instance
    std::shared_ptr<Epetra_Import> p_dof_import(
        new Epetra_Import(dof_update_map, rDx.Map()));
    this->mpDofImport.swap(p_dof_import);

    this->mImportIsInitialized = true;
}

template <>
void RelaxedDofUpdater<TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>>::Clear()
{
    this->mpDofImport.reset();
    this->mImportIsInitialized = false;
}

template <>
void RelaxedDofUpdater<TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>>::UpdateDofs(
    RelaxedDofUpdater<TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>>::DofsArrayType& rDofSet,
    const RelaxedDofUpdater<TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>>::SystemVectorType& rDx,
    const double RelaxationFactor)
{
    KRATOS_TRY;

    if (!this->mImportIsInitialized)
        this->Initialize(rDofSet, rDx);

    int system_size = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>::Size(rDx);

    // defining a temporary vector to gather all of the values needed
    Epetra_Vector local_dx(this->mpDofImport->TargetMap());

    // importing in the new temp vector the values
    int ierr = local_dx.Import(rDx, *this->mpDofImport, Insert);
    KRATOS_ERROR_IF(ierr != 0)
        << "Epetra failure found while trying to import Dx." << std::endl;

    int num_dof = rDofSet.size();

    for (int i = 0; i < num_dof; ++i)
    {
        auto it_dof = rDofSet.begin() + i;

        if (it_dof->IsFree())
        {
            int global_id = it_dof->EquationId();
            if (global_id < system_size)
            {
                double dx_i = local_dx[this->mpDofImport->TargetMap().LID(global_id)];
                it_dof->GetSolutionStepValue() += RelaxationFactor * dx_i;
            }
        }
    }

    KRATOS_CATCH("");
}

template <>
void RelaxedDofUpdater<TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>>::AssignDofs(
    RelaxedDofUpdater<TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>>::DofsArrayType& rDofSet,
    const RelaxedDofUpdater<TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>>::SystemVectorType& rX)
{
    KRATOS_TRY;

    if (!this->mImportIsInitialized)
        this->Initialize(rDofSet, rX);

    int system_size = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>::Size(rX);

    // defining a temporary vector to gather all of the values needed
    Epetra_Vector local_dx(this->mpDofImport->TargetMap());

    // importing in the new temp vector the values
    int ierr = local_dx.Import(rX, *this->mpDofImport, Insert);
    KRATOS_ERROR_IF(ierr != 0)
        << "Epetra failure found while trying to import Dx." << std::endl;

    int num_dof = rDofSet.size();

    for (int i = 0; i < num_dof; ++i)
    {
        auto it_dof = rDofSet.begin() + i;

        if (it_dof->IsFree())
        {
            int global_id = it_dof->EquationId();
            if (global_id < system_size)
            {
                double dx_i = local_dx[this->mpDofImport->TargetMap().LID(global_id)];
                it_dof->GetSolutionStepValue() = dx_i;
            }
        }
    }

    KRATOS_CATCH("");
}

template <>
std::string RelaxedDofUpdater<TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>>::Info() const
{
    std::stringstream buffer;
    buffer << "RelaxedDofUpdater - TrilinosSpace";
    return buffer.str();
}
template class RelaxedDofUpdater<TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>>;
} // namespace Kratos