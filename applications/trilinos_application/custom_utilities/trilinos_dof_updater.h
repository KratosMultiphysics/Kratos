//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#if !defined(KRATOS_TRILINOS_DOF_UPDATER_H_INCLUDED )
#define  KRATOS_TRILINOS_DOF_UPDATER_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/dof_updater.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
template< class TSparseSpace >
class TrilinosDofUpdater : public DofUpdater<TSparseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TrilinosDofUpdater
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosDofUpdater);

    using BaseType = DofUpdater<TSparseSpace>;
    using DofsArrayType = typename BaseType::DofsArrayType;
    using SystemVectorType = typename BaseType::SystemVectorType;
    using UniquePointer = typename BaseType::UniquePointer;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TrilinosDofUpdater():
        DofUpdater<TSparseSpace>()
    {}

    /// Deleted copy constructor
    TrilinosDofUpdater(TrilinosDofUpdater const& rOther) = delete;

    /// Destructor.
    ~TrilinosDofUpdater() override {}

    /// Deleted assignment operator
    TrilinosDofUpdater& operator=(TrilinosDofUpdater const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    UniquePointer Create() const override
    {
        return UniquePointer(new TrilinosDofUpdater());
    }

    void Initialize(
        DofsArrayType& rDofSet,
        SystemVectorType& rDx) override
    {

        int system_size = TSparseSpace::Size(rDx);
        int number_of_dofs = rDofSet.size();
        std::vector< int > index_array(number_of_dofs);

        // filling the array with the global ids
        int counter = 0;
        for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
            int id = i_dof->EquationId();
            if( id < system_size )
            {
                index_array[counter++] = id;
            }
        }

        std::sort(index_array.begin(),index_array.end());
        std::vector<int>::iterator new_end = std::unique(index_array.begin(),index_array.end());
        index_array.resize(new_end - index_array.begin());

        int check_size = -1;
        int tot_update_dofs = index_array.size();
        rDx.Comm().SumAll(&tot_update_dofs,&check_size,1);
        if ( (check_size < system_size) && (rDx.Comm().MyPID() == 0) )
        {
            std::stringstream msg;
            msg << "Dof count is not correct. There are less dofs then expected." << std::endl;
            msg << "Expected number of active dofs: " << system_size << ", dofs found: " << check_size << std::endl;
            KRATOS_ERROR << msg.str();
        }

        // defining a map as needed
        Epetra_Map dof_update_map(-1,index_array.size(), &(*(index_array.begin())),0,rDx.Comm() );

        // defining the import instance
        std::unique_ptr<Epetra_Import> p_dof_import(new Epetra_Import(dof_update_map,rDx.Map()));
        mpDofImport.swap(p_dof_import);

        mImportIsInitialized = true;
    }

    void Clear() override
    {
        mpDofImport.reset();
        mImportIsInitialized = false;
    }

    void UpdateDOF(
        DofsArrayType& rDofSet,
        SystemVectorType& rDx) override
    {
        KRATOS_TRY;

        if (!mImportIsInitialized)
            this->Initialize(rDofSet,rDx);

        int system_size = TSparseSpace::Size(rDx);

        // defining a temporary vector to gather all of the values needed
        Epetra_Vector local_dx( mpDofImport->TargetMap() );

        // importing in the new temp vector the values
        int ierr = local_dx.Import(rDx,*mpDofImport,Insert) ;
        KRATOS_ERROR_IF(ierr != 0) << "Epetra failure found while trying to import Dx." << std::endl;

        rDx.Comm().Barrier();

        // performing the update
        for (typename DofsArrayType::iterator it_dof = rDofSet.begin(); it_dof != rDofSet.end(); ++it_dof)
        {
            if (it_dof->IsFree()) {
                int global_id = it_dof->EquationId();
                if(global_id < system_size) {
                    double dx_i = local_dx[mpDofImport->TargetMap().LID(global_id)];
                    it_dof->GetSolutionStepValue() += dx_i;
                }
            }
        }

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "TrilinosDofUpdater" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << this->Info() << std::endl;
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << this->Info() << std::endl;
    }

    ///@}

private:

    //@name Member Variables
    ///@{

    bool mImportIsInitialized;

    std::unique_ptr<Epetra_Import> mpDofImport;

    ///@}

}; // Class TrilinosDofUpdater

///@}
///@name Input and output
///@{

/// input stream function
template< class TSparseSpace >
inline std::istream& operator >> (
    std::istream& rIStream,
    TrilinosDofUpdater<TSparseSpace>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TSparseSpace >
inline std::ostream& operator << (
    std::ostream& rOStream,
    const TrilinosDofUpdater<TSparseSpace>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_TRILINOS_DOF_UPDATER_H_INCLUDED  defined
