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

/// Utility class to update the values of degree of freedom (Dof) variables after solving the system.
/** This class encapsulates the operation of updating nodal degrees of freedom after a system solution.
 *  In pseudo-code, the operation to be performed is
 *  for each dof: dof.variable += dx[dof.equation_id]
 *  This operation is a simple loop in shared memory, but requires additional infrastructure in MPI,
 *  to obtain out-of-process update data. TrilinosDofUpdater manages both the update operation and the
 *  auxiliary infrastructure.
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

    /// Create a new instance of this class.
    /** This function is used by the SparseSpace class to create new
     *  DofUpdater instances of the appropriate type.
     *  @return a std::unique_pointer to the new instance.
     *  Note that the pointer is actually a pointer to the base class type.
     *  @see UblasSpace::CreateDofUpdater(), TrilinosSpace::CreateDofUpdater().
     */
    typename BaseType::UniquePointer Create() const override
    {
        return Kratos::make_unique<TrilinosDofUpdater>();
    }

    /// Initialize the DofUpdater in preparation for a subsequent UpdateDofs call.
    /** @param[in] rDofSet The list of degrees of freedom.
     *  @param[in] rDx The update vector.
     *  The DofUpdater needs to be initialized only if the dofset changes.
     *  If the problem does not require creating/destroying nodes or changing the
     *  mesh graph, it is in general enough to intialize this tool once at the
     *  begining of the problem.
     *  If the dofset only changes under certain conditions (for example because
     *  the domain is remeshed every N iterations), it is enough to call the
     *  Clear method to let this class know that its auxiliary data has to be re-generated
     *  and Initialize will be called as part of the next UpdateDofs call.
     */
    void Initialize(
        const DofsArrayType& rDofSet,
        const SystemVectorType& rDx) override
    {
        int system_size = TSparseSpace::Size(rDx);
        int number_of_dofs = rDofSet.size();
        std::vector< int > index_array(number_of_dofs);

        // filling the array with the global ids
        unsigned int counter = 0;
        for(typename DofsArrayType::const_iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
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

    /// Free internal storage to reset the instance and/or optimize memory consumption.
    void Clear() override
    {
        mpDofImport.reset();
        mImportIsInitialized = false;
    }

    /// Calculate new values for the problem's degrees of freedom using the update vector rDx.
    /** For each Dof in rDofSet, this function calculates the updated value for the corresponding
     *  variable as value += rDx[dof.EquationId()].
     *  @param[in/out] rDofSet The list of degrees of freedom.
     *  @param[in] rDx The update vector.
     *  This method will check if Initialize() was called before and call it if necessary.
     */
    void UpdateDofs(
        DofsArrayType& rDofSet,
        const SystemVectorType& rDx) override
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

        int num_dof = rDofSet.size();

        // performing the update
        #pragma omp parallel for
        for(int i = 0;  i < num_dof; ++i) {
            auto it_dof = rDofSet.begin() + i;

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

    /// This lets the class control if Initialize() was properly called.
    bool mImportIsInitialized = false;

    /// Auxiliary trilinos data structure to import out-of-process data in the update vector.
    std::unique_ptr<Epetra_Import> mpDofImport = nullptr;

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
