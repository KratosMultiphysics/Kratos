//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

/* Trilinos includes */
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Import.hpp>
#include <Teuchos_RCP.hpp>

// Project includes
#include "includes/define.h"
#include "utilities/dof_updater.h"

namespace Kratos
{
///@addtogroup TrilinosApplication
///@{

///@name Kratos Classes
///@{

/// Utility class to update the values of degree of freedom (Dof) variables after solving the system (Tpetra version).
/** This class encapsulates the operation of updating nodal degrees of freedom after a system solution.
 *  In pseudo-code, the operation to be performed is
 *  for each dof: dof.variable += dx[dof.equation_id]
 *  This operation is a simple loop in shared memory, but requires additional infrastructure in MPI,
 *  to obtain out-of-process update data. TrilinosDofUpdaterExperimental manages both the update operation
 *  and the auxiliary infrastructure, using Tpetra import objects created by the experimental space.
 *  @see TrilinosDofUpdater for the Epetra counterpart.
 */
template< class TSparseSpace >
class TrilinosDofUpdaterExperimental : public DofUpdater<TSparseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TrilinosDofUpdaterExperimental
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosDofUpdaterExperimental);

    using BaseType = DofUpdater<TSparseSpace>;
    using DofsArrayType = typename BaseType::DofsArrayType;
    using SystemVectorType = typename BaseType::SystemVectorType;

    /// Tpetra scalar/ordinal definitions (taken from the space)
    using ST = typename TSparseSpace::ST;
    using LO = typename TSparseSpace::LO;
    using GO = typename TSparseSpace::GO;
    using NT = typename TSparseSpace::NT;

    /// Import definitions
    using ImportPointerType = typename TSparseSpace::ImportPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TrilinosDofUpdaterExperimental():
        DofUpdater<TSparseSpace>()
    {}

    /// Deleted copy constructor
    TrilinosDofUpdaterExperimental(TrilinosDofUpdaterExperimental const& rOther) = delete;

    /// Destructor.
    ~TrilinosDofUpdaterExperimental() override {}

    /// Deleted assignment operator
    TrilinosDofUpdaterExperimental& operator=(TrilinosDofUpdaterExperimental const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// Create a new instance of this class.
    /** This function is used by the SparseSpace class to create new
     *  DofUpdater instances of the appropriate type.
     *  @return a std::unique_pointer to the new instance.
     *  Note that the pointer is actually a pointer to the base class type.
     *  @see UblasSpace::CreateDofUpdater(), TrilinosSpaceExperimental::CreateDofUpdater().
     */
    typename BaseType::UniquePointer Create() const override
    {
        return Kratos::make_unique<TrilinosDofUpdaterExperimental>();
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
        KRATOS_TRY;

        // CreateImport already performs the global DOF-count consistency check
        mpDofImport = TSparseSpace::CreateImport(rDofSet, rDx);
        mImportIsInitialized = true;

        KRATOS_CATCH("");
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
            this->Initialize(rDofSet, rDx);

        const std::size_t system_size = TSparseSpace::Size(rDx);

        // Defining a temporary vector to gather all of the values needed
        using MultiVectorType = Tpetra::MultiVector<ST, LO, GO, NT>;
        const auto p_target_map = mpDofImport->getTargetMap();
        MultiVectorType local_dx(p_target_map, 1);

        // Importing in the temporary vector the values
        local_dx.doImport(rDx, *mpDofImport, Tpetra::INSERT);

        // Performing the update
        // NOTE: Plain serial loop on purpose: mixing Kratos OpenMP parallel utilities with
        // Kokkos-backed Tpetra host views causes conflicts (see TrilinosBlockBuilderAndSolver)
        const auto data = local_dx.getData(0);
        for (auto it_dof = rDofSet.begin(); it_dof != rDofSet.end(); ++it_dof) {
            if (it_dof->IsFree()) {
                const std::size_t global_id = it_dof->EquationId();
                if (global_id < system_size) {
                    const LO local_id = p_target_map->getLocalElement(static_cast<GO>(global_id));
                    it_dof->GetSolutionStepValue() += data[local_id];
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
        buffer << "TrilinosDofUpdaterExperimental" ;
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
    ImportPointerType mpDofImport = nullptr;

    ///@}

}; // Class TrilinosDofUpdaterExperimental

///@}
///@name Input and output
///@{

/// input stream function
template< class TSparseSpace >
inline std::istream& operator >> (
    std::istream& rIStream,
    TrilinosDofUpdaterExperimental<TSparseSpace>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TSparseSpace >
inline std::ostream& operator << (
    std::ostream& rOStream,
    const TrilinosDofUpdaterExperimental<TSparseSpace>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

}  // namespace Kratos.
