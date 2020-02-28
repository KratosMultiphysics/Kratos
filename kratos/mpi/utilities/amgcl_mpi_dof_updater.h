//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher), based on the work of Jordi Cotela
//

#ifndef KRATOS_AMGCL_MPI_DOF_UPDATER_H_INCLUDED
#define KRATOS_AMGCL_MPI_DOF_UPDATER_H_INCLUDED

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
 *  to obtain out-of-process update data. AmgclMPIDofUpdater manages both the update operation and the
 *  auxiliary infrastructure.
 */
template< class TSparseSpace >
class AmgclMPIDofUpdater : public DofUpdater<TSparseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AmgclMPIDofUpdater
    KRATOS_CLASS_POINTER_DEFINITION(AmgclMPIDofUpdater);

    using BaseType = DofUpdater<TSparseSpace>;
    using DofsArrayType = typename BaseType::DofsArrayType;
    using SystemVectorType = typename BaseType::SystemVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AmgclMPIDofUpdater():
        DofUpdater<TSparseSpace>()
    {}

    /// Deleted copy constructor
    AmgclMPIDofUpdater(AmgclMPIDofUpdater const& rOther) = delete;

    /// Destructor.
    ~AmgclMPIDofUpdater() override {}

    /// Deleted assignment operator
    AmgclMPIDofUpdater& operator=(AmgclMPIDofUpdater const& rOther) = delete;

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
        return Kratos::make_unique<AmgclMPIDofUpdater>();
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
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
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
        KRATOS_ERROR << "THIS FUNCTION IS NOT YET IMPLEMENTED!" << std::endl;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "AmgclMPIDofUpdater" ;
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

}; // Class AmgclMPIDofUpdater

///@}
///@name Input and output
///@{

/// input stream function
template< class TSparseSpace >
inline std::istream& operator >> (
    std::istream& rIStream,
    AmgclMPIDofUpdater<TSparseSpace>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TSparseSpace >
inline std::ostream& operator << (
    std::ostream& rOStream,
    const AmgclMPIDofUpdater<TSparseSpace>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_AMGCL_MPI_DOF_UPDATER_H_INCLUDED  defined
