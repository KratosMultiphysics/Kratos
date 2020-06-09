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

#if !defined(KRATOS_RELAXED_DOF_UPDATER_H_INCLUDED )
#define  KRATOS_RELAXED_DOF_UPDATER_H_INCLUDED

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
class Epetra_Import;
#endif



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
 *  to obtain out-of-process update data. RelaxedDofUpdater takes care of both the operation and the eventual
 *  auxiliary infrastructure.
  */
template< class TSparseSpace >
class KRATOS_API(RANS_APPLICATION) RelaxedDofUpdater
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RelaxedDofUpdater
    KRATOS_CLASS_POINTER_DEFINITION(RelaxedDofUpdater);

	using DofType = Dof<typename TSparseSpace::DataType>;
	using DofsArrayType = PointerVectorSet<
		DofType,
		SetIdentityFunction<DofType>,
		std::less<typename SetIdentityFunction<DofType>::result_type>,
		std::equal_to<typename SetIdentityFunction<DofType>::result_type>,
		DofType* >;

    using SystemVectorType = typename TSparseSpace::VectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RelaxedDofUpdater(){}

    /// Deleted copy constructor
    RelaxedDofUpdater(RelaxedDofUpdater const& rOther) = delete;

    /// Destructor.
    virtual ~RelaxedDofUpdater() = default;

    /// Deleted assignment operator
    RelaxedDofUpdater& operator=(RelaxedDofUpdater const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// Create a new instance of this class.
    /** This function is used by the SparseSpace class to create new
     *  RelaxedDofUpdater instances of the appropriate type.
     *  @return a std::unique_pointer to the new instance.
     *  @see UblasSpace::CreateRelaxedDofUpdater(), TrilinosSpace::CreateRelaxedDofUpdater().
     */
    virtual typename RelaxedDofUpdater::UniquePointer Create() const
    {
        return Kratos::make_unique<RelaxedDofUpdater>();
    }

    /// Initialize the RelaxedDofUpdater in preparation for a subsequent UpdateDofs call.
    /** Note that the base RelaxedDofUpdater does not have internal data, so this does nothing.
     *  @param[in] rDofSet The list of degrees of freedom.
     *  @param[in] rDx The update vector.
     */
    virtual void Initialize(
        const DofsArrayType& rDofSet,
        const SystemVectorType& rDx);

    /// Free internal storage to reset the instance and/or optimize memory consumption.
    /** Note that the base RelaxedDofUpdater does not have internal data, so this does nothing.
     */
    virtual void Clear();

    /// Calculate new values for the problem's degrees of freedom using the update vector rDx.
    /** For each Dof in rDofSet, this function calculates the updated value for the corresponding
     *  variable as value += rDx[dof.EquationId()].
     *  @param[in/out] rDofSet The list of degrees of freedom.
     *  @param[in] rDx The update vector.
     *  @param[in] RelaxationFactor Relaxation factor
     *  This method will check if Initialize() was called before and call it if necessary.
     */
    virtual void UpdateDofs(
        DofsArrayType& rDofSet,
        const SystemVectorType& rDx,
        const double RelaxationFactor);

    /// Assign new values for the problem's degrees of freedom using the vector rX.
    /** For each Dof in rDofSet, this function assigns the value for the corresponding
     *  variable as value = rX[dof.EquationId()].
     *  @param[in/out] rDofSet The list of degrees of freedom.
     *  @param[in] rX The solution vector.
     *  This method will check if Initialize() was called before and call it if necessary.
     */
    virtual void AssignDofs(DofsArrayType& rDofSet, const SystemVectorType& rX);

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << this->Info() << std::endl;
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << this->Info() << std::endl;
    }

    ///@}
private:

    //@name Member Variables
    ///@{

    /// This lets the class control if Initialize() was properly called.
    bool mImportIsInitialized = false;

    #ifdef KRATOS_USING_MPI // mpi-parallel compilation
    /// Auxiliary trilinos data structure to import out-of-process data in the update vector.
    std::shared_ptr<Epetra_Import> mpDofImport = nullptr;
    #endif

    ///@}


}; // Class RelaxedDofUpdater


///@}
///@name Input and output
///@{

/// input stream function
template< class TSparseSpace >
inline std::istream& operator >> (
    std::istream& rIStream,
    RelaxedDofUpdater<TSparseSpace>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TSparseSpace >
inline std::ostream& operator << (
    std::ostream& rOStream,
    const RelaxedDofUpdater<TSparseSpace>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_RELAXED_DOF_UPDATER_H_INCLUDED  defined
