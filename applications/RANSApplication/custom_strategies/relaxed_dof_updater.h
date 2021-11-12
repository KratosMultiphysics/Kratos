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
#include "utilities/dof_updater.h"

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
 *  for each dof: dof.variable += dx[dof.equation_id] * Relaxation
 *  This operation is a simple loop in shared memory, but requires additional infrastructure in MPI,
 *  to obtain out-of-process update data. RelaxedDofUpdater takes care of both the operation and the eventual
 *  auxiliary infrastructure.
  */
template< class TSparseSpace >
class KRATOS_API(RANS_APPLICATION) RelaxedDofUpdater
    :public DofUpdater<TSparseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RelaxedDofUpdater
    KRATOS_CLASS_POINTER_DEFINITION(RelaxedDofUpdater);

    using BaseType = DofUpdater<TSparseSpace>;

    using DofType = typename BaseType::DofType;

    using DofsArrayType = typename BaseType::DofsArrayType;

    using SystemVectorType = typename TSparseSpace::VectorType;



    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RelaxedDofUpdater(const double RelaxationFactor)
        : mRelaxationFactor(RelaxationFactor),
          mInitialRelaxationFactor(0.1),
          mMinRelaxationFactor(1e-5),
          mMaxRelaxationFactor(1.0)
    {}

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
    typename BaseType::UniquePointer Create() const override
    {
        return Kratos::make_unique<RelaxedDofUpdater>(mRelaxationFactor);
    }

    /// Initialize the RelaxedDofUpdater in preparation for a subsequent UpdateDofs call.
    /** Note that the base RelaxedDofUpdater does not have internal data, so this does nothing.
     *  @param[in] rDofSet The list of degrees of freedom.
     *  @param[in] rDx The update vector.
     */
    void Initialize(
        const DofsArrayType& rDofSet,
        const SystemVectorType& rDx) override;

    void InitializeAitken() {
        mOldRelaxationFactor = 0.0;
    }

    /// Free internal storage to reset the instance and/or optimize memory consumption.
    /** Note that the base RelaxedDofUpdater does not have internal data, so this does nothing.
     */
    void Clear() override;

    /// Calculate new values for the problem's degrees of freedom using the update vector rDx.
    /** For each Dof in rDofSet, this function calculates the updated value for the corresponding
     *  variable as value += rDx[dof.EquationId()].
     *  @param[in/out] rDofSet The list of degrees of freedom.
     *  @param[in] rDx The update vector.
     *  @param[in] RelaxationFactor Relaxation factor
     *  This method will check if Initialize() was called before and call it if necessary.
     */
    void UpdateDofs(
        DofsArrayType& rDofSet,
        const SystemVectorType& rDx) override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

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
    const double mRelaxationFactor;

    const double mInitialRelaxationFactor;
    const double mMinRelaxationFactor;
    const double mMaxRelaxationFactor;

    double mOldRelaxationFactor;
    Vector mOldDx;
    Vector mDiffDx;

    #ifdef KRATOS_USING_MPI // mpi-parallel compilation
    /// Auxiliary trilinos data structure to import out-of-process data in the update vector.
    std::shared_ptr<Epetra_Import> mpDofImport = nullptr;
    #else

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
