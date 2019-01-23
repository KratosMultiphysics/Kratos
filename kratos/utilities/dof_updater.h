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

#if !defined(KRATOS_DOF_UPDATER_H_INCLUDED )
#define  KRATOS_DOF_UPDATER_H_INCLUDED

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

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
 *  to obtain out-of-process update data. DofUpdater takes care of both the operation and the eventual
 *  auxiliary infrastructure.
 *  @see TrilinosDofUpdater for the trilinos version.
 */
template< class TSparseSpace >
class DofUpdater
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DofUpdater
    KRATOS_CLASS_POINTER_DEFINITION(DofUpdater);

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
    DofUpdater(){}

    /// Deleted copy constructor
    DofUpdater(DofUpdater const& rOther) = delete;

    /// Destructor.
    virtual ~DofUpdater(){}

    /// Deleted assignment operator
    DofUpdater& operator=(DofUpdater const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// Create a new instance of this class.
    /** This function is used by the SparseSpace class to create new
     *  DofUpdater instances of the appropriate type.
     *  @return a std::unique_pointer to the new instance.
     *  @see UblasSpace::CreateDofUpdater(), TrilinosSpace::CreateDofUpdater().
     */
    virtual typename DofUpdater::UniquePointer Create() const
    {
        return Kratos::make_unique<DofUpdater>();
    }

    /// Initialize the DofUpdater in preparation for a subsequent UpdateDofs call.
    /** Note that the base DofUpdater does not have internal data, so this does nothing.
     *  @param[in] rDofSet The list of degrees of freedom.
     *  @param[in] rDx The update vector.
     */
    virtual void Initialize(
        const DofsArrayType& rDofSet,
        const SystemVectorType& rDx)
    {}

    /// Free internal storage to reset the instance and/or optimize memory consumption.
    /** Note that the base DofUpdater does not have internal data, so this does nothing.
     */
    virtual void Clear() {}

    /// Calculate new values for the problem's degrees of freedom using the update vector rDx.
    /** For each Dof in rDofSet, this function calculates the updated value for the corresponding
     *  variable as value += rDx[dof.EquationId()].
     *  @param[in/out] rDofSet The list of degrees of freedom.
     *  @param[in] rDx The update vector.
     *  This method will check if Initialize() was called before and call it if necessary.
     */
    virtual void UpdateDofs(
        DofsArrayType& rDofSet,
        const SystemVectorType& rDx)
    {
        const int num_dof = static_cast<int>(rDofSet.size());

        #pragma omp parallel for
        for(int i = 0;  i < num_dof; ++i) {
            auto it_dof = rDofSet.begin() + i;

			if (it_dof->IsFree())
                it_dof->GetSolutionStepValue() += TSparseSpace::GetValue(rDx,it_dof->EquationId());
        }
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "DofUpdater" ;
        return buffer.str();
    }

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

}; // Class DofUpdater


///@}
///@name Input and output
///@{

/// input stream function
template< class TSparseSpace >
inline std::istream& operator >> (
    std::istream& rIStream,
    DofUpdater<TSparseSpace>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TSparseSpace >
inline std::ostream& operator << (
    std::ostream& rOStream,
    const DofUpdater<TSparseSpace>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_DOF_UPDATER_H_INCLUDED  defined
