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


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

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
class DofUpdater
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DofUpdater
    KRATOS_CLASS_POINTER_DEFINITION(DofUpdater);

    using DofsArrayType = ModelPart::DofsArrayType;
    using SystemVectorType = typename TSparseSpace::VectorType;

    using UniquePointer = std::unique_ptr<DofUpdater>;

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

    virtual UniquePointer Create() const
    {
        return UniquePointer(new DofUpdater());
    }

    virtual void Initialize(
        DofsArrayType& rDofSet,
        SystemVectorType& rDx)
    {}

    virtual void Clear() {}

    virtual void UpdateDOF(
        DofsArrayType& rDofSet,
        SystemVectorType& rDx)
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
