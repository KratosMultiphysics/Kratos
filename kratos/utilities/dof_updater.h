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
class DOFUpdater
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DOFUpdater
    KRATOS_CLASS_POINTER_DEFINITION(DOFUpdater);

    using DofsArrayType = ModelPart::DofsArrayType;
    using SystemVectorType = typename TSparseSpace::VectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DOFUpdater(){}

    /// Deleted copy constructor
    DOFUpdater(DOFUpdater const& rOther) = delete;

    /// Destructor.
    virtual ~DOFUpdater(){}

    /// Deleted assignment operator
    DOFUpdater& operator=(DOFUpdater const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    virtual std::unique_ptr<DOFUpdater> Create()
    {
        return std::unique_ptr<DOFUpdater>(new DOFUpdater());
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
        buffer << "DOFUpdater" ;
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

}; // Class DOFUpdater

///@}
///@name Input and output
///@{

/// input stream function
template< class TSparseSpace >
inline std::istream& operator >> (
    std::istream& rIStream,
    DOFUpdater<TSparseSpace>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TSparseSpace >
inline std::ostream& operator << (
    std::ostream& rOStream,
    const DOFUpdater<TSparseSpace>& rThis)
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
