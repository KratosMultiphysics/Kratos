//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_RESIDUALBASED_INCREMENTALUPDATE_DRYWETTING_SCHEME_H_INCLUDED
#define KRATOS_RESIDUALBASED_INCREMENTALUPDATE_DRYWETTING_SCHEME_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"


namespace Kratos
{
///@addtogroup ShallowWaterApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
template<class TSparseSpace, class TDenseSpace >
class ResidualBasedIncrementalUpdateDryWettingScheme : public ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ResidualBasedIncrementalUpdateDryWettingScheme
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedIncrementalUpdateDryWettingScheme);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ResidualBasedIncrementalUpdateDryWettingScheme(){}

    /// Destructor.
    virtual ~ResidualBasedIncrementalUpdateDryWettingScheme(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ResidualBasedIncrementalUpdateDryWettingScheme";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class ResidualBasedIncrementalUpdateDryWettingScheme

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_RESIDUALBASED_INCREMENTALUPDATE_DRYWETTING_SCHEME_H_INCLUDED  defined


