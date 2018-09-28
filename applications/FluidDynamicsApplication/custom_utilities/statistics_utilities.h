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
//                   Author2 Fullname
//

#ifndef KRATOS_STATISTICAL_UTILITIES_H_INCLUDED
#define KRATOS_STATISTICAL_UTILITIES_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
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
template< class ValueContainterType >
class StatisticsUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of StatisticsUtilities
    KRATOS_CLASS_POINTER_DEFINITION(StatisticsUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    StatisticsUtilities() {}

    /// Destructor.
    virtual ~StatisticsUtilities() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Update(const ValueContainerType& rMeasurement) {}

    void Combine(const StatisticsUtilities& rOther) {}

    ValueContainerType Output() {}

    void LoadFromFile() {}

    void SaveToFile() {}

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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "StatisticsUtilities";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const { rOStream << "StatisticsUtilities"; }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const {}

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

    /// Assignment operator.
    StatisticsUtilities &operator=(StatisticsUtilities const &rOther) {}

    /// Copy constructor.
    StatisticsUtilities(StatisticsUtilities const &rOther) {}

    ///@}

}; // Class StatisticsUtilities

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                StatisticsUtilities &rThis) {}

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const StatisticsUtilities &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_STATISTICAL_UTILITIES_H_INCLUDED  defined
