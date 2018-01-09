//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                    
//



#if !defined(KRATOS_COUNTER_H_INCLUDED )
#define  KRATOS_COUNTER_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"


namespace Kratos
{

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
template<class TCountedType>
class Counter
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Counter
    KRATOS_CLASS_POINTER_DEFINITION(Counter);

    typedef unsigned int SizeType;

    ///@}
    ///@name Life Cycle
    ///@{


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{

    // return counter
    static SizeType GetCounter()
    {
        return Counter<TCountedType>::msCounter;
    }

    // return counter
    static SizeType Increment()
    {
        return ++Counter<TCountedType>::msCounter;
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return std::string("Counter");
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << Counter();
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

    /// Default constructor.
    Counter()
    {
    }

    /// Copy constructor.
    Counter(Counter<TCountedType> const& rOther)
    {
    }

    /// Destructor.
    virtual ~Counter()
    {
    }

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    /// Number of created objects.
    static SizeType msCounter;

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
    Counter& operator=(Counter const& rOther);

    ///@}

}; // Class Counter

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TCountedType>
inline std::istream& operator >> (std::istream& rIStream,
                                  Counter<TCountedType>& rThis);

/// output stream function
template<class TCountedType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Counter<TCountedType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : ";
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

// initialize counter with zero
template<typename TCountedType>
typename Counter<TCountedType>::SizeType Counter<TCountedType>::msCounter = 0;

}  // namespace Kratos.

#endif // KRATOS_COUNTER_H_INCLUDED  defined 


