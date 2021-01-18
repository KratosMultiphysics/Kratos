//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
// kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#if !defined(KRATOS_MEMORY_INFO_H_INCLUDED)
#define KRATOS_MEMORY_INFO_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos {
///@addtogroup Kratos Core
///@{

///@name Kratos Classes
///@{

/// MemoryInfo gives the OS information about the memory useage by Kratos.
/** This class provides the peak memory usage and current memory usage.
 *  The information is taken bu OS and may vary in each execution depending
 *  on page allocation and other factors.
 *  The supported platforms are Windows, Linux and OSX.
*/
class KRATOS_API(KRATOS_CORE) MemoryInfo {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MemoryInfo
    KRATOS_CLASS_POINTER_DEFINITION(MemoryInfo);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MemoryInfo() = default;

    /// Destructor.
    virtual ~MemoryInfo() = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * This function returns the peak memory used by this
     * process in bytes. It returns zero for not supported
     * operating systems. The main idea is from:
     * http://nadeausoftware.com/articles/2012/07/c_c_tip_how_get_process_resident_set_size_physical_memory_use
     * rewritten/reordered to have independent code for each platform and some
     *styling applied
    **/
    static std::size_t GetPeakMemoryUsage();

    /**
     * This function returns the physical memory used currently by this
     * process in bytes. It returns zero for not supported operating
     * systems. The main idea is from:
     * http://nadeausoftware.com/articles/2012/07/c_c_tip_how_get_process_resident_set_size_physical_memory_use
     * rewritten for c++ and with style change
     **/
    static std::size_t GetCurrentMemoryUsage();

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const { return "MemoryInfo"; }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const { rOStream << Info(); }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
        rOStream << "Current Memory Usage : "
                 << HumanReadableSize(GetCurrentMemoryUsage()) << std::endl;
        rOStream << "Peak Memory Usage    : "
                 << HumanReadableSize(GetPeakMemoryUsage()) << std::endl;
    }

  ///@}

private:
    ///@name Private Operations
    ///@{

    std::string HumanReadableSize(std::size_t InBytes) const;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    MemoryInfo &operator=(MemoryInfo const &rOther);

    /// Copy constructor.
    MemoryInfo(MemoryInfo const &rOther);

    ///@}

}; // Class MemoryInfo

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream, MemoryInfo &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const MemoryInfo &rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_MEMORY_INFO_H_INCLUDED  defined
