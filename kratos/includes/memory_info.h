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
#include <iostream>
#include <sstream>
#include <string>

// Project includes
#include "includes/define.h"

// to be included after define.h which defines the KRATOS_COMPILED_IN Macros
#if defined(KRATOS_COMPILED_IN_LINUX)
#include <sys/resource.h>
#include <unistd.h>

#elif defined(KRATOS_COMPILED_IN_OSX)
#include <mach/mach.h>

#elif defined(KRATOS_COMPILED_IN_WINDOWS)
#include <windows.h>
#include <psapi.h>
#endif

#include <fstream>
#include <iostream>

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
class MemoryInfo {
public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of MemoryInfo
    KRATOS_CLASS_POINTER_DEFINITION(MemoryInfo);

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  MemoryInfo() {}

  /// Destructor.
  virtual ~MemoryInfo() {}

  ///@}
  ///@name Operators
  ///@{

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
  static std::size_t GetPeakMemoryUsage() {

#if defined(KRATOS_COMPILED_IN_LINUX)
    struct rusage resource_usage;
    getrusage(RUSAGE_SELF, &resource_usage);

    return resource_usage.ru_maxrss * 1024;
#elif defined(KRATOS_COMPILED_IN_OSX)
    struct rusage resource_usage;
    getrusage(RUSAGE_SELF, &resource_usage);

    return resource_usage.ru_maxrss;
#elif defined(KRATOS_COMPILED_IN_WINDOWS)
  	PROCESS_MEMORY_COUNTERS memory_counter;
  	GetProcessMemoryInfo( GetCurrentProcess( ), &memory_counter, sizeof(memory_counter) );

  	return memory_counter.PeakWorkingSetSize;
#else
    return 0;
#endif
  }

  /**
   * This function returns the physical memory used currently by this
   * process in bytes. It returns zero for not supported operating
   * systems. The main idea is from:
   * http://nadeausoftware.com/articles/2012/07/c_c_tip_how_get_process_resident_set_size_physical_memory_use
   * rewritten for c++ and with style change
   **/
  static std::size_t GetCurrentMemoryUsage() {
#if defined(KRATOS_COMPILED_IN_LINUX)
    std::size_t program_size = 0;
    std::size_t resident_size = 0;

    std::ifstream process_file("/proc/self/statm");

    if (process_file.fail())
      return 0;

    process_file >> program_size;
    process_file >> resident_size;

    return resident_size * sysconf(_SC_PAGESIZE);
#elif defined(KRATOS_COMPILED_IN_OSX)
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;

    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info,
                  &infoCount) != KERN_SUCCESS)
      return 0;
    return info.resident_size;
#elif defined(KRATOS_COMPILED_IN_WINDOWS)
  	PROCESS_MEMORY_COUNTERS memory_counter;
  	GetProcessMemoryInfo( GetCurrentProcess( ), &memory_counter, sizeof(memory_counter) );

  	return memory_counter.WorkingSetSize;
#else
    return 0;
#endif
  }

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
  virtual std::string Info() const { return "MemoryInfo"; }

  /// Print information about this object.
  virtual void PrintInfo(std::ostream &rOStream) const { rOStream << Info(); }

  /// Print object's data.
  virtual void PrintData(std::ostream &rOStream) const {
    rOStream << "Current Memory Usage : "
             << HumanReadableSize(GetCurrentMemoryUsage()) << std::endl;
    rOStream << "Peak Memory Usage    : "
             << HumanReadableSize(GetPeakMemoryUsage()) << std::endl;
  }

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

  std::string HumanReadableSize(std::size_t InBytes) const {
    constexpr char extension[] = {'\0', 'K', 'M', 'G', 'T', 'P', 'E', 'E'};

    std::stringstream output;
    output.precision(4);

    double output_size = InBytes;
    int i = 0;
    for (; i < 7; i++) {
      if (output_size < 1024) {
          break;
      }
      output_size /= 1024;
    }
    output << output_size << " " << extension[i] << 'B';
    return output.str();
  }
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
