//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <iostream>

// Project includes
#include "includes/memory_info.h"

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

namespace Kratos
{
std::size_t MemoryInfo::GetPeakMemoryUsage() {
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

std::size_t MemoryInfo::GetCurrentMemoryUsage() {
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

std::string MemoryInfo::HumanReadableSize(std::size_t InBytes) const {
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

}  // namespace Kratos.
