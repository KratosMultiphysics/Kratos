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

#if !defined(KRATOS_MEMORY_MANAGER_H_INCLUDED )
#define  KRATOS_MEMORY_MANAGER_H_INCLUDED



// System includes
#include <string>
#include <iostream>



// External includes


// Project includes
#include "includes/define.h"

// to be included after define.h which defines the KRATOS_COMPILED_IN Macros
#if defined(KRATOS_COMPILED_IN_LINUX)
#include <unistd.h>
#include <sys/resource.h>

#elif defined(KRATOS_COMPILED_IN_APPLIE)
#include <mach/mach.h>
#endif

#include <iostream>
#include <fstream>

namespace Kratos
{
  ///@addtogroup ApplicationNameApplication
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
  class MemoryManager
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of MemoryManager
      KRATOS_CLASS_POINTER_DEFINITION(MemoryManager);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      MemoryManager();

      /// Destructor.
      virtual ~MemoryManager();


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
     * rewritten for c++ and with style change 
     **/ 
    std::size_t GetPeakMemoryUsage()
    {
    
    #if defined(KRATOS_COMPILED_IN_LINUX)
         struct rusage resource_usage;
        getrusage( RUSAGE_SELF, &resource_usage );
    
        return resource_usage.ru_maxrss * 1024;
    #elif defined(KRATOS_COMPILED_IN_APPLIE)
        struct rusage resource_usage;
        getrusage( RUSAGE_SELF, &resource_usage );
    
        return resource_usage.ru_maxrss;
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
    std::size_t GetCurrentMemoryUsage()
    {
    #if defined(KRATOS_COMPILED_IN_LINUX)
        std::size_t program_size = 0;
        std::size_t resident_size = 0;
     
        std::ifstream process_file("/proc/self/statm");
    
        if(process_file.fail())
            return 0;
    
        process_file >> program_size;
        process_file >> resident_size;
      
        return resident_size * sysconf( _SC_PAGESIZE);
    #elif defined(KRATOS_COMPILED_IN_APPLIE)
        struct mach_task_basic_info info;
        mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
        
        if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
            (task_info_t)&info, &infoCount ) != KERN_SUCCESS )
            return 0;		
        return info.resident_size;
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
      virtual std::string Info() const;

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const;

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const;


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
      MemoryManager& operator=(MemoryManager const& rOther);

      /// Copy constructor.
      MemoryManager(MemoryManager const& rOther);


      ///@}

    }; // Class MemoryManager

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    MemoryManager& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const MemoryManager& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MEMORY_MANAGER_H_INCLUDED  defined
