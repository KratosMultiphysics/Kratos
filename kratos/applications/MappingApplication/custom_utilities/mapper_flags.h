//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher

#if !defined(KRATOS_MAPPER_FLAGS_H_INCLUDED )
#define  KRATOS_MAPPER_FLAGS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "containers/flags.h"


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
  class MapperFlags
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of MapperFlags
      KRATOS_CLASS_POINTER_DEFINITION(MapperFlags);

      /// Local Flags
      KRATOS_DEFINE_LOCAL_FLAG( SWAP_SIGN );
      KRATOS_DEFINE_LOCAL_FLAG( ADD_VALUES );
      KRATOS_DEFINE_LOCAL_FLAG( CONSERVATIVE );
      KRATOS_DEFINE_LOCAL_FLAG( REMESHED );
      KRATOS_DEFINE_LOCAL_FLAG( NON_CONFORMING_INTERFACE );
      KRATOS_DEFINE_LOCAL_FLAG( INTERPOLATE_VALUES );

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      MapperFlags(){}

      /// Destructor.
      virtual ~MapperFlags(){}


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
      virtual std::string Info() const
      {
	       std::stringstream buffer;
         buffer << "MapperFlags" ;
         return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "MapperFlags";}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const {}


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
      MapperFlags& operator=(MapperFlags const& rOther);

    //   /// Copy constructor.
    //   MapperFlags(MapperFlags const& rOther){}


      ///@}

    }; // Class MapperFlags

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    MapperFlags& rThis)
    {
        return rIStream;
    }

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const MapperFlags& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MAPPER_FLAGS_H_INCLUDED  defined
