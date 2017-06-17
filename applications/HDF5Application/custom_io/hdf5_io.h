//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:    Jordi Cotela
//

#if !defined(KRATOS_HDF5_IO_H_INCLUDED )
#define  KRATOS_HDF5_IO_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/io.h"


namespace Kratos
{
  ///@addtogroup HDF5Application
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

  /// Base class for HDF5 I/O.
  class HDF5IO: public IO
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of HDF5IO
    KRATOS_CLASS_POINTER_DEFINITION(HDF5IO);

    ///@}

    ///@name Life Cycle
    ///@{

    /// Default constructor.
    HDF5IO(std::string FileName, Flags Options = IO::WRITE);

    /// Destructor.
    virtual ~HDF5IO() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    virtual void WriteModelPart(ModelPart& rModelPart) override;

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
    virtual std::string Info() const /*override*/;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const /*override*/;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const /*override*/;

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

    const std::string mFileName;

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
    HDF5IO& operator=(HDF5IO const& rOther);

    /// Copy constructor.
    HDF5IO(HDF5IO const& rOther);

    ///@}
    }; // Class HDF5IO

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    HDF5IO& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const HDF5IO& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_HDF5_IO_H_INCLUDED  defined
