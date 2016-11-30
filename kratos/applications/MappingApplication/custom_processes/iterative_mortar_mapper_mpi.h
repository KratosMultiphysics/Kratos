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


#if !defined(KRATOS_ITERATIVE_MORTAR_MAPPER_MPI_H_INCLUDED )
#define  KRATOS_ITERATIVE_MORTAR_MAPPER_MPI_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"


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
  class IterativeMortarMapperMPI : public IterativeMortarMapper
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of IterativeMortarMapperMPI
      KRATOS_CLASS_POINTER_DEFINITION(IterativeMortarMapperMPI);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      IterativeMortarMapperMPI();

      /// Destructor.
      virtual ~IterativeMortarMapperMPI();


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
      IterativeMortarMapperMPI& operator=(IterativeMortarMapperMPI const& rOther);

      /// Copy constructor.
      IterativeMortarMapperMPI(IterativeMortarMapperMPI const& rOther);


      ///@}

    }; // Class IterativeMortarMapperMPI

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    IterativeMortarMapperMPI& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const IterativeMortarMapperMPI& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ITERATIVE_MORTAR_MAPPER_MPI_H_INCLUDED  defined
