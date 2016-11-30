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
//


#if !defined(KRATOS_NEAREST_NEIGHBOR_MAPPER_MPI_H_INCLUDED )
#define  KRATOS_NEAREST_NEIGHBOR_MAPPER_MPI_H_INCLUDED



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
  class NearestNeighborMapperMPI : public NearestNeighborMapper
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of NearestNeighborMapperMPI
      KRATOS_CLASS_POINTER_DEFINITION(NearestNeighborMapperMPI);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      NearestNeighborMapperMPI();

      /// Destructor.
      virtual ~NearestNeighborMapperMPI();


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
      NearestNeighborMapperMPI& operator=(NearestNeighborMapperMPI const& rOther);

      /// Copy constructor.
      NearestNeighborMapperMPI(NearestNeighborMapperMPI const& rOther);


      ///@}

    }; // Class NearestNeighborMapperMPI

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    NearestNeighborMapperMPI& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const NearestNeighborMapperMPI& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NEAREST_NEIGHBOR_MAPPER_MPI_H_INCLUDED  defined
