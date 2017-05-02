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


#if !defined(KRATOS_CALCULATE_DISTANCE_TO_SKIN_PROCESS_H_INCLUDED )
#define  KRATOS_CALCULATE_DISTANCE_TO_SKIN_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"


namespace Kratos
{
  ///@addtogroup Kratos Core
  ///@{

  ///@name Kratos Classes
  ///@{

  /// Calculates the nodal distances using elemental discontinuous distances.
  /** This class calculates the nodal distances as a minimum elemental distances connected to it.
  */
  class CalculateDistanceToSkinProcess
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of CalculateDistanceToSkinProcess
      KRATOS_CLASS_POINTER_DEFINITION(CalculateDistanceToSkinProcess);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      CalculateDistanceToSkinProcess();

      /// Destructor.
      virtual ~CalculateDistanceToSkinProcess();


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
      CalculateDistanceToSkinProcess& operator=(CalculateDistanceToSkinProcess const& rOther);

      /// Copy constructor.
      CalculateDistanceToSkinProcess(CalculateDistanceToSkinProcess const& rOther);


      ///@}

    }; // Class CalculateDistanceToSkinProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    CalculateDistanceToSkinProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const CalculateDistanceToSkinProcess& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CALCULATE_DISTANCE_TO_SKIN_PROCESS_H_INCLUDED  defined
