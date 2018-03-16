//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (based on work of Pablo Becker)
//

#if !defined(KRATOS_POINT_LOCATOR_H_INCLUDED)
#define  KRATOS_POINT_LOCATOR_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"


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
  class PointLocator
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of PointLocator
      KRATOS_CLASS_POINTER_DEFINITION(PointLocator);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      PointLocator(ModelPart& rModelPart) : mrModelPart(rModelPart) {}

      /// Destructor.
      virtual ~PointLocator() {}


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{
      bool Find(const Point& rThePoint);

      void InterpolateValue(const Variable<double>& rVariable, double& rValue);
    //   void InterpolateValue(const Variable<T>& rVariable, T& rValue);


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

      ModelPart& mrModelPart;
      bool mIsInitalized = false;
      Element::Pointer mpFoundElement;


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
      PointLocator& operator=(PointLocator const& rOther);

      /// Copy constructor.
      PointLocator(PointLocator const& rOther);


      ///@}

    }; // Class PointLocator

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    PointLocator& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const PointLocator& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_POINT_LOCATOR_H_INCLUDED  defined
