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

      typedef Node<3> PointType;
      typedef typename PointType::CoordinatesArrayType CoordinatesArrayType;

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

      bool FindNode(const Point& rThePoint, int& rNodeId, double DistanceThreshold);

      bool FindElement(const Point& rThePoint, int& rObjectId, CoordinatesArrayType& rLocalCoordinates);

      bool FindCondition(const Point& rThePoint, int& rObjectId, CoordinatesArrayType& rLocalCoordinates);

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
          buffer << "PointLocator" ;
          return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "PointLocator";}

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

      ModelPart& mrModelPart;

      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{

      template<typename TObjectType>
      bool FindObject(const TObjectType& rObjects, const std::string& rObjectType,
                      const Point& rThePoint, int& rObjectId, CoordinatesArrayType& rLocalCoordinates);

      void CheckResults(const std::string& rObjectType,
                        const Point& rThePoint,
                        int GlobalObjectsFound);

      bool NodeIsCloseEnough(const Node<3>& rNode,
                             const Point& rThePoint,
                             double DistanceThreshold);

      ///@}
      ///@name Private  Access
      ///@{


      ///@}
      ///@name Private Inquiry
      ///@{


      ///@}
      ///@name Un accessible methods
      ///@{


      ///@}

    }; // Class PointLocator

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

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
