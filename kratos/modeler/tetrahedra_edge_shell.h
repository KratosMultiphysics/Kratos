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
//


#if !defined(KRATOS_TETRAHEDRA_EDGE_SHELL_H_INCLUDED )
#define  KRATOS_TETRAHEDRA_EDGE_SHELL_H_INCLUDED

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "geometries/geometry.h"


namespace Kratos
{
  ///@addtogroup KratosCore
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

  /// This class defines an edge shell for a mesh of tetrahedra.
  /** The definition of such a shell comes from the Paul L. George paper:
  Back to the edge flips in 3 dimentions
  http://www.imr.sandia.gov/papers/imr12/george03.pdf
  */
  class TetrahedraEdgeShell
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of TetrahedraEdgeShell
      KRATOS_CLASS_POINTER_DEFINITION(TetrahedraEdgeShell);

	  using PointType = Node<3>;
    using GeomertyType = Geometry<PointType>;

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor is deleted.
      TetrahedraEdgeShell() = delete;

	  
	  TetrahedraEdgeShell(PointType::Pointer pEdgePoint1, PointType::Pointer pEdgePoint2);

    TetrahedraEdgeShell(TetrahedraEdgeShell&& rOther) noexcept
		  :  mpPoint1(rOther.mpPoint1) , mpPoint2(rOther.mpPoint2), mShellPoints(rOther.mShellPoints), mTetrahedraEdge(rOther.mTetrahedraEdge), mIsClosed(rOther.mIsClosed)
	  {
	  }

      /// Destructor.
      virtual ~TetrahedraEdgeShell();


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      void AddElement(Element* pTheElement, char EdgeIndex);

	  void AddShellPoints();

	  double CalculateMinQuality(const GeomertyType::QualityCriteria QualityCriteria) const;


      ///@}
      ///@name Access
      ///@{

	  std::size_t GetNumberOfShellPoints() const {
		  return mShellPoints.size();
	  }

	  std::size_t GetNumberOfTetrahedra() const {
		  return mTetrahedraEdge.size();
	  }

	  Element* pGetElement(std::size_t ElementLocalIndex);

	  bool IsClosed() const { return mIsClosed; }

	  PointType::Pointer Point1() { return mpPoint1; }

	  PointType::Pointer Point2() { return mpPoint2; }

	  PointType::Pointer ShellPoint(std::size_t i) { return mShellPoints[i]; }



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

     private:
      ///@name Static Member Variables
      ///@{


      ///@}
      ///@name Member Variables
      ///@{

		 PointType::Pointer mpPoint1;
		 PointType::Pointer mpPoint2;

		 std::vector<PointType::Pointer> mShellPoints;
		 std::vector<std::pair<Element*, char> > mTetrahedraEdge;
		 bool mIsClosed;


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
      TetrahedraEdgeShell& operator=(TetrahedraEdgeShell const& rOther);

      /// Copy constructor.
      TetrahedraEdgeShell(TetrahedraEdgeShell const& rOther);


      ///@}

    }; // Class TetrahedraEdgeShell

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    TetrahedraEdgeShell& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const TetrahedraEdgeShell& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_TETRAHEDRA_EDGE_SHELL_H_INCLUDED  defined
