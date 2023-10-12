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

	  using PointType = Node;
    using GeomertyType = Geometry<PointType>;

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor is deleted.
      TetrahedraEdgeShell() = delete;

	  /// This constructor assumes that the neibour elements of the node are cacluclated
	  TetrahedraEdgeShell(PointType& EdgePoint1, PointType& EdgePoint2);

    TetrahedraEdgeShell(TetrahedraEdgeShell&& rOther) noexcept
		  :  mrPoint1(rOther.mrPoint1) , mrPoint2(rOther.mrPoint2), mShellPoints(rOther.mShellPoints), mTetrahedra(rOther.mTetrahedra)
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

      void AddTetrahedron(GeomertyType* TheTetrahedron);

      void AddShellPoints(PointType* pPoint1, PointType* pPoint2);


      ///@}
      ///@name Access
      ///@{

	  std::size_t GetNumberOfShellPoints() const {
		  return mShellPoints.size();
	  }

	  std::size_t GetNumberOfTetrahedra() const {
		  return mTetrahedra.size();
	  }


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

		 PointType const& mrPoint1;
		 PointType const& mrPoint2;

		 std::vector<std::pair<PointType*,PointType*> > mShellPoints;
		 std::vector<GeomertyType*> mTetrahedra;


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
