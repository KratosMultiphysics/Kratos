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

#if !defined(KRATOS_TETRAHEDRA_BALL_H_INCLUDED )
#define  KRATOS_TETRAHEDRA_BALL_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/node.h"
#include "geometries/tetrahedra_3d_4.h"



namespace Kratos
{
  ///@addtogroup KratosCore
  ///@{

  ///@name Kratos Classes
  ///@{

  /// Stores a ball of tetrahedra sourronding a node of mesh.
  /** This class contains the tetrahedra which share a node and make some ooperations over them.
  */
  class TetrahedraBall
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of TetrahedraBall
      KRATOS_CLASS_POINTER_DEFINITION(TetrahedraBall);

	  using NodeType = Node;

	  using TetrahedraContainerType = std::vector<Geometry<NodeType>* >;

      ///@}
      ///@name Life Cycle
      ///@{

	  TetrahedraBall() = delete;
	  
	  /// Copy constructor is deleted
	  TetrahedraBall(TetrahedraBall const& rOther) = delete;

	  /// Constructor which creates the ball for given node. 
	  /// This constructor uses the GetValue(NEIGHBOUR_ELEMENTS) of the node.
      TetrahedraBall(NodeType& rThisNode);

      /// Destructor.
	  virtual ~TetrahedraBall() {}


      ///@}
      ///@name Operators
      ///@{

	  /// Assignment operator is deleted
	  TetrahedraBall& operator=(TetrahedraBall const& rOther) = delete;

      ///@}
      ///@name Operations
      ///@{

	  std::size_t Size() const;

	  double CalculateMinQuality(const Geometry<NodeType>::QualityCriteria QualityCriteria) const;

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


    private:
      ///@name Static Member Variables
      ///@{


      ///@}
      ///@name Member Variables
      ///@{

		TetrahedraContainerType mTetrahedra;


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




      ///@}

    }; // Class TetrahedraBall

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    TetrahedraBall& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const TetrahedraBall& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_TETRAHEDRA_BALL_H_INCLUDED  defined
