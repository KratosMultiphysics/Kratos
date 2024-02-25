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

#if !defined(KRATOS_MESH_NODE_COLLAPSING_PROCESS_H_INCLUDED )
#define  KRATOS_MESH_NODE_COLLAPSING_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"


namespace Kratos
{
  ///@addtogroup Kratos Core
  ///@{

  ///@name Kratos Classes
  ///@{

  /// Remove the selected node from the mesh and collapse the connectivity around it.
  /** Detail class definition.
  */
  class MeshNodeCollapsingProcess : public Process
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of MeshNodeCollapsingProcess
      KRATOS_CLASS_POINTER_DEFINITION(MeshNodeCollapsingProcess);

	  ///@}
	  ///@name Flags
	  ///@{

	  KRATOS_DEFINE_LOCAL_FLAG(TO_COLLAPSE);

	  ///@}
      ///@name Life Cycle
      ///@{

	  ///Constructor to be used.
	  explicit MeshNodeCollapsingProcess(ModelPart& rModelPart);

	  /// Default constructor deleted.
	  MeshNodeCollapsingProcess() = delete;

	  /// Copy constructor deleted.
	  MeshNodeCollapsingProcess(MeshNodeCollapsingProcess const& rOther) = delete;

	  /// Destructor.
      ~MeshNodeCollapsingProcess() override;


      ///@}
      ///@name Operators
      ///@{

	  /// Assignment operator is deleted.
	  MeshNodeCollapsingProcess& operator=(MeshNodeCollapsingProcess const& rOther) = delete;

      ///@}
      ///@name Operations
      ///@{

	  void Execute() override;


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
      std::string Info() const override;

      /// Print information about this object.
      void PrintInfo(std::ostream& rOStream) const override;

      /// Print object's data.
      void PrintData(std::ostream& rOStream) const override;


      ///@}
      ///@name Friends
      ///@{


      ///@}

	protected:
		///@name Member Variables
		///@{
		ModelPart& mrModelPart;
		///@}

    private:
      ///@name Static Member Variables
      ///@{


      ///@}
      ///@name Member Variables
      ///@{


      ///@}
      ///@name Private Operations
      ///@{

		void CollapseNodes();

		void CollapseNode(Node& rThisNode);

		double CalculateQualityIfNodeCollapses(Node& rThisNode, Node const& rCoarseNode);

		double CalculateMinQualityOfNeighbourElements(Node& rThisNode, Node const& rCoarseNode);

		bool ElementHas(Element& rElement, Node const& rCoarseNode);

		void SwapElementNode(Element& rElement, Node const& rThisNode, Node::Pointer pCoarseNode);

      ///@}

    }; // Class MeshNodeCollapsingProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    MeshNodeCollapsingProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const MeshNodeCollapsingProcess& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MESH_NODE_COLLAPSING_PROCESS_H_INCLUDED  defined
