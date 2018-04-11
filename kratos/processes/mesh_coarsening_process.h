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

#if !defined(KRATOS_MESH_COARSENING_PROCESS_H_INCLUDED )
#define  KRATOS_MESH_COARSENING_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/model_part.h"
#include "processes/mesh_node_collapsing_process.h"


namespace Kratos
{
  ///@addtogroup KratosCore
  ///@{

  ///@name Kratos Classes
  ///@{

  /// Short class definition.
  /** Detail class definition.
  */
  class KRATOS_API(KRATOS_CORE) MeshCoarseningProcess : public MeshNodeCollapsingProcess
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of MeshCoarseningProcess
      KRATOS_CLASS_POINTER_DEFINITION(MeshCoarseningProcess);

	  ///@}
	  ///@name Flags
	  ///@{

	  KRATOS_DEFINE_LOCAL_FLAG(COARSE_MESH_NODE);

	  ///@}
      ///@name Life Cycle
      ///@{

      /// Is not default constructable.
      MeshCoarseningProcess() = delete;

	  /// Is not copyable.
	  MeshCoarseningProcess(MeshCoarseningProcess const& rOther) = delete;

	  /// The constructor to be called
	  MeshCoarseningProcess(ModelPart& rModelPart);

      /// Destructor.
      ~MeshCoarseningProcess() override;


      ///@}
      ///@name Operators
      ///@{

	  /// Is not assignable.
	  MeshCoarseningProcess& operator=(MeshCoarseningProcess const& rOther) = delete;


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

		void SelectCoarseMeshNodes();


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

    }; // Class MeshCoarseningProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    MeshCoarseningProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const MeshCoarseningProcess& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MESH_COARSENING_PROCESS_H_INCLUDED  defined
