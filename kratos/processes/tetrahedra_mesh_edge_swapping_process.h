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

#if !defined(KRATOS_TETRAHEDRA_MESH_EDGE_SWAPPING_PROCESS_H_INCLUDED )
#define  KRATOS_TETRAHEDRA_MESH_EDGE_SWAPPING_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <unordered_map>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "modeler/tetrahedra_edge_shell.h"

namespace Kratos
{
  ///@addtogroup Kratos Core
  ///@{

  ///@name Kratos Classes
  ///@{

  /// Short class definition.
  /** Detail class definition.
  */
  class KRATOS_API(KRATOS_CORE) TetrahedraMeshEdgeSwappingProcess : public Process
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of TetrahedraMeshEdgeSwappingProcess
      KRATOS_CLASS_POINTER_DEFINITION(TetrahedraMeshEdgeSwappingProcess);

      using PointType=Node;

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      TetrahedraMeshEdgeSwappingProcess() = delete;

      TetrahedraMeshEdgeSwappingProcess(ModelPart & rModelPart);

      /// Destructor.
      ~TetrahedraMeshEdgeSwappingProcess() override;


      ///@}
      ///@name Operators
      ///@{


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


      class Edge
      {
        PointType* mpPoint1;
        PointType* mpPoint2;
      public:
        Edge():mpPoint1(nullptr), mpPoint2(nullptr){}
        Edge(PointType* pPoint1, PointType* pPoint2): mpPoint1(pPoint1), mpPoint2(pPoint2){}

        bool operator == (Edge const& rOther) const {
          return ((mpPoint1 == rOther.mpPoint1) && (mpPoint2 == rOther.mpPoint2))
              || ((mpPoint1 == rOther.mpPoint2) && (mpPoint2 == rOther.mpPoint1));
        }

        std::size_t operator()(Edge const& TheEdge) const
          {
              std::size_t h1 = std::hash<PointType*>{}(TheEdge.mpPoint1);
              std::size_t h2 = std::hash<PointType*>{}(TheEdge.mpPoint2);
              return h1 + h2;
          }

          PointType const * GetPoint1() const {return mpPoint1;}
          PointType const * GetPoint2() const {return mpPoint2;}
      };

      using EdgesContainerType = std::unordered_map<Edge, TetrahedraEdgeShell, Edge>;


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
      TetrahedraMeshEdgeSwappingProcess& operator=(TetrahedraMeshEdgeSwappingProcess const& rOther);

      /// Copy constructor.
      TetrahedraMeshEdgeSwappingProcess(TetrahedraMeshEdgeSwappingProcess const& rOther);


      ///@}

    }; // Class TetrahedraMeshEdgeSwappingProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    TetrahedraMeshEdgeSwappingProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const TetrahedraMeshEdgeSwappingProcess& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_TETRAHEDRA_MESH_EDGE_SWAPPING_PROCESS_H_INCLUDED  defined
