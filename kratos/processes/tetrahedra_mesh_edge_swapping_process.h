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
  class TetrahedraMeshEdgeSwappingProcess : public Process
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of TetrahedraMeshEdgeSwappingProcess
      KRATOS_CLASS_POINTER_DEFINITION(TetrahedraMeshEdgeSwappingProcess);

      using PointType=Node<3>;

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      TetrahedraMeshEdgeSwappingProcess() = delete;

      TetrahedraMeshEdgeSwappingProcess(ModelPart & rModelPart);

      /// Destructor.
      virtual ~TetrahedraMeshEdgeSwappingProcess();


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
      virtual std::string Info() const override;

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const override;

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const override;


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

	  EdgesContainerType mEdges;

      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{

	  template<typename TEdgeSwappingCasesType>
	  void EdgeSwapping(TetrahedraEdgeShell & EdgeShell) {
		  TEdgeSwappingCasesType SwappingCases;
		  auto quality_criteria = Geometry<Node<3> >::QualityCriteria::VOLUME_TO_EDGE_LENGTH;
		  double original_min_quality = EdgeShell.CalculateMinQuality(quality_criteria);
		  Tetrahedra3D4<Node<3>> tetrahedra_1 = EdgeShell.pGetElement(0)->GetGeometry(); // This initialization is to avoid creating a dummy
		  Tetrahedra3D4<Node<3>> tetrahedra_2 = EdgeShell.pGetElement(0)->GetGeometry(); // It will be reinitialized afterward
		  typename const TEdgeSwappingCasesType::EdgeSwappingCaseType* p_best_case = nullptr;
		  double max_cases_quality = original_min_quality;

		  for (auto i_case = SwappingCases.GetCases().begin(); i_case != SwappingCases.GetCases().end(); i_case++) {
			  double case_min_quality = std::numeric_limits<double>::max();
			  if (i_case->GetMinQuality() > original_min_quality) {// There are no previously calculated tetrahedra with worse quality
				  for (std::size_t i = 0; i < SwappingCases.NumberOfTrianglesPerCase(); i++) {
					  SwappingCases.SetTetrahedraForCase(*i_case, i, EdgeShell, tetrahedra_1, tetrahedra_2);
					  double min_quality = std::min(tetrahedra_1.Quality(quality_criteria), tetrahedra_2.Quality(quality_criteria));
					  if (tetrahedra_1.Volume() < std::numeric_limits<double>::epsilon()) // workaround for issue #126
						  min_quality = -1000.00;
					  if (tetrahedra_2.Volume() < std::numeric_limits<double>::epsilon()) // workaround for issue #126
						  min_quality = -1000.00;
					  if (min_quality < case_min_quality) {
						  case_min_quality = min_quality;
					  }
				  }
				  if (case_min_quality > max_cases_quality) {
					  if (mrModelPart.NumberOfElements() == 5)
					  {
						  KRATOS_WATCH(tetrahedra_1);
						  KRATOS_WATCH(tetrahedra_2);
						  KRATOS_WATCH(tetrahedra_1.Quality(quality_criteria));
						  KRATOS_WATCH(tetrahedra_2.Quality(quality_criteria));
						  KRATOS_WATCH(case_min_quality);
						  KRATOS_WATCH(max_cases_quality);
					  }
					  p_best_case = &(*i_case);
					  max_cases_quality = case_min_quality;
					  // Todo: break if apt quality reached.
				  }
			  }
		  }
		  if (max_cases_quality > original_min_quality + std::numeric_limits<double>::epsilon()) {
			  
			  for (std::size_t i = 0; i < SwappingCases.NumberOfTrianglesPerCase(); i++) {
				  SwappingCases.SetTetrahedraForCase(*p_best_case, i, EdgeShell, tetrahedra_1, tetrahedra_2);
				  if (2 * i < EdgeShell.GetNumberOfTetrahedra()) {
					  if(i == 3)
						KRATOS_WATCH(tetrahedra_1);
					  EdgeShell.pGetElement(2 * i)->GetGeometry() = tetrahedra_1;
					  EdgeShell.pGetElement(2 * i)->Set(MODIFIED);
				  }
				  else {
					  KRATOS_WATCH(mrModelPart.NumberOfElements() + 1);
					  KRATOS_WATCH(max_cases_quality);
					  KRATOS_WATCH(original_min_quality);
					  KRATOS_WATCH(tetrahedra_1);
					  mrModelPart.AddElement(EdgeShell.pGetElement(0)->Clone(mrModelPart.NumberOfElements() + 1, tetrahedra_1));
					  KRATOS_WATCH(mrModelPart.GetElement(mrModelPart.NumberOfElements()));
				  }
				  if ((2 * i) + 1 < EdgeShell.GetNumberOfTetrahedra()){
					EdgeShell.pGetElement((2 * i) + 1)->GetGeometry() = tetrahedra_2;
					EdgeShell.pGetElement((2 * i) + 1)->Set(MODIFIED);
				  }
				  else {
					  mrModelPart.AddElement(EdgeShell.pGetElement(0)->Clone(mrModelPart.NumberOfElements() + 1, tetrahedra_2));
				  }
			  }
		  }
	  }

	  void EdgeSwapping3(TetrahedraEdgeShell & EdgeShell);
	  void EdgeSwapping4(TetrahedraEdgeShell & EdgeShell);

	  double TetrahedraMeshEdgeSwappingProcess::CalculateMinimumQualityForAShellTriangle(TetrahedraEdgeShell & EdgeShell, std::array<std::size_t, 3> const& TriangleConnectivity);

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
