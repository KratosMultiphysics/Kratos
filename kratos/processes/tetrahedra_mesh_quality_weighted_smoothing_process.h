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

#if !defined(KRATOS_TETRAHEDRA_MESH_QUALITY_WEIGHTED_SMOOTHING_PROCESS_H_INCLUDED )
#define  KRATOS_TETRAHEDRA_MESH_QUALITY_WEIGHTED_SMOOTHING_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "processes/tetrahedra_mesh_worst_element_smoothing_process.h"
#include "includes/model_part.h"



namespace Kratos
{
  ///@addtogroup Kratos Core
  ///@{

  ///@name Kratos Classes
  ///@{

  /// Short class definition.
  /** Detail class definition.
  */
  class KRATOS_API(KRATOS_CORE) TetrahedraMeshQualityWeightedSmoothingProcess : public TetrahedraMeshWorstElementSmoothingProcess
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of TetrahedraMeshQualityWeightedSmoothingProcess
      KRATOS_CLASS_POINTER_DEFINITION(TetrahedraMeshQualityWeightedSmoothingProcess);

	  typedef Node NodeType;

	  typedef GlobalPointersVector< Node > NeighboursVectorType;

	  typedef std::vector<Point > PointsVectorType;

	  ///@}
	  ///@name Flags
	  ///@{

	  ///@}
	  ///@name Life Cycle
	  ///@{

	  /// Constructor takes the modelpart to apply smoothing to its mesh 0.
      TetrahedraMeshQualityWeightedSmoothingProcess(ModelPart& rModelPart, double AptQuality = 0.75, std::size_t IterationsNumber=1);

      /// Destructor.
      ~TetrahedraMeshQualityWeightedSmoothingProcess() override;


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
      std::string Info() const override;


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


		void FindOptimumPositionsAndWeights(NodeType& rNode, PointsVectorType& rOptimumPoints, Vector& rWeights) override;

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
      TetrahedraMeshQualityWeightedSmoothingProcess& operator=(TetrahedraMeshQualityWeightedSmoothingProcess const& rOther);

      /// Copy constructor.
      TetrahedraMeshQualityWeightedSmoothingProcess(TetrahedraMeshQualityWeightedSmoothingProcess const& rOther);


      ///@}

    }; // Class TetrahedraMeshQualityWeightedSmoothingProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    TetrahedraMeshQualityWeightedSmoothingProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const TetrahedraMeshQualityWeightedSmoothingProcess& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_TETRAHEDRA_MESH_QUALITY_WEIGHTED_SMOOTHING_PROCESS_H_INCLUDED  defined
