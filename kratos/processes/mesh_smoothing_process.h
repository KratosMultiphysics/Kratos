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
	           
#if !defined(KRATOS_MESH_SMOOTHING_PROCESS_H_INCLUDED )
#define  KRATOS_MESH_SMOOTHING_PROCESS_H_INCLUDED



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
  ///@addtogroup Kratos
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
  class KRATOS_API(KRATOS_CORE) MeshSmoothingProcess : public Process
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of MeshSmoothingProcess
      KRATOS_CLASS_POINTER_DEFINITION(MeshSmoothingProcess);

	  typedef Node<3> NodeType;

	  typedef WeakPointerVector< Node<3> > NeighboursVectorType;

	  typedef std::vector<Point > PointsVectorType;
  
	  ///@}
	  ///@name Flags 
	  ///@{ 

	  KRATOS_DEFINE_LOCAL_FLAG(LAPLACIAN_SMOOTHING);
	  KRATOS_DEFINE_LOCAL_FLAG(EDGE_LENGTH_SMOOTHING);
	  KRATOS_DEFINE_LOCAL_FLAG(MOVEMENT_SMOOTHING);
	  KRATOS_DEFINE_LOCAL_FLAG(MULTI_LEVEL_SMOOTHING);
	  KRATOS_DEFINE_LOCAL_FLAG(COARSE_MESH_NODE);

	  ///@}
	  ///@name Life Cycle 
	  ///@{ 

	  /// Constructor takes the modelpart to apply smoothing to its mesh 0.
      MeshSmoothingProcess(ModelPart& rModelPart, Flags Options = LAPLACIAN_SMOOTHING, std::size_t IterationsNumber=10);

      /// Destructor.
      virtual ~MeshSmoothingProcess();
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      
	  virtual void Execute();
      
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

		ModelPart& mrModelPart;

		Flags mOptions;
		std::size_t mMaxIterationsNumber;
        
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{

		void PerformSmoothing();

		void PerformCoarsening(ModelPart::MeshType& rOriginalMesh, ModelPart::MeshType& rCoarsMesh);

		void SelectCoarseMeshNodes(ModelPart::MeshType& rOriginalMesh);

		void CollapseNodes(ModelPart::MeshType& rOriginalMesh);

		void CollapseToNearestCoarseNode(NodeType& rThisNode);

		void CreateCoarseMeshElements(ModelPart::MeshType& rOriginalMesh, ModelPart::MeshType& rCoarsMesh);

		void ChangeElmentCoarseNodes(Element& ThisElement);

		bool IsNotCollapsed(Element& ThisElement);

		void MoveNode(NodeType& rNode, NeighboursVectorType const& rNeighbours, PointsVectorType const& rOptimumPoints, Vector const& rWeights);

		void MeshSmoothingProcess::FindOptimumPositionsAndWeights(NodeType& rNode, NeighboursVectorType const& rNeighbours, PointsVectorType& rOptimumPoints, Vector& rWeights);

		void LaplacianSmoothingPositionsAndWeights(NodeType& rNode, NeighboursVectorType const& rNeighbours, PointsVectorType& rOptimumPoints, Vector& rWeights);

		void EdgeLengthSmoothingPositionsAndWeights(NodeType& rNode, NeighboursVectorType const& rNeighbours, PointsVectorType& rOptimumPoints, Vector& rWeights);

		void MovementSmoothingPositionsAndWeights(NodeType& rNode, NeighboursVectorType const& rNeighbours, PointsVectorType& rOptimumPoints, Vector& rWeights);

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
      MeshSmoothingProcess& operator=(MeshSmoothingProcess const& rOther);

      /// Copy constructor.
      MeshSmoothingProcess(MeshSmoothingProcess const& rOther);

        
      ///@}    
        
    }; // Class MeshSmoothingProcess 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    MeshSmoothingProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const MeshSmoothingProcess& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MESH_SMOOTHING_PROCESS_H_INCLUDED  defined 


