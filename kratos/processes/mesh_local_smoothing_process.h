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
	           
#if !defined(KRATOS_MESH_LOCAL_SMOOTHING_PROCESS_H_INCLUDED )
#define  KRATOS_MESH_LOCAL_SMOOTHING_PROCESS_H_INCLUDED



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
  
  /// The base class for local smoothing processes providing a laplacian smoothing.
  /** This class asks for optimum position of each node given by its neighbour elements
      and their corresponding weights to calculate the node position as described in the book
	  "Mesh Generation" of Pascal Frey and Paul-Louis George
  */
  class KRATOS_API(KRATOS_CORE) MeshLocalSmoothingProcess : public Process
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of MeshLocalSmoothingProcess
      KRATOS_CLASS_POINTER_DEFINITION(MeshLocalSmoothingProcess);

	  typedef Node<3> NodeType;

	  typedef WeakPointerVector< Node<3> > NeighboursVectorType;

	  typedef std::vector<Point > PointsVectorType;
  
	  ///@}
	  ///@name Flags 
	  ///@{ 

	  ///@}
	  ///@name Life Cycle 
	  ///@{ 

	  /// Constructor takes the modelpart to apply smoothing to its mesh 0.
	  MeshLocalSmoothingProcess(
          ModelPart &rModelPart,
          double AptQuality = 0.5,
          std::size_t MaxIterationsNumber = 10,
          const Flags &rBoundaryFlag = BOUNDARY);

      /// Destructor.
      ~MeshLocalSmoothingProcess() override;
      

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

		virtual void FindOptimumPositionsAndWeights(NodeType& rNode, PointsVectorType& rOptimumPoints, Vector& rWeights);

        
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

		std::size_t mMaxIterationsNumber;

		double mAptQuality;

		std::size_t mNumberOfLowQualityElements;

		double mMeshMinQuality;

		double mMeshQualityNorm;

        const Flags &mrBoundaryFlag;
        
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{

		void SelectLowQualityElementNodes();

		void PerformSmoothing();

		void InterpolateNodeOptimumPosition(PointsVectorType const& rOptimumPoints, Vector const& rWeights, Point& OptimumPosition);

		void MoveNodeIfImprovesMinimumQuality(NodeType& rNode, Point const& OptimumPosition);


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
      MeshLocalSmoothingProcess& operator=(MeshLocalSmoothingProcess const& rOther);

      /// Copy constructor.
      MeshLocalSmoothingProcess(MeshLocalSmoothingProcess const& rOther);

        
      ///@}    
        
    }; // Class MeshLocalSmoothingProcess 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    MeshLocalSmoothingProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const MeshLocalSmoothingProcess& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MESH_LOCAL_SMOOTHING_PROCESS_H_INCLUDED  defined 


