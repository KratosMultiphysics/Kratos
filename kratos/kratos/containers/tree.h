//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2008-03-12 14:20:36 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_TREE_H_INCLUDED )
#define  KRATOS_TREE_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <cmath>


// External includes 


// Project includes
//#include "includes/define.h"


namespace Kratos
{

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
	template<
		class TPointType,
		class TIteratorType,
		class TDistanceIteratorType>
	class TreeNode
	{
	public:
	    typedef std::size_t SizeType;
		typedef TPointType PointType;
		typedef double CoordinateType;
		typedef TIteratorType IteratorType;
		typedef TDistanceIteratorType DistanceIteratorType;

		virtual void PrintData(std::ostream& rOStream, std::string const& Perfix = std::string()) const{}
		TreeNode(){}
		virtual ~TreeNode(){}
		virtual IteratorType SearchNearestPoint(PointType const& ThisPoint, 
		  CoordinateType& rResultDistance)
	 {
		 // must be implemented in derived classes.
		 return msNull;
	 }

		virtual SizeType SearchInRadius(PointType const& ThisPoint, double Radius, IteratorType Results, 
			DistanceIteratorType ResultsDistances, SizeType MaximumNumberOfResults)
	 {
		 // must be implemented in derived classes.
		 return 0;
	 }

		virtual SizeType SearchInRadius(PointType const& ThisPoint, double Radius, double Radius2, IteratorType Results, 
			DistanceIteratorType ResultsDistances, SizeType MaximumNumberOfResults)
	 {
		 // must be implemented in derived classes.
		 return 0;
	 }

		virtual SizeType SearchInRadius(PointType const& ThisPoint, double Radius, IteratorType Results, 
						SizeType MaximumNumberOfResults)
	 {
		 // must be implemented in derived classes.
		 return 0;
	 }

		virtual SizeType SearchInRadius(PointType const& ThisPoint, double Radius, double Radius2, IteratorType Results, 
						SizeType MaximumNumberOfResults)
	 {
		 // must be implemented in derived classes.
		 return 0;
	 }

	 static IteratorType& NullIterator()
	 {
		 return msNull;
	 }

		//virtual void Search(PointType const& ThisPoint, SizeType NumberOfNearestPoints, 
		//  PointsContainerType& rResultPoints,  std::vector<double>& rResultDistances)
	 //{
		////  // must be implemented in derived classes.
	 //}
	private:
		static IteratorType msNull;


	};

	template<class TPointType, class TIteratorType, class TDistanceIteratorType> 
		typename TreeNode<TPointType, TIteratorType, TDistanceIteratorType>::IteratorType
		TreeNode<TPointType, TIteratorType, TDistanceIteratorType>::msNull;

  /// Short class definition.
  /** Detail class definition.
  */
  template<std::size_t TDimension, 
		class TPointType,
		class TPartitionType,
		class TLeafType,
		class TIteratorType,
		class TDistanceIteratorType>
  class Tree
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Tree
      //KRATOS_CLASS_POINTER_DEFINITION(Tree);

	  typedef TreeNode<TPointType, TIteratorType, TDistanceIteratorType> NodeType;
	
	  typedef TLeafType LeafType;

	  typedef TPartitionType PartitionType;
	  
	  typedef TPointType PointType;

	  typedef double CoordinateType;

	  typedef TIteratorType IteratorType;

	  typedef TDistanceIteratorType DistanceIteratorType;

	  typedef std::size_t SizeType;

	  typedef std::size_t IndexType;

  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Constructor.
	  Tree(IteratorType PointsBegin, IteratorType PointsEnd, SizeType BucketSize = 1) 
		  : mBucketSize(BucketSize), mPointsBegin(PointsBegin), mPointsEnd(PointsEnd)
	  {

		if(mPointsBegin == mPointsEnd)
		{
			mRoot = &msEmptyLeaf;
			return;
		}

		PointType max_point = **mPointsBegin;
		PointType min_point = **mPointsBegin;
		for(SizeType i = 0 ; i < TDimension ; i++)
			for(IteratorType point_iterator = mPointsBegin ;
				point_iterator != mPointsEnd ; point_iterator++)
			{
				if((**point_iterator)[i] > max_point[i]) 
					max_point[i] = (**point_iterator)[i];
				else if((**point_iterator)[i] < min_point[i]) 
					min_point[i] = (**point_iterator)[i];

			}
		mRoot = Construct(mPointsBegin, mPointsEnd, max_point, min_point);
	  }

      /// Destructor.
	  virtual ~Tree()
	  {
		  if(mRoot != &msEmptyLeaf)
		  	delete mRoot;
	  }
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{

	  IteratorType SearchNearestPoint(PointType const& ThisPoint, 
		CoordinateType& rResultDistance)
	  {
		 // searching the tree
		  IteratorType result = mRoot->SearchNearestPoint(ThisPoint,rResultDistance);
		
		  // Calculating the distance
		  rResultDistance = sqrt(rResultDistance);

		  return result;
	  }

	  IteratorType SearchNearestPoint(PointType const& ThisPoint)
	  {
			CoordinateType rResultDistance;

		 // searching the tree
		  IteratorType result = mRoot->SearchNearestPoint(ThisPoint,rResultDistance);

		  return result;
	  }

	  SizeType SearchInRadius(PointType const& ThisPoint, double Radius, IteratorType Results, 
			DistanceIteratorType ResultsDistances, SizeType MaximumNumberOfResults)
	 {
		 // Using the square of radius for avoiding square root calculation during search
		 double Radius2 = Radius * Radius;

		 // searching the tree
		 SizeType n = mRoot->SearchInRadius(ThisPoint, Radius, Radius2, Results, ResultsDistances, MaximumNumberOfResults);

		  // Calculating the distances
		 //for(SizeType i = 0 ; i < n ; i++)
		 //{
			 //*(ResultsDistances) = sqrt(*ResultsDistances);
			 //ResultsDistances++;
		 //}
		 return n;
	 }

	  SizeType SearchInRadius(PointType const& ThisPoint, double Radius, IteratorType Results, 
				  SizeType MaximumNumberOfResults)
	 {
		 // Using the square of radius for avoiding square root calculation during search
		 double Radius2 = Radius * Radius;

		 // searching the tree
		 return  mRoot->SearchInRadius(ThisPoint, Radius, Radius2, Results, MaximumNumberOfResults);
	 }


	  //void Search(PointType const& ThisPoint, SizeType NumberOfNearestPoints, 
		 // PointsContainerType& rResultPoints,  std::vector<double>& rResultDistances)
	  //{
		 // // mRoot->Search(ThisPoint,NumberOfNearestPoints,rResultPoints,rResultDistances)
	  //}
      
      
      ///@}
      ///@name Access
      ///@{

	  TPointType& BoundingBoxLowPoint()
	  {
		  return mBoundingBoxLowPoint;
	  }
      
	  TPointType& BoundingBoxHighPoint()
	  {
		  return mBoundingBoxHighPoint;
	  }
      
      
      ///@}
      ///@name Inquiry
      ///@{
      
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
	  virtual std::string Info() const
	  {
		  return "Tree";
	  }
      
      /// Print information about this object.
	  virtual void PrintInfo(std::ostream& rOStream) const
	  {
		  rOStream << "Tree";
	  }

      /// Print object's data.
	  virtual void PrintData(std::ostream& rOStream, std::string const& Perfix = std::string()) const
	  {
		  //for(typename IndicesType::const_iterator i = mIndices.begin() ; i != mIndices.end() ; i++)
			 // rOStream << "    " << (*mpPoints)[*i] << std::endl;
		  mRoot->PrintData(rOStream, "  ");
	  }
      
            
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
        
	  static LeafType msEmptyLeaf;	
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
        
	  SizeType mBucketSize;
        
      TPointType mBoundingBoxLowPoint;
	  TPointType mBoundingBoxHighPoint;

	  IteratorType mPointsBegin;
	  IteratorType mPointsEnd;

	  NodeType* mRoot;
        
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
        
	  NodeType* Construct(IteratorType PointsBegin, 
			IteratorType PointsEnd, 
			PointType HighPoint, 
			PointType LowPoint)
	  {
		  SizeType number_of_points = std::distance(PointsBegin,PointsEnd);
		  if (number_of_points == 0)
			  return &msEmptyLeaf;
		  else if (number_of_points <= mBucketSize)
			  return new LeafType(PointsBegin, PointsEnd); 
		  else 
		  {
			IndexType cutting_dimension;
			CoordinateType cutting_value;

			IteratorType partition = PartitionType::Partition(PointsBegin, PointsEnd, 
				cutting_dimension, cutting_value);

			PointType partition_high_point = HighPoint;
			PointType partition_low_point = LowPoint;

			partition_high_point[cutting_dimension] = cutting_value;
			partition_low_point[cutting_dimension] = cutting_value;

			return new PartitionType(cutting_dimension, cutting_value, 
				HighPoint[cutting_dimension], LowPoint[cutting_dimension],
				Construct(PointsBegin, partition, partition_high_point, LowPoint), 
				Construct(partition, PointsEnd, HighPoint, partition_low_point));						

		  }
	  }
        
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
      Tree& operator=(Tree const& rOther);

      /// Copy constructor.
      Tree(Tree const& rOther);

        
      ///@}    
        
    }; // Class Tree 
 
  template<std::size_t TDimension, 
		class TPointType,
		class TPartitionType,
		class TLeafType,
		class TIteratorType, 
		class TDistanceIteratorType>
		typename Tree<TDimension, TPointType, TPartitionType, TLeafType, TIteratorType, TDistanceIteratorType>::LeafType 
		Tree<TDimension, TPointType, TPartitionType, TLeafType, TIteratorType, TDistanceIteratorType>::msEmptyLeaf;

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template<std::size_t TDimension, 
		class TPointType,
		class TPartitionType,
		class TLeafType,
		class TIteratorType, 
		class TDistanceIteratorType>
  inline std::istream& operator >> (std::istream& rIStream, 
				    Tree<TDimension, TPointType, TPartitionType, TLeafType, TIteratorType, TDistanceIteratorType>& rThis);

  /// output stream function
  template<std::size_t TDimension, 
		class TPointType,
		class TPartitionType,
		class TLeafType,
		class TIteratorType, 
		class TDistanceIteratorType>
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Tree<TDimension, TPointType, TPartitionType, TLeafType, TIteratorType, TDistanceIteratorType>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_TREE_H_INCLUDED  defined 


