//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: clabra $
//   Date:                $Date: 2007-03-29 19:37:47 $
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
#include "search_structure.h"

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
        std::size_t TDimension,
		class TPointType,
        class TPointerType,
		class TIteratorType,
		class TDistanceIteratorType,
        class TIteratorIteratorType = typename std::vector<TIteratorType>::iterator
          >
	class TreeNode
	{
	public:

       // Global definitions
	    typedef std::size_t SizeType;
        typedef std::size_t IndexType;
		typedef double CoordinateType;

		typedef TPointType PointType;
        typedef TPointerType PointerType;
		typedef TIteratorType IteratorType;
		typedef TDistanceIteratorType DistanceIteratorType;

        typedef TreeNode<TDimension,TPointType,TPointerType,TIteratorType,TDistanceIteratorType> TreeNodeType;

        typedef typename std::vector<IteratorType>::iterator IteratorIteratorType;

        typedef SearchStructure<IndexType,SizeType,CoordinateType,TIteratorType,IteratorIteratorType,TDimension> SearchStructureType;

		virtual void PrintData(std::ostream& rOStream, std::string const& Perfix = std::string()) const{}
		
        TreeNode(){}
		
        virtual ~TreeNode(){}
		
        virtual void SearchNearestPoint(PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance) {}

		virtual void SearchNearestPoint(PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance,
            SearchStructureType& Auxiliar) {}

		virtual void SearchInRadius(PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results, 
			DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults)
		{
		  // must be implemented in derived classes.
		  return;
		}
		
        virtual void SearchInRadius(PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results, 
			DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType& Auxiliar)
		{
		  // must be implemented in derived classes.
		  return;
		}

		virtual void SearchInRadius(PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results, 
			SizeType& NumberOfResults, SizeType const& MaxNumberOfResults)
		{
		  // must be implemented in derived classes.
		  return;
		}

		virtual void SearchInRadius(PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results, 
			SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType& Auxiliar)
		{
		  // must be implemented in derived classes.
		  return;
		}

        virtual void SearchInBox(PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& Results, SizeType& NumberOfResults,
            SizeType const& MaxNumberOfResults )
        {
		  // must be implemented in derived classes.
		  return;
		}

		static IteratorType& NullIterator()
		{
		  return msNull;
		}

		static PointerType& NullPointer()
		{
		  return msNullPointer;
		}

		static TreeNode& NullLeaf()
		{
		  return msNullLeaf;
		}

	private:
		static IteratorType msNull;
		static PointerType msNullPointer;
		static TreeNode msNullLeaf;


	};

	template<std::size_t TDimension, class TPointType, class TPointerType, class TIteratorType, class TDistanceIteratorType, class TIteratorIteratorType> 
	   typename TreeNode<TDimension, TPointType, TPointerType, TIteratorType, TDistanceIteratorType, TIteratorIteratorType>::IteratorType
	   TreeNode<TDimension, TPointType, TPointerType, TIteratorType, TDistanceIteratorType, TIteratorIteratorType>::msNull;

	template<std::size_t TDimension, class TPointType, class TPointerType, class TIteratorType, class TDistanceIteratorType, class TIteratorIteratorType> 
	   typename TreeNode<TDimension, TPointType, TPointerType, TIteratorType, TDistanceIteratorType, TIteratorIteratorType>::PointerType
	   TreeNode<TDimension, TPointType, TPointerType, TIteratorType, TDistanceIteratorType, TIteratorIteratorType>::msNullPointer;


	template<std::size_t TDimension, class TPointType, class TPointerType, class TIteratorType, class TDistanceIteratorType, class TIteratorIteratorType> 
	   TreeNode<TDimension, TPointType, TPointerType, TIteratorType, TDistanceIteratorType, TIteratorIteratorType>
	   TreeNode<TDimension, TPointType, TPointerType, TIteratorType, TDistanceIteratorType, TIteratorIteratorType>::msNullLeaf;


  /// Short class definition.
  /** Detail class definition.
  */
  template< class TPartitionType >
  class Tree
    {
    public:
      ///@name Type Definitions
      ///@{

      class Partitions
      {
        public:
          Partitions( const std::size_t NumPartitions ) : mNumPartitions(NumPartitions) {}
          ~Partitions(){};
          std::size_t mNumPartitions;
      };

      /// Pointer definition of Tree
      KRATOS_CLASS_POINTER_DEFINITION(Tree);

	  typedef TPartitionType PartitionType;

	  typedef typename PartitionType::LeafType LeafType;

	  typedef typename PartitionType::PointType PointType;

	  typedef typename PartitionType::IteratorType IteratorType;

	  typedef typename PartitionType::DistanceIteratorType DistanceIteratorType;

      typedef typename PartitionType::PointerType PointerType;
      
	  typedef typename PartitionType::DistanceFunction DistanceFunction;

	  enum { Dimension = PartitionType::Dimension };

	  typedef TreeNode<Dimension,PointType,PointerType,IteratorType,DistanceIteratorType> NodeType;
	  
      typedef typename NodeType::CoordinateType CoordinateType;

      typedef typename NodeType::SizeType       SizeType;

      typedef typename NodeType::IndexType      IndexType;

      //typedef typename NodeType::SearchStructureType SearchStructureType;
      typedef typename PartitionType::SearchStructureType SearchStructureType;

      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Constructor.
	  Tree(IteratorType PointsBegin, IteratorType PointsEnd, SizeType BucketSize = 1) 
		  : mBucketSize(BucketSize), mPointsBegin(PointsBegin), mPointsEnd(PointsEnd)
	  {

		if(mPointsBegin == mPointsEnd)
			return;

		PointType max_point = **mPointsBegin;
		PointType min_point = **mPointsBegin;
        for(IteratorType point_iterator = mPointsBegin ; point_iterator != mPointsEnd ; point_iterator++)
          for(SizeType i = 0 ; i < Dimension ; i++)
          {
            if((**point_iterator)[i] > max_point[i]) 
              max_point[i] = (**point_iterator)[i];
            else if((**point_iterator)[i] < min_point[i]) 
              min_point[i] = (**point_iterator)[i];
          }

		mRoot = TPartitionType::Construct(mPointsBegin, mPointsEnd, max_point, min_point, mBucketSize);
	  }
	  
      Tree(IteratorType PointsBegin, IteratorType PointsEnd, Partitions Parts ) 
		  : mPointsBegin(PointsBegin), mPointsEnd(PointsEnd)
	  {

		if(mPointsBegin == mPointsEnd)
			return;

        SizeType NumPoints = SearchUtils::PointerDistance(mPointsBegin,mPointsEnd);
        mBucketSize = static_cast<std::size_t>( (double) NumPoints / (double) Parts.mNumPartitions ) + 1;

		PointType max_point = **mPointsBegin;
		PointType min_point = **mPointsBegin;
		for(IteratorType point_iterator = mPointsBegin ; point_iterator != mPointsEnd ; point_iterator++)
		   for(SizeType i = 0 ; i < Dimension ; i++)
		   {
			  if((**point_iterator)[i] > max_point[i]) 
				 max_point[i] = (**point_iterator)[i];
			  else if((**point_iterator)[i] < min_point[i]) 
				 min_point[i] = (**point_iterator)[i];
		   }

		mRoot = TPartitionType::Construct(mPointsBegin, mPointsEnd, max_point, min_point, mBucketSize);
	  }

      /// Destructor.
	  virtual ~Tree()
	  {
		  delete mRoot;
	  }
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{

      PointerType ExistPoint( PointerType const& ThisPoint, CoordinateType const Tolerance = static_cast<CoordinateType>(10.0*DBL_EPSILON) )
      {
        PointerType Result = *mPointsBegin;
        CoordinateType ResultDistance = static_cast<CoordinateType>(DBL_MAX);
		// searching the tree
		mRoot->SearchNearestPoint(ThisPoint,Result,ResultDistance);
        if (ResultDistance<Tolerance*Tolerance)
          return Result;
        return this->NullPointer();
      }

	  PointerType SearchNearestPoint(PointType const& ThisPoint, CoordinateType& rResultDistance)
	  {
		PointerType Result = *mPointsBegin;
		rResultDistance = static_cast<CoordinateType>(DBL_MAX); // DistanceFunction()(ThisPoint,**mPointsBegin);

		// searching the tree
		mRoot->SearchNearestPoint(ThisPoint,Result,rResultDistance);

		return Result;
	  }

	  PointerType SearchNearestPoint(PointType const& ThisPoint)
	  {
		PointerType Result = *mPointsBegin; // NULL ??
		CoordinateType rResultDistance = static_cast<CoordinateType>(DBL_MAX); // DistanceFunction()(ThisPoint,**mPointsBegin);

		// searching the tree
		mRoot->SearchNearestPoint(ThisPoint,Result,rResultDistance);

		return Result;
	  }

	  SizeType SearchInRadius(PointType const& ThisPoint, CoordinateType Radius, IteratorType Results, 
							  DistanceIteratorType ResultsDistances, SizeType MaxNumberOfResults)
	  {
		// Using the square of radius for avoiding square root calculation during search
		CoordinateType Radius2 = Radius * Radius;

		// searching the tree
        SizeType NumberOfResults = 0;
		mRoot->SearchInRadius(ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults);

		return NumberOfResults;
	  }

	  SizeType SearchInRadius(PointType const& ThisPoint, CoordinateType Radius, IteratorType Results, SizeType MaxNumberOfResults)
	  {
		// Using the square of radius for avoiding square root calculation during search
		CoordinateType Radius2 = Radius * Radius;

		// searching the tree
        SizeType NumberOfResults = 0;
		mRoot->SearchInRadius(ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults);
        return NumberOfResults;
	  }
      
      SizeType SearchInBox(PointType const& MinPointBox, PointType const& MaxPointBox, IteratorType Results, SizeType MaxNumberOfResults )
      {
        SizeType NumberOfResults = 0;
        mRoot->SearchInBox(MinPointBox,MaxPointBox,Results,NumberOfResults,MaxNumberOfResults);
        return NumberOfResults;
      }
     

      ///@}
      ///@name Access
      ///@{

	  PointType& BoundingBoxLowPoint()
	  {
		  return mBoundingBoxLowPoint;
	  }
      
	  PointType& BoundingBoxHighPoint()
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
        
      PointType mBoundingBoxLowPoint;
	  PointType mBoundingBoxHighPoint;

	  IteratorType mPointsBegin;
	  IteratorType mPointsEnd;

	  NodeType* mRoot;
        
      ///@} 
      ///@name Private Operators:
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
      Tree& operator=(Tree const& rOther);

      /// Copy constructor.
      Tree(Tree const& rOther);

        
      ///@}    
        
    }; // Class Tree 
 
  template< class TPartitionType >
	typename Tree<TPartitionType>::LeafType Tree<TPartitionType>::msEmptyLeaf;

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template<class TPartitionType>
  inline std::istream& operator >> (std::istream& rIStream, Tree<TPartitionType>& rThis);

  /// output stream function
  template<class TPartitionType>
  inline std::ostream& operator << (std::ostream& rOStream, const Tree<TPartitionType>& rThis)
  {
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
  }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_TREE_H_INCLUDED  defined 


