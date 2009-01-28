//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: antonia $
//   Date:                $Date: 2008-10-10 14:04:56 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_OCT_TREE_H_INCLUDED )
#define  KRATOS_OCT_TREE_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <cstddef>
#include <vector>


// External includes 


// Project includes
#include "tree.h"

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
	template< class TLeafType >
	class OCTreePartitionBase : public TreeNode< TLeafType::Dimension,
                                                 typename TLeafType::PointType,
	                                             typename TLeafType::PointerType,
											     typename TLeafType::IteratorType,
											     typename TLeafType::DistanceIteratorType >
	{
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of KDTree
      KRATOS_CLASS_POINTER_DEFINITION(OCTreePartionBase);

	  typedef TLeafType LeafType;

	  typedef typename LeafType::PointType PointType;

      typedef typename LeafType::ContainerType ContainerType;

  	  typedef typename LeafType::IteratorType IteratorType;

	  typedef typename LeafType::DistanceIteratorType DistanceIteratorType;

      typedef typename LeafType::PointerType PointerType;

      typedef typename LeafType::DistanceFunction DistanceFunction;

	  enum { Dimension = LeafType::Dimension };
	  
	  typedef TreeNode<Dimension, PointType, PointerType, IteratorType, DistanceIteratorType> TreeNodeType;

	  typedef typename TreeNodeType::CoordinateType CoordinateType;

	  typedef typename TreeNodeType::SizeType SizeType;

	  typedef typename TreeNodeType::IndexType IndexType;

	  static const SizeType number_of_childs = 1 << Dimension;

      typedef typename TreeNodeType::SearchStructure SearchStructure;

      ///@}
      ///@name Life Cycle 
      ///@{ 


      /// Constructor.
      // Defined in derived class
      OCTreePartitionBase(){}

      /// Destructor.
	  virtual ~OCTreePartitionBase()
	  {
		/* 		SizeType number_of_childs = 8; */
		for(SizeType i = 0 ; i < number_of_childs ; i++)
		  delete mpChilds[i];
	  }
	  
      virtual void PrintData(std::ostream& rOStream, std::string const& Perfix = std::string()) const
	  {
	      rOStream << Perfix << "Partition at point (" << mPosition[0];
		  for(IndexType j = 0 ; j < Dimension - 1 ; j++)
			 rOStream << "," << mPosition[j];
	      rOStream << std::endl;

		  for(SizeType j = 0 ; j < number_of_childs ; j++)
			 mpChilds[j]->PrintData(rOStream, Perfix + "  ");

	  }
      
      ///@}
      ///@name Operations
      ///@{
	  
	    void SearchNearestPoint(PointType const& rThisPoint, PointerType Result)
	    {
		   CoordinateType rResultDistance;
		   return SearchNearestPoint(rThisPoint, Result, rResultDistance);
		}
		
		void SearchNearestPoint(PointType const& ThisPoint, PointerType& Result, CoordinateType& rResultDistance)
		{
		   CoordinateType distances_to_partitions[number_of_childs];
		   CoordinateType ResidualDistance = 0.0;
		   
		   SizeType child_index = GetChildIndex(ThisPoint);

		   mpChilds[child_index]->SearchNearestPoint(ThisPoint, Result, rResultDistance, ResidualDistance);
		   
		   DistanceToPartitions(child_index, ThisPoint, distances_to_partitions, ResidualDistance);
		   
		   for(SizeType i = 0 ; i < number_of_childs ; i++)
			  if((i != child_index) && (distances_to_partitions[i] < rResultDistance)){
				 mpChilds[i]->SearchNearestPoint(ThisPoint, Result, rResultDistance, distances_to_partitions[i]);
			  }
		
	    }

		void SearchNearestPoint(PointType const& ThisPoint, PointerType& Result, CoordinateType& rResultDistance, CoordinateType& ResidualDistance)
		{
		  CoordinateType distances_to_partitions[number_of_childs];

		  SizeType child_index = GetChildIndex(ThisPoint);

		  mpChilds[child_index]->SearchNearestPoint(ThisPoint, Result, rResultDistance, ResidualDistance);

		  DistanceToPartitions(child_index, ThisPoint, distances_to_partitions, ResidualDistance);

		  for(SizeType i = 0 ; i < number_of_childs ; i++)
			if((i != child_index) && (distances_to_partitions[i] < rResultDistance)){
			  mpChilds[i]->SearchNearestPoint(ThisPoint, Result, rResultDistance, distances_to_partitions[i]);
			}

		}

		SizeType SearchInRadius(PointType const& ThisPoint, CoordinateType Radius, CoordinateType Radius2, IteratorType Results, DistanceIteratorType ResultsDistances, SizeType MaxNumberOfResults)
		{
		   //CoordinateType ResidualDistance = 0.0;
           SearchStructure Auxiliar;
		   return SearchInRadius(ThisPoint, Radius, Radius2, Results, ResultsDistances, MaxNumberOfResults, Auxiliar );
		}

		//SizeType SearchInRadius(PointType const& ThisPoint, CoordinateType Radius, CoordinateType Radius2, IteratorType Results, DistanceIteratorType ResultsDistances, SizeType MaxNumberOfResults, CoordinateType& ResidualDistance )
		SizeType SearchInRadius(PointType const& ThisPoint, CoordinateType Radius, CoordinateType Radius2, IteratorType Results, DistanceIteratorType ResultsDistances, SizeType MaxNumberOfResults, SearchStructure& Auxiliar )
		{

		   SizeType child_index = GetChildIndex(ThisPoint);
		   
		   // n is number of points found
		   SizeType n = mpChilds[child_index]->SearchInRadius(ThisPoint, Radius, Radius2, Results, ResultsDistances, MaxNumberOfResults, Auxiliar);

		   CoordinateType distances_to_partitions[number_of_childs];

		   DistanceToPartitions(child_index, ThisPoint, distances_to_partitions, Auxiliar.residual_distance[0]);

		   for(IndexType i = 0 ; i < number_of_childs ; i++)
			  if((i != child_index) && (distances_to_partitions[i] < Radius2))
				 n += mpChilds[i]->SearchInRadius(ThisPoint, Radius, Radius2, Results+n, ResultsDistances+n, MaxNumberOfResults-n, Auxiliar);

		   return n;
		}
      
		SizeType SearchInRadius(PointType const& ThisPoint, CoordinateType Radius, CoordinateType Radius2, IteratorType Results, SizeType MaxNumberOfResults)
		{
		   CoordinateType ResidualDistance = 0.0;
		   return SearchInRadius(ThisPoint, Radius, Radius2, Results, MaxNumberOfResults, ResidualDistance);
		}

		SizeType SearchInRadius(PointType const& ThisPoint, CoordinateType Radius, CoordinateType Radius2, IteratorType Results, SizeType MaxNumberOfResults, CoordinateType& ResidualDistance)
		{

		   SizeType child_index = GetChildIndex(ThisPoint);

		   SizeType n = mpChilds[child_index]->SearchInRadius(ThisPoint, Radius, Radius2, Results, MaxNumberOfResults, ResidualDistance); // n is number of points found
		   CoordinateType distances_to_partitions[number_of_childs];

		   DistanceToPartitions(child_index, ThisPoint, distances_to_partitions, ResidualDistance);

		   for(IndexType i = 0 ; i < number_of_childs ; i++)
			  if((i != child_index) && (distances_to_partitions[i] < Radius2))
				 n += mpChilds[i]->SearchInRadius(ThisPoint, Radius, Radius2, Results+n, MaxNumberOfResults-n, distances_to_partitions[i]);

		   return n;
		}


	private:

	    IndexType GetChildIndex(PointType const& rThisPoint) const
	    {
		   // Calculating the cell index 
		   IndexType child_index = 0;
		   const IndexType dim_mask[] = { 1, 2, 4 };

		   for( IndexType i = 0 ; i < Dimension ; i++)
 			  if( rThisPoint[i] >= mPosition[i] )
                child_index += dim_mask[i];

		   return child_index;
		}

	    void DistanceToPartitions(IndexType ContainingChildIndex, PointType const& rThisPoint, CoordinateType rDistances[], CoordinateType& ResidualDistance) const
		{
		    const IndexType coordinate_mask[] = { 1, 2, 4 };

		    CoordinateType offset_from_postition[Dimension];
		
		    for(IndexType j = 0 ; j < Dimension ; j++)
		    {
			   CoordinateType temp = rThisPoint[j] - mPosition[j];
			   offset_from_postition[j] = temp*temp;
			}
			
			for(IndexType i = 0 ; i < number_of_childs ; i++)
			{
			   rDistances[i] = ResidualDistance;
			   
			   IndexType partitions = ContainingChildIndex^i;
			   
               for(IndexType j = 0 ; j < Dimension ; j++)
                 if(coordinate_mask[j] & partitions)
                   rDistances[i] += offset_from_postition[j];
			}
		    
		}
	   
	private:

		IndexType mCutingDimension;
		PointType mPosition; // Position of partition
		TreeNodeType* mpChilds[number_of_childs];  // 8 is number of childs

	public:
/*        
        static TreeNodeType* Construct(IteratorType PointsBegin, IteratorType PointsEnd, 
            PointType HighPoint, PointType LowPoint, SizeType BucketSize)
        {
          SizeType number_of_points = std::distance(PointsBegin,PointsEnd);
          if (number_of_points == 0)
            return NULL;
          else if (number_of_points <= BucketSize)
          {
            return new LeafType(PointsBegin, PointsEnd); 
          }
          else 
          {
            return new OCTreePartition(PointsBegin, PointsEnd, LowPoint, HighPoint, BucketSize);
          }
        }
*/
	};

  
    /// Short class definition.
  /** Detail class definition.
  */

    template< class TLeafType >
      class OCTreePartitionAverageSplit : public OCTreePartitionBase<TLeafType>
    {
      public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of KDTree
        KRATOS_CLASS_POINTER_DEFINITION(OCTreePartionAverageSplit);

        //typedef OCTreePartitionBase<TLeafType>  BaseType;
        //typedef typename BaseType::IteratorType IteratorType;

        typedef TLeafType LeafType;
        typedef typename LeafType::PointType PointType;
        typedef typename LeafType::ContainerType ContainerType;
        typedef typename LeafType::IteratorType IteratorType;
        typedef typename LeafType::DistanceIteratorType DistanceIteratorType;
        typedef typename LeafType::PointerType PointerType;
        typedef typename LeafType::DistanceFunction DistanceFunction;
        enum { Dimension = LeafType::Dimension };
        typedef TreeNode<Dimension, PointType, PointerType, IteratorType, DistanceIteratorType> TreeNodeType;
        typedef typename TreeNodeType::CoordinateType CoordinateType;
        typedef typename TreeNodeType::SizeType SizeType;
        typedef typename TreeNodeType::IndexType IndexType;

        static const SizeType number_of_childs = 1 << Dimension;

        /// Partition constructor.
        OCTreePartitionAverageSplit(IteratorType PointsBegin, IteratorType PointsEnd, 
            PointType const& MinPoint, PointType const& MaxPoint,  SizeType BucketSize = 1)
        {
          PointType mid_cell_lenght;

          SizeType TempSize = std::distance(PointsBegin,PointsEnd);
          PointerType* Temp = new PointerType[ TempSize ];

          // Template definition of SplitMode
          Split(this->mPosition,MinPoint,MaxPoint,PointsBegin,PointsEnd);

          SizeType cell_sizes[number_of_childs];
          IteratorType cell_position[number_of_childs];

          for(IndexType i = 0 ; i < number_of_childs ; i++)
            cell_sizes[i] = 0;

          PointerType* i_temp = Temp;
          for(IteratorType i_point = PointsBegin ; i_point < PointsEnd ; i_point++, i_temp++)
          {
            *i_temp = *i_point;
            IndexType child_index = GetChildIndex(**i_point);
            cell_sizes[child_index]++;
          }

          cell_position[0] = PointsBegin;
          for(IndexType i = 1 ; i < number_of_childs ; i++)
          {
            cell_position[i] = cell_position[i-1];
            std::advance(cell_position[i], cell_sizes[i-1]);
          }

          // 	      for(IteratorType i_point = PointsBegin ; i_point < PointsEnd ; i_point++)
          // 	      {
          // 			 IndexType child_index = GetChildIndex(**i_point);
          // 			 *(cell_position[child_index]++) = *i_point;
          // 			 //swap(*i_point, *(cell_position[child_index]++));
          // 		  }
          //  The same in Bin_Static ( use temporal array for reorder the list )
          for(PointerType* i_point = Temp ; i_point < Temp+TempSize ; i_point++)
          {
            IndexType child_index = GetChildIndex(**i_point);
            *(cell_position[child_index]++) = *i_point;
          }

          // Creatig the first child
          if(cell_sizes[0] > BucketSize)
            this->mpChilds[0]= new OCTreePartitionAverageSplit(PointsBegin, cell_position[0], MinPoint, this->mPosition, BucketSize);
          else
            this->mpChilds[0]= new LeafType(PointsBegin, cell_position[0]);

          PointType new_min_point;
          PointType new_max_point;

          // Creating the rest of the childs
          for(IndexType i = 1 ; i < number_of_childs ; i++)
          {
            if(cell_sizes[i] > BucketSize)
            {
              IndexType factor = 1;
              for(IndexType j = 0 ; j < Dimension ; j++)
              {
                if(i & factor){
                  new_min_point[j] = this->mPosition[j];
                  new_max_point[j] = MaxPoint[j];
                } else {
                  new_min_point[j] = MinPoint[j];
                  new_max_point[j] = this->mPosition[j];
                }
                factor <<= 1;
              }

              this->mpChilds[i]= new OCTreePartitionAverageSplit(cell_position[i-1], cell_position[i], new_min_point, new_max_point, BucketSize);
            }
            else
              this->mpChilds[i]= new LeafType(cell_position[i-1], cell_position[i]);
          }

        }

        static TreeNodeType* Construct(IteratorType PointsBegin, IteratorType PointsEnd, 
            PointType HighPoint, PointType LowPoint, SizeType BucketSize)
        {
          SizeType number_of_points = std::distance(PointsBegin,PointsEnd);
          if (number_of_points == 0)
            return NULL;
          else if (number_of_points <= BucketSize)
            return new LeafType(PointsBegin, PointsEnd); 
          else 
            return new OCTreePartitionAverageSplit(PointsBegin, PointsEnd, LowPoint, HighPoint, BucketSize);
        }

      private:

        // Average Point Split
        void Split( PointType& Position, PointType const& MinPoint, PointType const& MaxPoint, IteratorType const& PointsBegin, IteratorType const& PointsEnd )
        {
          std::size_t i;
          for(i = 0 ; i < Dimension ; i++)
            Position[i] = 0.0;
          for(IteratorType i_point = PointsBegin ; i_point < PointsEnd ; i_point++)
            for(i = 0 ; i < Dimension ; i++)
              Position[i] += (**i_point)[i];
          for(i = 0 ; i < Dimension ; i++)
            Position[i] /= static_cast<CoordinateType>(std::distance(PointsBegin,PointsEnd));
        }
    };

    /// Short class definition.
  /** Detail class definition.
  */


    template< class TLeafType >
  class OCTreePartitionMidPointSplit : public OCTreePartitionBase<TLeafType>
  {
	public:

        ///@name Type Definitions
        ///@{

        /// Pointer definition of KDTree
        KRATOS_CLASS_POINTER_DEFINITION(OCTreePartionMidPointSplit);

        typedef TLeafType LeafType;
        typedef typename LeafType::PointType PointType;
        typedef typename LeafType::ContainerType ContainerType;
        typedef typename LeafType::IteratorType IteratorType;
        typedef typename LeafType::DistanceIteratorType DistanceIteratorType;
        typedef typename LeafType::PointerType PointerType;
        typedef typename LeafType::DistanceFunction DistanceFunction;
        enum { Dimension = LeafType::Dimension };
        typedef TreeNode<Dimension, PointType, PointerType, IteratorType, DistanceIteratorType> TreeNodeType;
        typedef typename TreeNodeType::CoordinateType CoordinateType;
        typedef typename TreeNodeType::SizeType SizeType;
        typedef typename TreeNodeType::IndexType IndexType;
        
        static const SizeType number_of_childs = 1 << Dimension;


        /// Partition constructor.
        OCTreePartitionMidPointSplit(IteratorType PointsBegin, IteratorType PointsEnd, 
            PointType const& MinPoint, PointType const& MaxPoint,  SizeType BucketSize = 1)
        {
          PointType mid_cell_lenght;

          SizeType TempSize = std::distance(PointsBegin,PointsEnd);
          PointerType* Temp = new PointerType[ TempSize ];

          // Template definition of SplitMode
          Split(this->mPosition,MinPoint,MaxPoint,PointsBegin,PointsEnd);

          SizeType cell_sizes[number_of_childs];
          IteratorType cell_position[number_of_childs];

          for(IndexType i = 0 ; i < number_of_childs ; i++)
            cell_sizes[i] = 0;

          PointerType* i_temp = Temp;
          for(IteratorType i_point = PointsBegin ; i_point < PointsEnd ; i_point++, i_temp++)
          {
            *i_temp = *i_point;
            IndexType child_index = GetChildIndex(**i_point);
            cell_sizes[child_index]++;
          }

          cell_position[0] = PointsBegin;
          for(IndexType i = 1 ; i < number_of_childs ; i++)
          {
            cell_position[i] = cell_position[i-1];
            std::advance(cell_position[i], cell_sizes[i-1]);
          }

          // 	      for(IteratorType i_point = PointsBegin ; i_point < PointsEnd ; i_point++)
          // 	      {
          // 			 IndexType child_index = GetChildIndex(**i_point);
          // 			 *(cell_position[child_index]++) = *i_point;
          // 			 //swap(*i_point, *(cell_position[child_index]++));
          // 		  }
          //  The same in Bin_Static ( use temporal array for reorder the list )
          for(PointerType* i_point = Temp ; i_point < Temp+TempSize ; i_point++)
          {
            IndexType child_index = GetChildIndex(**i_point);
            *(cell_position[child_index]++) = *i_point;
          }

          // Creatig the first child
          if(cell_sizes[0] > BucketSize)
            this->mpChilds[0]= new OCTreePartitionMidPointSplit(PointsBegin, cell_position[0], MinPoint, this->mPosition, BucketSize);
          else
            this->mpChilds[0]= new LeafType(PointsBegin, cell_position[0]);

          PointType new_min_point;
          PointType new_max_point;

          // Creating the rest of the childs
          for(IndexType i = 1 ; i < number_of_childs ; i++)
          {
            if(cell_sizes[i] > BucketSize)
            {
              IndexType factor = 1;
              for(IndexType j = 0 ; j < Dimension ; j++)
              {
                if(i & factor){
                  new_min_point[j] = this->mPosition[j];
                  new_max_point[j] = MaxPoint[j];
                } else {
                  new_min_point[j] = MinPoint[j];
                  new_max_point[j] = this->mPosition[j];
                }
                factor <<= 1;
              }

              this->mpChilds[i]= new OCTreePartitionMidPointSplit(cell_position[i-1], cell_position[i], new_min_point, new_max_point, BucketSize);
            }
            else
              this->mpChilds[i]= new LeafType(cell_position[i-1], cell_position[i]);

          }

        }

        static TreeNodeType* Construct(IteratorType PointsBegin, IteratorType PointsEnd, 
            PointType HighPoint, PointType LowPoint, SizeType BucketSize)
        {
          SizeType number_of_points = std::distance(PointsBegin,PointsEnd);
          if (number_of_points == 0)
            return NULL;
          else if (number_of_points <= BucketSize)
            return new LeafType(PointsBegin, PointsEnd); 
          else 
            return new OCTreePartitionMidPointSplit(PointsBegin, PointsEnd, LowPoint, HighPoint, BucketSize);
        }

    private:
	  void Split( PointType& Position, PointType const& MinPoint, PointType const& MaxPoint, IteratorType const& PointsBegin, IteratorType const& PointsEnd )
	  {
		// calculating the partition postion as midpoint of box
		for(std::size_t i = 0 ; i < Dimension ; i++)
		  Position[i] = (MaxPoint[i] + MinPoint[i]) * 0.500;
	  }
  };


  
}  // namespace Kratos.

#endif // KRATOS_OCT_TREE_H_INCLUDED   defined 


