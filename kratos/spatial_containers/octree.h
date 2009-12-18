//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: clabra $
//   Date:                $Date: 2007-03-27 17:02:19 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_OCT_TREE_H_INCLUDED )
#define  KRATOS_OCT_TREE_H_INCLUDED


// System includes
#include <cstddef>

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
  
  template< std::size_t Dimension, class PointType, class IteratorType, class CoordinateType >
  class OcTreeAverageSplit
  {
	public:
	  void operator()( PointType& Position, PointType const& MinPoint, PointType const& MaxPoint, IteratorType const& PointsBegin, IteratorType const& PointsEnd )
	  {
		for(std::size_t i = 0 ; i < Dimension ; i++)
		  Position[i] = 0.0;
		for(IteratorType i_point = PointsBegin ; i_point < PointsEnd ; i_point++)
		  for(std::size_t i = 0 ; i < Dimension ; i++)
			Position[i] += (**i_point)[i];
		for(std::size_t i = 0 ; i < Dimension ; i++)
		  Position[i] /= static_cast<CoordinateType>(SearchUtils::PointerDistance(PointsBegin,PointsEnd));
	  }
  };


  template< std::size_t Dimension, class PointType, class IteratorType, class CoordinateType >
  class OcTreeMidPointSplit
  {
	public:
	  void operator()( PointType& Position, PointType const& MinPoint, PointType const& MaxPoint, IteratorType const& PointsBegin, IteratorType const& PointsEnd )
	  {
		// calculating the partition postion as midpoint of box
		for(std::size_t i = 0 ; i < Dimension ; i++)
		  Position[i] = (MaxPoint[i] + MinPoint[i]) * 0.500;
	  }
  };



  /// Short class definition.
  /** Detail class definition.
  */
	template< class TLeafType >
	class OCTreePartition : public TreeNode< TLeafType::Dimension, 
                                             typename TLeafType::PointType,
	                                         typename TLeafType::PointerType,
											 typename TLeafType::IteratorType,
											 typename TLeafType::DistanceIteratorType >
	{
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of KDTree
      KRATOS_CLASS_POINTER_DEFINITION(OCTreePartition);

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

	  typedef OcTreeAverageSplit<Dimension,PointType,IteratorType,CoordinateType>  AverageSplit;
	  typedef OcTreeMidPointSplit<Dimension,PointType,IteratorType,CoordinateType> MidPointSplit;
	 
      typedef typename LeafType::SearchStructureType SearchStructureType;

	  static const SizeType number_of_childs = 1 << Dimension;


      ///@}
      ///@name Life Cycle 
      ///@{ 


      /// Partition constructor.
	  OCTreePartition(IteratorType PointsBegin, IteratorType PointsEnd, 
			PointType const& MinPoint, PointType const& MaxPoint,  SizeType BucketSize = 1)
	  {
		/* 	      const SizeType number_of_childs = 8; */
	      PointType mid_cell_lenght;
		  
		  SizeType TempSize = SearchUtils::PointerDistance(PointsBegin,PointsEnd);
		  PointerType* Temp = new PointerType[ TempSize ];

		  // Template definition of SplitMode
		  // SplitMode()(mPostion,MinPoint,MaxPoint,mPointBegin,mPointEnd);
		  AverageSplit()(mPosition,MinPoint,MaxPoint,PointsBegin,PointsEnd);
	      
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
			 mpChilds[0]= new OCTreePartition(PointsBegin, cell_position[0], MinPoint, mPosition, BucketSize);
	      else
			 mpChilds[0]= new LeafType(PointsBegin, cell_position[0]);
		  
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
                  new_min_point[j] = mPosition[j];
                  new_max_point[j] = MaxPoint[j];
                } else {
                  new_min_point[j] = MinPoint[j];
                  new_max_point[j] = mPosition[j];
                }
                factor <<= 1;
              }

              mpChilds[i]= new OCTreePartition(cell_position[i-1], cell_position[i], 
                  new_min_point , new_max_point,  BucketSize);
            }
            else
              mpChilds[i]= new LeafType(cell_position[i-1], cell_position[i]);

          }
		  
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

      /// Destructor.
	  virtual ~OCTreePartition()
	  {
		/* 		SizeType number_of_childs = 8; */
		for(SizeType i = 0 ; i < number_of_childs ; i++)
		  delete mpChilds[i];
	  }
      
      ///@}
      ///@name Operations
      ///@{

		void SearchNearestPoint( PointType const& ThisPoint, PointerType& Result, CoordinateType& ResultDistance )
		{
          SearchStructureType Auxiliar;
          SearchNearestPoint(ThisPoint,Result,ResultDistance,Auxiliar);
	    }

		void SearchNearestPoint(PointType const& ThisPoint, PointerType& Result, CoordinateType& rResultDistance,
            SearchStructureType& Auxiliar )
		{
		  CoordinateType distances_to_partitions[number_of_childs];

		  SizeType child_index = GetChildIndex(ThisPoint);

		  mpChilds[child_index]->SearchNearestPoint(ThisPoint, Result, rResultDistance, Auxiliar );

		  DistanceToPartitions(child_index, ThisPoint, distances_to_partitions, Auxiliar.residual_distance[0]);

		  for(SizeType i = 0 ; i < number_of_childs ; i++)
			if((i != child_index) && (distances_to_partitions[i] < rResultDistance)){
              Auxiliar.residual_distance[0] = distances_to_partitions[i];
			  mpChilds[i]->SearchNearestPoint(ThisPoint, Result, rResultDistance, Auxiliar);
			}

		}

		void SearchInRadius(PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
            DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults)
		{
           SearchStructureType Auxiliar;
		   SearchInRadius(ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Auxiliar );
		}

		void SearchInRadius(PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
            DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType& Auxiliar )
		{
		   SizeType child_index = GetChildIndex(ThisPoint);
		   
		   // n is number of points found
		   mpChilds[child_index]->SearchInRadius(ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Auxiliar);

		   CoordinateType distances_to_partitions[number_of_childs];

		   DistanceToPartitions(child_index, ThisPoint, distances_to_partitions, Auxiliar.residual_distance[0]);

		   for(IndexType i = 0 ; i < number_of_childs ; i++)
			  if((i != child_index) && (distances_to_partitions[i] < Radius2))
				 mpChilds[i]->SearchInRadius(ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Auxiliar);
		}
      
		void SearchInRadius(PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
            SizeType& NumberOfResults, SizeType const& MaxNumberOfResults)
		{
           SearchStructureType Auxiliar;
		   SearchInRadius(ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Auxiliar );
		}

		void SearchInRadius(PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
            SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType& Auxiliar )
		{
		   SizeType child_index = GetChildIndex(ThisPoint);
		   
		   // n is number of points found
		   mpChilds[child_index]->SearchInRadius(ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Auxiliar);

		   CoordinateType distances_to_partitions[number_of_childs];

		   DistanceToPartitions(child_index, ThisPoint, distances_to_partitions, Auxiliar.residual_distance[0]);

		   for(IndexType i = 0 ; i < number_of_childs ; i++)
			  if((i != child_index) && (distances_to_partitions[i] < Radius2))
				 mpChilds[i]->SearchInRadius(ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Auxiliar);
		}


	private:

	    IndexType GetChildIndex(PointType const& rThisPoint) const
	    {
		   // Calculating the cell index 
		   IndexType child_index = 0;
		   const IndexType dim_mask[] = { 1, 2, 4 };

		   for( IndexType i = 0 ; i < Dimension ; i++)
 			  if( rThisPoint[i] >= mPosition[i] ) child_index += dim_mask[i];

//		   for(IndexType i = 0 ; i < Dimension ; i++)
//			  child_index += IndexType(rThisPoint[i] >= mPosition[i]) << i;

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
			   
			   //IndexType coordinate_mask = 1;
			   //IndexType coordinate_mask = 0;
			   for(IndexType j = 0 ; j < Dimension ; j++)
			   {
				  //rDistances[i] += ( coordinate_mask & partitions ) * offset_from_postition[j];
				  //coordinate_mask <<= 1;
				  if(coordinate_mask[j] & partitions) rDistances[i] += offset_from_postition[j];
			   }
			}
		    
		}
	   
	private:

		IndexType mCutingDimension;
		PointType mPosition; // Position of partition
    
		TreeNodeType* mpChilds[8];  // 8 is number of childs



	public:
		static TreeNodeType* Construct(IteratorType PointsBegin, 
									   IteratorType PointsEnd, 
									   PointType HighPoint, 
									   PointType LowPoint,
									   SizeType BucketSize)
		{
		  SizeType number_of_points = SearchUtils::PointerDistance(PointsBegin,PointsEnd);
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

	};
  
  
}  // namespace Kratos.

#endif // KRATOS_OCT_TREE_H_INCLUDED   defined 


