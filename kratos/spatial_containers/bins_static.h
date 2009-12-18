//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: clabra $
//   Date:                $Date: 2007-03-27 17:02:19 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_BINS_CONTAINER_H_INCLUDE)
#define KRATOS_BINS_CONTAINER_H_INCLUDE

#include "tree.h"


namespace Kratos {



template<  std::size_t TDimension,
           class TPointType,
		   class TContainerType,
		   class TPointerType = typename TContainerType::value_type,
		   class TIteratorType = typename TContainerType::iterator,
		   class TDistanceIteratorType = typename std::vector<double>::iterator,
		   class TDistanceFunction = Kratos::SearchUtils::SquaredDistanceFunction<TDimension,TPointType> >
class Bins : public TreeNode<TDimension,TPointType, TPointerType, TIteratorType, TDistanceIteratorType> 
{


   public:

      /// Pointer definition of Bins
      KRATOS_CLASS_POINTER_DEFINITION(Bins);

     typedef TPointType                         PointType;
	 typedef TContainerType                     ContainerType;
	 typedef TIteratorType                      IteratorType;
	 typedef TDistanceIteratorType              DistanceIteratorType;
	 typedef TPointerType                       PointerType;
	 typedef TDistanceFunction                  DistanceFunction;
	 enum { Dimension = TDimension };

	 typedef TreeNode<Dimension,PointType,PointerType,IteratorType,DistanceIteratorType> TreeNodeType;

	 typedef typename TreeNodeType::SizeType        SizeType;
	 typedef typename TreeNodeType::IndexType       IndexType;
	 typedef typename TreeNodeType::CoordinateType  CoordinateType;

	 typedef TreeNodeType LeafType;

	 // not always PointVector == ContainerType ( if ContainerType = C array )
	 typedef std::vector<PointerType>        PointVector;
	 typedef typename PointVector::iterator  PointIterator;

	 typedef std::vector<IteratorType>          IteratorVector;
	 typedef typename IteratorVector::iterator  IteratorIterator;
	 typedef typename IteratorVector::const_iterator IteratorConstIterator;

	 typedef Tvector<IndexType,TDimension>   CellType;
    
     typedef typename TreeNodeType::IteratorIteratorType IteratorIteratorType;
     typedef typename TreeNodeType::SearchStructureType SearchStructureType;

     typedef Kratos::SearchUtils::SearchNearestInRange<PointType,PointerType,IteratorType,DistanceFunction,CoordinateType> SearchNearestInRange;
     typedef Kratos::SearchUtils::SearchRadiusInRange<PointType,IteratorType,DistanceIteratorType,DistanceFunction,SizeType,CoordinateType> SearchRadiusInRange;
     typedef Kratos::SearchUtils::SearchBoxInRange<PointType,IteratorType,SizeType,TDimension> SearchBoxInRange;


   public:


	 //************************************************************************

	 // constructor 1
	 Bins() : mPointBegin(this->NullIterator()), mPointEnd(this->NullIterator()) {};

	 //************************************************************************

	 Bins( IteratorType const& PointBegin, IteratorType const& PointEnd, SizeType BucketSize = 1 )
	   : mPointBegin(PointBegin), mPointEnd(PointEnd)
	 {
	   if(mPointBegin==mPointEnd)
		 return;
	   CalculateBoundingBox();
	   CalculateCellSize();
	   AllocateCellsContainer();
	   GenerateBins();
	 }

	 //************************************************************************

	 Bins( IteratorType const& PointBegin, IteratorType const& PointEnd, PointType const& MinPoint, PointType const& MaxPoint, SizeType BucketSize = 1 )
	   : mPointBegin(PointBegin), mPointEnd(PointEnd)
	 {
	   if(mPointBegin==mPointEnd)
		 return;

	   for(SizeType i = 0 ; i < TDimension ; i++)
	   {
		 mMinPoint[i] = MinPoint[i];
		 mMaxPoint[i] = MaxPoint[i];
	   }

	   CalculateCellSize();
	   AllocateCellsContainer();
	   GenerateBins();
	 }

	 //************************************************************************

	 Bins( IteratorType const& PointBegin, IteratorType const& PointEnd, CoordinateType BoxSize, SizeType BucketSize = 1 )
	   : mPointBegin(PointBegin), mPointEnd(PointEnd)
	 {
	   if(mPointBegin==mPointEnd)
		 return;
	   CalculateBoundingBox();
	   CalculateCellSize(BoxSize);
	   AllocateCellsContainer();
	   GenerateBins();
	 }

	 //************************************************************************

	 // destructor
	 virtual ~Bins(){ }

	 //************************************************************************

	 IteratorType Begin() { return mPointBegin; }

	 //************************************************************************

	 IteratorType End() { return mPointBegin; }

	 //************************************************************************

	 CoordinateType CellSize( SizeType const& iDim ) { return mCellSize[iDim]; }

	 //************************************************************************

	 SizeType NumCell( SizeType const& iDim ) { return mN[iDim]; }

	 //************************************************************************
	 
   private:

     //************************************************************************

	 void CalculateBoundingBox() {
	   for(SizeType i = 0 ; i < TDimension ; i++){
		 mMinPoint[i] = (**mPointBegin)[i];
		 mMaxPoint[i] = (**mPointBegin)[i];
	   }
	   for(IteratorType Point = mPointBegin ; Point != mPointEnd ; Point++)
		 for(SizeType i = 0 ; i < TDimension ; i++){
		   if( (**Point)[i] < mMinPoint[i] ) mMinPoint[i] = (**Point)[i];
		   if( (**Point)[i] > mMaxPoint[i] ) mMaxPoint[i] = (**Point)[i];
		 }
	 }

	 //************************************************************************

	 void CalculateCellSize() {

       CoordinateType delta[TDimension];
       CoordinateType alpha[TDimension];
       CoordinateType mult_delta = 1.00;
       SizeType index = 0;
       for(SizeType i = 0 ; i < TDimension ; i++) {
         delta[i] = mMaxPoint[i] - mMinPoint[i];
         if ( delta[i] > delta[index] )
           index = i;
         delta[i] = (delta[i] == 0.00) ? 1.00 : delta[i];
       }

       for(SizeType i = 0 ; i < TDimension ; i++){
         alpha[i] = delta[i] / delta[index];
         mult_delta *= alpha[i];
       }

       mN[index] = static_cast<SizeType>( pow(static_cast<CoordinateType>(SearchUtils::PointerDistance(mPointBegin,mPointEnd)/mult_delta), 1.00/TDimension)+1 );

       for(SizeType i = 0 ; i < TDimension ; i++){
         if(i!=index) {
           mN[i] = static_cast<SizeType>(alpha[i] * mN[index]);
           mN[i] = ( mN[i] == 0 ) ? 1 : mN[i];
         }
       }

       for(SizeType i = 0 ; i < TDimension ; i++){
         mCellSize[i] = delta[i] / mN[i];
         mInvCellSize[i] = 1.00 / mCellSize[i];
       }
	 }

	 //************************************************************************

	 void CalculateCellSize( CoordinateType BoxSize ) {
	   for(SizeType i = 0 ; i < TDimension ; i++){
		 mCellSize[i] = BoxSize;
		 mInvCellSize[i] = 1.00 / mCellSize[i];
		 mN[i] = static_cast<SizeType>( (mMaxPoint[i]-mMinPoint[i]) / mCellSize[i]) + 1;
	   }
	 }

	 //************************************************************************

	 void AllocateCellsContainer() {
	   SizeType Size = 1;
	   for(SizeType i = 0 ; i < TDimension ; i++)
		 Size *= mN[i];
	   mIndexCell.resize(Size+1);
	   mIndexCellBegin = mIndexCell.begin();
	   mIndexCellEnd   = mIndexCell.end();
	 }

	 //************************************************************************

	 void GenerateBins( ){

	   PointVector TempPoint(mPointBegin,mPointEnd);

	   // Reset index vector
	   for( IteratorIterator Iter = mIndexCell.begin(); Iter != mIndexCell.end(); Iter++)
		 *Iter = mPointBegin;

	   // Update storage counter, storing ahead
	   for( IteratorType Point = mPointBegin ; Point != mPointEnd ; Point++)
		 mIndexCell[ CalculateIndex(**Point) + 1 ]++;

	   // Storage/reshufing pass 1

	   // Update storage counter and store
	   for( IteratorIterator Iter = mIndexCell.begin()+1 ; Iter != mIndexCell.end() ; Iter++)
		 *Iter = *(Iter-1) + SearchUtils::PointerDistance(mPointBegin,*Iter);

	   // Point pass 2
	   // Store the points in lbin1

	   // Update storage counter, storing in lbin1
	   for( PointIterator Point = TempPoint.begin() ; Point != TempPoint.end() ; Point++)
		 *(mIndexCell[CalculateIndex(**Point)]++) = *Point;
/* 
	   // TEST  !!! OJO -> No aumenta el contador del IteratorIterator !!!!
	   for( PointIterator Point = TempPoint.begin() ; Point != TempPoint.end() ; Point++)
	   {
	     Iter = mIndexCell[CalculateIndex(**Point)];
	     while( Iter != Point )
		 {
	       Iter2 = mIndexCell[CalculateIndex(**Iter)];
	       std::swap(*Iter,*Iter2);
	     }
	   }
*/

	   // Storage/reshuffing pass 2

	   // Loop over bins, in reverse order
	   for(IteratorIterator Iter = mIndexCell.end()-1; Iter != mIndexCell.begin(); Iter--)
		 *Iter = *(Iter-1);
	   mIndexCell[0] = mPointBegin;

	 }

	 //************************************************************************

	 IndexType CalculatePosition( CoordinateType const& ThisCoord, SizeType ThisDimension )
	 {
	   CoordinateType d_index = (ThisCoord - mMinPoint[ThisDimension]) * mInvCellSize[ThisDimension];
	   IndexType index = static_cast<SizeType>( (d_index < 0.00) ? 0.00 : d_index );
	   return  (index > mN[ThisDimension]-1) ? mN[ThisDimension]-1 : index;
	 }

	 //************************************************************************

	 IndexType CalculateIndex( PointType const& ThisPoint )
	 {
	   IndexType Index = 0;
	   for(SizeType iDim = TDimension-1 ; iDim > 0 ; iDim--){
		 Index += CalculatePosition(ThisPoint[iDim],iDim);
		 Index *= mN[iDim-1];
	   }
	   Index += CalculatePosition(ThisPoint[0],0);
	   return Index;
	 }

	 //************************************************************************

	 IndexType CalculateIndex( CellType const& ThisIndex )
	 {
	   IndexType Index = 0;
	   for(SizeType iDim = TDimension-1 ; iDim > 0 ; iDim--){
		 Index += ThisIndex[iDim];
		 Index *= mN[iDim-1];
	   }
	   Index += ThisIndex[0];
	   return Index;
	 }

	 //************************************************************************

	 CellType CalculateCell( PointType const& ThisPoint ){
	   CellType Cell;
	   for(SizeType i = 0 ; i < TDimension ; i++)
		 Cell[i] = CalculatePosition(ThisPoint[i],i);
	   return Cell;
	 }

	 CellType CalculateCell( PointType const& ThisPoint, CoordinateType Radius ){
	   CellType Cell;
	   for(SizeType i = 0 ; i < TDimension ; i++)
		 Cell[i] = CalculatePosition(ThisPoint[i]+Radius,i);
	   return Cell;
	 }

	 //************************************************************************

   public:

	 //************************************************************************
     //************************************************************************

     PointerType ExistPoint( PointerType const& ThisPoint, CoordinateType const Tolerance = static_cast<CoordinateType>(10.0*DBL_EPSILON) )
     {
       PointerType Nearest;
       CoordinateType Distance = static_cast<CoordinateType>(DBL_MAX);
       bool Found;
       SearchStructureType Box( CalculateCell(*ThisPoint,-Tolerance), CalculateCell(*ThisPoint,Tolerance), mN );
       SearchNearestInBox( *ThisPoint, Nearest, Distance, Box, Found );
       if(Found)
         return Nearest;
       return this->NullPointer();
     }

     //************************************************************************
     //************************************************************************

	 PointerType SearchNearestPoint( PointType const& ThisPoint )
	 {
	   PointerType Result            = *mPointBegin;                           //static_cast<PointerType>(NULL);
	   CoordinateType ResultDistance = static_cast<CoordinateType>(DBL_MAX);
       SearchStructureType Box( CalculateCell(ThisPoint), mN, mIndexCellBegin );
	   SearchNearestPointLocal( ThisPoint, Result, ResultDistance, Box );
	   return Result;
	 }

	 //************************************************************************

	 PointerType SearchNearestPoint( PointType const& ThisPoint, CoordinateType& rResultDistance )
	 {
	   PointerType Result = *mPointBegin;                           //static_cast<PointerType>(NULL);
	   rResultDistance    = static_cast<CoordinateType>(DBL_MAX);
       SearchStructureType Box( CalculateCell(ThisPoint), mN, mIndexCellBegin );
	   SearchNearestPointLocal( ThisPoint, Result, rResultDistance, Box);
	   return Result;
	 }

	 //************************************************************************
	 
     // New Thread Safe!!!
	 PointerType SearchNearestPoint( PointType const& ThisPoint, CoordinateType& rResultDistance, SearchStructureType& Box )
	 {
	   PointerType Result            = *mPointBegin;                           //static_cast<PointerType>(NULL);
	   CoordinateType ResultDistance = static_cast<CoordinateType>(DBL_MAX);
       Box.Set( CalculateCell(ThisPoint), mN, mIndexCellBegin );
	   SearchNearestPointLocal( ThisPoint, Result, rResultDistance, Box);
	   return Result;
	 }
     
	 //************************************************************************
	 //************************************************************************

     void SearchNearestPoint( PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance )
     {
       SearchStructureType Box( CalculateCell(ThisPoint), mN, mIndexCellBegin );
       SearchNearestPointLocal(ThisPoint,rResult,rResultDistance,Box);
     }

     //************************************************************************

	 void SearchNearestPoint( PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance, SearchStructureType& Box )
	 {
	   // This case is when BinStatic is a LeafType in Other Spacial Structure
	   // Then, it is possible a better Result before this search
       Box.Set( CalculateCell(ThisPoint), mN, mIndexCellBegin );
	   SearchNearestPointLocal( ThisPoint, rResult, rResultDistance, Box );
	 }

	 //************************************************************************
	 //************************************************************************

     // **** THREAD SAFE  -> The user pass the SearchStructure (BinBox)
	 void SearchNearestPointLocal( PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance, SearchStructureType& Box ) 
	 {
	   if( mPointBegin == mPointEnd )
		 return;

       bool Found;

	   // initial search
	   ++Box;
	   SearchNearestInBox( ThisPoint, rResult, rResultDistance, Box, Found );
	   // increase mBox and try again
	   while(!Found){
		 ++Box;
		 SearchNearestInBox( ThisPoint, rResult, rResultDistance, Box, Found );
	   }
	 
     }

	 //************************************************************************
	 //************************************************************************

	 SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, IteratorType Results,
         DistanceIteratorType ResultsDistances, SizeType const& MaxNumberOfResults )
	 {
	   CoordinateType Radius2 = Radius * Radius;
       SizeType NumberOfResults = 0;
       SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
	   SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box );
       return NumberOfResults;
	 }

	 //************************************************************************
	 
     SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, IteratorType Results,
         DistanceIteratorType ResultsDistances, SizeType const& MaxNumberOfResults, SearchStructureType& Box )
	 {
	   CoordinateType Radius2 = Radius * Radius;
       SizeType NumberOfResults = 0;
	   Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
	   SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box );
       return NumberOfResults;
	 }

     //************************************************************************

	 void SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
         DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults )
	 {
       SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
	   SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box);
	 }

     //************************************************************************
	 
     void SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
         DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType& Box )
	 {
       Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
	   SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box);
	 }
     
     
     //************************************************************************

     // **** THREAD SAFE
	
     // Dimension = 1
     void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
         DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
	 {
       SearchRadiusInRange()(*(Box.RowBegin),*(Box.RowEnd),ThisPoint,Radius2,Results,ResultsDistances,NumberOfResults,MaxNumberOfResults);
	 }
	
     // Dimension = 2
     void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
         DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
	 {
       for(IndexType I = Box.Axis[1].Begin() ; I <= Box.Axis[1].End() ; I += Box.Axis[1].Block )
         SearchRadiusInRange()(Box.RowBegin[I],Box.RowEnd[I],ThisPoint,Radius2,Results,ResultsDistances,NumberOfResults,MaxNumberOfResults);
	 }
	
     // Dimension = 3
     void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
         DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
	 {
       for(IndexType II = Box.Axis[2].Begin() ; II <= Box.Axis[2].End() ; II += Box.Axis[2].Block )
         for(IndexType I = II + Box.Axis[1].Begin() ; I <= II + Box.Axis[1].End() ; I += Box.Axis[1].Block )
           SearchRadiusInRange()(Box.RowBegin[I],Box.RowEnd[I],ThisPoint,Radius2,Results,ResultsDistances,NumberOfResults,MaxNumberOfResults);
	 }

	 //************************************************************************
	 //************************************************************************

	 SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType Radius, IteratorType Results, SizeType MaxNumberOfResults )
	 {
	   CoordinateType Radius2 = Radius * Radius;
       SizeType NumberOfResults = 0;
       SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
	   SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Box );
       return NumberOfResults;
	 }

	 //************************************************************************

	 SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType Radius, IteratorType Results,
         SizeType MaxNumberOfResults, SearchStructureType& Box )
	 {
	   CoordinateType Radius2 = Radius * Radius;
       SizeType NumberOfResults = 0;
	   Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
	   SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Box );
       return NumberOfResults;
	 }

	 //************************************************************************

	 void SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
         SizeType& NumberOfResults, SizeType const& MaxNumberOfResults )
	 {
       SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
	   SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Box );
	 }

	 //************************************************************************

	 void SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
         SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType& Box )
	 {
       Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
	   SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Box );
	 }

	 //************************************************************************

     // Dimension = 1
	 void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
         SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, 
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
	 {
       SearchRadiusInRange()(*(Box.RowBegin),*(Box.RowEnd),ThisPoint,Radius2,Results,NumberOfResults,MaxNumberOfResults);
	 }
	 
     // Dimension = 2
     void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
         SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, 
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
	 {
       for(IndexType I = Box.Axis[1].Begin() ; I <= Box.Axis[1].End() ; I += Box.Axis[1].Block )
         SearchRadiusInRange()(Box.RowBegin[I],Box.RowEnd[I],ThisPoint,Radius2,Results,NumberOfResults,MaxNumberOfResults);
	 }
     
     // Dimension = 3
     void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
         SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, 
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
	 {
       for(IndexType II = Box.Axis[2].Begin() ; II <= Box.Axis[2].End() ; II += Box.Axis[2].Block )
         for(IndexType I = II + Box.Axis[1].Begin() ; I <= II + Box.Axis[1].End() ; I += Box.Axis[1].Block )
           SearchRadiusInRange()(Box.RowBegin[I],Box.RowEnd[I],ThisPoint,Radius2,Results,NumberOfResults,MaxNumberOfResults);
	 }

	 //************************************************************************
	 //************************************************************************

     // Dimension = 1
	 void SearchNearestInBox( PointType const& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance, 
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box, bool& Found )
     {
       Found = false;
       SearchNearestInRange()( *(Box.RowBegin), *(Box.RowEnd), ThisPoint, ResultPoint, ResultDistance, Found );
     }

     // Dimension = 2
	 void SearchNearestInBox( PointType const& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance, 
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box, bool& Found )
     {
       Found = false;
       for(IndexType I = Box.Axis[1].Begin() ; I <= Box.Axis[1].End() ; I += Box.Axis[1].Block )
         SearchNearestInRange()( Box.RowBegin[I], Box.RowEnd[I], ThisPoint, ResultPoint, ResultDistance, Found );
     }
	 
     // Dimension = 3
     void SearchNearestInBox( PointType const& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance, 
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box, bool& Found )
     {
       Found = false;
       for(IndexType II = Box.Axis[2].Begin() ; II <= Box.Axis[2].End() ; II += Box.Axis[2].Block )
         for(IndexType I = II + Box.Axis[1].Begin() ; I <= II + Box.Axis[1].End() ; I += Box.Axis[1].Block )
           SearchNearestInRange()( Box.RowBegin[I], Box.RowEnd[I], ThisPoint, ResultPoint, ResultDistance, Found );
     }

	 //************************************************************************
	 //************************************************************************
	 
	 SizeType SearchInBox( PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType Results,
         SizeType MaxNumberOfResults )
     {
       SizeType NumberOfResults = 0;
       SearchStructureType Box( CalculateCell(SearchMinPoint), CalculateCell(SearchMaxPoint), mN, mIndexCellBegin );
       SearchInBoxLocal( SearchMinPoint, SearchMaxPoint, Results, NumberOfResults, MaxNumberOfResults, Box );
       return NumberOfResults;
     }
     
     //************************************************************************
	 
     void SearchInBox(PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& Results, SizeType& NumberOfResults,
            SizeType const& MaxNumberOfResults )
     {
       NumberOfResults = 0;
       SearchStructureType Box( CalculateCell(SearchMinPoint), CalculateCell(SearchMaxPoint), mN, mIndexCellBegin );
       SearchInBoxLocal( SearchMinPoint, SearchMaxPoint, Results, NumberOfResults, MaxNumberOfResults, Box );
     }
     
     //************************************************************************

     // Dimension = 1
	 void SearchInBoxLocal( PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& ResultsPoint,
         SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
     {
       SearchBoxInRange()(SearchMinPoint,SearchMaxPoint,*(Box.RowBegin),*(Box.RowEnd),ResultsPoint,NumberOfResults,MaxNumberOfResults);
     }

     // Dimension = 2
	 void SearchInBoxLocal( PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& ResultsPoint,
         SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
     {
       for(IndexType I = Box.Axis[1].Begin() ; I <= Box.Axis[1].End() ; I += Box.Axis[1].Block )
         SearchBoxInRange()(SearchMinPoint,SearchMaxPoint,Box.RowBegin[I],Box.RowEnd[I],ResultsPoint,NumberOfResults,MaxNumberOfResults);
     }

     // Dimension = 3
	 void SearchInBoxLocal( PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& ResultsPoint,
         SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
     {
       for(IndexType II = Box.Axis[2].Begin() ; II <= Box.Axis[2].End() ; II += Box.Axis[2].Block )
         for(IndexType I = II + Box.Axis[1].Begin() ; I <= II + Box.Axis[1].End() ; I += Box.Axis[1].Block )
           SearchBoxInRange()(SearchMinPoint,SearchMaxPoint,Box.RowBegin[I],Box.RowEnd[I],ResultsPoint,NumberOfResults,MaxNumberOfResults);
     }

	 //************************************************************************
	 //************************************************************************

	 /// Turn back information as a string.
	 virtual std::string Info() const
	 {
	   return "BinsContainer";
	 }

	 /// Print information about this object.
	 virtual void PrintInfo(std::ostream& rOStream) const
	 {
	   rOStream << "BinsContainer";
	 }

	 /// Print object's data.
	 virtual void PrintData(std::ostream& rOStream, std::string const& Perfix = std::string()) const
	 {
	   rOStream << Perfix << "Bin[" << SearchUtils::PointerDistance(mPointBegin, mPointEnd) << "] : " << std::endl;
	   for(IteratorConstIterator i_cell = mIndexCell.begin() ; i_cell != mIndexCell.end()-1 ; i_cell++)
	   {
		 rOStream << Perfix << "[ " ;
		 for(IteratorType i_point = *i_cell ; i_point != *(i_cell+1) ; i_point++)
		   rOStream << **i_point << " ";
		 rOStream << "]" << std::endl;
	   }
	   rOStream << std::endl;
	 }

	 /// Print Size of Container
	 void PrintSize( std::ostream& rout ){
	   rout << " BinsSize: ";
	   for(SizeType i = 0 ; i < TDimension ; i++)
		 rout << "[" << mN[i] << "]";
	   rout << std::endl;
	 }

	 /// Print Limits Points of the Container
	 void PrintBox( std::ostream& rout ){
	   rout << " BinsBox: Min [";  mMinPoint.Print(rout);
	   rout <<       "];  Max [";  mMaxPoint.Print(rout);
	   rout <<       "];  Size ["; mCellSize.Print(rout);
	   rout << "]" << std::endl;
	 }

	 /// Assignment operator.
	 Bins& operator=(Bins const& rOther);

	 /// Copy constructor.
	 Bins(Bins const& rOther);

	 TPointType GetMinPoint()
	 {
	   TPointType point;
	   for(SizeType i = 0 ; i < TDimension ; i++)
		 point[i] = mMinPoint[i];
	   return point;
	 }

	 TPointType GetMaxPoint()
	 {
	   TPointType point;
	   for(SizeType i = 0 ; i < TDimension ; i++)
		 point[i] = mMaxPoint[i];
	   return point;
	 }


   private:

     // Point Access Iterators (vector reordered!!)
	 IteratorType     mPointBegin;
	 IteratorType     mPointEnd;

     // Bin Parameters (Sizes,BoundingBox,...)
	 PointType  mMinPoint;
	 PointType  mMaxPoint;
	 PointType  mCellSize;
	 PointType  mInvCellSize;
	 Tvector<SizeType,TDimension>  mN;

	 // Bins Access Vector ( vector<Iterator> )
	 IteratorVector           mIndexCell;
	 IteratorIterator         mIndexCellBegin;
	 IteratorIterator         mIndexCellEnd;

	 // Work Variables ( For non-copy of Search Variables )
	 //SearchStructureType mBox;

   public:

	 static TreeNodeType* Construct(IteratorType PointsBegin, IteratorType PointsEnd, PointType MaxPoint, PointType MinPoint, SizeType BucketSize)
	 {
	   SizeType number_of_points = SearchUtils::PointerDistance(PointsBegin,PointsEnd);
	   if (number_of_points == 0)
		 return NULL;
	   else 
	   {
		 return new Bins( PointsBegin, PointsEnd, MinPoint, MaxPoint, BucketSize );
	   }
	 }


};

template< std::size_t TDimension, class TPointType, class TContainerType, class TPointerType,
          class TIteratorType, class TDistanceIteratorType, class TDistanceFunction >
std::ostream & operator<<( std::ostream& rOStream, Bins<TDimension,TPointType,TContainerType,TPointerType,TIteratorType,TDistanceIteratorType,TDistanceFunction>& rThis)
{
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintSize(rOStream);
	rThis.PrintData(rOStream);
	return rOStream;
}



}

#endif // KRATOS_BINS_CONTAINER_H_INCLUDE
