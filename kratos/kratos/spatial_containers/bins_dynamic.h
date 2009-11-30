//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2007-03-27 17:02:19 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_BINS_DYNAMIC_CONTAINER_H_INCLUDE)
#define KRATOS_BINS_DYNAMIC_CONTAINER_H_INCLUDE

#include <vector>
#include <cfloat>

#include "tree.h"

namespace Kratos {

   template<
      std::size_t TDimension,
      class TPointType,
      class TContainerType,
      class TPointerType = typename TContainerType::value_type,
      class TIteratorType = typename TContainerType::iterator,
      class TDistanceIteratorType = typename std::vector<double>::iterator,
      class TDistanceFunction = SpacialSearchSquaredDistanceFunction<TDimension,TPointType>
      >
   class BinsDynamic : public TreeNode<TDimension,TPointType, TPointerType, TIteratorType, TDistanceIteratorType>
	  {
		
	  public:

        /// Pointer definition of BinsDynamic
        KRATOS_CLASS_POINTER_DEFINITION(BinsDynamic);

		typedef TreeNode<TDimension,TPointType,TPointerType,TIteratorType,TDistanceIteratorType> TreeNodeType;
		typedef TPointType                         PointType;
		typedef TContainerType                     ContainerType;
		typedef TIteratorType                      IteratorType;
		typedef TDistanceIteratorType              DistanceIteratorType;
		typedef TPointerType                       PointerType;
		typedef TDistanceFunction                  DistanceFunction;
        enum { Dimension = TDimension };

		typedef typename TreeNodeType::CoordinateType  CoordinateType;  // double
		typedef typename TreeNodeType::SizeType        SizeType;        // std::size_t
		typedef typename TreeNodeType::IndexType       IndexType;       // std::size_t

		typedef TreeNodeType LeafType;

        typedef typename TreeNodeType::SearchStructureType SearchStructureType;

		// Local Container ( PointPointer Container per Cell )
		// can be different to ContainerType
        // not always PointVector == ContainerType ( if ContainerType = C array )
		typedef std::vector<PointerType>       PointVector;
		typedef typename PointVector::iterator PointIterator;
        //typedef std::vector<IteratorType>          IteratorVector;
        //typedef typename IteratorVector::iterator  IteratorIterator;
        //typedef typename IteratorVector::const_iterator IteratorConstIterator;

		// Global Container
		typedef std::vector<PointVector>        CellsContainerType; 

		typedef Tvector<IndexType,TDimension>   CellType;

	  public:

		//************************************************************************

		// constructor 1
		BinsDynamic() : mPointBegin(this->NullIterator()), mPointEnd(this->NullIterator()), mNumPoints(0)
		{};

        //************************************************************************
         
        BinsDynamic( IteratorType const& PointBegin, IteratorType const& PointEnd, SizeType BucketSize = 1 )
        : mPointBegin(PointBegin), mPointEnd(PointEnd)
        {
           if(mPointBegin==mPointEnd)
              return;
           mNumPoints = std::distance(mPointBegin,mPointEnd);
           CalculateBoundingBox();
           CalculateCellSize();
           AllocateCellsContainer();
           GenerateBins();
        }
        
        //************************************************************************
         
        BinsDynamic( IteratorType const& PointBegin, IteratorType const& PointEnd, PointType const& MinPoint, PointType const& MaxPoint, SizeType BucketSize = 1 )
        : mPointBegin(PointBegin), mPointEnd(PointEnd)
        {
           if(mPointBegin==mPointEnd)
              return;

           mNumPoints = std::distance(mPointBegin,mPointEnd);
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
         
        BinsDynamic( PointType const& MinPoint, PointType const& MaxPoint, SizeType BucketSize )
          : mNumPoints(0)
        {
		   for(SizeType i = 0 ; i < TDimension ; i++)
		   {
			 mMinPoint[i] = MinPoint[i];
			 mMaxPoint[i] = MaxPoint[i];
		   }
           CalculateCellSize(BucketSize);
           AllocateCellsContainer();
        }
        
        //************************************************************************
        
        BinsDynamic( IteratorType const& PointBegin, IteratorType const& PointEnd, CoordinateType BoxSize, SizeType BucketSize = 1 )
        : mPointBegin(PointBegin), mPointEnd(PointEnd)
        {
           if(mPointBegin==mPointEnd)
              return;
           mNumPoints = std::distance(mPointBegin,mPointEnd);
           CalculateBoundingBox();
           CalculateCellSize(BoxSize);
           AllocateCellsContainer();
           GenerateBins();
        }
        
        //************************************************************************

        // destructor
        virtual ~BinsDynamic(){ }
        
        //************************************************************************
        
        IteratorType Begin() { return mPointBegin; }

        //************************************************************************
        
        IteratorType End() { return mPointBegin; }

        //************************************************************************
        
        CoordinateType CellSize( SizeType const& iDim ) { return mCellSize[iDim]; }
        
        //************************************************************************
        
        SizeType NumCell( SizeType const& iDim ) { return mN[iDim]; }
        
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
            for(SizeType i = 0 ; i < TDimension ; i++) {
                delta[i] = mMaxPoint[i] - mMinPoint[i];
                delta[i] = (delta[i] == 0.00) ? 1.00 : delta[i];
            }

            alpha[0] = 1.00;
            for(SizeType i = 1 ; i < TDimension ; i++){
                alpha[i] = delta[i] / delta[0];
                mult_delta *= alpha[i];
            }

            mN[0] = static_cast<SizeType>( pow(static_cast<CoordinateType>(PointerDistance(mPointBegin,mPointEnd)/mult_delta), 1.00/TDimension)+1 );

            for(SizeType i = 1 ; i < TDimension ; i++)
                mN[i] = static_cast<SizeType>(alpha[i] * mN[0]);

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
            // Resize Global Container
            mPoints.resize(Size);
        }
        
        //************************************************************************
         
        void GenerateBins(){
           
           for(IteratorType i_point = mPointBegin ; i_point != mPointEnd ; i_point++)
              mPoints[CalculateIndex(**i_point)].push_back(*i_point); 
            
        }
        
        //************************************************************************

        IndexType CalculatePosition( CoordinateType const& ThisCoord, SizeType ThisDimension )
        {
            CoordinateType d_index = (ThisCoord - mMinPoint[ThisDimension]) * mInvCellSize[ThisDimension];
            IndexType index = static_cast<IndexType>( (d_index < 0.00) ? 0.00 : d_index );
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
        
         /*          void AddPoint( PointType const& ThisPoint ){ */
         /*             mPoints[CalculateIndex(ThisPoint)].push_back(&ThisPoint); */
         /*          } */
        
        //************************************************************************
        
         void AddPoint( PointerType const& ThisPoint ){
            mPoints[CalculateIndex(*ThisPoint)].push_back(ThisPoint);
            mNumPoints++;
         }
        
        //************************************************************************
        
        PointerType ExistPoint( PointType const& ThisPoint, CoordinateType const Tolerance = static_cast<CoordinateType>(10.0*DBL_EPSILON) )
        {
         CoordinateType Tolerance2 = Tolerance*Tolerance;
         PointerType nearest = SearchNearestPoint(ThisPoint);
         if( nearest != this->NullPointer() )
           if( TDistanceFunction()(*nearest,*ThisPoint) < Tolerance2 )
             return nearest;
         return this->NullPointer();
        }
        
        //************************************************************************
        
        PointerType ExistPoint( PointerType const& ThisPoint, CoordinateType const Tolerance = static_cast<CoordinateType>(10.0*DBL_EPSILON) )
        {
         CoordinateType Tolerance2 = Tolerance*Tolerance;
         PointerType nearest;
         SearchNearestPoint(*ThisPoint,nearest,Tolerance2);
         if( nearest != this->NullPointer() )
           if( TDistanceFunction()(*nearest,*ThisPoint) < Tolerance2 )
             return nearest;
         return this->NullPointer();
        }
        
        //************************************************************************

        class BinBox
        {

           public:
                CellType Min;
                CellType Max;
                Tvector<SizeType,TDimension> Size;
                Tvector<SizeType,TDimension> MaxSize;
                Tvector<SizeType,TDimension> Block;

           public:

                BinBox() : Min(0), Max(0), Size(0), MaxSize(0) , Block(1) {}

                BinBox( CellType const& MinCell, CellType const& MaxCell, Tvector<SizeType,TDimension> const& _MaxSize )
                : Min(MinCell), Max(MaxCell), MaxSize(_MaxSize)
                {
                   ConstrainSize();
                   SetSize();
                   SetBlock();
                }
                
                void Set( CellType const& MinCell, CellType const& MaxCell, Tvector<SizeType,TDimension> const& _MaxSize )
                {
                   for(SizeType i = 0; i < TDimension; i++){
                      Min[i] = MinCell[i];
                      Max[i] = MaxCell[i];
                      MaxSize[i] = _MaxSize[i];
                   }
                   ConstrainSize();
                   SetSize();
                   SetBlock();
                }

                ~BinBox(){};

                BinBox const& operator++()
                {
                   for(SizeType i = 0; i < TDimension; i++){
                      if( Min[i] > static_cast<IndexType>(0) ) Min[i]--;
                      if( Max[i] < MaxSize[i]-1 ) Max[i]++;
                   }
                   SetSize();
                   return *this;
                }
                
                BinBox const& operator--()
                {
                   for(SizeType i = 0; i < TDimension; i++){
                      Min[i]++;
                      Max[i]--;
                   }
                   SetSize();
                   return *this;
                }

                void ConstrainSize()
                {
                   for( SizeType i = 0; i < TDimension; i++){
                      Min[i] = std::max( Min[i], static_cast<IndexType>(0) );
                      Max[i] = std::min( Max[i], MaxSize[i]-1 );
                   }
                }

                void SetSize()
                {
                   for(SizeType i = 0 ; i < TDimension ; i++)
                      Size[i] = static_cast<SizeType>(Max[i]-Min[i]);
                }

                void SetBlock()
                {
                   Block[0] = 1;
                   for(SizeType i = 1; i < TDimension; i++)
                      Block[i] = Block[i-1] * MaxSize[i-1];
                }

                IndexType RowBegin( SizeType const& iDim )
                {
                   return Min[iDim] * Block[iDim];
                }
                
                IndexType RowEnd( SizeType const& iDim )
                {
                   return Max[iDim] * Block[iDim];
                }

        };

        //************************************************************************

        void SearchBoxCalculate( PointType const& Point ){
           CellType Cell = CalculateCell(Point);
           SearchBox.Set( Cell, Cell, mN );
        }
        void SearchBoxCalculate( PointType const& MinBoxPoint, PointType const& MaxBoxPoint ){
           SearchBox.Set( CalculateCell(MinBoxPoint), CalculateCell(MaxBoxPoint), mN );
        }
        void SearchBoxCalculate( PointType const& ThisPoint, CoordinateType Radius ){
           SearchBox.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
        }


        //************************************************************************
        
        void SearchInRadiusRow( PointIterator const& StoreBegin, PointIterator const& StoreEnd, PointType const& ThisPoint,
              CoordinateType const& Radius, IteratorType& Result, DistanceIteratorType& Distance, SizeType const& MaxResult, SizeType& nRes )
        {
           for(PointIterator Point = StoreBegin ; (Point != StoreEnd) && (nRes < MaxResult) ; Point++){
              CoordinateType dist = TDistanceFunction()(**Point,ThisPoint);
              if( dist < Radius ){
                 *Result++   = *Point;
                 *Distance++ = dist;
                 nRes++;
              }
           }
        }
         
         //************************************************************************
        
        void SearchInRadiusRow( PointIterator const& StoreBegin, PointIterator const& StoreEnd, PointType const& ThisPoint,
              CoordinateType const& Radius, IteratorType& Result, SizeType const& MaxResult, SizeType& nRes )
        {
           for(PointIterator Point = StoreBegin ; (Point != StoreEnd) && (nRes < MaxResult) ; Point++){
              if( TDistanceFunction()(**Point,ThisPoint) < Radius ){
                 *Result++   = *Point;
                 nRes++;
              }
           }
        }
         
         //************************************************************************
         
         void SearchNearestInRow( PointIterator const& mRowBegin, PointIterator const& mRowEnd, PointType const& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance ){
                  
            for(PointIterator Point = mRowBegin ; Point != mRowEnd ; Point++){
               CoordinateType Distance = TDistanceFunction()(**Point,ThisPoint);
               if( Distance < ResultDistance ){
                  ResultPoint = *Point;
                  ResultDistance = Distance;
               }
            }

         }

         //************************************************************************
         
         void CopyPointInRow( PointIterator const& mRowBegin, PointIterator const& mRowEnd, IteratorType& ResultsPoint, SizeType const& MaxResults, SizeType& NumResults )
         {
            PointIterator Point = mRowBegin;
            while( Point != mRowEnd && NumResults < MaxResults )
               *ResultsPoint++ = *Point++;
         }
         
         //************************************************************************
         
         PointerType SearchNearestPoint( PointType const& ThisPoint )
         {
            CoordinateType ResultDistance = static_cast<CoordinateType>(1.0/DBL_EPSILON);
            PointerType Result;
            SearchNearestPoint( ThisPoint, Result, ResultDistance );
            return Result;
         }
         
         //************************************************************************
         
         PointerType SearchNearestPoint( PointType const& ThisPoint, CoordinateType ResultDistance )
         {
            ResultDistance = static_cast<CoordinateType>(1.0/DBL_EPSILON);
            PointerType Result;
            SearchNearestPoint( ThisPoint, Result, ResultDistance );
            return Result;
         }
         
         //************************************************************************
/*
         PointerType SearchNearestPoint( PointType const& ThisPoint, CoordinateType& rResultDistance = static_cast<CoordinateType>(1.0/DBL_EPSILON) )
         {
            PointerType Result;
            SearchNearestPoint( ThisPoint, Result, rResultDistance );
            return Result;
         }
*/
         //************************************************************************

         void SearchNearestPoint( PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance = static_cast<CoordinateType>(1.0/DBL_EPSILON) )
         {
            // Verification in Tree constructor
            rResult = this->NullPointer();
            if( mNumPoints == 0 )
              return;

            SearchBoxCalculate(ThisPoint,ThisPoint);

            ++SearchBox;
            SearchNearestInBox( SearchBox, ThisPoint, rResult, rResultDistance );

            while( rResult == this->NullPointer() ){
               ++SearchBox;
               SearchNearestInBox( SearchBox, ThisPoint, rResult, rResultDistance );
            }

         }

         //************************************************************************
         
         SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType Radius, IteratorType Result, DistanceIteratorType Distance, SizeType MaxResult )
         {
            CoordinateType Radius2 = Radius * Radius;
            return SearchInRadius( ThisPoint, Radius, Radius2, Result, Distance, MaxResult );
         }
         
         //************************************************************************
         
         SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType Radius, CoordinateType Radius2, IteratorType Result, DistanceIteratorType Distance, SizeType MaxResult )
         {
            SizeType nRes = 0;
            SearchBoxCalculate(ThisPoint,Radius);

            // search in sub_box
            if( TDimension == 3 )
            {
               for(IndexType III = SearchBox.RowBegin(2) ; III <= SearchBox.RowEnd(2) ; III += SearchBox.Block[2] )
                  for(IndexType II = III + SearchBox.RowBegin(1) ; II <= III + SearchBox.RowEnd(1) ; II += SearchBox.Block[1] )
                     for(IndexType I = II + SearchBox.RowBegin(0) ; I <= II + SearchBox.RowEnd(0) ; I++ )
                        SearchInRadiusRow(mPoints[I].begin(),mPoints[I].end(),ThisPoint,Radius2,Result,Distance,MaxResult,nRes);
            }
            else if( TDimension == 2 )
            {
               for(IndexType II = SearchBox.RowBegin(1) ; II <= SearchBox.RowEnd(1) ; II += SearchBox.Block[1] )
                  for(IndexType I = II + SearchBox.RowBegin(0) ; I <= II + SearchBox.RowEnd(0) ; I++ )
                     SearchInRadiusRow(mPoints[I].begin(),mPoints[I].end(),ThisPoint,Radius2,Result,Distance,MaxResult,nRes);
            }
            else if( TDimension == 1 )
            {
               for(IndexType I = SearchBox.RowBegin(0) ; I <= SearchBox.RowEnd(0) ; I++ )
                  SearchInRadiusRow(mPoints[I].begin(),mPoints[I].end(),ThisPoint,Radius2,Result,Distance,MaxResult,nRes);
            }
            else
            {
               std::cout << Info() << " Dimension " << TDimension << " Not Implemented" << std::endl;
            }

            return nRes;
         }
         
         //************************************************************************
         
         SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType Radius, IteratorType Result, SizeType MaxResult )
         {
            CoordinateType Radius2 = Radius * Radius;
            return SearchInRadius(ThisPoint,Radius,Radius2,Result,MaxResult);
         }
         
         //************************************************************************
         
         SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType Radius, CoordinateType Radius2, IteratorType Result, SizeType MaxResult )
         {
            SizeType nRes = 0;
            SearchBoxCalculate(ThisPoint,Radius);

            // search in sub_box
            if( TDimension == 3 )
            {
               for(IndexType III = SearchBox.RowBegin(2) ; III <= SearchBox.RowEnd(2) ; III += SearchBox.Block[2] )
                  for(IndexType II = III + SearchBox.RowBegin(1) ; II <= III + SearchBox.RowEnd(1) ; II += SearchBox.Block[1] )
                     for(IndexType I = II + SearchBox.RowBegin(0) ; I <= II + SearchBox.RowEnd(0) ; I++ )
                        SearchInRadiusRow(mPoints[I].begin(),mPoints[I].end(),ThisPoint,Radius2,Result,MaxResult,nRes);
            }
            else if( TDimension == 2 )
            {
               for(IndexType II = SearchBox.RowBegin(1) ; II <= SearchBox.RowEnd(1) ; II += SearchBox.Block[1] )
                  for(IndexType I = II + SearchBox.RowBegin(0) ; I <= II + SearchBox.RowEnd(0) ; I++ )
                     SearchInRadiusRow(mPoints[I].begin(),mPoints[I].end(),ThisPoint,Radius2,Result,MaxResult,nRes);
            }
            else if( TDimension == 1 )
            {
               for(IndexType I = SearchBox.RowBegin(0) ; I <= SearchBox.RowEnd(0) ; I++ )
                  SearchInRadiusRow(mPoints[I].begin(),mPoints[I].end(),ThisPoint,Radius2,Result,MaxResult,nRes);
            }
            else
            {
               std::cout << Info() << " Dimension " << TDimension << " Not Implemented" << std::endl;
            }


            return nRes;
         }
         
         //************************************************************************
        
         void SearchNearestInBox( BinBox const& Box, PointType const& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance )
         {

            // search in sub_box
            if( TDimension == 3 )
            {
               for(IndexType III = SearchBox.RowBegin(2) ; III <= SearchBox.RowEnd(2) ; III += SearchBox.Block[2] )
                  for(IndexType II = III + SearchBox.RowBegin(1) ; II <= III + SearchBox.RowEnd(1) ; II += SearchBox.Block[1] )
                     for(IndexType I = II + SearchBox.RowBegin(0) ; I <= II + SearchBox.RowEnd(0) ; I++ )
                        SearchNearestInRow( mPoints[I].begin(), mPoints[I].end(), ThisPoint, ResultPoint, ResultDistance );
            }
            else if( TDimension == 2 )
            {
               for(IndexType II = SearchBox.RowBegin(1) ; II <= SearchBox.RowEnd(1) ; II += SearchBox.Block[1] )
                  for(IndexType I = II + SearchBox.RowBegin(0) ; I <= II + SearchBox.RowEnd(0) ; I++ )
                     SearchNearestInRow( mPoints[I].begin(), mPoints[I].end(), ThisPoint, ResultPoint, ResultDistance );
            }
            else if( TDimension == 1 )
            {
               for(IndexType I = SearchBox.RowBegin(0) ; I <= SearchBox.RowEnd(0) ; I++ )
                  SearchNearestInRow( mPoints[I].begin(), mPoints[I].end(), ThisPoint, ResultPoint, ResultDistance );
            }
            else
            {
               std::cout << Info() << " Dimension " << TDimension << " Not Implemented" << std::endl;
            }


         }

         //************************************************************************
         
         SizeType SearchInBox( PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& ResultsPoint, SizeType MaxResults )
         {
            SizeType nRes = 0;
            SearchBoxCalculate(SearchMinPoint,SearchMaxPoint);

            // search in Box
            if( TDimension == 3 )
            {
               for(IndexType III = SearchBox.RowBegin(2) ; III <= SearchBox.RowEnd(2) ; III += SearchBox.Block[2] )
                  for(IndexType II = III + SearchBox.RowBegin(1) ; II <= III + SearchBox.RowEnd(1) ; II += SearchBox.Block[1] )
                     for(IndexType I = II + SearchBox.RowBegin(0) ; I <= II + SearchBox.RowEnd(0) ; I++ )
                        CopyPointInRow( mPoints[I].begin(), mPoints[I].end(), ResultsPoint, MaxResults, nRes );
            }
            else if( TDimension == 2 )
            {
               for(IndexType II = SearchBox.RowBegin(1) ; II <= SearchBox.RowEnd(1) ; II += SearchBox.Block[1] )
                  for(IndexType I = II + SearchBox.RowBegin(0) ; I <= II + SearchBox.RowEnd(0) ; I++ )
                     CopyPointInRow( mPoints[I].begin(), mPoints[I].end(), ResultsPoint, MaxResults, nRes );
            }
            else if( TDimension == 1 )
            {
               for(IndexType I = SearchBox.RowBegin(0) ; I <= SearchBox.RowEnd(0) ; I++ )
                  CopyPointInRow( mPoints[I].begin(), mPoints[I].end(), ResultsPoint, MaxResults, nRes );
            }
            else
            {
               std::cout << Info() << " Dimension " << TDimension << " Not Implemented" << std::endl;
            }

            return nRes;
         }
         
         //************************************************************************

         /// Turn back information as a string.
         virtual std::string Info() const
         {
            return "BinsDynamic";
         }
      
         /// Print information about this object.
         virtual void PrintInfo(std::ostream& rOStream) const
         {
            rOStream << "BinsDynamic";
         }

         /// Print object's data.
         virtual void PrintData(std::ostream& rOStream, std::string const& Perfix = std::string()) const
         {
			rOStream << Perfix << "Bin[" << PointerDistance(mPointBegin, mPointEnd) << "] : " << std::endl;
			for(typename CellsContainerType::const_iterator i_cell = mPoints.begin() ; i_cell != mPoints.end() ; i_cell++)
			{
			    rOStream << Perfix << "[ " ;
			    for(typename PointVector::const_iterator i_point = i_cell->begin() ; i_point != i_cell->end() ; i_point++)
				rOStream << **i_point << "    ";
			    rOStream << " ]" << std::endl;
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
         BinsDynamic& operator=(BinsDynamic const& rOther);

         /// Copy constructor.
         BinsDynamic(BinsDynamic const& rOther);

    private:

         IteratorType     mPointBegin;
         IteratorType     mPointEnd;
         
         Tvector<CoordinateType,TDimension>  mMinPoint;
         Tvector<CoordinateType,TDimension>  mMaxPoint;
         Tvector<CoordinateType,TDimension>  mCellSize;
         Tvector<CoordinateType,TDimension>  mInvCellSize;
         Tvector<SizeType,TDimension>        mN;
         SizeType                            mNumPoints;

         // Bins Access Vector ( vector<Iterator> )
         CellsContainerType mPoints;

         // Work Variables ( For non-copy of Search Variables )
         BinBox SearchBox;

	public:
	  static TreeNodeType* Construct(IteratorType PointsBegin, IteratorType PointsEnd, PointType MaxPoint, PointType MinPoint, SizeType BucketSize)
	  {
		 
		SizeType number_of_points = PointerDistance(PointsBegin,PointsEnd);
		if (number_of_points == 0)
		  return NULL;
		else 
		{
		  return new BinsDynamic( PointsBegin, PointsEnd, MinPoint, MaxPoint, BucketSize );
		}

	  }

};

   template<
      std::size_t TDimension,
      class TPointType,
      class TContainerType,
      class TPointerType,
      class TIteratorType,
      class TDistanceIteratorType,
      class TDistanceFunction >
std::ostream & operator<<( std::ostream& rOStream,
      BinsDynamic<TDimension,TPointType,TContainerType,TPointerType,TIteratorType,TDistanceIteratorType,TDistanceFunction>& rThis)
{
   rThis.PrintInfo(rOStream);
   rOStream << std::endl;
   rThis.PrintSize(rOStream);
   rThis.PrintData(rOStream);
   return rOStream;
};



};

#endif // KRATOS_BINS_DYNAMIC_CONTAINER_H_INCLUD
