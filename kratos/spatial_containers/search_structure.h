//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: clabra $
//   Date:                $Date: 2007-03-27 17:02:19 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_SEARCH_STRUCTURE_H_INCLUDED )
#define  KRATOS_SEARCH_STRUCTURE_H_INCLUDED
// System includes
#include <vector>
// External includes 

// Project includes

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
  template< class T, std::size_t TDimension >
  class Tvector {
    private:
      T data[TDimension];
    public:
      Tvector(){
        for(std::size_t i = 0; i < TDimension; i++)
          data[i] = T(0);
      }
      Tvector( T const& value ){
        for(std::size_t i = 0; i < TDimension; i++)
          data[i] = value;
      }
      Tvector( Tvector<T,TDimension> const& Other ){
        for(std::size_t i = 0; i < TDimension; i++)
          data[i] = Other[i];
      }
      ~Tvector(){}
      T& operator[](std::size_t index) { return data[index]; }
      T const& operator[](std::size_t index) const { return data[index]; }
  };


  /// SEARCH UTILS ////

   namespace SearchUtils
   {

    
     template< class TPointerType >
       std::size_t PointerDistance( TPointerType const& PointerBegin, TPointerType const& PointerEnd )
       {
         //return std::distance(PointerBegin,PointerEnd);
         return ( PointerEnd - PointerBegin );  // required for SUN Compiler
       };
  

     template< std::size_t TDimension, class TPointType >
     bool PointInBox( TPointType const& BoxMinPoint, TPointType const& BoxMaxPoint, TPointType const& ThisPoint )
      {
        for(std::size_t i = 0 ; i < TDimension ; i++)
          if( ThisPoint[i] < BoxMinPoint[i] || ThisPoint[i] > BoxMaxPoint[i] )
            return false;
        return true;
      }
    
     template< std::size_t TDimension, class TPointType >
       class SquaredDistanceFunction
       {
         public:
           double operator()( TPointType const& p1, TPointType const& p2 ){
             double tmp = p1[0] - p2[0];
             double dist = tmp*tmp;
             for( std::size_t i = 1 ; i < TDimension ; i++){
               tmp = p1[i] - p2[i];
               dist += tmp*tmp;
             }
             return dist;
           }
       };

     template< class TIteratorType, class TSizeType >
       class CopyRange
       {
         public:
           void operator()( TIteratorType const& RangeBegin, TIteratorType const& RangeEnd, TIteratorType& Results,
               TSizeType& NumberOfResults, TSizeType const& MaxNumberOfResults )
           {
             for( TIteratorType Point = RangeBegin; (Point != RangeEnd)&&(NumberOfResults < MaxNumberOfResults); Point++)
             {
               *Results++ = *Point;
               NumberOfResults++;
             }
           }
       };

     template< class TPointType, class TIteratorType, class TSizeType, std::size_t TDimension >
       class SearchBoxInRange
       {
         public:
           void operator()( TPointType const& MinBoxPoint, TPointType const& MaxBoxPoint, TIteratorType const& RangeBegin, TIteratorType const& RangeEnd,
                            TIteratorType& Results, TSizeType& NumberOfResults, TSizeType const& MaxNumberOfResults )
           {
             for( TIteratorType Point = RangeBegin; (Point != RangeEnd)&&(NumberOfResults < MaxNumberOfResults); Point++)
               if ( PointInBox<TDimension,TPointType>(MinBoxPoint,MaxBoxPoint,**Point) )
               {
                 *Results++ = *Point;
                 NumberOfResults++;
               }
           }
       };


     template< class TPointType, class TPointerType, class TIteratorType, class TDistanceFunction, class TCoordinateType >
       class SearchNearestInRange
       {
         public:
           void operator()( const TIteratorType& RangeBegin, const TIteratorType& RangeEnd, const TPointType& ThisPoint, TPointerType& Result, TCoordinateType& Distance )
           {
             TCoordinateType NewDistance;
             for(TIteratorType Point = RangeBegin ; Point != RangeEnd ; Point++){
               NewDistance = TDistanceFunction()(**Point,ThisPoint);
               if( NewDistance < Distance ){
                 Result = *Point;
                 Distance = NewDistance;
               }
             }
           }
       };


     template< class TPointType, class TIteratorType, class TDistanceIteratorType, class TDistanceFunction, class TSizeType, class TCoordinateType >
       class SearchRadiusInRange
       {
         public:

           void operator()( TIteratorType const& RangeBegin, TIteratorType const& RangeEnd, TPointType const& ThisPoint, TCoordinateType const& Radius,
               TIteratorType& Results, TSizeType& NumberOfResults, TSizeType const& MaxNumberOfResults )
           {
             TCoordinateType distance;
             for(TIteratorType Point = RangeBegin ; (Point != RangeEnd) && (NumberOfResults < MaxNumberOfResults) ; Point++){
               distance = TDistanceFunction()(**Point,ThisPoint); // squared distance function
               if( distance < Radius ){
                 *Results++   = *Point;
                 NumberOfResults++;
               }
             }
           }

           void operator()( TIteratorType const& RangeBegin, TIteratorType const& RangeEnd, TPointType const& ThisPoint, TCoordinateType const& Radius,
               TIteratorType& Results, TDistanceIteratorType& Distances, TSizeType& NumberOfResults, TSizeType const& MaxNumberOfResults )
           {
             TCoordinateType distance;
             for(TIteratorType Point = RangeBegin ; (Point != RangeEnd) && (NumberOfResults < MaxNumberOfResults) ; Point++){
               distance = TDistanceFunction()(**Point,ThisPoint); // squared distance function
               if( distance < Radius ){
                 *Results++   = *Point;
                 *Distances++ = distance;
                 NumberOfResults++;
               }
             }
           }
       };


   };



  /// TOOLS UTILS ///

  template< class IndexType, class SizeType>
  class SubBinAxis
  {
    public:
      IndexType Min;
      IndexType Max;
      IndexType MaxSize;
      IndexType Block;

      SubBinAxis() : Min(0), Max(0), MaxSize(0), Block(1) {}

      SubBinAxis(IndexType const& Min_, IndexType const& Max_, IndexType const& MaxSize_, IndexType const& Block_)
      {
        Set(Min_,Max_,MaxSize_,Block_);
      }
      ~SubBinAxis() {}
      void Set(IndexType const& iCell, IndexType const& MaxSize_, IndexType const& Block_)
      {
        Set(iCell,iCell,MaxSize_,Block_);
      }
      void Set(IndexType const& Min_, IndexType const& Max_, IndexType const& MaxSize_, IndexType const& Block_)
      {
        MaxSize = MaxSize_;
        Block = Block_;
        Min = std::max( Min_, static_cast<IndexType>(0) );
        Max = std::min( Max_, MaxSize-1 );
      }
      IndexType Begin() { return Min*Block; }
      IndexType End()   { return Max*Block; }
      SizeType Size()   { return static_cast<SizeType>(Max-Min); }
      SubBinAxis const& operator++()
      {
        if( Min > static_cast<IndexType>(0) ) Min--;
        if( Max < MaxSize-1 ) Max++;
        return *this;
      }
      SubBinAxis const& operator--()
      {
        Min++; Max--;
        return *this;
      }
  };


  template< 
    class IndexType,
    class SizeType,
    class CoordinateType,
    class IteratorType,
    class IteratorIteratorType,
    std::size_t Dimension >
  class SearchStructure
  {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of SearchStructure
      KRATOS_CLASS_POINTER_DEFINITION(SearchStructure);

      typedef Tvector<IndexType,Dimension> IndexVector;
      typedef Tvector<SizeType,Dimension> SizeVector;

      typedef SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,Dimension> ThisType;

    public:
      // Bin
      //SubBinAxis<IndexType,SizeType> Axis[Dimension];
      SubBinAxis<IndexType,SizeType> Axis[3];
      IteratorIteratorType RowBegin;
      IteratorIteratorType RowEnd;  
      // KDTree
      CoordinateType distance_to_partition;
      CoordinateType distance_to_partition2;
      CoordinateType residual_distance[Dimension];
      SizeType BucketCounter;

    public:

      SearchStructure(){}

      SearchStructure( IndexVector const& Min_, IndexVector const& Max_, SizeVector const& MaxSize_, IteratorIteratorType const& IteratorBegin ) {
        Set(Min_,Max_,MaxSize_,IteratorBegin);
      }

      ~SearchStructure(){}

      void Set( IndexVector const& IndexCell, SizeVector const& _MaxSize, IteratorIteratorType const& IteratorBegin ){
        Set( IndexCell, IndexCell, _MaxSize, IteratorBegin );
      }

      void Set( IndexVector const& Min_, IndexVector const& Max_, SizeVector const& MaxSize_, IteratorIteratorType const& IteratorBegin ){
        IndexType Block = 1;
        Axis[0].Set(Min_[0],Max_[0],MaxSize_[0],Block);
        for(SizeType i = 1; i < Dimension; i++){
          Block *= MaxSize_[i-1];
          Axis[i].Set(Min_[i],Max_[i],MaxSize_[i],Block);
        }
        RowBegin = IteratorBegin + Axis[0].Min;
        RowEnd   = IteratorBegin + Axis[0].Max + 1;
      }

      SearchStructure const& operator++(){
        for(SizeType i = 0; i < Dimension; i++)
          ++(Axis[i]);
        return *this;
      }
      SearchStructure const& operator--(){
        for(SizeType i = 0; i < Dimension; i++)
          (Axis[i])--;
        return *this;
      }

  };

}  // namespace Kratos.

#endif
