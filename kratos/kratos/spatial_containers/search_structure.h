//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: antonia $
//   Date:                $Date: 2008-10-10 14:04:56 $
//   Revision:            $Revision: 1.1 $
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
    std::size_t Dimension,
    class IndexType,
    class SizeType,
    class CoordinateType,
    class IteratorType,
    class IteratorIteratorType = typename std::vector<IteratorType>::iterator >
  class AuxiliarSearchStructure
  {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of AuxiliarSearchStructure
      KRATOS_CLASS_POINTER_DEFINITION(AuxiliarSearchStructure);

      typedef Tvector<IndexType,Dimension> IndexVector;
      typedef Tvector<SizeType,Dimension> SizeVector;

    public:
      // Bin
      SubBinAxis<IndexType,SizeType> Axis[Dimension];
      IteratorIteratorType RowBegin;
      IteratorIteratorType RowEnd;  
      // KDTree
      CoordinateType distance_to_partition;
      CoordinateType distance_to_partition2;
      CoordinateType residual_distance[Dimension];
      SizeType BucketCounter;

    public:

      AuxiliarSearchStructure(){}

      AuxiliarSearchStructure( IndexVector const& Min_, IndexVector const& Max_, SizeVector const& MaxSize_, IteratorIteratorType const& IteratorBegin ) {
        Set(Min_,Max_,MaxSize_,IteratorBegin);
      }

      ~AuxiliarSearchStructure(){}

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

      AuxiliarSearchStructure const& operator++(){
        for(SizeType i = 0; i < Dimension; i++)
          ++(Axis[i]);
        return *this;
      }
      AuxiliarSearchStructure const& operator--(){
        for(SizeType i = 0; i < Dimension; i++)
          (Axis[i])--;
        return *this;
      }

  };
  
}  // namespace Kratos.

#endif
