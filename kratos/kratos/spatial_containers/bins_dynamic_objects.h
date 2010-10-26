/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/


//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Nelson $
//   Date:                $Date: 29-09-2010 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_BINS_DYNAMIC_OBJECTS_CONTAINER_H_INCLUDED)
#define  KRATOS_BINS_DYNAMIC_OBJECTS_CONTAINER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <cmath>
#include <algorithm>
#include <time.h>

// External includes 
//#include "boost/smart_ptr.hpp"


// Project includes
#include "tree.h"
#include "cell.h"
#include "bounding_box.h"

#ifdef _OPENMP
#include <omp.h>
#endif


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
    class TContainerType,
    class TCalculateBoundingBoxFunction,  /// funcion que crea las cajas
    class TIntersectionFunction,          /// funcion que calcula la interseccion de objetos y celdas
    class TDistanceFunction,              /// Verifica si los objetos contenidos en la celda estan en contacto  
    class TDistanceIteratorType = typename std::vector<double>::iterator,
    class TPointerType          = typename TContainerType::value_type,
    class TIteratorType         = typename TContainerType::iterator
    > 
  class BinsObjectDynamic  
    {
    public:
      ///@name Type Definitions
      ///@{
      typedef Cell<TPointerType, TContainerType> CellType;     
      typedef TPointType                         PointType;
      typedef TContainerType                     ContainerType;
      typedef TIteratorType                      IteratorType;
      typedef TDistanceIteratorType              DistanceIteratorType;
      typedef TPointerType                       PointerType;
      typedef TDistanceFunction                  DistanceFunction;
      typedef TCalculateBoundingBoxFunction      CalculateBoundingBoxFunction;  
      
      enum { Dimension = TDimension };
           
      typedef std::vector<CellType> CellContainerType;
      typedef typename CellContainerType::iterator CellContainerIterator;
      
      typedef BoundingBox< TPointType,TPointerType> BoundingBoxType;
      typedef std::vector< BoundingBoxType > BoundingBoxContainerType;
      
      
      typedef typename BoundingBoxContainerType::iterator BoundingBoxIterator; 
      
      typedef TreeNode<TDimension,TPointType,TPointerType,TIteratorType,TDistanceIteratorType> TreeNodeType;
      typedef typename TreeNodeType::CoordinateType  CoordinateType;  // double
      typedef typename TreeNodeType::SizeType        SizeType;        // std::size_t
      typedef typename TreeNodeType::IndexType       IndexType;       // std::size_t
            
      
      //typedef TreeNodeType LeafType;    
       typedef typename TreeNodeType::IteratorIteratorType IteratorIteratorType;
       typedef typename TreeNodeType::SearchStructureType  SearchStructureType;
       
      	
      /// Pointer definition of BinsObjectDynamic
      KRATOS_CLASS_POINTER_DEFINITION(BinsObjectDynamic);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      BinsObjectDynamic() {}     
      /// Constructor de bins a bounding box
      
      BinsObjectDynamic (TIteratorType const& ObjectsBegin, TIteratorType const& ObjectsEnd) 
           : mObjectsBegin(ObjectsBegin), mObjectsEnd(ObjectsEnd)
      {
	
	 CalculateBoundingBox();
	 CalculateCellSize();
	 AllocateCellsContainer();
	 GenerateBins();
	 	 
	 
      }
     

      /// Destructor.
      virtual ~BinsObjectDynamic(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      
      
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
      virtual std::string Info() const{ return "Bins " ; }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const{}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const{}
      
      CellContainerType& GetCellContainer()
       {
         return mCells; 
       }
       
      Tvector<SizeType,TDimension>& GetDivisions()
       {
	 return mN;
       }
       
      Tvector<CoordinateType,TDimension>& GetCellSize()  
      {
       return mCellSize;
      }
      

//************************************************************************   
//************************************************************************

             
      
      
            
      ///@}      
      ///@name Friends
      ///@{
      
            
      ///@}
      
      ///@name Protected static Member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Protected member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Protected Operators
      ///@{ 
      

	
	
    /// Computa los boxes de cada uno de los elementos del model part		
    void CalculateBoundingBox()
    {
      KRATOS_TRY 
      
      TPointType Low, High;  
      CalculateBoundingBoxFunction()(*mObjectsBegin,mMinPoint,mMaxPoint);
      
      std::size_t size = std::distance(mObjectsBegin, mObjectsEnd);     
      #ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
      #else
      int number_of_threads = 1;
      #endif

      vector<unsigned int> node_partition;
      CreatePartition(number_of_threads, size, node_partition);
    
      std::vector<TPointType> Max(number_of_threads);
      std::vector<TPointType> Min(number_of_threads);
       
      for(int k=0; k<number_of_threads; k++ )
      {
	 Max[k] = mMaxPoint;   
	 Min[k] = mMinPoint;
      }
      
       
      #ifdef _OPENMP
      double start_prod = omp_get_wtime();
      #endif
      
      #pragma omp parallel for  private(High, Low)  // if(size>1000000) shared(mMaxPoint, mMinPoint)
      for(int k=0; k<number_of_threads; k++)
 	{
 	  TIteratorType i_begin = mObjectsBegin + node_partition[k];
 	  TIteratorType i_end   = mObjectsBegin + node_partition[k+1];
      
          for ( TIteratorType i_object  = i_begin ; i_object != i_end ; i_object++ )
          //for (TIteratorType i_object = mObjectsBegin ; i_object != mObjectsEnd ; i_object++)
	    { 
	         CalculateBoundingBoxFunction()(*i_object, Low, High);
		 for(std::size_t i = 0 ; i < TDimension ; i++)
		 {
		   
		    Max[k][i]      =   (Max[k][i]  < High[i]) ? High[i] : Max[k][i];
		    Min[k][i]      =   (Min[k][i]  > Low[i])  ? Low[i]  : Min[k][i];
		   //mMaxPoint[i]  =   (mMaxPoint[i]  < High[i]) ? High[i] : mMaxPoint[i];
		   //mMinPoint[i]  =   (mMinPoint[i]  > Low[i])  ? Low[i]  : mMinPoint[i]; 
		 }
	    }
	}	

          for(int k=0; k<number_of_threads; k++)
 	    {
               for(std::size_t i = 0 ; i < TDimension ; i++)
		 {
		    mMaxPoint[i]  =   (mMaxPoint[i]  < Max[k][i]) ? Max[k][i] : mMaxPoint[i];
		    mMinPoint[i]  =   (mMinPoint[i]  > Min[k][i]) ? Min[k][i] : mMinPoint[i]; 
		 }
	    }
	
	
       #ifdef _OPENMP
       double stop_prod = omp_get_wtime();
       std::cout << "Time Calculating Bounding Boxes  = " << stop_prod - start_prod << std::endl;
       #endif
	 //std::cout<< mMaxPoint[0] << " " << mMaxPoint[1] << std::endl;
	 //std::cout<< mMinPoint[0] << " " << mMinPoint[1] << std::endl;  
      
      KRATOS_CATCH("")    
     }
     
    
//************************************************************************   
//************************************************************************   
    
    void CalculateCellSize()
    {
      
      KRATOS_TRY 
      
      double delta[TDimension];
      double alpha[TDimension];
      double mult_delta = 1.00;
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

      
      mN[index] = static_cast<SizeType>( pow(static_cast<CoordinateType>(SearchUtils::PointerDistance(mObjectsBegin,mObjectsEnd)/mult_delta), 1.00/TDimension) +1 );
      
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
      
      KRATOS_CATCH("")  
      
    }
    
    
//************************************************************************   
//************************************************************************
   
   void GenerateBins()
   {    
      TPointType Low, High;      
      SearchStructureType Box;
      /// LLenando las celdas con los objectos
      for(IteratorType i_object = mObjectsBegin ; i_object != mObjectsEnd ; i_object++)
           {
	     CalculateBoundingBoxFunction()(*i_object, Low, High);
             Box.Set( CalculateCell(Low), CalculateCell(High), mN );
             FillObject(Box, i_object);
           }
   }    
      
            
            
            
//************************************************************************   
//************************************************************************

    SizeType SearchObjects(PointerType& ThisObject, IteratorType& Results,  const SizeType& MaxNumberOfResults )
    {
      TPointType Low, High;      
      SearchStructureType Box;
      SizeType NumberOfResults = 0;
      CalculateBoundingBoxFunction()(ThisObject, Low, High);
      Box.Set( CalculateCell(Low), CalculateCell(High), mN );
      SearchInBoxLocal(ThisObject, Results, NumberOfResults, MaxNumberOfResults, Box );
      return NumberOfResults;
     } 
     
     

//************************************************************************   
//************************************************************************
            
                  
    // **** THREAD SAFE

    // Dimension = 1
//     void SearchInBoxLocal(const PointerType& ThisObject, IteratorType& Results, DistanceIteratorType& ResultsDistances, 
// 			  SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
//                           SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
//      {
//        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block )
//           SearchBoxInRange(mCells[I].begin(),mCells[I].end(),ThisObject,Results, ResultsDistances,NumberOfResults,MaxNumberOfResults);
//      }

    // Dimension = 2
    void SearchInBoxLocal(PointerType& ThisObject, IteratorType& Results, 
			  SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
     for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block )
        for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block ) 
           SearchBoxInRange(mCells[I].Begin(), mCells[I].End(), ThisObject, Results, NumberOfResults, MaxNumberOfResults);
    }

    // Dimension = 3
//     void SearchInBoxLocal(const PointerType& ThisObject, IteratorType& Results, DistanceIteratorType& ResultsDistances, 
// 			  SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
//                           SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
//      {
//      for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block )
//        for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block )
//  	  for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
//  	      SearchBoxInRange(mPoints[I].begin(),mPoints[I].end(),ThisObject,Results, ResultsDistances,NumberOfResults,MaxNumberOfResults);
//      }
	
	
//************************************************************************   
//************************************************************************

 void SearchBoxInRange(const IteratorType& RangeBegin, const IteratorType& RangeEnd, PointerType& ThisObject, 
		       IteratorType& Results, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults )
   {
      for(IteratorType i_object = RangeBegin ; (i_object != RangeEnd) && (NumberOfResults < MaxNumberOfResults) ; i_object++)
      {
           if(DistanceFunction()(ThisObject, *i_object))
	   {
               Results   = i_object;
               Results++;
               NumberOfResults++;
            }
      }
   }
   
   
//************************************************************************   
//************************************************************************
        
     Tvector<IndexType,TDimension>  CalculateCell( const TPointType& ThisPoint )
     {
        Tvector<IndexType,TDimension>  Cell;
        for(SizeType i = 0 ; i < TDimension ; i++)
            Cell[i] = CalculatePosition(ThisPoint[i],i);
        return Cell; 
     }
     
     
//************************************************************************   
//************************************************************************

    // Dimension = 1 
    void FillObject( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box, IteratorType i_object)
    {
        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block )
            mCells[I].Add(*i_object);
    }

   
//************************************************************************   
//************************************************************************

    // Dimension = 2
    void FillObject( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box, const IteratorType& i_object)
    {
       
       TPointType  LowPointCell;
       TPointType  HighPointCell;
       SizeType number;
       const IndexType columns = Box.Axis[1].Block;   
       for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block ) 
          {
           for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block ) 
	       {
	        number           = round(static_cast<SizeType>(II / columns));;
 		LowPointCell[0]  = (I-II) * mCellSize[0];
 		LowPointCell[1]  = number * mCellSize[1];
 		HighPointCell[0] = (I-II + 1  ) * mCellSize[0];
		HighPointCell[1] = (number + 1) * mCellSize[1];		
		if (TIntersectionFunction()(*i_object, LowPointCell, HighPointCell))
		  {
	           mCells[I].Add(*i_object);
	          }
               }
           }
    }
    
    
//************************************************************************   
//************************************************************************

    // Dimension = 2
    void FillObjectNew( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box, IteratorType i_object)
    {
       for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block )
           for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block ) 
	        mCells[I].Add(*i_object);
    }   

//************************************************************************   
//************************************************************************

    // Dimension = 3
    void FillObject( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box, IteratorType i_object)
         {
           for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block )
             for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block )
               for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
                   mCells[I].Add(*i_object);
         }

//************************************************************************
//************************************************************************


    IndexType CalculatePosition( CoordinateType const& ThisCoord, SizeType& ThisDimension )
        {
            CoordinateType d_index = (ThisCoord - mMinPoint[ThisDimension]) * mInvCellSize[ThisDimension];
            IndexType index = static_cast<IndexType>( (d_index < 0.00) ? 0.00 : d_index );
            return  (index > mN[ThisDimension]-1) ? mN[ThisDimension]-1 : index;
	    
        }



//************************************************************************
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
//************************************************************************ 

void AllocateCellsContainer() 
       {
	  SizeType Size = 1;
	  for(SizeType i = 0 ; i < TDimension ; i++)
	      Size *= mN[i];
	  ///Resize Global Container
	  mCells.resize(Size);	    
        }
             
        
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
        
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
        
      TPointType    mMinPoint;
      TPointType    mMaxPoint;  

      TIteratorType mObjectsBegin;
      TIteratorType mObjectsEnd;

      Tvector<CoordinateType,TDimension>  mCellSize;
      Tvector<CoordinateType,TDimension>  mInvCellSize;
      Tvector<SizeType,TDimension>  mN;
      
      CellContainerType mCells;  ///The bin    
      

      
 
	
	
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
      
    inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, vector<unsigned int>& partitions)
    {
      partitions.resize(number_of_threads+1);
      int partition_size = number_of_rows / number_of_threads;
      partitions[0] = 0;
      partitions[number_of_threads] = number_of_rows;
      for(unsigned int i = 1; i<number_of_threads; i++)
      partitions[i] = partitions[i-1] + partition_size ;
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
      BinsObjectDynamic& operator=(BinsObjectDynamic const& rOther){}

      /// Copy constructor.
      BinsObjectDynamic(BinsObjectDynamic const& rOther){}

        
      ///@}    
        
    }; // Class BinsObjectDynamic 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  /*
  template<
    std::size_t TDimension,
    class TPointType,
    class TContainerType,
    class TCalculateBoundingBoxFunction,
    class TIntersectionFunction,
    class TDistanceFunction, 
    class TDistanceIteratorType, = typename std::vector<double>::iterator,
    class TPointerType,          = typename TContainerType::value_type,
    class TIteratorType,         = typename TContainerType::iterator
    > 
  inline std::istream& operator >> (std::istream& rIStream, 
				    BinsObjectDynamic<TDimension, TPointerType, TPointType, TGeometryType, TBoundingBoxFunction, 
				    TIntersectionFunction, TIteratorType>& rThis);
				    

  /// output stream function
  template<
    std::size_t TDimension,
    class TPointType,
    class TContainerType,
    class TCalculateBoundingBoxFunction,
    class TIntersectionFunction,
    class TDistanceFunction, 
    class TDistanceIteratorType = typename std::vector<double>::iterator,
    class TPointerType          = typename TContainerType::value_type,
    class TIteratorType         = typename TContainerType::iterator
    > 
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const BinsObjectDynamic<TDimension, TPointerType, TPointType, TGeometryType,TBoundingBoxFunction, 
				    TIntersectionFunction, TIteratorType>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  */
  
}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined 


