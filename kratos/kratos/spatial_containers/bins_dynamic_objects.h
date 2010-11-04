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
//   Last Modified by:    $Author: Nelson Lafontaine $
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
#include <ctime>
#include <vector>

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
  template<class TConfigure> 
  class BinsObjectDynamic  
    {
    public:
      ///@name Type Definitions
      ///@{
       
      enum { Dimension = TConfigure::Dimension };   
	
      typedef TConfigure                                   Configure; 
      typedef typename TConfigure::PointType               PointType;
      typedef typename TConfigure::PointerType             PointerType;
      typedef typename TConfigure::ContainerType           ContainerType;
      typedef typename TConfigure::IteratorType            IteratorType;
      typedef typename TConfigure::ResultContainerType     ResultContainerType; // 
      typedef typename TConfigure::ResultIteratorType      ResultIteratorType;
        
      typedef Cell<Configure> CellType;  
      typedef std::vector<CellType> CellContainerType;
      typedef typename CellContainerType::iterator CellContainerIterator;

      
      typedef BoundingBox< PointType,PointerType> BoundingBoxType;
      typedef std::vector< BoundingBoxType > BoundingBoxContainerType;    
      typedef typename BoundingBoxContainerType::iterator BoundingBoxIterator; 
      
      typedef TreeNode<Dimension, PointType, PointerType, IteratorType,  typename TConfigure::DistanceIteratorType> TreeNodeType;
      typedef typename TreeNodeType::CoordinateType  CoordinateType;  // double
      typedef typename TreeNodeType::SizeType        SizeType;        // std::size_t
      typedef typename TreeNodeType::IndexType       IndexType;       // std::size_t
           
      
       ///Contact Pair
       typedef typename TConfigure::ContainerContactType  ContainerContactType;
       typedef typename TConfigure::IteratorContactType   IteratorContactType;
      
      ///typedef TreeNodeType LeafType;    
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
      
      BinsObjectDynamic (IteratorType const& ObjectsBegin, IteratorType const& ObjectsEnd) 
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
       
      Tvector<SizeType, Dimension>& GetDivisions()
       {
	 return mN;
       }
       
      Tvector<CoordinateType, Dimension>& GetCellSize()  
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
  
      PointType Low, High;  
      TConfigure::CalculateBoundingBox(*mObjectsBegin,mMinPoint,mMaxPoint);
      
      std::size_t size = std::distance(mObjectsBegin, mObjectsEnd);     
      #ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
      #else
      int number_of_threads = 1;
      #endif

      std::vector<unsigned int> node_partition;
      CreatePartition(number_of_threads, size, node_partition);
    
      std::vector<PointType> Max(number_of_threads);
      std::vector<PointType> Min(number_of_threads);
       
      for(int k=0; k<number_of_threads; k++ )
      {
	 Max[k] = mMaxPoint;   
	 Min[k] = mMinPoint;
      }
      
       
      #ifdef _OPENMP
      double start_prod = omp_get_wtime();
      #endif
      
      #pragma omp parallel for  private(High, Low) 
      for(int k=0; k<number_of_threads; k++)
 	{
 	  IteratorType i_begin = mObjectsBegin + node_partition[k];
 	  IteratorType i_end   = mObjectsBegin + node_partition[k+1];
      
          for ( IteratorType i_object  = i_begin ; i_object != i_end ; i_object++ )
	    { 
	         TConfigure::CalculateBoundingBox(*i_object, Low, High);
		 for(std::size_t i = 0 ; i < Dimension ; i++)
		 {
		   
		    Max[k][i]      =   (Max[k][i]  < High[i]) ? High[i] : Max[k][i];
		    Min[k][i]      =   (Min[k][i]  > Low[i])  ? Low[i]  : Min[k][i];
		 }
	    }
	}	

          for(int k=0; k<number_of_threads; k++)
 	    {
               for(std::size_t i = 0 ; i < Dimension ; i++)
		 {
		    mMaxPoint[i]  =   (mMaxPoint[i]  < Max[k][i]) ? Max[k][i] : mMaxPoint[i];
		    mMinPoint[i]  =   (mMinPoint[i]  > Min[k][i]) ? Min[k][i] : mMinPoint[i]; 
		 }
	    }
	
       #ifdef _OPENMP
       double stop_prod = omp_get_wtime();
       std::cout << "Time Calculating Bounding Boxes  = " << stop_prod - start_prod << std::endl;
       #endif
         
     }
     
    
//************************************************************************   
//************************************************************************   
    
    void CalculateCellSize()
    {
      
      
      double delta[Dimension];
      double alpha[Dimension];
      double mult_delta = 1.00;
      SizeType index = 0;
      for(SizeType i = 0 ; i < Dimension ; i++) {
         delta[i] = mMaxPoint[i] - mMinPoint[i];
         if ( delta[i] > delta[index] )
            index = i;
         delta[i] = (delta[i] == 0.00) ? 1.00 : delta[i];
      }

      for(SizeType i = 0 ; i < Dimension ; i++){
         alpha[i] = delta[i] / delta[index];
         mult_delta *= alpha[i];
      }

      
      mN[index] = static_cast<SizeType>( pow(static_cast<CoordinateType>(SearchUtils::PointerDistance(mObjectsBegin,mObjectsEnd)/mult_delta), 1.00/Dimension) +1 );
      
      for(SizeType i = 0 ; i < Dimension ; i++){
      if(i!=index) {
      mN[i] = static_cast<SizeType>(alpha[i] * mN[index]);
      mN[i] = ( mN[i] == 0 ) ? 1 : mN[i];
      }
      }

      for(SizeType i = 0 ; i < Dimension ; i++){
      mCellSize[i] = delta[i] / mN[i];
      mInvCellSize[i] = 1.00 / mCellSize[i];
      }
      
      
    }
    
    
//************************************************************************   
//************************************************************************
   
   void GenerateBins()
   {    
      PointType Low, High;      
      SearchStructureType Box;
      /// LLenando las celdas con los objectos
      for(IteratorType i_object = mObjectsBegin ; i_object != mObjectsEnd ; i_object++)
           {
	     TConfigure::CalculateBoundingBox(*i_object, Low, High);
             Box.Set( CalculateCell(Low), CalculateCell(High), mN );
             FillObject(Box, i_object);
           }
   }    
      
            
            
            
//************************************************************************   
//************************************************************************

    SizeType SearchObjects(PointerType& ThisObject, ResultIteratorType& Result,  const SizeType& MaxNumberOfResults )
    {
      PointType Low, High;      
      SearchStructureType Box;
      SizeType NumberOfResults = 0;
      TConfigure::CalculateBoundingBox(ThisObject, Low, High);
      Box.Set( CalculateCell(Low), CalculateCell(High), mN );
      SearchInBoxLocal(ThisObject, Result, NumberOfResults, MaxNumberOfResults, Box );
      return NumberOfResults;
     } 
       

//************************************************************************   
//************************************************************************
    
    SizeType SearchObjects(PointerType& ThisObject, ResultContainerType& Result)
    {
      PointType Low, High;      
      SearchStructureType Box;
      TConfigure::CalculateBoundingBox(ThisObject, Low, High);
      Box.Set( CalculateCell(Low), CalculateCell(High), mN );
      SearchInBoxLocal(ThisObject, Result, Box );
      return Result.size();
     }          
            
 //************************************************************************   
//************************************************************************          
          
    void SearchContact(ContainerContactType& Result)
    {
       for (CellContainerIterator icell = mCells.begin() ; icell!= mCells.end(); icell++)
        {
           icell->SearchContact(Result);
        }
       return; 
     }          
            
 //************************************************************************   
//************************************************************************       
          
      SizeType SearchContact(IteratorContactType& Result, const SizeType& MaxNumberOfResults )
    {
       SizeType NumberOfResults = 0;
       for (CellContainerIterator icell = mCells.begin() ; icell!= mCells.end(); icell++)
        {
           icell->SearchContact(Result, NumberOfResults, MaxNumberOfResults);
        }
       return NumberOfResults; 
     }          
            
 //************************************************************************   
//************************************************************************         
          
          
                  
    // **** THREAD SAFE
    
    // Dimension = 1
    void SearchInBoxLocal(PointerType& ThisObject, ResultIteratorType& Result, 
			  SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
     {
       for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block )
           mCells[I].SearchObjects(ThisObject, Result, NumberOfResults, MaxNumberOfResults);
     }

    // Dimension = 2
    void SearchInBoxLocal(PointerType& ThisObject, ResultIteratorType& Result, 
			  SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
     for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block )
        for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block ) 
	   mCells[I].SearchObjects(ThisObject, Result, NumberOfResults, MaxNumberOfResults);
    } 

    // Dimension = 3
    void SearchInBoxLocal(PointerType& ThisObject, ResultIteratorType& Result, 
			  SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
     {
     for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block )
       for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block )
 	  for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
 	       mCells[I].SearchObjects(ThisObject, Result, NumberOfResults, MaxNumberOfResults);
     }

 //************************************************************************   
//************************************************************************          
            
                  
    // **** THREAD SAFE
    
    // Dimension = 1
    void SearchInBoxLocal(PointerType& ThisObject, ResultContainerType& Result, 
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
     {
       for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block )
           mCells[I].SearchObjects(ThisObject, Result);
     }

    // Dimension = 2
    void SearchInBoxLocal(PointerType& ThisObject, ResultContainerType& Result, 
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
     for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block )
        for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block ) 
	   mCells[I].SearchObjects(ThisObject, Result);
    } 

    // Dimension = 3
    void SearchInBoxLocal(PointerType& ThisObject, ResultContainerType& Result, 
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
     {
     for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block )
       for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block )
 	  for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
 	       mCells[I].SearchObjects(ThisObject, Result);
     }



//************************************************************************   
//************************************************************************
        
     Tvector<IndexType,Dimension>  CalculateCell( const PointType& ThisPoint )
     {
        Tvector<IndexType,Dimension>  Cell;
        for(SizeType i = 0 ; i < Dimension ; i++)
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
       
       PointType  LowPointCell;
       PointType  HighPointCell;
       SizeType number;
       const IndexType columns = Box.Axis[1].Block;   
       for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block ) 
          {
           for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block ) 
	       {
	        number           =  static_cast<SizeType>(II / columns);;
 		LowPointCell[0]  = (I-II) * mCellSize[0];
 		LowPointCell[1]  = number * mCellSize[1];
 		HighPointCell[0] = (I-II + 1  ) * mCellSize[0];
		HighPointCell[1] = (number + 1) * mCellSize[1];		
		if (TConfigure::IntersectionBox(*i_object, LowPointCell, HighPointCell))
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
            for(SizeType iDim = Dimension-1 ; iDim > 0 ; iDim--){
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
            for(SizeType iDim = Dimension-1 ; iDim > 0 ; iDim--){
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
	  for(SizeType i = 0 ; i < Dimension ; i++)
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
        
      PointType    mMinPoint;
      PointType    mMaxPoint;  

      IteratorType mObjectsBegin;
      IteratorType mObjectsEnd;

      Tvector<CoordinateType,Dimension>  mCellSize;
      Tvector<CoordinateType,Dimension>  mInvCellSize;
      Tvector<SizeType,Dimension>  mN;
      
      CellContainerType mCells;  ///The bin    
      

      
 
	
	
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
      
    inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, std::vector<unsigned int>& partitions)
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
 template<class TConfigure> 
 inline std::istream& operator >> (std::istream& rIStream, 
				    BinsObjectDynamic<TConfigure>& rThis){return rIStream;  }
				    

  /// output stream function
  template<class TConfigure>  
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const BinsObjectDynamic<TConfigure> & rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
 
  
}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined 


