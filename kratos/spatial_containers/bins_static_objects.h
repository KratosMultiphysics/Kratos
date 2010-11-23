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


#if !defined(KRATOS_BINS_STATIC_OBJECTS_CONTAINER_H_INCLUDED)
#define  KRATOS_BINS_DYNAMIC_OBJECTS_CONTAINER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <cmath>
#include <algorithm>
//#include <time.h>

// Project includes
#include "tree.h"
//#include "cell.h"

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
  class BinsObjectStatic  
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
        
      //typedef Cell<Configure> CellType;  
      //typedef std::vector<CellType> CellContainerType;
      //typedef typename CellContainerType::iterator CellContainerIterator;

      
      /*       typedef BoundingBox< PointType,PointerType> BoundingBoxType; */
      /*       typedef std::vector< BoundingBoxType > BoundingBoxContainerType;     */
      /*       typedef typename BoundingBoxContainerType::iterator BoundingBoxIterator;  */
      
      typedef TreeNode<Dimension, PointType, PointerType, IteratorType,  typename TConfigure::DistanceIteratorType> TreeNodeType;
      typedef typename TreeNodeType::CoordinateType  CoordinateType;  // double
      typedef typename TreeNodeType::SizeType        SizeType;        // std::size_t
      typedef typename TreeNodeType::IndexType       IndexType;       // std::size_t
           
  
	 typedef std::vector<PointerType>        PointVector;
	 typedef typename PointVector::iterator  PointIterator;
  
     typedef std::vector<PointerType>        LocalContainerType;
     typedef typename LocalContainerType::iterator LocalIteratorType;

       ///Contact Pair
       typedef typename TConfigure::ContainerContactType  ContainerContactType;
       typedef typename TConfigure::IteratorContactType IteratorContactType;
      
      ///typedef TreeNodeType LeafType;    
       typedef typename TreeNodeType::IteratorIteratorType IteratorIteratorType;
       typedef typename TreeNodeType::SearchStructureType  SearchStructureType;
  
	 typedef std::vector<IteratorType>          IteratorVector;
	 typedef typename IteratorVector::iterator  IteratorIterator;
	 typedef typename IteratorVector::const_iterator IteratorConstIterator;
       
      typedef std::vector<IndexType> IndexContainer;
	 typedef typename IndexContainer::iterator  IndexIterator;

      /// Pointer definition of BinsObjectStatic
      KRATOS_CLASS_POINTER_DEFINITION(BinsObjectStatic);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      BinsObjectStatic() {}     
      /// Constructor de bins a bounding box
      
      BinsObjectStatic (IteratorType const& ObjectsBegin, IteratorType const& ObjectsEnd) 
           : mObjectsBegin(ObjectsBegin), mObjectsEnd(ObjectsEnd)
      {
	
	 CalculateBoundingBox();
	 CalculateCellSize();
	 //AllocateCellsContainer();
	 GenerateBins();
	 	 
	 
      }
     

      /// Destructor.
      virtual ~BinsObjectStatic(){}
      

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
      virtual std::string Info() const
      {
        return "BinsObjectStatic : ";
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
      {
        rOStream << Info();
      }

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
      {
        // Container Size
        rOStream << "   Container Size: ";
        for(SizeType i = 0 ; i < Dimension ; i++)
          rOStream << "[" << mN[i] << "]";
        rOStream << std::endl;
        // CellSize
        rOStream << "   Cell Size: ";
        for(SizeType i = 0 ; i < Dimension ; i++)
          rOStream << "[" << mCellSize[i] << "]";
        rOStream << std::endl;
        rOStream << "   Contained Objects: " << SearchUtils::PointerDistance(mObjectsBegin,mObjectsEnd) << std::endl;
        rOStream << "   Total Object Storaged: " << mObjectList.size() << std::endl;
      }
         
      /// Print Size of Container
      void PrintSize( std::ostream& rout ){
        rout << "  Container Size: ";
        for(SizeType i = 0 ; i < Dimension ; i++)
          rout << "[" << this->mN[i] << "]";
        rout << std::endl;
      }

      /// Print Limits Points of the Container
      void PrintBox( std::ostream& rout ){
        rout << " BinsBox: Min [";  mMinPoint.Print(rout);
        rout <<       "];  Max [";  mMaxPoint.Print(rout);
        rout <<       "];  Size ["; mCellSize.Print(rout);
        rout << "]" << std::endl;
      }


/*      
      CellContainerType& GetCellContainer()
       {
         return mCells; 
       }
*/       
      Tvector<SizeType, Dimension>& GetDivisions()
       {
	 return mN;
       }
/*       
      Tvector<CoordinateType, Dimension>& GetCellSize()  
      {
       return mCellSize;
      }
*/      

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
      
      std::size_t size = SearchUtils::PointerDistance(mObjectsBegin,mObjectsEnd);  // std::distance(mObjectsBegin, mObjectsEnd);     
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

      SizeType n_objects = SearchUtils::PointerDistance(mObjectsBegin,mObjectsEnd);

      // Allocate CellAcess

      SizeType Size = 1;
      for(SizeType i = 0 ; i < Dimension ; i++)
        Size *= mN[i];
      mObjectsAccess.resize(Size+1,0);

/*
	   PointVector TempObjects(n_objects); //(mPointBegin,mPointEnd);
       IteratorType iter_object = mObjectsBegin;
       for(PointIterator iter = TempObjects.begin(); iter != TempObjects.end(); iter++, iter_object++)
         *iter = *iter_object;

	   // Reset index vector
	   for( IteratorIterator Iter = mIndexCell.begin(); Iter != mIndexCell.end(); Iter++)
		 *Iter = mObjectsBegin;
*/
	   // Update storage counter, storing ahead
	   for( IteratorType i_object = mObjectsBegin ; i_object != mObjectsEnd ; i_object++)
       {
         TConfigure::CalculateBoundingBox(*i_object,Low,High);
         Box.Set( CalculateCell(Low), CalculateCell(High), mN );
         CountObject(Box,*i_object);
       }

	   // Storage/reshufing pass 1

	   // Update storage counter and store
	   for( IndexIterator cell = mObjectsAccess.begin()+1 ; cell != mObjectsAccess.end() ; cell++)
		 *cell += *(cell-1);

       mObjectList.resize(mObjectsAccess[Size]);

	   // Point pass 2
	   // Store the points in lbin1

	   // Update storage counter, storing in lbin1
	   for( IteratorType i_object = mObjectsBegin ; i_object != mObjectsEnd ; i_object++)
       {
         TConfigure::CalculateBoundingBox(*i_object,Low,High);
         Box.Set( CalculateCell(Low), CalculateCell(High), mN );
         FillObject(Box,*i_object);
       }

	   // Storage/reshuffing pass 2

	   // Loop over bins, in reverse order
	   for(IndexIterator Iter = mObjectsAccess.end()-1; Iter != mObjectsAccess.begin(); Iter--)
		 *Iter = *(Iter-1);
	   mObjectsAccess[0] = 0;

   }    
      
            
            
            
//************************************************************************   
//************************************************************************

    // Dimension = 1
    void CountObject( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box, PointerType object)
    {
      PointType  MinCell, MaxCell;

      MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  // 
      MaxCell[0] = MinCell[0] + mCellSize[0];

      for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] )
        if(TConfigure::IntersectionBox(object,MinCell,MaxCell))
          mObjectsAccess[I+1]++;
    }

//************************************************************************   
//************************************************************************

    // Dimension = 2
    void CountObject( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box, PointerType object )
    {
      PointType  MinCell, MaxCell;
      PointType  MinBox, MaxBox;

      for(SizeType i = 0; i < 2; i++){
        MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  // 
        MaxBox[i] = MinBox[i] + mCellSize[i];
      }

      MinCell[1] = MinBox[1];
      MaxCell[1] = MaxBox[1];
      for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1]+=mCellSize[1], MaxCell[1]+=mCellSize[1] )
      {
        MinCell[0] = MinBox[0];
        MaxCell[0] = MaxBox[0];
        for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] )
          if(TConfigure::IntersectionBox(object,MinCell,MaxCell))
            mObjectsAccess[I+1]++;
      }
    }

//************************************************************************   
//************************************************************************

    // Dimension = 3
    void CountObject( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box, PointerType object )
    {
      PointType  MinCell, MaxCell;
      PointType  MinBox, MaxBox;

      for(SizeType i = 0; i < 3; i++){
        MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  // 
        MaxBox[i] = MinBox[i] + mCellSize[i];
      }

      MinCell[2] = MinBox[2];
      MaxCell[2] = MaxBox[2];
      for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block, MinCell[2]+=mCellSize[2], MaxCell[2]+=mCellSize[2] )
      {
        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];
        for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1]+=mCellSize[1], MaxCell[1]+=mCellSize[1] )
        {
          MinCell[0] = MinBox[0];
          MaxCell[0] = MaxBox[0];
          for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] )
            if(TConfigure::IntersectionBox(object,MinCell,MaxCell))
              mObjectsAccess[I+1]++;
        }
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
      SearchInBoxLocal_(ThisObject, Result, NumberOfResults, MaxNumberOfResults, Box );
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
/*          
    void SearchContact(ContainerContactType& Results)
    {
      for (CellContainerIterator icell = mCells.begin() ; icell!= mCells.end(); icell++)
        if(icell->Size()>1)
          icell->SearchContact(Results);
    }          
*/            
 //************************************************************************   
//************************************************************************       
/*          
    SizeType SearchContact(IteratorContactType& Result, const SizeType& MaxNumberOfResults )
    {
      SizeType NumberOfResults = 0;
      for (CellContainerIterator icell = mCells.begin() ; icell!= mCells.end(); icell++)
        if(icell->Size()>1)
          icell->SearchContact(Result, NumberOfResults, MaxNumberOfResults);
      return NumberOfResults; 
    }          
*/            
 //************************************************************************   
//************************************************************************         
      
    void SearchObjectRow(PointerType& ThisObject, LocalIteratorType RowBegin, LocalIteratorType RowEnd, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults)
    {
      for(LocalIteratorType iter = RowBegin ; iter != RowEnd && NumberOfResults < MaxNumberOfResults ; iter++)
      {
        if(TConfigure::Intersection(ThisObject,*iter))
        {
          if( std::find(Result-NumberOfResults, Result, *iter) == Result ) 
          {
            *Result = *iter;
            Result++;
            NumberOfResults++;
          }
        }
      }
    }  
                  
    // **** THREAD SAFE
    
    // Dimension = 1
    void SearchInBoxLocal(PointerType& ThisObject, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
        SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
     {
       SearchObjectRow(ThisObject,mObjectList.begin()+mObjectsAccess[Box.Axis[0].Begin()],mObjectList.begin()+mObjectsAccess[Box.Axis[0].End()+1],Result,NumberOfResults,MaxNumberOfResults);
     }

    // Dimension = 2
    void SearchInBoxLocal(PointerType& ThisObject, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
     for(IndexType I = Box.Axis[1].Begin() ; I <= Box.Axis[1].End() ; I += Box.Axis[1].Block )
       SearchObjectRow(ThisObject,mObjectList.begin()+mObjectsAccess[I+Box.Axis[0].Begin()],mObjectList.begin()+mObjectsAccess[I+Box.Axis[0].End()+1],Result,NumberOfResults,MaxNumberOfResults);
    } 

    // Dimension = 3
    void SearchInBoxLocal(PointerType& ThisObject, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {
      for(IndexType II = Box.Axis[2].Begin() ; II <= Box.Axis[2].End() ; II += Box.Axis[2].Block )
        for(IndexType I = II + Box.Axis[1].Begin() ; I <= II + Box.Axis[1].End() ; I += Box.Axis[1].Block )
          SearchObjectRow(ThisObject,mObjectList.begin()+mObjectsAccess[I+Box.Axis[0].Begin()],mObjectList.begin()+mObjectsAccess[I+Box.Axis[0].End()+1],Result,NumberOfResults,MaxNumberOfResults);
    }

    // Dimension = 3
    void SearchInBoxLocal_(PointerType& ThisObject, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {
      PointType  MinCell, MaxCell;
      PointType  MinBox, MaxBox;
      IndexType  objects_begin, objects_end;

      for(SizeType i = 0; i < 3; i++){
        MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  // 
        MaxBox[i] = MinBox[i] + mCellSize[i];
      }
      CoordinateType MaxBox_ = static_cast<CoordinateType>(Box.Axis[0].Min+1) * mCellSize[0] + mMinPoint[0];  // 

      MinCell[2] = MinBox[2];
      MaxCell[2] = MaxBox[2];
      for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block, MinCell[2]+=mCellSize[2], MaxCell[2]+=mCellSize[2] )
      {
        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];
        for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1]+=mCellSize[1], MaxCell[1]+=mCellSize[1] )
        {
          MinCell[0] = MinBox[0];
          MaxCell[0] = MaxBox[0];
          objects_begin = mObjectsAccess[II + Box.Axis[0].Begin()];
          for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] ){
            if(TConfigure::IntersectionBox(ThisObject,MinCell,MaxCell)){
              objects_begin = mObjectsAccess[I];
              break;
            }
          }
          MinCell[0] = MaxBox_-mCellSize[0];
          MaxCell[0] = MaxBox_;
          objects_end = mObjectsAccess[II+Box.Axis[0].End()+1];
          for(IndexType I = II + Box.Axis[0].End() ; I >= II + Box.Axis[0].Begin() ; I -= Box.Axis[0].Block, MinCell[0]-=mCellSize[0], MaxCell[0]-=mCellSize[0] ){
            if(TConfigure::IntersectionBox(ThisObject,MinCell,MaxCell)){
              objects_end = mObjectsAccess[I+1];
              break;
            }
          }
          SearchObjectRow(ThisObject,mObjectList.begin()+objects_begin,mObjectList.begin()+objects_end,Result,NumberOfResults,MaxNumberOfResults);
        }
      }
    }

 //************************************************************************   
//************************************************************************          
            
      
    void SearchObjectRow(PointerType& ThisObject, LocalIteratorType RowBegin, LocalIteratorType RowEnd, ResultContainerType& Results )
    {
      for(LocalIteratorType iter = RowBegin ; iter != RowEnd ; iter++)
      {
        if(TConfigure::Intersection(ThisObject,*iter))
          if( std::find(Results.begin(), Results.end(), *iter) == Results.end() ) 
            Results.push_back(*iter);
      }
    }  
                  
    // **** THREAD SAFE
    
    // Dimension = 1
    void SearchInBoxLocal(PointerType& ThisObject, ResultContainerType& Results, 
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
    {
      SearchObjectRow(ThisObject,mObjectList.begin()+mObjectsAccess[Box.Axis[0].Begin()],mObjectList.begin()+mObjectsAccess[Box.Axis[0].End()+1],Results);
    }

    // Dimension = 2
    void SearchInBoxLocal(PointerType& ThisObject, ResultContainerType& Results, 
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
     for(IndexType I = Box.Axis[1].Begin() ; I <= Box.Axis[1].End() ; I += Box.Axis[1].Block )
       SearchObjectRow(ThisObject,mObjectList.begin()+mObjectsAccess[I+Box.Axis[0].Begin()],mObjectList.begin()+mObjectsAccess[I+Box.Axis[0].End()+1],Results);
    } 

    // Dimension = 3
    void SearchInBoxLocal(PointerType& ThisObject, ResultContainerType& Results, 
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {
      for(IndexType II = Box.Axis[2].Begin() ; II <= Box.Axis[2].End() ; II += Box.Axis[2].Block )
        for(IndexType I = II + Box.Axis[1].Begin() ; I <= II + Box.Axis[1].End() ; I += Box.Axis[1].Block )
          SearchObjectRow(ThisObject,mObjectList.begin()+mObjectsAccess[I+Box.Axis[0].Begin()],mObjectList.begin()+mObjectsAccess[I+Box.Axis[0].End()+1],Results);
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
    void FillObject( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box, const PointerType& object)
    {
      PointType  MinCell, MaxCell;

      MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  // 
      MaxCell[0] = MinCell[0] + mCellSize[0];
      for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] )
        if(TConfigure::IntersectionBox(object,MinCell,MaxCell))
          mObjectList[mObjectsAccess[I]++] = object;
    }
   

    // Dimension = 2
    void FillObject( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box, const PointerType& object)
    {
      PointType  MinCell, MaxCell;
      PointType  MinBox, MaxBox;

      for(SizeType i = 0; i < 2; i++){
        MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  // 
        MaxBox[i] = MinBox[i] + mCellSize[i];
      }

      MinCell[1] = MinBox[1];
      MaxCell[1] = MaxBox[1];
      for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1]+=mCellSize[1], MaxCell[1]+=mCellSize[1] )
      {
        MinCell[0] = MinBox[0];
        MaxCell[0] = MaxBox[0];
        for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] )
          if(TConfigure::IntersectionBox(object,MinCell,MaxCell))
            mObjectList[mObjectsAccess[I]++] = object;
      }
    }


    // Dimension = 3
    void FillObject( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box, const PointerType&  object)
    {
      PointType  MinCell, MaxCell;
      PointType  MinBox, MaxBox;

      for(SizeType i = 0; i < 3; i++){
        MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  // 
        MaxBox[i] = MinBox[i] + mCellSize[i];
      }

      MinCell[2] = MinBox[2];
      MaxCell[2] = MaxBox[2];
      for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block, MinCell[2]+=mCellSize[2], MaxCell[2]+=mCellSize[2] )
      {
        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];
        for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1]+=mCellSize[1], MaxCell[1]+=mCellSize[1] )
        {
          MinCell[0] = MinBox[0];
          MaxCell[0] = MaxBox[0];
          for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] )
            if(TConfigure::IntersectionBox(object,MinCell,MaxCell))
              mObjectList[mObjectsAccess[I]++] = object;
        }
      }
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
      
      //CellContainerType mCells;  ///The bin    
      
	 // Bins Access Vector ( vector<Iterator> )
	 //IteratorVector           mIndexCell;
	 //IteratorIterator         mIndexCellBegin;
	 //IteratorIterator         mIndexCellEnd;
     LocalContainerType       mObjectList;

     IndexContainer          mObjectsAccess;
      
 
	
	
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
      BinsObjectStatic& operator=(BinsObjectStatic const& rOther){}

      /// Copy constructor.
      BinsObjectStatic(BinsObjectStatic const& rOther){}

        
      ///@}    
      
     
    }; // Class BinsObjectStatic 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
 template<class TConfigure> 
 inline std::istream& operator >> (std::istream& rIStream,BinsObjectStatic<TConfigure>& rThis)
 {
   return rIStream;
 }
				    

  /// output stream function
  template<class TConfigure>  
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const BinsObjectStatic<TConfigure> & rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
 
  
};  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined 


