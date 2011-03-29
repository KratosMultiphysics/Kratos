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
#include <time.h>

// Project includes
#include "tree.h"
#include "cell.h"

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
      typedef typename TConfigure::ResultContainerType     ResultContainerType; 
      typedef typename TConfigure::ResultIteratorType      ResultIteratorType;
        
      typedef Cell<Configure> CellType;  
      typedef std::vector<CellType> CellContainerType;
      typedef typename CellContainerType::iterator CellContainerIterator;
      
      typedef TreeNode<Dimension, PointType, PointerType, IteratorType,  typename TConfigure::DistanceIteratorType> TreeNodeType;
      typedef typename TreeNodeType::CoordinateType  CoordinateType;  // double
      typedef typename TreeNodeType::SizeType        SizeType;        // std::size_t
      typedef typename TreeNodeType::IndexType       IndexType;       // std::size_t
           

      typedef Tvector<IndexType,Dimension>      IndexArray;
      typedef Tvector<SizeType,Dimension>       SizeArray;
      typedef Tvector<CoordinateType,Dimension> CoordinateArray;
      
       ///Contact Pair
      typedef typename TConfigure::ContainerContactType  ContainerContactType;
      typedef typename TConfigure::IteratorContactType IteratorContactType;
      
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
        mObjectsSize = SearchUtils::PointerDistance(mObjectsBegin,mObjectsEnd);
        CalculateBoundingBox();           // Calculate mMinPoint, mMaxPoint
        CalculateCellSize(mObjectsSize);  // Calculate number of Cells
        AllocateContainer();              // Allocate cell list
        GenerateBins();                   // Fill Cells with objects

      }
     
      BinsObjectDynamic (const PointType& MinPoint, const PointType& MaxPoint, CoordinateType CellSize)
        : mObjectsSize(0), mObjectsBegin(0), mObjectsEnd(0)
      {

        for(SizeType i = 0 ; i < Dimension ; i++)
        {
          mMinPoint[i] = MinPoint[i];
          mMaxPoint[i] = MaxPoint[i];
        }
        CalculateCellSize(CellSize);
        AllocateContainer();
      }
     
      BinsObjectDynamic (const PointType& MinPoint, const PointType& MaxPoint, SizeType NumPoints) 
        : mObjectsSize(0), mObjectsBegin(0), mObjectsEnd(0)
      {

        for(SizeType i = 0 ; i < Dimension ; i++)
        {
          mMinPoint[i] = MinPoint[i];
          mMaxPoint[i] = MaxPoint[i];
        }
        CalculateCellSize(NumPoints);
        AllocateContainer();
      }

      
      /// Destructor.
      virtual ~BinsObjectDynamic(){}
      
     
     
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
     
 SizeType SearchObjectsInner(PointerType& ThisObject, ResultIteratorType& Result,  const SizeType& MaxNumberOfResults )
    {
      PointType Low, High;      
      SearchStructureType Box;
      SizeType NumberOfResults = 0;
      TConfigure::CalculateBoundingBox(ThisObject, Low, High);
      Box.Set( CalculateCell(Low), CalculateCell(High), mN );
      SearchObjectLocalInner(ThisObject, Result, NumberOfResults, MaxNumberOfResults, Box );
      return NumberOfResults;
     } 
       

//************************************************************************   
//************************************************************************
    
 SizeType SearchObjectsInner(PointerType& ThisObject, ResultContainerType& Result)
    {
      PointType Low, High;      
      SearchStructureType Box;
      TConfigure::CalculateBoundingBox(ThisObject, Low, High);
      Box.Set( CalculateCell(Low), CalculateCell(High), mN );
      SearchObjectLocalInner(ThisObject, Result, Box );
      return Result.size();
     }          
            
//************************************************************************   
//************************************************************************          
          
    void SearchContact(ContainerContactType& Result)
    {
      for (CellContainerIterator icell = mCells.begin() ; icell!= mCells.end(); icell++)
        icell->SearchContact(Result);
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

    void AddObject(const PointerType& ThisObject)
    {
      PointType Low, High;      
      SearchStructureType Box;
      TConfigure::CalculateBoundingBox(ThisObject, Low, High);
      Box.Set( CalculateCell(Low), CalculateCell(High), mN );
      FillObject(Box,ThisObject);
      mObjectsSize++;
    }  

    void RemoveObject(const PointerType& ThisObject)
    {
      PointType Low, High;      
      SearchStructureType Box;
      TConfigure::CalculateBoundingBox(ThisObject, Low, High);
      Box.Set( CalculateCell(Low), CalculateCell(High), mN );
      RemoveObjectLocal(Box,ThisObject);
      mObjectsSize--;
    }  
         
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
      virtual std::string Info() const{
	return "BinsObjectDynamic" ; }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
      {
        rOStream << Info();
      }

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream, std::string const& Perfix = std::string()) const
      {	
        rOStream << " BinsSize: ";
        for(SizeType i = 0 ; i < Dimension ; i++)
          rOStream << "[" << mN[i] << "]";
        rOStream << std::endl;
        rOStream << "  CellSize: ";
        for(SizeType i = 0 ; i < Dimension ; i++)
          rOStream << "[" << mCellSize[i] << "]";
        rOStream << std::endl;
        SizeType nn = 0;
        for(SizeType i = 0 ; i < mCells.size(); i++)
          nn += mCells[i].Size();
        rOStream << "NumPointers: " << nn << std::endl;
      }
      
     /// Print Size of Container
    void PrintSize( std::ostream& rout ){
      rout << " BinsSize: ";
      for(SizeType i = 0 ; i < Dimension ; i++)
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
      
      
    CellContainerType& GetCellContainer()
    {
      return mCells; 
    }

    SizeArray& GetDivisions()
    {
      return mN;
    }

    CoordinateArray& GetCellSize()  
    {
      return mCellSize;
    }
      
      

//************************************************************************   
//************************************************************************

             
      private:
      
            
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

      #ifdef _OPENMP
      SizeType number_of_threads = omp_get_max_threads();
      #else
      SizeType number_of_threads = 1;
      #endif

      std::vector<SizeType> node_partition;
      CreatePartition(number_of_threads, mObjectsSize, node_partition);

      std::vector<PointType> Max(number_of_threads);
      std::vector<PointType> Min(number_of_threads);

      for(SizeType k=0; k<number_of_threads; k++ )
      {
        Max[k] = mMaxPoint;   
        Min[k] = mMinPoint;
      }

      #pragma omp parallel for  private(High, Low) 
      for(SizeType k=0; k<number_of_threads; k++)
      {
        IteratorType i_begin = mObjectsBegin + node_partition[k];
        IteratorType i_end   = mObjectsBegin + node_partition[k+1];

        for (IteratorType i_object  = i_begin ; i_object != i_end ; i_object++ )
        { 
          TConfigure::CalculateBoundingBox(*i_object, Low, High);
          for(SizeType i = 0 ; i < Dimension ; i++)
          {
            Max[k][i] = (Max[k][i]  < High[i]) ? High[i] : Max[k][i];
            Min[k][i] = (Min[k][i]  > Low[i])  ? Low[i]  : Min[k][i];
          }
        }
      }	

      for(SizeType k=0; k<number_of_threads; k++)
      {
        for(SizeType i = 0 ; i < Dimension ; i++)
        {
          mMaxPoint[i]  = (mMaxPoint[i]  < Max[k][i]) ? Max[k][i] : mMaxPoint[i];
          mMinPoint[i]  = (mMinPoint[i]  > Min[k][i]) ? Min[k][i] : mMinPoint[i]; 
        }
      }
      //for(std::size_t i = 0 ; i < Dimension ; i++){
        //mMaxPoint[i]  += 10.0 * DBL_EPSILON;
        //mMinPoint[i]  -= 10.0 * DBL_EPSILON;
      //}
         
     }
     
    
//************************************************************************   
//************************************************************************   
    
    void CalculateCellSize()
    {
      
      CoordinateType delta[Dimension];
      CoordinateType alpha[Dimension];
      CoordinateType mult_delta = 1.00;
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

      
      mN[index] = static_cast<SizeType>( pow(static_cast<CoordinateType>(mObjectsSize/mult_delta), 1.00/Dimension) +1 );
      
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
   
    void CalculateCellSize(const CoordinateType& CellSize)
    {
      for(SizeType i = 0 ; i < Dimension ; i++){
        mCellSize[i] = CellSize;
        mInvCellSize[i] = 1.00 / mCellSize[i];
        mN[i] = static_cast<SizeType>( (mMaxPoint[i]-mMinPoint[i]) / mCellSize[i]) + 1;
      }
    }
    
    void CalculateCellSize( const SizeType& NumPoints )
    {

      CoordinateType delta[Dimension];
      CoordinateType alpha[Dimension];
      CoordinateType mult_delta = 1.00;
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

      mN[index] = static_cast<SizeType>( pow(static_cast<CoordinateType>(NumPoints/mult_delta), 1.00/Dimension) +1 );

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
      /// Fill container with objects
      for(IteratorType i_object = mObjectsBegin ; i_object != mObjectsEnd ; i_object++)
      {
        TConfigure::CalculateBoundingBox(*i_object, Low, High);
        Box.Set( CalculateCell(Low), CalculateCell(High), mN );
        FillObject(Box, *i_object);
      }
   }    
      
            
//************************************************************************   
//************************************************************************
                  
// **** THREAD SAFE
// Dimension = 1
   void SearchInBoxLocal(PointerType& ThisObject, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
       SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
   {
     PointType  MinCell, MaxCell;
     PointType  MinBox, MaxBox;
     MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  // 
     MaxCell[0] = MinCell[0] + mCellSize[0];
     for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0])
       if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
         mCells[I].SearchObjects(ThisObject, Result, NumberOfResults, MaxNumberOfResults);
   }

   // Dimension = 2
   void SearchInBoxLocal(PointerType& ThisObject, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
       SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
   {
     PointType  MinCell, MaxCell;
     PointType  MinBox, MaxBox;

     for(SizeType i = 0; i < 2; i++){
       MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];   
       MaxBox[i] = MinBox[i] + mCellSize[i];
     }

     MinCell[1] = MinBox[1];
     MaxCell[1] = MaxBox[1];
     for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] )  {
       MinCell[0] = MinBox[0];
       MaxCell[0] = MaxBox[0];
       for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] ) {
         if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
           mCells[I].SearchObjects(ThisObject, Result, NumberOfResults, MaxNumberOfResults);
       }
     }
   } 

   // Dimension = 3
   void SearchInBoxLocal(PointerType& ThisObject, ResultIteratorType& Result, 
       SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
       SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
   {

     PointType  MinCell, MaxCell;
     PointType  MinBox, MaxBox;

     for(SizeType i = 0; i < 3; i++){
       MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  // 
       MaxBox[i] = MinBox[i] + mCellSize[i];
     }

     MinCell[2] = MinBox[2];
     MaxCell[2] = MaxBox[2];  
     for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block, MinCell[2] += mCellSize[2], MaxCell[2] += mCellSize[2] ){
       MinCell[1] = MinBox[1];
       MaxCell[1] = MaxBox[1];
       for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] ) {
         MinCell[0] = MinBox[0];
         MaxCell[0] = MaxBox[0]; 
         for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] ){
           if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell)){
             mCells[I].SearchObjects(ThisObject, Result, NumberOfResults, MaxNumberOfResults);
           }
         }
       }
     }
   }

 //************************************************************************   
//************************************************************************
                  
// **** THREAD SAFE
// Dimension = 1
   void SearchObjectLocalInner(PointerType& ThisObject, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
       SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
   {
     PointType  MinCell, MaxCell;
     PointType  MinBox, MaxBox;
     MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  // 
     MaxCell[0] = MinCell[0] + mCellSize[0];
     for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0])
       if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
         mCells[I].SearchObjectsInner(ThisObject, Result, NumberOfResults, MaxNumberOfResults);
   }

   // Dimension = 2
   void SearchObjectLocalInner(PointerType& ThisObject, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
       SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
   {
     PointType  MinCell, MaxCell;
     PointType  MinBox, MaxBox;

     for(SizeType i = 0; i < 2; i++){
       MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];   
       MaxBox[i] = MinBox[i] + mCellSize[i];
     }

     MinCell[1] = MinBox[1];
     MaxCell[1] = MaxBox[1];
     for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] )  {
       MinCell[0] = MinBox[0];
       MaxCell[0] = MaxBox[0];
       for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] ) {
         if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
           mCells[I].SearchObjectsInner(ThisObject, Result, NumberOfResults, MaxNumberOfResults);
       }
     }
   } 

   // Dimension = 3
   void SearchObjectLocalInner(PointerType& ThisObject, ResultIteratorType& Result, 
       SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
       SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
   {

     PointType  MinCell, MaxCell;
     PointType  MinBox, MaxBox;

     for(SizeType i = 0; i < 3; i++){
       MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  // 
       MaxBox[i] = MinBox[i] + mCellSize[i];
     }

     MinCell[2] = MinBox[2];
     MaxCell[2] = MaxBox[2];  
     for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block, MinCell[2] += mCellSize[2], MaxCell[2] += mCellSize[2] ){
       MinCell[1] = MinBox[1];
       MaxCell[1] = MaxBox[1];
       for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] ) {
         MinCell[0] = MinBox[0];
         MaxCell[0] = MaxBox[0]; 
         for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] ){
           if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell)){
             mCells[I].SearchObjectsInner(ThisObject, Result, NumberOfResults, MaxNumberOfResults);
           }
         }
       }
     }
   }

 //************************************************************************   
//************************************************************************          
            
                  
    // **** THREAD SAFE
    
    // Dimension = 1
    void SearchInBoxLocal(PointerType& ThisObject, ResultContainerType& Result, 
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
     {
       
       PointType  MinCell, MaxCell;
       PointType  MinBox, MaxBox;
       MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  // 
       MaxCell[0] = MinCell[0] + mCellSize[0];
       for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block ) {
	   if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
              mCells[I].SearchObjects(ThisObject, Result);
       }
     }

    // Dimension = 2
    void SearchInBoxLocal(PointerType& ThisObject, ResultContainerType& Result, 
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
      PointType  MinCell, MaxCell;
      PointType  MinBox, MaxBox;

      for(SizeType i = 0; i < 2; i++){
        MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];   
        MaxBox[i] = MinBox[i] + mCellSize[i];
      }

      MinCell[1] = MinBox[1];
      MaxCell[1] = MaxBox[1];

      for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block ){
        MinCell[0] = MinBox[0];
        MaxCell[0] = MaxBox[0];
        for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block ) 
        {
          if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
            mCells[I].SearchObjects(ThisObject, Result);
        }
      }
    } 

    // Dimension = 3
    void SearchInBoxLocal(PointerType& ThisObject, ResultContainerType& Result, 
        SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {
      PointType  MinCell, MaxCell;
      PointType  MinBox, MaxBox;

      for(SizeType i = 0; i < 3; i++){
        MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  // 
        MaxBox[i] = MinBox[i] + mCellSize[i];
      }

      MinCell[2] = MinBox[2];
      MaxCell[2] = MaxBox[2];

      for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block ){
        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];
        for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block ){
          MinCell[0] = MinBox[0];
          MaxCell[0] = MaxBox[0]; 
          for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block ){
            if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell)){
              mCells[I].SearchObjects(ThisObject, Result);
            }
          }
        }
      }
    }



//************************************************************************   
//************************************************************************          
            
                  
    // **** THREAD SAFE
    
    // Dimension = 1
    void SearchObjectLocalInner(PointerType& ThisObject, ResultContainerType& Result, 
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
     {
       PointType  MinCell, MaxCell;
       PointType  MinBox, MaxBox;
       MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  // 
       MaxCell[0] = MinCell[0] + mCellSize[0];
       for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block )
         if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
           mCells[I].SearchObjectsInner(ThisObject, Result);

     }

    // Dimension = 2
    void SearchObjectLocalInner(PointerType& ThisObject, ResultContainerType& Result, 
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
       PointType  MinCell, MaxCell;
       PointType  MinBox, MaxBox;

       for(SizeType i = 0; i < 2; i++){
        MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];   
        MaxBox[i] = MinBox[i] + mCellSize[i];
       }

      MinCell[1] = MinBox[1];
      MaxCell[1] = MaxBox[1];
      
      for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block ){
        MinCell[0] = MinBox[0];
        MaxCell[0] = MaxBox[0];
        for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block ) 
	{
	   if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
	       mCells[I].SearchObjectsInner(ThisObject, Result);
	}
     }
    } 

    // Dimension = 3
    void SearchObjectLocalInner(PointerType& ThisObject, ResultContainerType& Result, 
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {
      PointType  MinCell, MaxCell;
      PointType  MinBox, MaxBox;

      for(SizeType i = 0; i < 3; i++){
        MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  // 
        MaxBox[i] = MinBox[i] + mCellSize[i];
      }

      MinCell[2] = MinBox[2];
      MaxCell[2] = MaxBox[2];

      for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block ){
        MinCell[2] = MinBox[2];
        MaxCell[2] = MaxBox[2];
        for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block ){
          MinCell[0] = MinBox[0];
          MaxCell[0] = MaxBox[0]; 
          for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block ){
            if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell)){
              mCells[I].SearchObjectsInner(ThisObject, Result);
            }
          }
        }
      }
    }



//************************************************************************   
//************************************************************************
        
     IndexArray  CalculateCell( const PointType& ThisPoint )
     {
        IndexArray IndexCell;
        for(SizeType i = 0 ; i < Dimension ; i++)
            IndexCell[i] = CalculatePosition(ThisPoint[i],i);
        return IndexCell; 
     }
     
     
//************************************************************************   
//************************************************************************

    // Dimension = 1 
    void FillObject( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box, const PointerType& i_object)
    {
      PointType  MinCell, MaxCell;

      MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  // 
      MaxCell[0] = MinCell[0] + mCellSize[0];
      for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] )
      {
        if(TConfigure::IntersectionBox(i_object, MinCell, MaxCell))
          mCells[I].Add(i_object);
      }
    }

   
//************************************************************************   
//************************************************************************

    // Dimension = 2
    void FillObject( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box, const PointerType& i_object)
    {
      PointType  MinCell, MaxCell;
      PointType  MinBox, MaxBox;

      for(SizeType i = 0; i < 2; i++){
        MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];   
        MaxBox[i] = MinBox[i] + mCellSize[i];
      }

      MinCell[1] = MinBox[1];
      MaxCell[1] = MaxBox[1];
      for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1]+=mCellSize[1], MaxCell[1]+=mCellSize[1] )
      {
        MinCell[0] = MinBox[0];
        MaxCell[0] = MaxBox[0];
        for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] )
        {
          if(TConfigure::IntersectionBox(i_object,MinCell,MaxCell))
            mCells[I].Add(i_object);
        }
      }
    }
    

//************************************************************************   
//************************************************************************

    // Dimension = 3
    void FillObject( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box, const PointerType& i_object)
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
          {
            if(TConfigure::IntersectionBox(i_object,MinCell,MaxCell))
              mCells[I].Add(i_object);
          }
        }
      }
    }

//************************************************************************
//************************************************************************

    // Dimension = 1 
    void RemoveObjectLocal( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box, const PointerType& i_object)
    {
      PointType  MinCell, MaxCell;

      MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  // 
      MaxCell[0] = MinCell[0] + mCellSize[0];
      for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] )
      {
        if(TConfigure::IntersectionBox(i_object, MinCell, MaxCell))
          mCells[I].Remove(i_object);
      }
    }

   
//************************************************************************   
//************************************************************************

    // Dimension = 2
    void RemoveObjectLocal( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box, const PointerType& i_object)
    {
      PointType  MinCell, MaxCell;
      PointType  MinBox, MaxBox;

      for(SizeType i = 0; i < 2; i++){
        MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];   
        MaxBox[i] = MinBox[i] + mCellSize[i];
      }

      MinCell[1] = MinBox[1];
      MaxCell[1] = MaxBox[1];
      for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1]+=mCellSize[1], MaxCell[1]+=mCellSize[1] )
      {
        MinCell[0] = MinBox[0];
        MaxCell[0] = MaxBox[0];
        for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] )
        {
          if(TConfigure::IntersectionBox(i_object,MinCell,MaxCell))
            mCells[I].Remove(i_object);
        }
      }
    }
    

//************************************************************************   
//************************************************************************

    // Dimension = 3
    void RemoveObjectLocal( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box, const PointerType& i_object)
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
          {
            if(TConfigure::IntersectionBox(i_object,MinCell,MaxCell))
              mCells[I].Remove(i_object);
          }
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

    void AllocateContainer() 
    {
      SizeType Size = mN[0];
      for(SizeType i = 1 ; i < Dimension ; i++)
        Size *= mN[i];
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
      SizeType     mObjectsSize;

      CoordinateArray  mCellSize;
      CoordinateArray  mInvCellSize;
      SizeArray        mN;
      
      CellContainerType mCells;  ///The bin    
	
	
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
      
    inline void CreatePartition(SizeType number_of_threads, const SizeType number_of_rows, std::vector<SizeType>& partitions)
    {
      partitions.resize(number_of_threads+1);
      SizeType partition_size = number_of_rows / number_of_threads;
      partitions[0] = 0;
      partitions[number_of_threads] = number_of_rows;
      for(SizeType i = 1; i<number_of_threads; i++)
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
      ///@} 
      
      
      public:
      /// Assignment operator.
      BinsObjectDynamic<TConfigure> & operator=(const BinsObjectDynamic<TConfigure> & rOther)
      {
	mMinPoint            = rOther.mMinPoint;     
	mMaxPoint            = rOther.mMaxPoint;  
	mObjectsBegin        = rOther.mObjectsBegin; 
 	mObjectsEnd          = rOther.mObjectsEnd; 
	mObjectsSize         = rOther.mObjectsSize;  
	mCellSize            = rOther.mCellSize;
	mInvCellSize         = rOther.mInvCellSize; 
	mN                   = rOther.mN; 
	mCells               = rOther.mCells; 
	return *this;
      }

      /// Copy constructor.
      BinsObjectDynamic(const BinsObjectDynamic& rOther)
      {
        *this =  rOther;
      }

        
        
      
     
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


