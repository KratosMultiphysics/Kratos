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


#if !defined(KRATOS_BINS_CONTAINER_H_INCLUDED)
#define  KRATOS_BINS_CONTAINER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <cmath>
#include <algorithm>

// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "tree.h"
#include "cell.h"
#include "bounding_box.h"


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
    class TGeometryType,
    //class TBoundingBoxFunction,
    class TDistanceIteratorType = typename std::vector<double>::iterator,
    class TPointerType          = typename TContainerType::value_type,
    class TIteratorType         = typename TContainerType::iterator,
    class TDistanceFunction     = Kratos::SearchUtils::SquaredDistanceFunction<TDimension,TPointType>  
    > 
  class Bins  
    {
    public:
      ///@name Type Definitions
      ///@{
      typedef Cell<TDimension, TPointerType> CellType;     
      typedef TPointType                     PointType;
      typedef TContainerType                 ContainerType;
      typedef TIteratorType                  IteratorType;
      typedef TDistanceIteratorType          DistanceIteratorType;
      typedef TPointerType                   PointerType;
      typedef TDistanceFunction              DistanceFunction;
      
      enum { Dimension = TDimension };
           
      typedef std::vector< CellType> CellContainerType;
      typedef typename CellContainerType::iterator CellContainerIterator;
      
      typedef BoundingBox< TPointType,TGeometryType > BoundingBoxType;
      typedef std::vector< BoundingBoxType > BoundingBoxContainerType;
      typedef typename BoundingBoxContainerType::iterator BoundingBoxIterator; 
      
      typedef TreeNode<TDimension,TPointType,TPointerType,TIteratorType,TDistanceIteratorType> TreeNodeType;
      typedef typename TreeNodeType::CoordinateType  CoordinateType;  // double
      typedef typename TreeNodeType::SizeType        SizeType;        // std::size_t
      typedef typename TreeNodeType::IndexType       IndexType;       // std::size_t
            
      
      //typedef TreeNodeType LeafType;    
       typedef typename TreeNodeType::IteratorIteratorType IteratorIteratorType;
       typedef typename TreeNodeType::SearchStructureType SearchStructureType;
       
      	
      /// Pointer definition of Bins
      KRATOS_CLASS_POINTER_DEFINITION(Bins);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Bins() {}     
      /// Constructor de bins a bounding box
      
      Bins (TIteratorType const& ObjectsBegin, TIteratorType const& ObjectsEnd) 
           : mObjectsBegin(ObjectsBegin), mObjectsEnd(ObjectsEnd)
      {
	
	 
	 CalculateBoundingBoxes();
	 CalculateBoundingBox();
         CalculateCellSize();
         AllocateCellsContainer();
         GenerateBins();
	 
      }
     
     

      /// Destructor.
      virtual ~Bins(){}
      

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
      }
      
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
      
      
      
            
      ///@}      
      ///@name Friends
      ///@{
      
            
      ///@}
      
    protected:
      ///@name Protected static Member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Protected member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Protected Operators
      ///@{ 
      

	
	
    /// Computa los boxes de cada uno de los elementos del model part		
    void CalculateBoundingBoxes()
    {
      KRATOS_TRY 	
      TPointType Low, High;
      
      for (TIteratorType i_object = mObjectsBegin ; i_object != mObjectsEnd ; i_object++)
	 {
	        i_object->GetGeometry().Bounding_Box(High, Low); 
		 BoundingBoxType rThisBoundingBox(High, Low, &i_object->GetGeometry());
	         mBoxes.push_back(rThisBoundingBox ); 
	 }
	 
      KRATOS_CATCH("")    
     }
     
    
    
              
    //************************************************************************          
    /// Computa el punto minimo y maximo de todos los bounding box          
    
    void CalculateBoundingBox()
    { 
     KRATOS_TRY 
     BoundingBoxIterator it; 
     for(size_t i = 0 ; i < TDimension ; i++){
            mMinPoint[i] = (mBoxes[0].LowPoint())[i];
            mMaxPoint[i] = (mBoxes[0].HighPoint())[i];  
       } 	    
      for(it=mBoxes.begin(); it!=mBoxes.end() ; it++)
       {
 	 for(size_t i = 0 ; i < TDimension ; i++)
 	    {
 	      mMaxPoint[i]  =   (mMaxPoint[i]  < ((*it).HighPoint())[i] ) ? ((*it).HighPoint())[i] : mMaxPoint[i];
 	      mMinPoint[i]  =   (mMinPoint[i]  > ((*it).LowPoint())[i] )  ? ((*it).LowPoint())[i]  : mMinPoint[i]; 
 	    }
       }
      KRATOS_CATCH("")      
    } 
    
    
    
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
      Tvector<CoordinateType,TDimension>  MinPoint;
      Tvector<CoordinateType,TDimension>  MaxPoint;
      //Tvector<IndexType,TDimension> Cell;
      
      SearchStructureType Box;
      for(IteratorType i_object = mObjectsBegin ; i_object != mObjectsEnd ; i_object++)
           {
	     i_object->GetGeometry().Bounding_Box(High, Low);
	     for(SizeType i = 0 ; i < TDimension ; i++)
	       { MinPoint[i] =  Low[i]; MaxPoint[i] = High[i];}
	     
             Box.Set( CalculateCell(Low), CalculateCell(High), mN );
             FillObject(Box,i_object);
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
    void FillObject( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box, IteratorType i_object)
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

    IndexType CalculatePosition( CoordinateType const& ThisCoord, SizeType ThisDimension )
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

void AllocateCellsContainer() {
            SizeType Size = 1;
            for(SizeType i = 0 ; i < TDimension ; i++)
                Size *= mN[i];
            // Resize Global Container
            mCells.resize(Size);
	    
	    AllocateCell() ;
        }
        
//************************************************************************
//************************************************************************ 

void AllocateCell() 
       {
                 
            const std::size_t Size = 10;/// mN[0];  WARNING = Calcular el maximo  
            for(CellContainerIterator icell = GetCellContainer().begin() ; icell != GetCellContainer().end() ; icell++)
	    {     
               icell->AllocateCell(Size);
	    }
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
      
       CellContainerType mCells; // The bin    
       BoundingBoxContainerType mBoxes;
      

      
 
	
	
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
        
        
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
      Bins& operator=(Bins const& rOther){}

      /// Copy constructor.
      Bins(Bins const& rOther){}

        
      ///@}    
        
    }; // Class Bins 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  
  template<
  std::size_t TDimension,
  class TPointType,
  class TContainerType,
  class TGeometryType,
  class TBoundingBoxFunction,
  class TDistanceIteratorType, 
  class TPointerType,
  class TIteratorType,
  class TDistanceFunction
  >
  inline std::istream& operator >> (std::istream& rIStream, 
				    Bins<TDimension, TPointerType, TPointType, TGeometryType, TBoundingBoxFunction, TIteratorType, TDistanceFunction>& rThis);
				    

  /// output stream function
  template<
  std::size_t TDimension,
  class TPointType,
  class TContainerType,
  class TGeometryType,
  class TBoundingBoxFunction,
  class TDistanceIteratorType, 
  class TPointerType,
  class TIteratorType,
  class TDistanceFunction
  >
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Bins<TDimension, TPointerType, TPointType, TGeometryType,TBoundingBoxFunction, TIteratorType, TDistanceFunction>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined 


