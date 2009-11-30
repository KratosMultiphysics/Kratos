//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: clabra $
//   Date:                $Date: 2007-03-27 17:02:19 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_BUCKET_H_INCLUDED )
#define  KRATOS_BUCKET_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <vector>



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
  
  /// Short class definition.
  /** Detail class definition.
  */
	template<
       std::size_t TDimension,
       class TPointType,
       class TContainerType,
       class TPointerType = typename TContainerType::value_type,
       class TIteratorType = typename TContainerType::iterator,
       class TDistanceIteratorType = typename std::vector<double>::iterator,
       class TDistanceFunction = SpacialSearchSquaredDistanceFunction<TDimension,TPointType>
       >
  class Bucket  : public TreeNode<TDimension,TPointType,TPointerType,TIteratorType,TDistanceIteratorType>
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Bucket
      KRATOS_CLASS_POINTER_DEFINITION(Bucket);

        typedef TreeNode<TDimension, TPointType, TPointerType, TIteratorType, TDistanceIteratorType> BaseType;

		typedef TPointType PointType;

        typedef TContainerType ContainerType;

		typedef TIteratorType IteratorType;

		typedef TDistanceIteratorType DistanceIteratorType;

        typedef TPointerType PointerType;
        
        typedef TDistanceFunction  DistanceFunction;

  		typedef typename BaseType::SizeType       SizeType;

		typedef typename BaseType::IndexType      IndexType;

		typedef typename BaseType::CoordinateType CoordinateType;

		enum { Dimension = TDimension };

        typedef typename BaseType::SearchStructureType SearchStructureType;


      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Bucket()
			: mPointsBegin(this->NullIterator()), mPointsEnd(this->NullIterator())
		{}

      Bucket(IteratorType PointsBegin,IteratorType PointsEnd)
			: mPointsBegin(PointsBegin), mPointsEnd(PointsEnd)
		{}

      /// Destructor.
		virtual ~Bucket(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{

		/* 		SizeType const Dimension() const { return TDimension; } */

		IteratorType Begin()
		{ 
          return mPointsBegin; 
        }

        IteratorType End()
        { 
          return mPointsEnd; 
        }

        void SearchNearestPoint(PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance )
        {
          if(mPointsBegin == mPointsEnd) 
            return;

          // Loop over all points in bucket
          for(IteratorType iPoint = mPointsBegin ; iPoint != mPointsEnd ; iPoint++)
          {	
            CoordinateType distance = TDistanceFunction()(ThisPoint,**iPoint);
            if(distance < rResultDistance)
            {
              rResult = *iPoint;
              rResultDistance = distance;
            }
          }

        }
        
        void SearchNearestPoint(PointType const& ThisPoint, PointerType& Result, CoordinateType& ResultDistance, SearchStructureType& Auxiliar )
        {
		  SearchNearestPoint(ThisPoint,Result,ResultDistance);
		}


        void SearchInRadius(PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results, 
              DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults)
        {
          
          if(mPointsBegin == mPointsEnd)
            return;

          CoordinateType distance2; // For holding distance to the point for 
          //Auxiliar.BucketCounter++;  // DEBUG

          // Loop over all points in bucket
          for(IteratorType i = mPointsBegin ; (i != mPointsEnd) && (NumberOfResults < MaxNumberOfResults) ; i++)
          {	
            distance2 = TDistanceFunction()(ThisPoint,**i);
            if(distance2 < Radius2)
            {
              *(Results) = *i;
              *(ResultsDistances) = distance2;
              Results++;
              ResultsDistances++;
              NumberOfResults++;
            }
          }

        }
        
		void SearchInRadius(PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results, 
			DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType& Auxiliar)
        {
           SearchInRadius(ThisPoint,Radius,Radius2,Results,ResultsDistances,NumberOfResults,MaxNumberOfResults);
        }


        void SearchInRadius(PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results, 
              SizeType& NumberOfResults, SizeType const& MaxNumberOfResults)
        {
          
          if(mPointsBegin == mPointsEnd)
            return;

          CoordinateType distance2; // For holding distance to the point for 
          //Auxiliar.BucketCounter++;  // DEBUG

          // Loop over all points in bucket
          for(IteratorType i = mPointsBegin ; (i != mPointsEnd) && (NumberOfResults < MaxNumberOfResults) ; i++)
          {	
            distance2 = TDistanceFunction()(ThisPoint,**i);
            if(distance2 < Radius2)
            {
              *(Results) = *i;
              Results++;
              NumberOfResults++;
            }
          }

        }
        
        void SearchInRadius(PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results, 
              SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType& Auxiliar)
        {
           SearchInRadius(ThisPoint,Radius,Radius2,Results,NumberOfResults,MaxNumberOfResults);
        }
     
        void SearchInBox(PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& Results, SizeType& NumberOfResults,
            SizeType const& MaxNumberOfResults )
        {
          for(IteratorType i = mPointsBegin ; (i != mPointsEnd) && (NumberOfResults < MaxNumberOfResults) ; i++)
          {
            if( SpatialSearchPointInBox<Dimension,PointType>(SearchMinPoint,SearchMaxPoint,**i) )
            {
              *(Results) = *i;
              Results++;
              NumberOfResults++;
            }
          }
        }


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
		  return "Bucket";
	  }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	  {
		  rOStream << "Bucket";
	  }

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream, std::string const& Perfix = std::string()) const
      {
         rOStream << Perfix << "Leaf[" << PointerDistance(mPointsBegin, mPointsEnd) << "] : ";
         for(IteratorType i = mPointsBegin ; i != mPointsEnd ; i++)
            rOStream << **i << "    ";
         rOStream << std::endl;
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
       
      /*       bool PointInBox(PointType const& BoxMinPoint, PointType const& BoxMaxPoint, PointType const& ThisPoint) */
      /*       { */
      /*         for(SizeType i = 0 ; i < Dimension ; i++) */
      /*           if( ThisPoint[i] < BoxMinPoint[i] || ThisPoint[i] > BoxMaxPoint[i] ) */
      /*             return false; */
      /*         return true; */
      /*       } */
        
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
		IteratorType mPointsBegin;
		IteratorType mPointsEnd;
        
        
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
      Bucket& operator=(Bucket const& rOther);

      /// Copy constructor.
      Bucket(Bucket const& rOther);

        
      ///@}    
        
    }; // Class Bucket 

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
    class TPointerType,
    class TIteratorType,
    class TDistanceIteratorType,
    class TDistanceFunction>
 inline std::istream& operator >> (std::istream& rIStream, 
       Bucket<TDimension, TPointType, TContainerType, TPointerType, TIteratorType, TDistanceIteratorType, TDistanceFunction>& rThis);

  /// output stream function
  template<
    std::size_t TDimension,
    class TPointType,
    class TContainerType,
    class TPointerType,
    class TIteratorType,
    class TDistanceIteratorType,
    class TDistanceFunction>
  inline std::ostream& operator << (std::ostream& rOStream, 
       const Bucket<TDimension, TPointType, TContainerType, TPointerType, TIteratorType, TDistanceIteratorType, TDistanceFunction>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined 


