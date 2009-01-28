//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2008-03-12 14:20:36 $
//   Revision:            $Revision: 1.2 $
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
//#include "includes/define.h"
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
	template<std::size_t TDimension, class TPointType, class TIteratorType, class TDistanceIteratorType, class TDistanceFunction>
  class Bucket  : public TreeNode<TPointType, TIteratorType, TDistanceIteratorType>
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Bucket
      //KRATOS_CLASS_POINTER_DEFINITION(Bucket);

  		typedef std::size_t SizeType;

		typedef std::size_t IndexType;

		typedef TPointType PointType;

		typedef double CoordinateType;

		typedef TIteratorType IteratorType;

		typedef TDistanceIteratorType DistanceIteratorType;

		//typedef std::vector<IndexType> IndicesType;


      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Bucket()
			: mPointsBegin(this->NullIterator()), mPointsEnd(this->NullIterator())
		{}

      Bucket(IteratorType PointsBegin, 
			IteratorType PointsEnd)
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
		IteratorType Begin()
		{ 
			return mPointsBegin; 
		}

		IteratorType End()
		{ 
			return mPointsEnd; 
		}
 
		IteratorType SearchNearestPoint(PointType const& ThisPoint, 
			CoordinateType& rResultDistance)
		{
			CoordinateType distance; // For holding distance to the point for 
			//CoordinateType max_distance = -1;

			if(mPointsBegin == mPointsEnd)
				return this->NullIterator();

      /* 			rResultDistance = 0; */

			// initializing rResultDistance with distance to first point
            rResultDistance = TDistanceFunction()(ThisPoint,**mPointsBegin);
			IteratorType result = mPointsBegin;

			// Loop over all points in bucket
			for(IteratorType i = mPointsBegin + 1 ; i != mPointsEnd ; i++)
			{	
               distance = TDistanceFunction()(ThisPoint,**i);
               if(distance < rResultDistance)
               {
                  result = i;
                  rResultDistance = distance;
               }
			}

			return result;
		}

        SizeType SearchInRadius(PointType const& ThisPoint, double Radius, IteratorType Results, 
              DistanceIteratorType ResultsDistances, SizeType MaximumNumberOfResults)
        {
           double Radius2 = Radius * Radius;
           return SearchInRadius(ThisPoint,Radius,Radius2,Results,ResultsDistances,MaximumNumberOfResults);

        }

        SizeType SearchInRadius(PointType const& ThisPoint, double Radius, double Radius2, IteratorType Results, 
              DistanceIteratorType ResultsDistances, SizeType MaximumNumberOfResults)
        {
           CoordinateType distance; // For holding distance to the point for 

           if(mPointsBegin == mPointsEnd)
              return 0;

			SizeType n = 0; // Number of founded points

			// Loop over all points in bucket
			for(IteratorType i = mPointsBegin ; (i != mPointsEnd) && (n < MaximumNumberOfResults) ; i++)
			{	
               distance = TDistanceFunction()(ThisPoint,**i);
               if(distance < Radius2)
               {
                  *(Results) = *i;
                  *(ResultsDistances) = distance;
                  Results++;
                  ResultsDistances++;
                  n++;
               }
            }

			return n;
        }

        SizeType SearchInRadius(PointType const& ThisPoint, double Radius, IteratorType Results, 
              SizeType MaximumNumberOfResults)
        {
           double Radius2 = Radius * Radius;
           return SearchInRadius(ThisPoint,Radius,Radius2,Results,MaximumNumberOfResults);

        }

        SizeType SearchInRadius(PointType const& ThisPoint, double Radius, double Radius2, IteratorType Results, 
              SizeType MaximumNumberOfResults)
        {
           CoordinateType distance; // For holding distance to the point for 
           
           if(mPointsBegin == mPointsEnd)
              return 0;
           
           SizeType n = 0; // Number of founded points
           
           // Loop over all points in bucket
           for(IteratorType i = mPointsBegin ; (i != mPointsEnd) && (n < MaximumNumberOfResults) ; i++)
           {	
              distance = TDistanceFunction()(ThisPoint,**i);
              
              if(distance < Radius2)
              {
                 *(Results) = *i;
                 Results++;
                 n++;
              }
           }
           
           return n;
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
         rOStream << Perfix << "Leaf[" << std::distance(mPointsBegin, mPointsEnd) << "] : ";
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
  template<std::size_t TDimension, class TPointType, class TIteratorType, class TDistanceIteratorType, class TDistanceFunction>
 inline std::istream& operator >> (std::istream& rIStream, 
				    Bucket<TDimension, TPointType, TIteratorType, TDistanceIteratorType, TDistanceFunction>& rThis);

  /// output stream function
  template<std::size_t TDimension, class TPointType, class TIteratorType, class TDistanceIteratorType, class TDistanceFunction>
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Bucket<TDimension, TPointType, TIteratorType, TDistanceIteratorType, TDistanceFunction>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined 


