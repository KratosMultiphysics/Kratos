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
//   Date:                $Date: 2008-10-23 11:30:49 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_BOUNDING_BOX_H_INCLUDED )
#define  KRATOS_BOUNDING_BOX_H_INCLUDED 



// System includes
#include <string>
#include <iostream> 
#include <sstream>
#include <cstddef>



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
  template<class TPointType,  class TPointerType>   
  class BoundingBox
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of BoundingBox
      KRATOS_CLASS_POINTER_DEFINITION(BoundingBox);
      
      typedef TPointerType  PointerType;
      typedef BoundingBox<TPointType, PointerType > BoundingBoxType;   
        
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      BoundingBox() : mHighPoint(), mLowPoint()
      {
	std::cout<< "Calling empty constructor" <<std::endl;
      }
      
      BoundingBox(const TPointType& Point) :  mLowPoint(Point), mHighPoint(Point)
       {
       }

      BoundingBox(const TPointType& LowPoint, const TPointType& HighPoint ) :  mLowPoint(LowPoint), mHighPoint(HighPoint)
       {
       } 

       BoundingBox(const PointerType Object, const TPointType& LowPoint, const TPointType& HighPoint) :  mLowPoint(LowPoint), mHighPoint(HighPoint), mObject(Object)
       {
       }    

      
      /// Destructor.
      virtual ~BoundingBox(){};
            

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      
      
      ///@}
      ///@name Access
      ///@{ 
      
      void Set(const TPointerType Object, const TPointType& LowPoint, const TPointType& HighPoint)   
      { 
	mLowPoint  = LowPoint;
	mHighPoint = HighPoint; 
	mObject     = Object;
	
      }
      
      TPointType const& HighPoint() const
        {
           return mHighPoint; 
        }
      

      TPointType& HighPoint() 
        {
           return mHighPoint; 
        }
      
      TPointType const& LowPoint() const
        {
           return mLowPoint; 
        }
      

      TPointType& LowPoint() 
        {
           return mLowPoint; 
        }


      /// Assignment operator.
      BoundingBox& operator=(BoundingBox const& rOther)
       {
	  mHighPoint = rOther.mHighPoint;
	  mLowPoint  = rOther.mLowPoint;
	  mObject     = rOther.mObject;
	  return *this;
        }
       

      /// Copy constructor.
      BoundingBox(BoundingBox const& rOther) : 
          mHighPoint(rOther.mHighPoint), 
          mLowPoint(rOther.mLowPoint),
          mObject(rOther.mObject)     
       {             
       }



      ///@}
      ///@name Inquiry
      ///@{
      
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
       {
         return "BoundingBox";
       }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
      {
	  rOStream << Info();
      }

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
      {
	  rOStream << mHighPoint << " , " << mLowPoint;
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
      TPointType mHighPoint;        
      TPointType mLowPoint;           
      TPointerType  mObject;

        
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
      

        
      ///@}    

        
    }; // Class BoundingBox 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template<class TPointType, class TPointerType> 
  inline std::istream& operator >> (std::istream& rIStream, 
				     BoundingBox<TPointType, TPointerType>& rThis);

  /// output stream function
  template<class TPointType, class TPointerType> 
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const BoundingBox<TPointType, TPointerType>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
 
}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined 


