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
//   Date:                $Date: 2010-10-08 16:07:33 $
//   Revision:            $Revision: 1.0$
//
//


#if !defined(KRATOS_CELL_H_INCLUDED)
#define  KRATOS_CELL_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <cmath>
#include <algorithm>




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
  
  template< class  TConfigure> 
  class Cell
    {
    public:
      ///@name Type Definitions
      ///@{
     
	
      /// configure types
      typedef std::size_t  SizeType;      
      typedef typename TConfigure::PointType               PointType;
      typedef typename TConfigure::PointerType             PointerType;
      typedef typename TConfigure::ContainerType           ContainerType;
      typedef typename TConfigure::IteratorType            IteratorType;
      typedef typename TConfigure::ResultContainerType     ResultContainerType; 
      typedef typename TConfigure::ResultIteratorType      ResultIteratorType;
      
      ///configure Contact Pair
      typedef typename TConfigure::ContainerContactType  ContainerContactType;
      typedef typename TConfigure::ContactPairType       ContactPairType;
      typedef typename TConfigure::IteratorContactType   IteratorContactType;
     
      typedef std::vector<PointerType>     LocalContainerType;
      typedef typename LocalContainerType::iterator LocalIteratorType;

      
      /// Pointer definition of Cell
      KRATOS_CLASS_POINTER_DEFINITION(Cell);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Cell()
      {
      }

      /// Destructor.
      virtual ~Cell(){}
      
      
      void Add(const PointerType& ThisObject)
      {
	mObjects.push_back(ThisObject);
      }
      
      void Remove(const PointerType& ThisObject)
      {
       mObjects.erase(std::remove(mObjects.begin(),mObjects.end(),ThisObject),mObjects.end());
      }
      
      void Clear()
      {
	mObjects.clear();
      }
      
      
      void AllocateCell(const std::size_t size)
      {
	mObjects.reserve(size);
      }
      
      
      /// Assignment operator.
      Cell& operator=(Cell const& rOther)
       {
	  mObjects   = rOther.mObjects;
	  return *this;
       }
       

      /// Copy constructor.
      Cell(Cell const& rOther) : 
          mObjects(rOther.mObjects)    
       {             
       }

//************************************************************************   
//************************************************************************       

      void SearchObjects(PointerType& rThisObject, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults)
	    {
	      for(LocalIteratorType i_object = Begin() ; i_object != End()  && NumberOfResults < MaxNumberOfResults ; i_object++){
		if(TConfigure::Intersection(rThisObject, *i_object)) {
		ResultIteratorType repeated_object = std::find(Result-NumberOfResults, Result, *i_object);	
		if(repeated_object==Result) 
		  {
		    *Result   = *i_object;
		    Result++;
		    NumberOfResults++;
		    }
		}
	      }
	    }  
      
//************************************************************************   
//************************************************************************       

      void SearchObjects(PointerType& rThisObject, ResultContainerType& Result)
	    {
	      for(LocalIteratorType i_object = Begin() ; i_object != End(); i_object++){
	      if(TConfigure::Intersection(rThisObject, *i_object))
	      {
		ResultIteratorType repeated_object = std::find(Result.begin(), Result.end(), *i_object);
		if(repeated_object==Result.end()) 
		  {   
		    Result.push_back(*i_object);
		  }
		}
	      }
	    }   
      
//************************************************************************   
//************************************************************************       

      void SearchObjectsInner(PointerType& rThisObject, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults)
      {
        for(LocalIteratorType i_object = Begin() ; i_object != End()  && NumberOfResults < MaxNumberOfResults ; i_object++){
          if( rThisObject != *i_object )
          {
            if(TConfigure::Intersection(rThisObject, *i_object)) {
              ResultIteratorType repeated_object = std::find(Result-NumberOfResults, Result, *i_object);	
              if(repeated_object==Result) 
              {
                *Result   = *i_object;
                Result++;
                NumberOfResults++;
              }
            }
          }
        }
      }  
      
//************************************************************************   
//************************************************************************       

      void SearchObjectsInner(PointerType& rThisObject, ResultContainerType& Result)
      {
        for(LocalIteratorType i_object = Begin() ; i_object != End(); i_object++){
          if( rThisObject != *i_object )
          {
            if(TConfigure::Intersection(rThisObject, *i_object))
            {
              ResultIteratorType repeated_object = std::find(Result.begin(), Result.end(), *i_object);
              if(repeated_object==Result.end()) 
              {   
                Result.push_back(*i_object);
              }
            }
          }
        }
      }   
      
//************************************************************************   
//************************************************************************ 
  
 
      void SearchContact(ContainerContactType& Result)
	    {
	      ContactPairType Pair;
	      for(LocalIteratorType i_object_1 = Begin(); i_object_1 != End(); i_object_1 ++)
		      { 
			Pair[0] = *i_object_1;
			for(LocalIteratorType i_object_2 = i_object_1 + 1 ; i_object_2!= End(); i_object_2++ )
			{    
			  if(TConfigure::Intersection(*i_object_1, *i_object_2)){
			    Pair[1] = *i_object_2;
				IteratorContactType repeated_par_1 = std::find(Result.begin(), Result.end(), Pair);
				if(repeated_par_1==Result.end())  
				    Result.push_back(Pair);}
			}
		  }
	    } 


//************************************************************************   
//************************************************************************ 
     
      void SearchContact(IteratorContactType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults )
	    {
		ContactPairType Pair;
		for(LocalIteratorType i_object_1 = Begin(); i_object_1 != End(); i_object_1 ++)
		{ 
		  Pair[0] = *i_object_1; 
		  for(LocalIteratorType i_object_2 = i_object_1 + 1 ; i_object_2!= End() && NumberOfResults < MaxNumberOfResults ; i_object_2++ )
		  {    
		    Pair[1] = *i_object_2;
		    if(TConfigure::Intersection(*i_object_1, *i_object_2)){
			  IteratorContactType repeated_par_1 = std::find(Result-NumberOfResults, Result, Pair);
			  if(repeated_par_1==Result){
			      *Result =  Pair;
			      NumberOfResults++;
			      Result++;
			    }
			}
		  }
		}
	    }

//************************************************************************   
//************************************************************************  

      LocalIteratorType Begin() 
	    {
	      return mObjects.begin();
	    }

      LocalIteratorType End()  
	    {
	      return mObjects.end();
	    }
      
      SizeType Size()
      {
        return mObjects.size();
      }
       
      LocalIteratorType Begin() const   
	      { 	 
		return mObjects.begin();
	      }
        
      LocalIteratorType End() const  
	    {
	      return mObjects.end();
	    }
      
      SizeType Size() const
      {
        return mObjects.size();
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
      virtual std::string Info() const
      {
	return "Cell Class "; 
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const{ return; }

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const{return;}

      
      

            
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
        

      std::vector<PointerType> mObjects;
        
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
        
    }; // Class Cell 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template< class  TConfigure> 
  inline std::istream& operator >> (std::istream& rIStream, 
				    Cell<TConfigure>& rThis){ return rIStream;}

  /// output stream function
  template< class  TConfigure> 
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Cell<TConfigure>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined 


