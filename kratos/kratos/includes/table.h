//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_TABLE_H_INCLUDED )
#define  KRATOS_TABLE_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"


namespace Kratos
{
  ///@addtogroup Kratos Core
  ///@{

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
  
  /// This class represents the value of its variable depending to other variable.
  /** Table class stores the value of its first variable respect to the value of its second variable.
  *   It also provides a piecewise linear interpolator/extrapolator for getting intermediate values.
  */
  template<class TArgumentType, class TResultType>
  class Table
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Table
      KRATOS_CLASS_POINTER_DEFINITION(Table);

	  typedef TArgumentType argument_type; // To be STL conformance.
	  typedef TResultType result_type; // To be STL conformance.

	  typedef std::vector<std::pair<TArgumentType,TResultType> > TableContainerType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Table()
	  {
	  }

      /// Destructor.
      virtual ~Table()
	  {
	  }
      

      ///@}
      ///@name Operators 
      ///@{
      
	  // This operator calculates the piecewise linear interpolation for 
	  // given argument
	  result_type operator()(argument_type const& X) const
	  {
		  return GetValue(X);
	  }

	  // This operator gives the result for the nearest value to argument found in table 
	  result_type const & operator[](argument_type const& X) const
	  {
	  }

	  // This operator gives the result for the nearest value to argument found in table 
	  result_type & operator[](argument_type& X)
	  {
	  }
      
      ///@}
      ///@name Operations
      ///@{
      
	  // Get the value for the given argument using piecewise linear
	  result_type GetValue(argument_type const& X) const
	  {
		  std::size_t size = mData.size();

		  if(size == 0)
			  KRATOS_ERROR(std::invalid_argument, "Get value from empty table", "");

		  if(size==1) // constant table. Returning the only value we have.
			  return mData.begin()->second;

		  result_type y1 = mData[0].second;
		  result_type y2 = mData[1].second;

		  result_type result;
		  if(x <= mData[0].first)
			  return Interpolate(x, mData[0].first, mData[0].second, mData[1].first, mData[1].second, result);

		  for(std::size_t i = 1 ; i < size ; i++)
			if(x <= mData[i].first)
			  return Interpolate(x, mData[i-1].first, mData[i-1].second, mData[i].first, mData[i].second, result);

		  // now the x is outside the table and we hae to extrapolate it using last two records of table.
		  return Interpolate(x, mData[size-2].first, mData[size-2].second, mData[size-1].first, mData[size-1].second, result);
	  }

	  result_type& Interpolate(argument_type const& X, argument_type const& X1, result_type const& Y1, argument_type const& X2, result_type const& Y2, Result)
	  {
		  argument_type dx = X2-X1;
		  result_type dy=Y2-Y1;

		  double dx_norm = NormOf(dx);
		  double scale = 0.00;

		  if(dx_norm)
			  scale = NormOf(X-X1) / dx_norm;

		  Result = Y1 + dy * scale;

		  return Result;

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
      virtual std::string Info() const;
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const;

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const;
      
            
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
        
		TableContainerType mData;
        
      ///@} 
      ///@name Private Operators
      ///@{ 

		double NormOf(double Value)
		{
			return Value;
		}
        
        
        
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
      Table& operator=(Table const& rOther);

      /// Copy constructor.
      Table(Table const& rOther);

        
      ///@}    
        
    }; // Class Table 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    Table& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Table& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_TABLE_H_INCLUDED  defined 


