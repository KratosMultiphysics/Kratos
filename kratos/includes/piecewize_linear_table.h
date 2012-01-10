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
#include <boost/array.hpp>  


// Project includes
#include "includes/define.h"
#include "containers/variable.h"


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
  /** Table class stores the value of its second variable respect to the value of its first variable.
  *   It also provides a piecewise linear interpolator/extrapolator for getting intermediate values.
  */
  template<class TArgumentType, class TResultType = TArgumentType, std::size_t TResultsColumns = 1>
  class Table
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Table
      KRATOS_CLASS_POINTER_DEFINITION(Table);

	  typedef TArgumentType argument_type; // To be STL conformance.
	  typedef TResultType result_type; // To be STL conformance.

	  typedef boost::array<TResultType, TResultsColumns>  result_array_type;

      typedef std::pair<argument_type, result_array_type> RecordType;

	  typedef std::vector<RecordType> TableContainerType;

	  typedef Variable<TArgumentType> XVariableType;
	  typedef Variable<TResultType> YVariableType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Table() : mData(), mpXVariable(NULL) , mpYVariable(NULL)
	  {
	  }

      /// Default constructor.
      Table(XVariableType const& XVariable, YVariableType const& YVariable) : mData(), mpXVariable(&XVariable) , mpYVariable(&YVariable)
	  {
	  }

      /// Destructor.
      virtual ~Table()
	  {
	  }
      
      
      /// Assignment operator.
      Table& operator=(Table const& rOther)
	  {
		  mData = rOther.mData;
		  return *this;
	  }


      ///@}
      ///@name Operators 
      ///@{

	  // I want to put operator(i,j) for accessing, operator(i) for first column and operator[i] for getting the complete row
      
	  // This operator calculates the piecewise linear interpolation for 
	  // given argument
	  result_type operator()(argument_type const& X) const
	  {
		  return GetValue(X);
	  }

	  // This operator gives the result for the nearest value to argument found in table 
	  result_type const & operator[](argument_type const& X) const
	  {
              GetNearestValue(X);
	  }

	  // This operator gives the result for the nearest value to argument found in table 
	  result_type & operator[](argument_type& X)
	  {
              GetNearestValue(X);
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

		  result_type result;
		  if(X <= mData[0].first)
			  return Interpolate(X, mData[0].first, mData[0].second, mData[1].first, mData[1].second, result);

		  for(std::size_t i = 1 ; i < size ; i++)
			if(X <= mData[i].first)
			  return Interpolate(X, mData[i-1].first, mData[i-1].second, mData[i].first, mData[i].second, result);

		  // now the x is outside the table and we hae to extrapolate it using last two records of table.
		  return Interpolate(X, mData[size-2].first, mData[size-2].second, mData[size-1].first, mData[size-1].second, result);
	  }

	  // Get the nesrest value for the given argument
	  result_array_type& GetNearestRow(argument_type const& X)
	  {
		  std::size_t size = mData.size();

		  if(size == 0)
			  KRATOS_ERROR(std::invalid_argument, "Get value from empty table", "");

		  if(size==1) // constant table. Returning the only value we have.
			  return mData.begin()->second;

		  if(X <= mData[0].first)
			  return mData[0].second;

		  for(std::size_t i = 1 ; i < size ; i++)
			if(X <= mData[i].first)
			  return ((X - mData[i-1].first) < (mData[i].first - X)) ? mData[i-1].second : mData[i].second;

		  // now the x is outside the table and we hae to extrapolate it using last two records of table.
		  return mData[size-1].second;
	  }

	  // Get the nesrest value for the given argument
	  result_type const& GetNearestValue(argument_type const& X)  const
	  {
		  std::size_t size = mData.size();

		  if(size == 0)
			  KRATOS_ERROR(std::invalid_argument, "Get value from empty table", "");

		  if(size==1) // constant table. Returning the only value we have.
			  return mData.begin()->second;

		  if(X <= mData[0].first)
			  return mData[0].second;

		  for(std::size_t i = 1 ; i < size ; i++)
			if(X <= mData[i].first)
			  return ((X - mData[i-1].first) < (mData[i].first - X)) ? mData[i-1].second : mData[i].second;

		  // now the x is outside the table and we hae to extrapolate it using last two records of table.
		  return mData[size-1].second;
        }

        result_type& Interpolate(argument_type const& X, argument_type const& X1, result_type const& Y1, argument_type const& X2, result_type const& Y2, result_type& Result)
        {
            double epsilon = 1e-12;

            double dx = X2 - X1;
            result_type dy = Y2 - Y1;

            double scale = 0.00;

            if (dx > epsilon)
                scale = (X - X1) / dx_norm;

            Result = Y1 + dy * scale;

            return Result;

        }

          void PushBack(argument_type const& X, result_type const& Y)
          {
			  result_array_type a = {{Y}};
              mData.push_back(RecordType(X,a));
          }
     
      ///@}
      ///@name Access
      ///@{ 
      

          TableContainerType& Data()
          {
              return mData;
          }

          TableContainerType const& Data() const
          {
              return mData;
          }

          XVariableType& GetXVariable()
          {
              return *mpXVariable;
          }

          YVariableType& GetYVariable()
          {
              return *mpYVariable;
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
		  return "Table";
	  }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	  {
		  rOStream << Info();
	  }

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
	  {
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
        
		TableContainerType mData;
                const XVariableType* mpXVariable;
                const YVariableType* mpYVariable;
        
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
  template<class TArgumentType, class TResultType>
  inline std::istream& operator >> (std::istream& rIStream, 
				    Table<TArgumentType, TResultType>& rThis);

  /// output stream function
  template<class TArgumentType, class TResultType>
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Table<TArgumentType, TResultType>& rThis)
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


