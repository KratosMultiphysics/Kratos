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
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:33 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_QUADRATURE_H_INCLUDED )
#define  KRATOS_QUADRATURE_H_INCLUDED



// System includes
#include <string>
#include <vector>
#include <iostream>
#include <cstddef>


// External includes 


// Project includes
#include "includes/define.h"
#include "integration/integration_point.h"


namespace Kratos
{
  template<std::size_t TOrder> 
    struct Pow
    {
      static inline std::size_t Value(std::size_t X)
      {
	return X * Pow<TOrder - 1>::Value(X);
      }
    };

  template<> 
    struct Pow<1>
    {
      static inline std::size_t Value(std::size_t X)
      {
	return X;
      }
    };
 
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
	template<class TQuadraturePointsType, std::size_t TDimension = TQuadraturePointsType::Dimension, class TIntegrationPointType = IntegrationPoint<TDimension> >
  class Quadrature
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Quadrature
      KRATOS_CLASS_POINTER_DEFINITION(Quadrature);
  
      typedef TIntegrationPointType IntegrationPointType;
 
	  typedef std::vector<IntegrationPointType> IntegrationPointsArrayType;

	  typedef typename TQuadraturePointsType::IntegrationPointsArrayType InitialIntegrationPointsArrayType;

      typedef typename IntegrationPointType::PointType PointType;
 
	  typedef std::size_t SizeType;
 
	  typedef std::size_t IndexType;
 
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Quadrature(){}

      /// Copy constructor.
      Quadrature(Quadrature const& rOther){}

       /// Destructor.
      virtual ~Quadrature(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      /// Assignment operator.
      Quadrature& operator=(Quadrature const& rOther){}

      
      ///@}
      ///@name Operations
      ///@{
      
       static SizeType IntegrationPointsNumber()
	{
	  return Pow<TDimension - TQuadraturePointsType::Dimension + 1>::Value(TQuadraturePointsType::IntegrationPointsNumber());
	}

     static const IntegrationPointsArrayType& IntegrationPoints()
       {
 	    return msIntegrationPointsVector; 
       }

       static IntegrationPointsArrayType GenerateIntegrationPoints()
 	{
 	  IntegrationPointsArrayType results;
 	  GenerateIntegrationPoints(results, Quadrature<TQuadraturePointsType, TDimension, TIntegrationPointType>());
 	  return results;
 	}

      template<class TIntegrationPointsArrayType>
      static void GenerateIntegrationPoints(TIntegrationPointsArrayType& Result,
					    Quadrature<TQuadraturePointsType, TDimension, TIntegrationPointType> const& Dummy)
	{
	  typename Quadrature<TQuadraturePointsType, TDimension - 1, TIntegrationPointType>::IntegrationPointsArrayType points;
	  Quadrature<TQuadraturePointsType, TDimension - 1, TIntegrationPointType>::GenerateIntegrationPoints(points, Quadrature<TQuadraturePointsType, TDimension - 1, TIntegrationPointType>());

	  for(SizeType i = 0 ; i < TQuadraturePointsType::IntegrationPointsNumber() ; i++)
	    for(SizeType j = 0 ; j < points.size() ; j++)
	      {
		IntegrationPointType temp(points[j]);
 		temp[TDimension -1] = TQuadraturePointsType::IntegrationPoints()[i].X();
		temp.Weight() *= TQuadraturePointsType::IntegrationPoints()[i].Weight();
		Result.push_back(temp);
	      }
	}

      static void GenerateIntegrationPoints(IntegrationPointsArrayType& Result,
					    Quadrature<TQuadraturePointsType, TQuadraturePointsType::Dimension, TIntegrationPointType> const& Dummy)
      {
		  for(SizeType i = 0 ; i < TQuadraturePointsType::IntegrationPointsNumber() ; i++)
	      {
			Result.push_back(TQuadraturePointsType::IntegrationPoints()[i]);
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
	  std::stringstream buffer;
	  buffer << TDimension << " dimensional quadrature with " << IntegrationPointsNumber() << " integration points";
	  return buffer.str();
      }
      
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
      {
	  rOStream << Info();
      }

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
      {
 	IndexType i;
 	for(i = 0 ; i < msIntegrationPointsVector.size() - 1 ; i++)
 	  rOStream << msIntegrationPointsVector[i] << " , " << std::endl;
 	rOStream << msIntegrationPointsVector[i];
      }
      
            
      ///@}      
      ///@name Friends
      ///@{

	  template<class TOtherQuadraturePointsType, std::size_t TOtherDimension, class TOtherIntegrationPointType> friend class Quadrature;
      
            
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

		

      static const IntegrationPointsArrayType msIntegrationPointsVector;
        
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
        
        
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
        
    }; // Class Quadrature 

  ///@}
#ifdef __SUNPRO_CC
  template<class TQuadraturePointsType, std::size_t TDimension, class TIntegrationPointType>
   const typename Quadrature<TQuadraturePointsType, TDimension, TIntegrationPointType>::IntegrationPointsArrayType 
   Quadrature<TQuadraturePointsType, TDimension, TIntegrationPointType>::msIntegrationPointsVector = 
     Quadrature<TQuadraturePointsType, TDimension, TIntegrationPointType>::GenerateIntegrationPoints();
#else
  template<class TQuadraturePointsType, std::size_t TDimension, class TIntegrationPointType> const typename Quadrature<TQuadraturePointsType, TDimension, TIntegrationPointType>::IntegrationPointsArrayType 		Quadrature<TQuadraturePointsType, TDimension, TIntegrationPointType>::msIntegrationPointsVector( Quadrature<TQuadraturePointsType, TDimension, TIntegrationPointType>::GenerateIntegrationPoints());
#endif
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template<class TQuadraturePointsType, std::size_t TDimension, class TIntegrationPointType>
  inline std::istream& operator >> (std::istream& rIStream, 
				    Quadrature<TQuadraturePointsType, TDimension, TIntegrationPointType>& rThis);

  /// output stream function
  template<class TQuadraturePointsType, std::size_t TDimension, class TIntegrationPointType>
  inline std::ostream& operator << (std::ostream& rOStream,
				    const Quadrature<TQuadraturePointsType, TDimension, TIntegrationPointType>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_QUADRATURE_H_INCLUDED  defined 


