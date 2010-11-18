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
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2008-11-17 17:46:04 $
//   Revision:            $Revision: 1.6 $
//
//


#if !defined(KRATOS_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <numeric>
#include <cstddef>


// External includes
#include <boost/array.hpp>  


// Project includes
#include "includes/define.h"
#include "integration/integration_point.h"


namespace Kratos
{
    class GaussLegendreIntegrationPoints1
    {
    public:
      	KRATOS_CLASS_POINTER_DEFINITION(GaussLegendreIntegrationPoints1);
  	typedef std::size_t SizeType;
	
	static const unsigned int Dimension = 1;

	typedef IntegrationPoint<1> IntegrationPointType;

      	typedef boost::array<IntegrationPointType, 1> IntegrationPointsArrayType;

      	typedef IntegrationPointType::PointType PointType;

      	static SizeType IntegrationPointsNumber(){  return 1; }
	
      	static IntegrationPointsArrayType& IntegrationPoints()
	{
	  // This is added to solve the problem of static initialization. Pooyan.
	  msIntegrationPoints[0] = IntegrationPointType(0.00, 2.00);
	  return msIntegrationPoints;
	}

     	std::string Info() const
      	{
	  std::stringstream buffer;
	  buffer << "GaussLegendre quadrature 1 ";
	  return buffer.str();
      	}
    protected:

    private:

    	static IntegrationPointsArrayType msIntegrationPoints;

    }; // Class QuadrilateralGaussianIntegrationPoints1


   class GaussLegendreIntegrationPoints2
    {
    public:
      	KRATOS_CLASS_POINTER_DEFINITION(GaussLegendreIntegrationPoints2);
  	typedef std::size_t SizeType;
	
	static const unsigned int Dimension = 1;

	typedef IntegrationPoint<1> IntegrationPointType;

      	typedef boost::array<IntegrationPointType, 2> IntegrationPointsArrayType;

      	typedef IntegrationPointType::PointType PointType;

      	static SizeType IntegrationPointsNumber(){  return 2; }
	
      	static IntegrationPointsArrayType& IntegrationPoints()
	{
	  // This is added to solve the problem of static initialization. Pooyan.
	  msIntegrationPoints[0] = IntegrationPointType(-std::sqrt(1.00 / 3.00), 1.00);
	  msIntegrationPoints[1] = IntegrationPointType( std::sqrt(1.00 / 3.00), 1.00);
	  return msIntegrationPoints;
	}

     	std::string Info() const
      	{
	  std::stringstream buffer;
	  buffer << "GaussLegendre quadrature 2 ";
	  return buffer.str();
      	}
    protected:

    private:

    	static IntegrationPointsArrayType msIntegrationPoints;

    }; // Class QuadrilateralGaussianIntegrationPoints2
    

   class GaussLegendreIntegrationPoints3
    {
    public:
      	KRATOS_CLASS_POINTER_DEFINITION(GaussLegendreIntegrationPoints3);
  	typedef std::size_t SizeType;
	
	static const unsigned int Dimension = 1;

	typedef IntegrationPoint<1> IntegrationPointType;

      	typedef boost::array<IntegrationPointType, 3> IntegrationPointsArrayType;

      	typedef IntegrationPointType::PointType PointType;

      	static SizeType IntegrationPointsNumber(){  return 3; }
	
      	static IntegrationPointsArrayType& IntegrationPoints()
	{
	  // This is added to solve the problem of static initialization. Pooyan.
	  msIntegrationPoints[0] = IntegrationPointType(-std::sqrt(3.00 / 5.00), 5.00 / 9.00);
	  msIntegrationPoints[1] = IntegrationPointType( 0.00                  , 8.00 / 9.00);
	  msIntegrationPoints[2] = IntegrationPointType( std::sqrt(3.00 / 5.00), 5.00 / 9.00);
	  return msIntegrationPoints;
	}

     	std::string Info() const
      	{
	  std::stringstream buffer;
	  buffer << "GaussLegendre quadrature 3 ";
	  return buffer.str();
      	}
    protected:

    private:

    	static IntegrationPointsArrayType msIntegrationPoints;

    }; // Class QuadrilateralGaussianIntegrationPoints3



   class GaussLegendreIntegrationPoints4
    {
    public:
      	KRATOS_CLASS_POINTER_DEFINITION(GaussLegendreIntegrationPoints4);
  	typedef std::size_t SizeType;
	
	static const unsigned int Dimension = 1;

	typedef IntegrationPoint<1> IntegrationPointType;

      	typedef boost::array<IntegrationPointType, 4> IntegrationPointsArrayType;

      	typedef IntegrationPointType::PointType PointType;

      	static SizeType IntegrationPointsNumber(){  return 4; }
	
      	static IntegrationPointsArrayType& IntegrationPoints()
	{
	  // This is added to solve the problem of static initialization. Pooyan.
	  msIntegrationPoints[0] = IntegrationPointType(-0.861136311594053, 0.347854845137454);
	  msIntegrationPoints[1] = IntegrationPointType(-0.339981043584856, 0.652145154862546);
	  msIntegrationPoints[2] = IntegrationPointType( 0.339981043584856, 0.652145154862546);
	  msIntegrationPoints[3] = IntegrationPointType( 0.861136311594053, 0.347854845137454);
	  return msIntegrationPoints;
	}

     	std::string Info() const
      	{
	  std::stringstream buffer;
	  buffer << "GaussLegendre quadrature 4 ";
	  return buffer.str();
      	}
    protected:

    private:

    	static IntegrationPointsArrayType msIntegrationPoints;

    }; // Class QuadrilateralGaussianIntegrationPoints4


   
   class GaussLegendreIntegrationPoints5
    {
    public:
      	KRATOS_CLASS_POINTER_DEFINITION(GaussLegendreIntegrationPoints5);
  	typedef std::size_t SizeType;
	
	static const unsigned int Dimension = 1;

	typedef IntegrationPoint<1> IntegrationPointType;

      	typedef boost::array<IntegrationPointType, 5> IntegrationPointsArrayType;

      	typedef IntegrationPointType::PointType PointType;

      	static SizeType IntegrationPointsNumber(){  return 5; }
	
      	static IntegrationPointsArrayType& IntegrationPoints()
	{
	  // This is added to solve the problem of static initialization. Pooyan.
	  msIntegrationPoints[0] = IntegrationPointType(-0.906179845938664, 0.236926885056189);
	  msIntegrationPoints[1] = IntegrationPointType(-0.538469310105683, 0.478628670499366);
	  msIntegrationPoints[2] = IntegrationPointType( 0.000000000000000, 0.568888888888889);
	  msIntegrationPoints[3] = IntegrationPointType( 0.538469310105683, 0.478628670499366);
	  msIntegrationPoints[4] = IntegrationPointType( 0.906179845938664, 0.236926885056189);
	  return msIntegrationPoints;
	}

     	std::string Info() const
      	{
	  std::stringstream buffer;
	  buffer << "GaussLegendre quadrature 5 ";
	  return buffer.str();
      	}
    protected:

    private:

    	static IntegrationPointsArrayType msIntegrationPoints;

    }; // Class QuadrilateralGaussianIntegrationPoints4









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
// 	template<std::size_t TOrder>
//     class GaussLegendreIntegrationPoints
//     {
//     public:
//      ///@name Type Definitions
//       ///@{
//       
//       /// Pointer definition of GaussLegendreIntegrationPoints
//       KRATOS_CLASS_POINTER_DEFINITION(GaussLegendreIntegrationPoints);
// 
// 	  typedef std::size_t SizeType;
//   
// 	  static const SizeType Dimension = 1;
// 
//       typedef IntegrationPoint<Dimension> IntegrationPointType;
// 
//       typedef boost::array<IntegrationPointType, TOrder> IntegrationPointsArrayType;
// 
//       typedef typename IntegrationPointType::PointType PointType;
 
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      ///@}
      ///@name Operators 
      ///@{
      
      ///@}
      ///@name Operations
      ///@{

//       static SizeType IntegrationPointsNumber()
// 	{
// 	  return TOrder;
// 	}
//       
//       static IntegrationPointsArrayType& IntegrationPoints()
// 	{
// 	  return msIntegrationPoints;
// 	}
      
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
//       std::string Info() const
//       {
// 	  std::stringstream buffer;
// 	  buffer << "Gauss-Legendre Integration points with order " << TOrder;
// 	  return buffer.str();
//       }
      
      
      ///@}      
      ///@name Friends
      ///@{

            
      ///@}
      
//     protected:
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
      
//     private:
      ///@name Static Member Variables 
      ///@{

//       static IntegrationPointsArrayType msIntegrationPoints;

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
        
//     }; // Class GaussLegendreIntegrationPoints 

  ///@}



//   template<> GaussLegendreIntegrationPoints<1>::IntegrationPointsArrayType GaussLegendreIntegrationPoints<1>::msIntegrationPoints =
//     {
// 		{
// 		IntegrationPointType(0.00, 2.00)
// 		}
//     };
  
//   template<> GaussLegendreIntegrationPoints<2>::IntegrationPointsArrayType GaussLegendreIntegrationPoints<2>::msIntegrationPoints =
//     {
// 		{
//       IntegrationPointType(-std::sqrt(1.00 / 3.00), 1.00),
//       IntegrationPointType( std::sqrt(1.00 / 3.00), 1.00)
// 		}
//     };
/*  
  template<> GaussLegendreIntegrationPoints<3>::IntegrationPointsArrayType GaussLegendreIntegrationPoints<3>::msIntegrationPoints =
    {
		{
      IntegrationPointType(-std::sqrt(3.00 / 5.00), 5.00 / 9.00),
      IntegrationPointType( 0.00                  , 8.00 / 9.00),
      IntegrationPointType( std::sqrt(3.00 / 5.00), 5.00 / 9.00)
		}
    };
*/
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED  defined 


