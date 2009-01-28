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
//   Revision:            $Revision: 1.4 $
//
//


#if !defined(KRATOS_TRIANGLE_GAUSSIAN_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_TRIANGLE_GAUSSIAN_INTEGRATION_POINTS_H_INCLUDED



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

    class TriangleGaussianIntegrationPoints1
    {
    public:
      	KRATOS_CLASS_POINTER_DEFINITION(TriangleGaussianIntegrationPoints1);
  	typedef std::size_t SizeType;
	
	static const unsigned int Dimension = 2;

	typedef IntegrationPoint<2> IntegrationPointType;

      	typedef boost::array<IntegrationPointType, 1> IntegrationPointsArrayType;

      	typedef IntegrationPointType::PointType PointType;

      	static SizeType IntegrationPointsNumber(){  return 1; }
	
      	static IntegrationPointsArrayType& IntegrationPoints()
	{
	  // This is added to solve the problem of static initialization. Pooyan.
	  msIntegrationPoints[0] = IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00 , 1.00 / 2.00 );
	  return msIntegrationPoints;
	}

     	std::string Info() const
      	{
	  std::stringstream buffer;
	  buffer << "Triangle gaussian quadrature 1 ";
	  return buffer.str();
      	}
    protected:

    private:

    	static IntegrationPointsArrayType msIntegrationPoints;

    }; // Class TriangleGaussianIntegrationPoints1


/*  TriangleGaussianIntegrationPoints1::IntegrationPointsArrayType TriangleGaussianIntegrationPoints1::msIntegrationPoints =
    {
		IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00 , 1.00 / 2.00 )
};*/

    
    
    
    
    
    
    
    class TriangleGaussianIntegrationPoints2
    {
	    public:
		    KRATOS_CLASS_POINTER_DEFINITION(TriangleGaussianIntegrationPoints2);
		    typedef std::size_t SizeType;
		    
		    static const unsigned int Dimension = 2;

		    typedef IntegrationPoint<2> IntegrationPointType;

		    typedef boost::array<IntegrationPointType, 3> IntegrationPointsArrayType;

		    typedef IntegrationPointType::PointType PointType;

		    static SizeType IntegrationPointsNumber()    {  return 3; }
		    
		    static IntegrationPointsArrayType& IntegrationPoints()
		    {
		      // This is added to solve the problem of static initialization. Pooyan.
		      msIntegrationPoints[0] = IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00 , 1.00 / 6.00 );
		      msIntegrationPoints[1] = IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00 , 1.00 / 6.00 );
		      msIntegrationPoints[2] = IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00 , 1.00 / 6.00 );
			    return msIntegrationPoints;
		    }

		    std::string Info() const
		    {
			    std::stringstream buffer;
			    buffer << "Triangle gaussian quadrature 2 ";
			    return buffer.str();
		    }
	    protected:

	    private:

		    static IntegrationPointsArrayType msIntegrationPoints;

    }; // Class TriangleGaussianIntegrationPoints2


//   TriangleGaussianIntegrationPoints2::IntegrationPointsArrayType TriangleGaussianIntegrationPoints2::msIntegrationPoints =
//     {
//  		IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00 , 1.00 / 6.00 ),
// 		IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00 , 1.00 / 6.00 ),
// 		IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00 , 1.00 / 6.00 )
//    };
  
  
   
   
   
   
   class TriangleGaussianIntegrationPoints3
   {
	   public:
		   KRATOS_CLASS_POINTER_DEFINITION(TriangleGaussianIntegrationPoints3);
		   typedef std::size_t SizeType;
		   
		   static const unsigned int Dimension = 2;

		   typedef IntegrationPoint<2> IntegrationPointType;

		   typedef boost::array<IntegrationPointType, 4> IntegrationPointsArrayType;

		   typedef IntegrationPointType::PointType PointType;

		   static SizeType IntegrationPointsNumber()    {  return 4; }
		    
		   static IntegrationPointsArrayType& IntegrationPoints()
		   {
			   return msIntegrationPoints;
		   }

		   std::string Info() const
		   {
			   std::stringstream buffer;
			   buffer << "Triangle gaussian quadrature 3 ";
			   return buffer.str();
		   }
	   protected:

	   private:

		   static IntegrationPointsArrayType msIntegrationPoints;

   }; // Class TriangleGaussianIntegrationPoints2
   
   
//    TriangleGaussianIntegrationPoints3::IntegrationPointsArrayType TriangleGaussianIntegrationPoints3::msIntegrationPoints =
//     {
//  		IntegrationPointType( 1.00 / 5.00 , 1.00 / 5.00 , 25.00 / 96.00 ),
//  		IntegrationPointType( 3.00 / 5.00 , 1.00 / 5.00 , 25.00 / 96.00 ),
//  		IntegrationPointType( 1.00 / 5.00 , 3.00 / 5.00 , 25.00 / 96.00 ),
// 		IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00 , -27.00 / 96.00 )
//    };
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_TRIANGLE_GAUSSIAN_INTEGRATION_POINTS_H_INCLUDED  defined 


