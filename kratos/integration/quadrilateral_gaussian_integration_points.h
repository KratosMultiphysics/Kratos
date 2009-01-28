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
//   Last Modified by:    $Author: janosch $
//   Date:                $Date: 2007-03-19 10:49:46 $
//   Revision:            $Revision: 1.4 $
//
//


#if !defined(KRATOS_QUADRILATERAL_GAUSSIAN_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_QUADRILATERAL_GAUSSIAN_INTEGRATION_POINTS_H_INCLUDED



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

    class QuadrilateralGaussianIntegrationPoints1
    {
    public:
      	KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralGaussianIntegrationPoints1);
  	typedef std::size_t SizeType;
	
	static const unsigned int Dimension = 2;

	typedef IntegrationPoint<2> IntegrationPointType;

      	typedef boost::array<IntegrationPointType, 1> IntegrationPointsArrayType;

      	typedef IntegrationPointType::PointType PointType;

      	static SizeType IntegrationPointsNumber(){  return 1; }
	
      	static IntegrationPointsArrayType& IntegrationPoints()
	{
	  return msIntegrationPoints;
	}

     	std::string Info() const
      	{
	  std::stringstream buffer;
	  buffer << "Quadrilateral gaussian quadrature 1 ";
	  return buffer.str();
      	}
    protected:

    private:

    	static IntegrationPointsArrayType msIntegrationPoints;

    }; // Class QuadrilateralGaussianIntegrationPoints1

    class QuadrilateralGaussianIntegrationPoints2
    {
	    public:
		    KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralGaussianIntegrationPoints2);
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
			    buffer << "Quadrilateral gaussian quadrature 2 ";
			    return buffer.str();
		    }
	    protected:

	    private:

		    static IntegrationPointsArrayType msIntegrationPoints;

    }; // Class QuadrilateralGaussianIntegrationPoints2


   class QuadrilateralGaussianIntegrationPoints3
   {
	   public:
		   KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralGaussianIntegrationPoints3);
		   typedef std::size_t SizeType;
		   
		   static const unsigned int Dimension = 2;

		   typedef IntegrationPoint<2> IntegrationPointType;

		   typedef boost::array<IntegrationPointType, 9> IntegrationPointsArrayType;

		   typedef IntegrationPointType::PointType PointType;

		   static SizeType IntegrationPointsNumber()    {  return 9; }
		    
		   static IntegrationPointsArrayType& IntegrationPoints()
		   {
			   return msIntegrationPoints;
		   }

		   std::string Info() const
		   {
			   std::stringstream buffer;
			   buffer << "Quadrilateral gaussian quadrature 3 ";
			   return buffer.str();
		   }
	   protected:

	   private:

		   static IntegrationPointsArrayType msIntegrationPoints;

   }; // Class QuadrilateralGaussianIntegrationPoints2

  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_TRIANGLE_GAUSSIAN_INTEGRATION_POINTS_H_INCLUDED  defined 


