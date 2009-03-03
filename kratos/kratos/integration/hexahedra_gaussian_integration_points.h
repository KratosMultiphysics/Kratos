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
//   Date:                $Date: 2007-03-19 10:49:03 $
//   Revision:            $Revision: 1.4 $
//
//


#if !defined(KRATOS_HEXAHEDRA_GAUSSIAN_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_HEXAHEDRA_GAUSSIAN_INTEGRATION_POINTS_H_INCLUDED



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

    class HexahedraGaussianIntegrationPoints1
    {
    public:
      	KRATOS_CLASS_POINTER_DEFINITION(HexahedraGaussianIntegrationPoints1);
  	typedef std::size_t SizeType;
	
	static const unsigned int Dimension = 3;

	typedef IntegrationPoint<3> IntegrationPointType;

      	typedef boost::array<IntegrationPointType, 1> IntegrationPointsArrayType;

      	typedef IntegrationPointType::PointType PointType;

      	static SizeType IntegrationPointsNumber(){  return 1; }
	
      	static IntegrationPointsArrayType& IntegrationPoints()
	{
	  // This is added to solve the problem of static initialization. Pooyan.
		 msIntegrationPoints[0] = IntegrationPointType( 0.00 , 0.00, 0.00 , 8.00 );
	  return msIntegrationPoints;
	}

     	std::string Info() const
      	{
	  std::stringstream buffer;
	  buffer << "Hexahedra gaussian quadrature 1 ";
	  return buffer.str();
      	}
    protected:

    private:

    	static IntegrationPointsArrayType msIntegrationPoints;

    }; // Class HexahedraGaussianIntegrationPoints1






    class HexahedraGaussianIntegrationPoints2
    {
	    public:
		    KRATOS_CLASS_POINTER_DEFINITION(HexahedraGaussianIntegrationPoints2);
		    typedef std::size_t SizeType;
		    
		    static const unsigned int Dimension = 3;

		    typedef IntegrationPoint<3> IntegrationPointType;

		    typedef boost::array<IntegrationPointType, 8> IntegrationPointsArrayType;

		    typedef IntegrationPointType::PointType PointType;

		    static SizeType IntegrationPointsNumber()    {  return 8; }
		    
		    static IntegrationPointsArrayType& IntegrationPoints()
		    {
				// This is added to solve the problem of static initialization. Pooyan.
        		 msIntegrationPoints[0] = IntegrationPointType( -1.00/std::sqrt(3.0) , -1.00/std::sqrt(3.0), -1.00/std::sqrt(3.0), 1.00 );
        		 msIntegrationPoints[1] = IntegrationPointType(  1.00/std::sqrt(3.0) , -1.00/std::sqrt(3.0), -1.00/std::sqrt(3.0), 1.00 );
        		 msIntegrationPoints[2] = IntegrationPointType(  1.00/std::sqrt(3.0) ,  1.00/std::sqrt(3.0), -1.00/std::sqrt(3.0), 1.00 );
        		 msIntegrationPoints[3] = IntegrationPointType( -1.00/std::sqrt(3.0) ,  1.00/std::sqrt(3.0), -1.00/std::sqrt(3.0), 1.00 );
        		 msIntegrationPoints[4] = IntegrationPointType( -1.00/std::sqrt(3.0) , -1.00/std::sqrt(3.0),  1.00/std::sqrt(3.0), 1.00 );
        		 msIntegrationPoints[5] = IntegrationPointType(  1.00/std::sqrt(3.0) , -1.00/std::sqrt(3.0),  1.00/std::sqrt(3.0), 1.00 );
        		 msIntegrationPoints[6] = IntegrationPointType(  1.00/std::sqrt(3.0) ,  1.00/std::sqrt(3.0),  1.00/std::sqrt(3.0), 1.00 );
        		 msIntegrationPoints[7] = IntegrationPointType( -1.00/std::sqrt(3.0) ,  1.00/std::sqrt(3.0),  1.00/std::sqrt(3.0), 1.00 );

			    return msIntegrationPoints;
		    }

		    std::string Info() const
		    {
			    std::stringstream buffer;
			    buffer << "Hexahedra gaussian quadrature 2 ";
			    return buffer.str();
		    }
	    protected:

	    private:

		    static IntegrationPointsArrayType msIntegrationPoints;

    }; // Class HexahedraGaussianIntegrationPoints2





   class HexahedraGaussianIntegrationPoints3
   {
	   public:
		   KRATOS_CLASS_POINTER_DEFINITION(HexahedraGaussianIntegrationPoints3);
		   typedef std::size_t SizeType;
		   
		   static const unsigned int Dimension = 3;

		   typedef IntegrationPoint<3> IntegrationPointType;

		   typedef boost::array<IntegrationPointType, 27> IntegrationPointsArrayType;

		   typedef IntegrationPointType::PointType PointType;

		   static SizeType IntegrationPointsNumber()    {  return 27; }
		    
		   static IntegrationPointsArrayType& IntegrationPoints()
		   {
			   return msIntegrationPoints;
		   }

		   std::string Info() const
		   {
			   std::stringstream buffer;
			   buffer << "Hexadra gaussian quadrature 3 ";
			   return buffer.str();
		   }
	   protected:

	   private:

		   static IntegrationPointsArrayType msIntegrationPoints;

   }; // Class HexahedraGaussianIntegrationPoints2
   
   
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_TETRAHEDRA_GAUSSIAN_INTEGRATION_POINTS_H_INCLUDED  defined 


