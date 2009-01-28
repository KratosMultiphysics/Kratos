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
//   Date:                $Date: 2007-12-13 14:22:10 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_PRISM_GAUSSIAN_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_PRISM_GAUSSIAN_INTEGRATION_POINTS_H_INCLUDED



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

    class PrismGaussianIntegrationPoints1
    {
    public:
      	KRATOS_CLASS_POINTER_DEFINITION(PrismGaussianIntegrationPoints1);
  	typedef std::size_t SizeType;
	
	static const unsigned int Dimension = 3;

	typedef IntegrationPoint<3> IntegrationPointType;

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
	  buffer << "Prism gaussian quadrature 1 ";
	  return buffer.str();
      	}
    protected:

    private:

    	static IntegrationPointsArrayType msIntegrationPoints;

    }; // Class PrismGaussianIntegrationPoints1






    class PrismGaussianIntegrationPoints2
    {
	    public:
		    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussianIntegrationPoints2);
		    typedef std::size_t SizeType;
		    
		    static const unsigned int Dimension = 3;

		    typedef IntegrationPoint<3> IntegrationPointType;

		    typedef boost::array<IntegrationPointType, 6> IntegrationPointsArrayType;

		    typedef IntegrationPointType::PointType PointType;

		    static SizeType IntegrationPointsNumber()    {  return 6; }
		    
		    static IntegrationPointsArrayType& IntegrationPoints()
		    {
			    return msIntegrationPoints;
		    }

		    std::string Info() const
		    {
			    std::stringstream buffer;
			    buffer << "Prism gaussian quadrature 2 ";
			    return buffer.str();
		    }
	    protected:

	    private:

		    static IntegrationPointsArrayType msIntegrationPoints;

    }; // Class PrismGaussianIntegrationPoints2
    
    class PrismGaussianIntegrationPoints3
    {
        public:
            KRATOS_CLASS_POINTER_DEFINITION(PrismGaussianIntegrationPoints3);
            typedef std::size_t SizeType;
		    
            static const unsigned int Dimension = 3;

            typedef IntegrationPoint<3> IntegrationPointType;

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
                buffer << "Prism gaussian quadrature 3 ";
                return buffer.str();
            }
        protected:

        private:

            static IntegrationPointsArrayType msIntegrationPoints;

    }; // Class PrismGaussianIntegrationPoints3

 
   
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_PRISM_GAUSSIAN_INTEGRATION_POINTS_H_INCLUDED  defined 


