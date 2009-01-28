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
//   Date:                $Date: 2008-11-17 14:50:46 $
//   Revision:            $Revision: 1.6 $
//
//


#if !defined(KRATOS_TETRAHEDRA_GAUSSIAN_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_TETRAHEDRA_GAUSSIAN_INTEGRATION_POINTS_H_INCLUDED



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

    class TetrahedraGaussianIntegrationPoints1
    {
    public:
      	KRATOS_CLASS_POINTER_DEFINITION(TetrahedraGaussianIntegrationPoints1);
  	typedef std::size_t SizeType;
	
	static const unsigned int Dimension = 3;

	typedef IntegrationPoint<3> IntegrationPointType;

      	typedef boost::array<IntegrationPointType, 1> IntegrationPointsArrayType;

      	typedef IntegrationPointType::PointType PointType;

      	static SizeType IntegrationPointsNumber(){  return 1; }
	
      	static IntegrationPointsArrayType& IntegrationPoints()
	{
	  // This is added to solve the problem of static initialization. Pooyan.
	  msIntegrationPoints[0] = IntegrationPointType( 0.25,0.25,0.25 , 1.00 / 6.00 );
	  return msIntegrationPoints;
	}

     	std::string Info() const
      	{
	  std::stringstream buffer;
	  buffer << "Tetrahedra gaussian quadrature 1 ";
	  return buffer.str();
      	}
    protected:

    private:

    	static IntegrationPointsArrayType msIntegrationPoints;

    }; // Class TetrahedraGaussianIntegrationPoints1






    class TetrahedraGaussianIntegrationPoints2
    {
	    public:
		    KRATOS_CLASS_POINTER_DEFINITION(TetrahedraGaussianIntegrationPoints2);
		    typedef std::size_t SizeType;
		    
		    static const unsigned int Dimension = 3;

		    typedef IntegrationPoint<3> IntegrationPointType;

		    typedef boost::array<IntegrationPointType, 4> IntegrationPointsArrayType;

		    typedef IntegrationPointType::PointType PointType;

		    static SizeType IntegrationPointsNumber()    {  return 4; }
		    
		    static IntegrationPointsArrayType& IntegrationPoints()
		    {
		      // This is added to solve the problem of static initialization. Pooyan.
		      msIntegrationPoints[0] = IntegrationPointType( 0.58541020,0.13819660,0.13819660 , 1.00 / 24.00 );
		      msIntegrationPoints[1] = IntegrationPointType( 0.13819660,0.58541020,0.13819660 , 1.00 / 24.00 );
		      msIntegrationPoints[2] = IntegrationPointType( 0.13819660,0.13819660,0.58541020 , 1.00 / 24.00 );
		      msIntegrationPoints[3] = IntegrationPointType( 0.13819660,0.13819660,0.13819660 , 1.00 / 24.00 );
			    return msIntegrationPoints;
		    }

		    std::string Info() const
		    {
			    std::stringstream buffer;
			    buffer << "Tetrahedra gaussian quadrature 2 ";
			    return buffer.str();
		    }
	    protected:

	    private:

		    static IntegrationPointsArrayType msIntegrationPoints;

    }; // Class TetrahedraGaussianIntegrationPoints2





   class TetrahedraGaussianIntegrationPoints3
   {
	   public:
		   KRATOS_CLASS_POINTER_DEFINITION(TetrahedraGaussianIntegrationPoints3);
		   typedef std::size_t SizeType;
		   
		   static const unsigned int Dimension = 3;

		   typedef IntegrationPoint<3> IntegrationPointType;

		   typedef boost::array<IntegrationPointType, 5> IntegrationPointsArrayType;

		   typedef IntegrationPointType::PointType PointType;

		   static SizeType IntegrationPointsNumber()    {  return 5; }
		    
		   static IntegrationPointsArrayType& IntegrationPoints()
		   {
			   return msIntegrationPoints;
		   }

		   std::string Info() const
		   {
			   std::stringstream buffer;
			   buffer << "Tetrahedra gaussian quadrature 3 ";
			   return buffer.str();
		   }
	   protected:

	   private:

		   static IntegrationPointsArrayType msIntegrationPoints;

   }; // Class TetrahedraGaussianIntegrationPoints3
   
   class TetrahedraGaussianIntegrationPoints4
   {
	   public:
		   KRATOS_CLASS_POINTER_DEFINITION(TetrahedraGaussianIntegrationPoints4);
		   typedef std::size_t SizeType;
		   
		   static const unsigned int Dimension = 3;
		   
		   typedef IntegrationPoint<3> IntegrationPointType;
		   
		   typedef boost::array<IntegrationPointType, 10> IntegrationPointsArrayType;
		   
		   typedef IntegrationPointType::PointType PointType;
		   
		   static SizeType IntegrationPointsNumber()	{ return 10; }
		   
		   static IntegrationPointsArrayType& IntegrationPoints()
		   {
			   return msIntegrationPoints;
		   }
		   
		   std::string Info() const
		   {
			   std::stringstream buffer;
			   buffer << "Tetrahedra gaussian quadrature 4 ";
			   return buffer.str();
		   }
		   
	   protected:
		   
	   private:
		   
		   static IntegrationPointsArrayType msIntegrationPoints;
		
   };
   
   class TetrahedraGaussianIntegrationPoints5
   {
	   public:
		   KRATOS_CLASS_POINTER_DEFINITION(TetrahedraGaussianIntegrationPoints4);
		   typedef std::size_t SizeType;
		   
		   static const unsigned int Dimension = 3;
		   
		   typedef IntegrationPoint<3> IntegrationPointType;
		   
		   typedef boost::array<IntegrationPointType, 11> IntegrationPointsArrayType;
		   
		   typedef IntegrationPointType::PointType PointType;
		   
		   static SizeType IntegrationPointsNumber()	{ return 11; }
		   
		   static IntegrationPointsArrayType& IntegrationPoints()
		   {
			   return msIntegrationPoints;
		   }
		   
		   std::string Info() const
		   {
			   std::stringstream buffer;
			   buffer << "Tetrahedra gaussian quadrature 5 ";
			   return buffer.str();
		   }
		   
	   protected:
		   
	   private:
		   
		   static IntegrationPointsArrayType msIntegrationPoints;
		
   };
   
   
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_TETRAHEDRA_GAUSSIAN_INTEGRATION_POINTS_H_INCLUDED  defined 


