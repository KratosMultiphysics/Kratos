/*
==============================================================================
KratosTestApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2010
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

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
//   Date:                $Date: 2007-03-06 10:30:31 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_COMBUSTION_UTILITIES_INCLUDED )
#define  KRATOS_COMBUSTION_UTILITIES_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
#include "geometries/tetrahedra_3d_4.h"
#include "incompressible_fluid_application.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/timer.h"
//#include "processes/node_erase_process.h"
//#include "utilities/binbased_fast_point_locator.h"

namespace Kratos
{


class CombustionUtilities
{

public:

	CombustionUtilities()
	{}

	/// Destructor
	~CombustionUtilities() {}

	//**********************************************************************************************
	//**********************************************************************************************

	void Mixture_Fraction(ModelPart& full_model_part)
	{
		KRATOS_TRY;

		//double yo,ypr, yf, Z;
		double r=4.0;
		double p=79.0/6.0;
		double q1=11.0/4.0;
		double q2=9.0/4.0;
		double gamma1=1.0;
		double gamma2=-1.0/(r+p);

		double ych4, yo2, yco2, yh2o, yn2, Z;

		for(ModelPart::NodesContainerType::iterator in = full_model_part.NodesBegin() ; in != full_model_part.NodesEnd() ; ++in)
		{
			Z =  in->FastGetSolutionStepValue(MIXTURE_FRACTION);
			if(Z>=0.055)
			{

				ych4=gamma2+(gamma1-gamma2)*Z;
				yo2=0.0;
				yco2= q1 * (r+1.0) / ((q1 + q2)* (r + p +1.0)) * ((1.0-gamma2)+(gamma2-gamma1)*Z);
				yh2o= q2 * (r+1.0) / ((q1+q2)*(r + p +1.0)) * ((1.0- gamma2) + (gamma2 - gamma1)*Z);
				yn2= p /(r + p +1.0) * ((1.0-gamma2)+(gamma2-gamma1)*Z);

				if(ych4>1.0)ych4=1.0;
				if(ych4<0.0)ych4=0.0;

				if(yo2>1.0)yo2=1.0;
				if(yo2<0.0)yo2=0.0;

				if(yco2>1.0)yco2=1.0;
				if(yco2<0.0)yco2=0.0;

				if(yh2o>1.0)yh2o=1.0;
				if(yh2o<0.0)yh2o=0.0;

				if(yn2>1.0)yn2=1.0;
				if(yn2<0.0)yn2=0.0;

				in->FastGetSolutionStepValue(YCH4)=ych4;
				in->FastGetSolutionStepValue(YO2)=yo2;
				in->FastGetSolutionStepValue(YCO2)=yco2;
				in->FastGetSolutionStepValue(YH2O)=yh2o;
				in->FastGetSolutionStepValue(YN2)=yn2;
			}
			else
			{

				ych4=0.0;
				yo2=-r * (gamma2+(gamma1-gamma2)*Z);
				yco2= q1 * (r+1.0) / ((q1 + q2)* (r + p +1.0)) * (1.0 + (r + p)*(gamma2+(gamma1-gamma2)*Z));
				yh2o= q2 / (q1 + q2) * (r+1.0) / (r + p +1.0) * (1.0+(r + p)* (gamma2 + (gamma1-gamma2)*Z));
				yn2= p /(r + p + 1.0) * ((1.0-gamma2)+(gamma2-gamma1)*Z);

				if(ych4>1.0) ych4=1.0;
				if(ych4<0.0) ych4=0.0;

				if(yo2>1.0) yo2=1.0;
				if(yo2<0.0) yo2=0.0;

				if(yco2>1.0) yco2=1.0;
				if(yco2<0.0) yco2=0.0;

				if(yh2o>1.0) yh2o=1.0;
				if(yh2o<0.0) yh2o=0.0;

				if(yn2>1.0) yn2=1.0;
				if(yn2<0.0) yn2=0.0;

				in->FastGetSolutionStepValue(YCH4)=ych4;
				in->FastGetSolutionStepValue(YO2)=yo2;
				in->FastGetSolutionStepValue(YCO2)=yco2;
				in->FastGetSolutionStepValue(YH2O)=yh2o;
				in->FastGetSolutionStepValue(YN2)=yn2;
			}
		}
		KRATOS_CATCH("");

	}




	void Enthalpy(ModelPart& full_model_part)
	{

		KRATOS_TRY;
//	  double c1, c2, c3, c4, c5/*, c*/, b ,a/*,h*/;
		double c1;

		double /*T, yf,*/ ych4; //, yo2, yco2, yh2o, yn2;


		for(ModelPart::NodesContainerType::iterator in = full_model_part.NodesBegin() ; in != full_model_part.NodesEnd() ; ++in)
		{
			if(in->FastGetSolutionStepValue(FLAG_VARIABLE)==1.0 || in->FastGetSolutionStepValue(FLAG_VARIABLE)==2.0)
			{

//	      T=in->FastGetSolutionStepValue(TEMPERATURE);
//	      yf=in->FastGetSolutionStepValue(Yfuel);
				ych4=in->FastGetSolutionStepValue(YCH4);
//				yo2=in->FastGetSolutionStepValue(YO2);
//				yco2=in->FastGetSolutionStepValue(YCO2);
//				yh2o=in->FastGetSolutionStepValue(YH2O);
//				yn2=in->FastGetSolutionStepValue(YN2);
//	      h=in->FastGetSolutionStepValue(ENTHALPY);

				c1= (-7.485 * 10000.0 - 20.5 * 298.0 - 5.25 * (0.01) * 298.0 *298.0/ 2.0) * ych4;
//	      c2= (0.0 -29.44 * 298.0 - 5.05 * (0.001) * 298.0 *298.0/2.0) * yo2;
//	      c3= (0.0 -29.68 * 298.0 - 3.47 * (0.001) * 298.0 *298.0 /2.0) * yn2;
//	      c4= (-3.935 * 100000.0 -28.60 * 298.0 - 2.88 * (0.01) * 298.0 *298.0 /2.0) * yco2;
//	      c5=(-2.418*100000.0 - 33.84 * 298.0 - 6.48*0.001 *298.0*298.0/2.0)* yh2o;
//
//	      c= c1 + c2 + c3 + c4 + c5;

//	      b= 20.05 * ych4 * 298.0  + 29.44 * yo2  * 298.0 + 29.68 * yn2 * 298.0 + 28.60 * yco2 * 298.0 + 33.84 * yh2o * 298.0;
//
//	      a= 5.25 * 0.01 * ych4 * 298.0 * 298.0 * 0.5 + 5.05 * 0.001 * yo2 * 298.0 * 298.0 * 0.5 + 3.47 * 0.001 * yn2 * 298.0 * 298.0 * 0.5 + 2.88 * 0.01 * yco2 * 298.0 * 298.0 *0.5 + 6.48 * 0.001 * yh2o * 298.0 * 298.0 * 0.5;

//	      h= c1 + c2 + c3 + c4 + c5+ b * 298.0 + a * 298.0 * 298.0;

				in->FastGetSolutionStepValue(ENTHALPY)=c1;
				(in)->Fix(ENTHALPY);
			}
		}

		KRATOS_CATCH("");
	}


	void Temperature(ModelPart& full_model_part)
	{
		KRATOS_TRY;
		double ych4,yo2,yco2,yh2o,yn2,h,c1,c2,c3, c4,c5,c,b,a,TEMP;

		for(ModelPart::NodesContainerType::iterator in = full_model_part.NodesBegin() ; in != full_model_part.NodesEnd() ; ++in)
		{

			if(in->FastGetSolutionStepValue(FLAG_VARIABLE)==0.0)
			{
				ych4=in->FastGetSolutionStepValue(YCH4);
				yo2=in->FastGetSolutionStepValue(YO2);
				yco2=in->FastGetSolutionStepValue(YCO2);
				yh2o=in->FastGetSolutionStepValue(YH2O);
				yn2=in->FastGetSolutionStepValue(YN2);
				h=in->FastGetSolutionStepValue(ENTHALPY);

				c1= (-7.485 * 10000.0 - 20.5 * 298.0 - 5.25 * (0.01) * 298.0 *298.0/ 2.0) * ych4;
				c2= (0.0 -29.44 * 298.0 - 5.05 * (0.001) * 298.0 *298.0/2.0) * yo2;
				c3= (0.0 -29.68 * 298.0 - 3.47 * (0.001) * 298.0 *298.0 /2.0) * yn2;
				c4= (-3.935 * 100000.0 -28.60 * 298.0 - 2.88 * (0.01) * 298.0 *298.0 /2.0) * yco2;
				c5=(-2.418*100000.0 - 33.84 * 298.0 - 6.48*0.001 *298.0*298.0/2.0)* yh2o;

				c= c1 + c2 + c3 + c4 + c5;

				b= 20.05 * ych4  + 29.44 * yo2 + 29.68 * yn2 + 28.60 * yco2 + 33.87 * yh2o;

				a= 5.25 * (0.01) * ych4 / 2.0 + 5.05 * (0.001) * yo2 /2.0  + 3.47 * (0.001) * yn2 / 2.0+ 2.88 * (0.01) * yco2 + 6.48 * 0.001 * yh2o;

				TEMP=(-b + sqrt(b*b - 4.0 * a * (c-h)))/(2.0*a);
				in->FastGetSolutionStepValue(TEMPERATURE)=TEMP;
			}
		}
		KRATOS_CATCH("");

	}
};


} // namespace Kratos.

#endif // KRATOS_LAGRANGIAN_PARTICLES_UTILITIES_INCLUDED  defined
