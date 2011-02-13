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
//   Date:                $Date: 2007-03-06 10:30:34 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_GEOMETRY_UTILITIES_INCLUDED )
#define  KRATOS_GEOMETRY_UTILITIES_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// External includes 


// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"


namespace Kratos
{
        ///this function provides basic routines for working with simplicial meshes. 
        ///It is faster than using Geometry as it is more specialized
	class GeometryUtils
	{
	public:

		/**this function is designed to compute the shape function derivatives, shape functions and volume in 3D
		 * @param geom it is the array of nodes. It is expected to be a tetrahedra
		 * @param a stack matrix of size 4*3 to store the shape function's derivatives
		 * @param an array_1d to store the shape functions at the barycenter
		 * @param the volume of the element
		 */
		static inline void CalculateGeometryData(
			Element::GeometryType& geom,
			boost::numeric::ublas::bounded_matrix<double,4,3>& DN_DX, 
			array_1d<double,4>& N, 
			double& Volume)
		{
			double x10 = geom[1].X() - geom[0].X();
			double y10 = geom[1].Y() - geom[0].Y();
			double z10 = geom[1].Z() - geom[0].Z();

			double x20 = geom[2].X() - geom[0].X();
			double y20 = geom[2].Y() - geom[0].Y();
			double z20 = geom[2].Z() - geom[0].Z();

			double x30 = geom[3].X() - geom[0].X();
			double y30 = geom[3].Y() - geom[0].Y();
			double z30 = geom[3].Z() - geom[0].Z();

			double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;

			DN_DX(0,0) = -y20 * z30 + y30 * z20 + y10 * z30 - z10 * y30 - y10 * z20 + z10 * y20;
			DN_DX(0,1) = -z20 * x30 + x20 * z30 - x10 * z30 + z10 * x30 + x10 * z20 - z10 * x20;
			DN_DX(0,2) = -x20 * y30 + y20 * x30 + x10 * y30 - y10 * x30 - x10 * y20 + y10 * x20;
			DN_DX(1,0) = y20 * z30 - y30 * z20;
			DN_DX(1,1) = z20 * x30 - x20 * z30;
			DN_DX(1,2) = x20 * y30 - y20 * x30;
			DN_DX(2,0) = -y10 * z30 + z10 * y30;
			DN_DX(2,1) = x10 * z30 - z10 * x30;
			DN_DX(2,2) = -x10 * y30 + y10 * x30;
			DN_DX(3,0) = y10 * z20 - z10 * y20;
			DN_DX(3,1) = -x10 * z20 + z10 * x20;
			DN_DX(3,2) = x10 * y20 - y10 * x20;

			DN_DX /= detJ;

			N[0] = 0.25;
			N[1] = 0.25;
			N[2] = 0.25;
			N[3] = 0.25;

			Volume = detJ*0.1666666666666666666667;
		}

		/**this function computes the element's volume (with sign)
		 * @param geom it is the array of nodes. It expects a tetrahedra
		 */
		static inline double CalculateVolume3D(
			Element::GeometryType& geom)
		{
			double x10 = geom[1].X() - geom[0].X();
			double y10 = geom[1].Y() - geom[0].Y();
			double z10 = geom[1].Z() - geom[0].Z();

			double x20 = geom[2].X() - geom[0].X();
			double y20 = geom[2].Y() - geom[0].Y();
			double z20 = geom[2].Z() - geom[0].Z();

			double x30 = geom[3].X() - geom[0].X();
			double y30 = geom[3].Y() - geom[0].Y();
			double z30 = geom[3].Z() - geom[0].Z();

			double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
			return  detJ*0.1666666666666666666667;
		}

		//********************************************************************************
		//********************************************************************************
		/**this function is designed to compute the shape function derivatives, shape functions and volume in 3D
		 * @param geom it is the array of nodes. It is expected to be a triangle
		 * @param a stack matrix of size 3*2 to store the shape function's derivatives
		 * @param an array_1d to store the shape functions at the barycenter
		 * @param the volume of the element
		 */
		static inline void CalculateGeometryData(
			Element::GeometryType& geom,
			boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, 
			array_1d<double,3>& N, 
			double& Area)
		{
		double x10 = geom[1].X() - geom[0].X();
		double y10 = geom[1].Y() - geom[0].Y();
		
		double x20 = geom[2].X() - geom[0].X();
		double y20 = geom[2].Y() - geom[0].Y();
		
		//Jacobian is calculated:
		//  |dx/dxi  dx/deta|	|x1-x0   x2-x0|
		//J=|				|=	|			  |
		//  |dy/dxi  dy/deta|	|y1-y0   y2-y0|


		double detJ = x10 * y20-y10 * x20;

		DN_DX(0,0) = -y20 + y10; DN_DX(0,1) = x20 - x10;
		DN_DX(1,0) =  y20	   ; DN_DX(1,1) = -x20     ;
		DN_DX(2,0) = -y10	   ; DN_DX(2,1) = x10	   ;

		DN_DX /= detJ;
		N[0] = 0.333333333333333;
		N[1] = 0.333333333333333;
		N[2] = 0.333333333333333;
		
		Area = 0.5*detJ;
		}


		//********************************************************************************
		//********************************************************************************
		/**this function computes the element's volume (with sign)
		 * @param geom it is the array of nodes. It expects a triangle
		 */
		static inline double CalculateVolume2D(
			Element::GeometryType& geom)
		{
			double x10 = geom[1].X() - geom[0].X();
			double y10 = geom[1].Y() - geom[0].Y();
			
			double x20 = geom[2].X() - geom[0].X();
			double y20 = geom[2].Y() - geom[0].Y();
			
			double detJ = x10 * y20-y10 * x20;
			return 0.5*detJ;
		}

		//********************************************************************************
		//********************************************************************************
		/** this function compute the maximum and minimum edge lenghts */
		static inline void SideLenghts2D(
			Element::GeometryType& geom,
			double& hmin, double& hmax)
		{
			double x10 = geom[1].X() - geom[0].X();
			double y10 = geom[1].Y() - geom[0].Y();
			
			double x20 = geom[2].X() - geom[0].X();
			double y20 = geom[2].Y() - geom[0].Y();

			double l = x20*x20 + y20*y20;
			hmax = l;
			hmin = l;

			if(l>hmax) hmax = l;
			else if(l<hmin) hmin = l;

			l = (x20-x10)*(x20-x10) + (y20-y10)*(y20-y10);
			if(l>hmax) hmax = l;
			else if(l<hmin) hmin = l;

			hmax = sqrt(hmax);
			hmin = sqrt(hmin);
		}




	};

}  // namespace Kratos.

#endif // KRATOS_GEOMETRY_UTILITIES_INCLUDED  defined 


